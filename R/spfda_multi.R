make_block <- function(x, mat){
  n_factor = ncol(mat)
  do.call(
    Matrix::cBind,
    lapply(seq_len(n_factor), function(ii) {
      do.call(
        Matrix::rBind,
        lapply(mat[, ii], function(v) {
          v * x
        }
        ))
    }
    ))
}


#' @importFrom splines splineDesign
#' @import Matrix
#' @export
spfda_multi <- function(X, Y, design, nknots = ncol(X) / 4, weights = 1,
                        alpha = 0.01, lambda = 2, kappa = 0, niter = 50,
                        collapse_method = mean, method = 'admm'){
  re = list()

  # Vars
  ntime_points = ncol(X)
  time = seq_len(ntime_points)


  Y = as.factor(Y)
  re$Y = Y
  n_signals = nrow(X)

  y_levels = unique(Y)

  y_count = sapply(y_levels, function(fa){sum(Y == fa)})
  weights = weights / mean(weights)
  weights = weights * rep(y_count / n_signals, each = ntime_points)


  n_factor = dim(design)[2]
  if(
    n_factor != length(y_levels)
  ){
    stop('Number of columns of design matrix should match with Y factors')
  }

  # Generate knots and Bsplines
  knots = seq(min(time), max(time), length.out = nknots)
  nknots = length(knots)
  B = Matrix(splineDesign(knots, x = time, outer.ok = T))
  B0 = Matrix(0, nrow = nrow(B), ncol = ncol(B))
  nsplines = ncol(B)

  if(length(alpha) == 1){
    alpha = rep(alpha, n_factor)
  }
  if(length(alpha) == n_factor){
    alpha = rep(alpha, each = ntime_points)
  }
  re$alpha = matrix(alpha, ncol = n_factor)[1,]

  if(length(lambda) == 1){
    lambda = rep(lambda, n_factor)
  }
  if(length(lambda) == n_factor){
    lambda = rep(lambda, each = ntime_points)
  }
  re$lambda = matrix(lambda, ncol = n_factor)[1,]
  lambda = lambda / nknots

  if(length(kappa) == 1){
    kappa = rep(kappa, n_factor)
  }
  if(length(kappa) == n_factor){
    kappa = rep(kappa, each = nsplines)
  }
  re$kappa = matrix(kappa, ncol = nsplines)[1,]





  # Init x, y, z for ADMM

  BB = make_block(B, design)
  I = Matrix::Diagonal(nsplines)

  A = do.call(Matrix::bdiag, replicate(n_factor, I, simplify = F)) # For simple X-Z

  S = - diag(1, ncol(B))
  lag = 1 #ceiling(ntime_points / nknots)

  row = seq_len(nrow(S) - lag)
  col = row + lag
  S[row + (col-1) * nrow(S)] = 1
  S = t(S) %*% S


  BTB = t(BB) %*% (weights * BB)
  BSB = do.call(Matrix::bdiag, replicate(n_factor, S, simplify = F))#kappa * t(S) %*% (S)
  ATA = t(A) %*% A
  BB1 = t(BB != 0)




  # Set rho for augmented lagrangian
  # start = log(quantile(abs(X), 0.1) / 10)
  # end = log(max(10 * range(X)^2, 100))
  # rhos = exp(seq(start, end, length.out = niter * 2 / 3))
  # rhos = c(rhos, rep(end, niter/3))
  v = eigen(BTB)$value
  step = min(v[v > 0]) / nsplines
  rhos = seq(step, by = step, length.out = niter)
  # Warm start: use ridge
  {
    rho = rhos[1]

    x_centers = lapply(y_levels, function(fa){
      sel = (Y == fa)
      apply(X[sel,], 2, collapse_method)
    })
    x_collapsed = do.call(c, x_centers)




    x = solve(BTB + 0.1 * A) %*% t(BB) %*% x_collapsed
    z = x
    y = rho * (A%*%x - z)

    mask = z[] * 0
  }

  for(rho in rhos){
    # ADMM Step 1: fix z,y, optimize x:
    # x = solve(BTB + rho * ATA + kappa * BSB) %*% (t(BB) %*% rowMeans(Y) + t(A) %*% (rho * z - y))
    if(method == 'admm'){
      x_new = solve(BTB + (BSB + rho * ATA) ) %*% (t(weights * BB) %*% x_collapsed + t(A) %*% (rho * z - y))
    }else{
      x_new = solve(BTB + (BSB) ) %*% (t(weights * BB) %*% x_collapsed + t(A) %*% (- y))
    }
    # x = (x + x_new) / 2
    x = x_new

    # x[abs(x) < 10^-5] = 0
    # ADMM Step 2: fix x, y, optimize z:

    ###### Method 1:
    tau = lambda^(1/(1-alpha)) * alpha^(alpha/(1-alpha)) * (1-alpha); tau

    theta = ((1-alpha) / (alpha * tau))^alpha * (colSums(abs(BB1 * as.vector(z))) ^ alpha)
    ctheta = theta^(1 - 1 / alpha)
    nas = sum(!is.finite(ctheta)) ; nas
    if(nas > 0){
      ctheta[!is.finite(ctheta)] = 10000000
    }
    ctheta = (BB1 %*% ctheta)
    # Another ADMM ???
    z_minus = (y + ctheta) / rho + A %*% x
    z_plus = (y - ctheta) / rho + A %*% x
    # z_mid = A %*% x
    # z_plus > 0 & z_minus < 0

    # Soft thresholding
    z = (z_plus < 0) * z_minus + (z_minus > 0) * z_plus #+ (z_plus > 0 & z_minus < 0) * (z_plus + z_minus) / 10
    # z = (z + x) / 2
    mask = mask | (z_plus > 0 & z_minus < 0)
    z[mask] = 0

    ########## Method 2


    #0.5 * z + (1-0.5) *

    # ADMM Step 3: Update dual y:
    y = y + rho * (A %*% x - z)
    # print(rho)
    # message(rho)
    # print(sqrt(mean((x_collapsed - as.vector(BB %*% x))^2)))
    # print(range(x-z))
  }

  re$X = X
  re$multiplier = list(x=x,y=y,z=z,rho=rho)
  re$W = z
  re$err = t(matrix(x_collapsed - as.vector(BB %*% z), ncol = n_factor))
  re$time = time
  re$niter = niter
  re$B = B
  re$BB = BB
  re$comps = lapply(seq_len(n_factor), function(ii){
    B = BB[length(time) * (ii - 1) + seq_along(time), ]
    center = x_centers[[ii]]
    x_centered = t(X[Y == y_levels[ii], ] - center)
    list(
      B = B,
      recall = B %*% z,
      factor = re$B %*% z[nsplines * (ii - 1) + seq_len(nsplines),],
      center,
      cov = x_centered %*% t(x_centered) / y_count[ii],
      n = y_count[ii]
    )
  })
  re$.plot = T
  re$knots = knots

  class(re) = c('spfda_multi', 'spfda_base')
  return(re)

}


#' @export
plot.spfda_multi <- function(re, xlab = 'Time', ylab = '', std = 2,
                             pal = c("#8DA0CB", "#66C2A5", "#FFD92F", "#FC8D62", "#E78AC3", "#A6D854", "#E5C494", "#B3B3B3"), ...){
  factors = unique(re$Y)
  par(mfrow = c(length(factors), 2))

  N = nrow(re$X)

  sapply(re$comps, function(comp){
    as.vector(comp$cov * comp$n)
  }) ->
    tmp
  tmp = rowSums(tmp) / N
  dim(tmp) = c(sqrt(length(tmp)),sqrt(length(tmp)))
  A = solve(t(re$BB) %*% re$BB / N + diag(re$multiplier$rho, ncol(re$BB))) %*% t(re$BB / N)

  cov = do.call(Matrix::bdiag, replicate(length(factors), tmp, simplify = F))
  jv = A %*% cov %*% t(A)
  # djv = sqrt(diag(re$B %*% jv %*% t(re$B))) * 2

  nsplines = ncol(re$B)


  lapply(seq_along(factors), function(ii){
    ind = (ii-1) * nsplines + seq_len(nsplines)
    recall = re$comps[[ii]]$recall
    x = re$X[re$Y == factors[ii], , drop = F]
    x_mean = colMeans(x)
    x_quant = apply(x, 2, quantile, c(0.25,0.5,0.75))
    y_range = range(range(x_mean), range(recall), range(x_quant))
    plot(range(re$time), y_range, type = 'n', xlab = xlab, ylab = ylab, main = sprintf('Signal %s', as.character(factors[ii])), ...)
    points(re$time, x_mean, type = 'l', lwd = 1, col = pal[2])
    matpoints(re$time, t(x_quant), type = 'l', lty = 2, col = pal[3])
    points(re$time, recall, type = 'l', lwd = 2, col = pal[1])
    abline(v = re$knots, lty = 3, lwd = 0.3, col = 'grey80')

    comp_fa = re$comps[[ii]]$factor
    # Functional CI
    djv = sqrt(diag(re$B %*% jv[ind,ind] %*% t(re$B))) * std
    plot(range(re$time), c(max(djv+comp_fa), min(comp_fa - djv)), xlab = xlab, ylab = '', type = 'n', main = sprintf('Component %d', ii), ...)
    abline(h = 0, col = pal[8], lty = 4)
    points(re$time, comp_fa, col = pal[4], type = 'l')

    points(re$time, comp_fa + djv, type = 'l', col = pal[6], lty = 2)
    points(re$time, comp_fa - djv, type = 'l', col = pal[6], lty = 2)

    abline(v = re$knots, lty = 3, lwd = 0.3, col = 'grey80')
  })
}
