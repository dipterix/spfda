#' @export
plot_coef = function(coef, se = NA, ylim, beta, time){
  layout(matrix(c(1,1,1,2,3,4), nrow = 2, byrow = T), heights = c(lcm(1.5), 1))
  par(mar=c(0,0,0,0))
  plot.new()
  mtext('Coefficients with 95% Joint CI', 1, line = -1, cex = 2)
  # mtext('95% Joint CI, W', 1, line = -1, cex = 2)

  par(mar=c(3.1,2.5,3.1,1.1))
  # coef = res$get_coef(); zz = 2
  # se = res$raw$f_sd
  if(missing(time)){
    time = seq(0,1,length.out = ncol(coef))
  }

  if(missing(ylim)){
    yat = pretty(range(coef, coef + se, coef - se, na.rm = TRUE))
    ylim = range(yat)
  }else{
    yat = pretty(ylim)
  }


  for(ii in 1:3){
    # coef = res$get_coef()

    rutabaga::plot_clean(
      time, ylim,
      main = eval(parse(text = sprintf('bquote(beta[%d])', ii))), cex.main = 1.8)
    abline(h = 0, lty=1, col = 'grey80')
    if(!missing(beta)){
      points(time, beta[ii,], type='l', lty = 2, col = 'grey80')
    }

    # points(time, flm$coefficients[ii,], type='l', lty = 2, col = 'grey80')
    rutabaga::ebar_polygon(time, coef[ii,], se[ii, ] * 2 , col = ii, lwd = 2, border = 'black', alpha = 50)

    rutabaga::ruta_axis(1, pretty(time))
    rutabaga::ruta_axis(2, yat)
  }
}

#' @export
simulate_data <- function(n = 1000, n_timepoints = 100, seed = 1){

  n_coef = 3
  beta1 = function(t){0}
  beta2 = function(t){sin(pi * t)}
  beta3 = function(t){
    if(t < 0.2 || t >= 0.8){
      return(0)
    }else if (t >= 0.4 && t < 0.6){
      return(1)
    }else if (t < 0.4){
      return(sin((5 * t - 1) * pi / 2))
    }else{
      return(sin((5 * t - 2) * pi / 2))
    }
  }

  # Generate discrete coefficients
  time = seq(0, 1, length.out = n_timepoints)
  beta = t(sapply(1:n_coef, function(ii){
    f = get(paste0('beta', ii))
    sapply(time, f)
  }))

  set.seed(seed)
  on.exit(set.seed(-1), add = TRUE, after = TRUE)
  # X: n x 3 rnorm data
  X = rnorm(n * n_coef)
  dim(X) = c(n, n_coef)

  # have some correlation
  junk = sample(n * n_coef, n*2)
  X[junk] = X[junk] + 1
  base_curve = ((time > 0.8)*3 + (time > 0.4)*2 + (time >= 0)) / 10

  err = t(replicate(n, {
    arima.sim(model=list(ar=c(.9)),n=n_timepoints) * base_curve
  }))

  Y = X %*% beta + err
  Y = Y + rnorm(length(Y))


  plot_input_data <- function(){
    # Plot coef
    layout(matrix(c(1,1,2), nrow = 1))
    par('mar' = c(5.1, 4.5, 4.1, 2.1))
    rutabaga::plot_clean(time, beta, xlab = 'Time', ylab = expression(beta),
                         main = bquote('Underlying '~beta['i, i=1,2,3']))
    matpoints(time, t(beta), type = 'l', lty=1, lwd = 3)

    rutabaga::ruta_axis(1, pretty(time))
    rutabaga::ruta_axis(2, pretty(beta))
    legend('topleft', c(
      expression(beta[1]),expression(beta[2]),expression(beta[3])
    ), lty=1,lwd=3,col=1:3, bty = 'n')

    rutabaga::plot_clean(time, c(0,base_curve), xlab = 'Time', main = expression(sigma(t)),
                         ylab = '')
    points(time[1:40], base_curve[1:40], type='S', lty=1, lwd=3)
    points(time[41:80], base_curve[41:80], type='S', lty=1, lwd=3)
    points(time[81:100], base_curve[81:100], type='S', lty=1, lwd=3)
    rutabaga::ruta_axis(1, pretty(time))
    rutabaga::ruta_axis(2, pretty(c(0,base_curve)))
  }


  return(environment())

}

#' @export
get_weights = function(X, Y, b = 0.1, parts = list(c(0, 0.4), c(0.4,0.8), c(0.8,1))){
  flm = lm(Y ~ X)
  E = flm$residuals
  time = seq(0,1,length.out = ncol(Y))
  # ols_res = spfda(Y = Y, X = X, time = time, nsp =K, ord = 4, lambda = 0, alpha = 1, W = NULL)
  # E = Y - res$predict(X)

  rands = matrix(0, nrow = nrow(X), ncol = length(time))
  sigma_t_sel = matrix(FALSE, nrow = length(time), ncol = length(time))

  for(ii in seq_along(parts)){
    part = parts[[ii]]
    if(ii == length(parts)){
      sel_pa = time >= part[1]
    }else{
      sel_pa = time >= part[1] & time < part[2]
    }
    sigma_t_sel[sel_pa,sel_pa] = TRUE
    theta = sapply(which(sel_pa), function(m){
      sel_pa[m] = FALSE
      delta_t = dnorm((time - time[m]) /b) / b
      tmt = diag(delta_t[sel_pa])
      zz = cbind(1, (time - time[m])[sel_pa] / b)
      ee = E[, sel_pa]
      rands = solve(t(zz) %*% tmt %*% zz) %*% t(zz) %*% tmt %*% t(ee)
      rands[1,]
    })
    rands[, sel_pa] = theta
  }
  # matplot(t(rands), type = 'l', col = 'grey80')

  sigma_e = apply(E-rands, 2, var) #colMeans((E-rands)^2)
  # sigma_e = var(as.vector(E-rands))
  sigma_t = cov(rands)
  sigma_t[!sigma_t_sel] = 0

  # assign('sigma_t', sigma_t, envir = globalenv())

  # plot(diag(sigma_t))
  # plot(apply(sigma_t, 2, var))
  # lines(sigma_e)

  eig = eigen(sigma_t)
  eig$values[eig$values < 0] = 0
  sigma_t = eig$vectors %*% diag(eig$values) %*% t(eig$vectors)

  diag(sigma_t) = sigma_e + diag(sigma_t)

  eig = eigen(sigma_t)
  eig$values[eig$values < 1e-4] = 1e-4
  eig$values = eig$values + max(eig$values) / 10

  dd = 1 / (sqrt(eig$values))
  dd = dd/max(dd)
  W = eig$vectors %*% diag(dd) %*% t(eig$vectors)
  W
}

# calculate log-likelihood
#' @export
logLik.spfda.model <- function (object, ...) {
  res = object$Y - object$predict(object$X)
  W = object$W
  p = ncol(object$X)
  N = nrow(res)
  nT = ncol(res)
  K = object$K
  coef = object$get_coef()

  if(!is.matrix(W)){
    det_w = 1
  }else{
    det_w = det(W)
    res = res %*% W
  }

  # assume each row of res is uncor
  val <- sum(apply(res, 1, function(x){
    - sum(x ^ 2) / 2 + nT * (log(det_w) - log(2 * pi))
  }))

  # calculate df

  df = sum(coef != 0)
  df_alt = sum(abs(coef) > 0.001 * max(abs(coef)))

  attr(val, "nall") <- N
  attr(val, "nobs") <- N
  attr(val, "df") <- df
  attr(val, "df_alt") <- df_alt
  attr(val, "pK") <- p * K
  class(val) <- "logLik"
  val
}


#' @export
BIC.spfda.model <- function(object, nu = 0.5, ...){
  loglik = logLik.spfda.model(object, ...)
  df = attr(loglik, "df")
  nobs = attr(loglik, "nobs")
  pK = attr(loglik, "pK")
  (-2 * loglik + df * log(nobs) + nu * df * log(pK)) / nobs
}
