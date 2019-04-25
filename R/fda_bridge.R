


fda_cv <- function(Y, X, ..., nfold = 5, method = 'fda_bridge'){
  n = nrow(Y)
  samp = sample(n, floor(n / nfold) * nfold, replace = F)
  dim(samp) = c(length(samp) / nfold, nfold)
  args = list(...)
  apply(samp, 2, function(test_idx){
    # train
    re = do.call(method, c(list(Y = Y[-test_idx, , drop = F], X = X[-test_idx, , drop = F]), args))
    # predict
    res = Y[test_idx, , drop = F] - X[test_idx, , drop = F] %*% re$eta %*% re$B
    sqrt(mean(res^2))
  })
}

fda_bridge_cb <- function(Y, X, ..., nsamp = 1000){
  n = nrow(Y)
  args = list(...)
  replicate(nsamp, {
    idx = sample(n, n, replace = T)
    Y_sub = Y[idx, ]
    X_sub = X[idx, ]
    re = do.call(fda_bridge, c(list(Y = Y_sub, X = X_sub), args))
    re$eta
  })
}

fda_bridge <- function(
  Y, X, time = seq(0, 1, length.out = ncol(Y)),
  nknots = 100, lambda = 0.1, alpha = 0.1,
  W = NULL, D = NULL, init = NULL, ord = 4,
  max_iter = 50, ...){

  # Initialize params
  n = nrow(Y)
  p = ncol(X)
  n_timepoints = length(time)
  W %?<-% diag(1, n_timepoints)

  time = sort(time)

  knots = c(rep(time[1], ord-1), seq(time[1], time[length(time)], length.out = nknots - ord), rep(time[length(time)], ord-1))
  B = t(
    # ns(x = time, df = nknots)
    splineDesign(knots, time, ord = ord)
  )
  class(B) = "matrix"
  # B = t(Matrix(array(B, dim = dim(B))))
  IB = B > 0

  # ------------ A1 ------------
  xtx = crossprod(X, X)
  ident_p = diag(1, p)

  xty = crossprod(X, Y)
  wwt = tcrossprod(W, W)
  bw = B %*% W
  bwwb = tcrossprod(bw, bw)

  init %?<-% (solve(xtx + ident_p) %*% xty %*% wwt %*% t(B) %*% solve(bwwb))
  gamma = init

  # ------------ A2 ------------
  svd_x = svd(X)
  svd_bw = svd(bw)


  # ------------ A3 ------------
  V = t(svd_x$v)
  U = svd_bw$u
  Z0 = as.matrix((svd_x$d)^2) %*% ((svd_bw$d)^2)
  rhos = seq(0.1, 60, length.out = max_iter) * n

  # ------------ A4 ------------
  last_mse = Inf
  for(ii in seq_len(max_iter)){
    rho = rhos[ii]
    Z = Z0 + rho

    # ------------ A4.1 ------------
    D = apply(abs(gamma), 1, function(g){
      tmp = (colSums(IB * g)^(alpha - 1)) * t(IB)
      tmp[!is.finite(tmp)] = 0
      colSums(tmp)
    })
    D = t(D)

    # ------------ A4.2 ------------
    eta = xi = gamma; theta = array(0, dim(gamma))
    for(m in 1:5){
      xi = t(V) %*% ((V %*% xty %*% wwt %*% t(B) %*% U - theta +
                       rho * V %*% eta %*% U) / Z) %*% t(U)
      eta = xi + 1/rho * t(V) %*% theta %*% t(U)
      eta_plus = eta - lambda/rho * D
      eta_minus = eta + lambda/rho * D
      eta = eta_plus * (eta_plus > 0) + eta_minus * (eta_minus < 0)

      theta = theta + V %*% (xi - eta) %*% U


    }

    gamma = xi

    mse = (sqrt(mean((Y - X %*% (eta %*% B)) ^ 2)))
    if(last_mse < mse * 0.8){
      break
    }
  }

  list(
    gamma = gamma,
    eta = eta,
    B = B,
    mse = (sqrt(mean((Y - X %*% (eta %*% B)) ^ 2))),
    niter = ii
  )
}
