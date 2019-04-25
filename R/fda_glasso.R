fda_glasso <- function(Y, X, time = seq(0, 1, length.out = ncol(Y)),
                       nknots = 100, lambda = 0.1,
                       W = NULL, D = NULL, init = NULL, ord = 4,
                       max_iter = 50, ..., ridge_factor = 1){
  alpha = 1

  # Initialize params
  n = nrow(Y)
  p = ncol(X)
  n_timepoints = length(time)
  W %?<-% diag(1, n_timepoints)

  # B = t(ns(x = time, df = nknots))
  # class(B) = "matrix"
  knots = c(rep(time[1], ord-1), seq(time[1], time[length(time)], length.out = nknots - ord), rep(time[length(time)], ord-1))
  B = t(
    # ns(x = time, df = nknots)
    splineDesign(knots, time, ord = ord)
  )
  class(B) = "matrix"
  IB = B > 0

  # ------------ A1 ------------
  xtx = crossprod(X, X)
  ident_p = diag(1, p)

  xty = crossprod(X, Y)
  wwt = tcrossprod(W, W)
  bw = B %*% W
  bwwb = tcrossprod(bw, bw)

  init %?<-% (solve(xtx + ridge_factor * ident_p) %*% xty %*% wwt %*% t(B) %*% solve(bwwb))
  gamma = init

  # ------------ A2 ------------
  svd_x = svd(X)
  svd_bw = svd(bw)


  # ------------ A3 ------------
  V = t(svd_x$v)
  U = svd_bw$u
  Z0 = as.matrix((svd_x$d)^2) %*% ((svd_bw$d)^2)
  rhos = seq(0.1, 60, length.out = max_iter) * n

  # ------------ B4 ------------
  D = matrix(1, nrow = p, ncol = n_timepoints) %*% t(IB)

  # ------------ B5 ------------
  eta = gamma; theta = array(0, dim(eta))
  for(ii in seq_len(max_iter)){
    rho = rhos[ii]
    Z = Z0 + rho

    gamma = t(V) %*% ((V %*% xty %*% wwt %*% t(B) %*% U - theta + rho * V %*% eta %*% U) / Z) %*% t(U)

    eta = gamma + 1/rho * t(V) %*% theta %*% t(U)
    eta_plus = eta - lambda/rho * D
    eta_minus = eta + lambda/rho * D
    eta = eta_plus * (eta_plus > 0) + eta_minus * (eta_minus < 0)

    theta = theta + V %*% (gamma - eta) %*% U
  }

  list(
    gamma = gamma,
    eta = eta,
    B = B,
    mse = (sqrt(mean((Y - X %*% (eta %*% B)) ^ 2)))
  )
}


fda_glasso_cb <- function(Y, X, ..., nsamp = 1000){
  n = nrow(Y)
  args = list(...)
  replicate(nsamp, {
    idx = sample(n, n, replace = T)
    Y_sub = Y[idx, ]
    X_sub = X[idx, ]
    re = do.call(fda_glasso, c(list(Y = Y_sub, X = X_sub), args))
    re$eta
  })
}
