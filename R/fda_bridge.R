


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
  future.apply::future_replicate(nsamp, {
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
  max_iter = 50, CI = FALSE, ...){

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
  ident_p = diag(10, p)

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
    rho = rhos[ii] # * (lambda * ii / max_iter)
    Z = Z0 + rho

    # ------------ A4.1 ------------
    D = apply(abs(gamma), 1, function(g){
      tmp = (colSums(IB * g))
      tmp[tmp > 0] = tmp[tmp > 0]^(alpha - 1)
      #tmp[!is.finite(tmp)] = 0
      tmp %*% t(IB)
    })
    D = alpha * t(D)

    # ------------ A4.2 ------------
    eta = xi = gamma; theta = array(0, dim(gamma))
    for(m in 1:ifelse(max_iter == ii, 100, 5)){
      xi = t(V) %*% ((V %*% xty %*% wwt %*% t(B) %*% U - theta +
                       rho * V %*% eta %*% U) / Z) %*% t(U)
      eta = xi + 1/rho * t(V) %*% theta %*% t(U)
      eta_plus = eta - lambda/rho * D
      eta_minus = eta + lambda/rho * D
      eta = eta_plus * (eta_plus > 0) + eta_minus * (eta_minus < 0)

      theta = theta + rho * V %*% (xi - eta) %*% U


    }

    gamma = eta

    mse = (sqrt(mean((Y - X %*% (eta %*% B)) ^ 2)))
    if(last_mse < mse * 0.8){
      break
    }
  }

  re_gamma = gamma
  if(CI){

    tmp = gamma
    tmp[tmp == 0] = runif(sum(tmp == 0), min = -1e-7, max = 1e-7)
    DD = U %*% t(theta) %*% V
    DD = t(DD) * sign(gamma) * lambda / tmp
    # DD[is.infinite(DD)] = 0
    residual = Y - X %*% solve(t(X) %*% X) %*% t(X) %*% Y %*% t(B) %*% solve(B %*% t(B)) %*% B
    # residual = Y - X %*% gamma %*% B
    sigma = cov(residual)

    p = dim(X)[2]
    K = dim(gamma)[2]

    # for covariance
    s = abs(re_gamma) > 0

    idjk = which(abs(tmp) > 0, arr.ind = T)

    P = apply(idjk, 1, function(idx){
      j = idx[1]; k = idx[2]
      p = xtx[,j, drop = F] %*% bwwb[k, , drop = F]
      d = DD[j,k]
      if( is.finite(d) && d > 0 ){
        p[j,k] = p[j,k] + d
      }

      # print(DD[j,k] / gamma[j, k])
      # if(p[j,k] > 1000){
      #   print(idx)
      # }
      apply(idjk, 1, function(idy){
        p[idy[1], idy[2]]
      })
      # as.vector(t(p))[as.vector(t(s))]
    });

    # image.plot(P)
    # t(P) %*% ? %*% P = Q
    P_solve = solve(P)


    wwtbt = W %*% t(B %*% W)
    Q = apply(idjk, 1, function(idx1){
      j1 = idx1[1]; k1 = idx1[2]
      q1 = wwtbt[, k1, drop = F]
      apply(idjk, 1, function(idx2){
        j2 = idx2[1]; k2 = idx2[2]
        q2 = wwtbt[, k2, drop = F]

        t(q1) %*% sigma %*% q2 * (xtx[j1, j2])
      })
    })

    avar_gamma = t(P_solve) %*% Q %*% (P_solve)
    idxx = idjk[,2] + (idjk[,1]-1) * K
    avar_g = matrix(0, nrow = p*K, ncol = p*K)
    # avar_g[as.vector(t(s)),as.vector(t(s))] = avar_gamma
    for(ii in seq_along(idxx)){
      avar_g[idxx[ii], idxx] = avar_gamma[ii, ]
    }

    idjk = which(!s, arr.ind = T)
    idxx = idjk[,2] + (idjk[,1]-1) * K
    # avar_g[as.vector(t(s)),as.vector(t(s))] = avar_gamma
    # ss = mean(re_gamma == 0); a = 0
    # ss = a * ss / (1 + ss * (a-1)); ss = sqrt(1 - ss)
    # avar_g[idxx, ] = avar_g[idxx, ] * ss
    # avar_g[, idxx] = avar_g[, idxx] * ss

    # s1 = as.vector(t(abs(gamma) > 0))
    # vars = eigen(avar_g[s1, s1])$value
    # vars[vars < 0] = 0
    # min_eigen = min(vars[cumsum(vars) / sum(vars) > 0.99999][1], 1)
    #
    #
    # ee = eigen(avar_g)
    # ee$values = ee$values - min(ee$values) + min_eigen
    #
    # avar_g = (ee$vectors) %*% diag(ee$values) %*% t(ee$vectors)


    # get var for coef
    f_sd = t(sapply(seq_len(p), function(j){
      idx = seq_len(K) + (j-1) * K
      cov = t(B) %*% avar_g[idx, idx] %*% B

      # ee = eigen(cov)
      # vars = ee$value
      # d = which(cumsum(vars) / sum(vars) > 0.999)[1]
      #
      # vars[-(1:d)] = vars[d]
      # cov = (ee$vectors) %*% diag(vars) %*% t(ee$vectors)
      #
      #
      # tmp = mvtnorm::rmvnorm(1000, sigma = cov)
      #
      # # calculate simultanous CI
      # sds = apply(tmp, 2, sd)
      #
      # quantile(t(tmp) / sds, 0.975) * sds

      v = diag(cov)
      # plot(sqrt(diag(cov)))
      sqrt(v)
    }))

    # matplot(t(f_sd), type = 'l')


  }else{
    f_sd = avar_g = NULL
  }



  list(
    gamma = re_gamma,
    eta = eta,
    theta = theta,
    rho = rho,
    D = D,
    Vt = V,
    Ut = U,
    B = B,
    mse = (sqrt(mean((Y - X %*% (eta %*% B)) ^ 2))),
    niter = ii,
    avar_g = avar_g,
    f_sd = f_sd
  )
}
