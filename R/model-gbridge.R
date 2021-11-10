
# Model file

fos_gp_bridge <- function(
  X, Y, time = seq(0, 1, length.out = ncol(Y)), nknots = 100, lambda = 0.1, alpha = 0.5,
  W = NULL, D = NULL, init = NULL, ord = 4, max_iter = 50, inner_iter = 5, CI = FALSE, ...){

  # Initialize params
  n <- nrow(Y)
  p <- ncol(X)
  n_timepoints <- length(time)

  if(!length(W) || !is.matrix(W)){
    W <- diag(1, n_timepoints)
  }
  WT = t(W)

  time <- sort(time)

  # B-spline
  knots <- c(rep(time[1], ord-1), seq(time[1], time[length(time)], length.out = nknots - ord), rep(time[length(time)], ord-1))
  B <- t(splineDesign(knots, time, ord = ord))
  class(B) <- "matrix"
  IB <- B > 0

  # ------------ A1: initializing ------------
  xtx <- crossprod(X, X)
  ident_p <- diag(10, p)

  xty <- crossprod(X, Y)
  bw <- B %*% W
  bwwb <- tcrossprod(bw, bw)

  if(!length(init) || !is.matrix(init)){
    init <- (((((solve(xtx + ident_p) %*% xty) %*% W) %*% WT) %*% t(B)) %*% solve(bwwb))
  }
  gamma <- init

  # ------------ A2 ------------
  svd_x <- svd(X)
  svd_bw <- svd(bw)


  # ------------ A3 schedule ------------
  V <- t(svd_x$v)
  U <- svd_bw$u
  Z0 <- as.matrix((svd_x$d)^2) %*% ((svd_bw$d)^2)
  # This is important, rho must be large enough, but large initial rho will result in
  # slowing down optimization. Gradually increase rho
  # 60 is just randomly chosen. might need to check whether 60 is proper ?
  rhos <- seq(0, 3.5 * log10(svd_x$d[[1]]), length.out = max_iter)
  rhos <- 10^rhos * n
    # seq(0.1, 60 * (svd_x$d[[1]])^2, length.out = max_iter) * n
  inner_iters <- ceiling(inner_iter / (exp(1) - 1) * exp(seq_len(max_iter) / max_iter))

  vxtywwttbu <- ((((V %*% xty) %*% W) %*% WT) %*% t(B)) %*% U
  force(vxtywwttbu)
  # ------------ A4 ------------
  # last_mse = Inf
  for(ii in seq_len(max_iter)){
    rho <- rhos[ii] # * (lambda * ii / max_iter)
    Z <- Z0 + rho

    # ------------ A4.1 ------------
    D <- apply(abs(gamma), 1, function(g){
      tmp <- (colSums(IB * g))
      tmp[tmp > 0] <- tmp[tmp > 0]^(alpha - 1)
      #tmp[!is.finite(tmp)] = 0
      tmp %*% t(IB)
    })
    D <- alpha * t(D)

    # ------------ A4.2 ------------
    eta <- xi <- gamma
    theta <- array(0, dim(gamma))
    for(m in seq_len(inner_iters[[ii]])){
      # xi <- (t(V) %*% ((vxtywwttbu - theta +
      #                   ((rho * V) %*% eta) %*% U) / Z)) %*% t(U)
      xi <- tcrossprod(
        crossprod( V,
                   (vxtywwttbu - theta + ((rho * V) %*% eta) %*% U) / Z),
        U)
      # eta <- xi + ((1/rho * t(V)) %*% theta) %*% t(U)
      eta <- xi + tcrossprod(crossprod(V /rho, theta), U)
      eta_plus <- eta - lambda/rho * D
      eta_minus <- eta + lambda/rho * D
      eta <- eta_plus * (eta_plus > 0) + eta_minus * (eta_minus < 0)

      theta <- theta + ((rho * V) %*% (xi - eta)) %*% U
    }

    gamma <- eta

    ##
    # No early termination because rho is increasing
    # early termination might result in sparsity not recovered

    # mse <- (sqrt(mean((Y - X %*% (eta %*% B)) ^ 2)))
    # print(mse)
    # if(ii > max_iter / 2 && last_mse * 0.9 < mse && mse < last_mse ){
    #   break
    # } else {
    #   last_mse <- mse
    # }
  }
  re_gamma <- gamma
  if(CI){

    tmp <- gamma
    tmp[tmp == 0] <- runif(sum(tmp == 0), min = -1e-7, max = 1e-7)
    DD <- U %*% t(theta) %*% V
    DD <- t(DD) * sign(gamma) * lambda / tmp

    residual <- Y - X %*% solve(t(X) %*% X) %*% t(X) %*% Y %*% t(B) %*% solve(B %*% t(B)) %*% B

    sigma <- cov(residual)

    p <- dim(X)[2]
    K <- dim(gamma)[2]

    # for covariance
    s <- abs(re_gamma) > 0

    idjk <- which(abs(tmp) > 0, arr.ind = TRUE)

    P <- apply(idjk, 1, function(idx){
      j <- idx[1]
      k <- idx[2]
      p <- xtx[,j, drop = FALSE] %*% bwwb[k, , drop = FALSE]
      d <- DD[j,k]
      if( is.finite(d) && d > 0 ){
        p[j,k] <- p[j,k] + d
      }

      apply(idjk, 1, function(idy){
        p[idy[1], idy[2]]
      })
    })

    # image.plot(P)
    # t(P) %*% ? %*% P <- Q
    P_solve <- solve(P)


    wwtbt <- W %*% t(B %*% W)
    Q <- apply(idjk, 1, function(idx1){
      j1 <- idx1[1]
      k1 <- idx1[2]
      q1 <- wwtbt[, k1, drop = FALSE]
      apply(idjk, 1, function(idx2){
        j2 <- idx2[1]
        k2 <- idx2[2]
        q2 <- wwtbt[, k2, drop = FALSE]

        t(q1) %*% sigma %*% q2 * (xtx[j1, j2])
      })
    })

    avar_gamma <- t(P_solve) %*% Q %*% (P_solve)
    idxx <- idjk[,2] + (idjk[,1]-1) * K
    avar_g <- matrix(0, nrow = p*K, ncol = p*K)
    for(ii in seq_along(idxx)){
      avar_g[idxx[ii], idxx] <- avar_gamma[ii, ]
    }

    idjk <- which(!s, arr.ind = TRUE)
    idxx <- idjk[,2] + (idjk[,1]-1) * K

    f_sd <- t(sapply(seq_len(p), function(j){
      idx <- seq_len(K) + (j-1) * K
      cov <- t(B) %*% avar_g[idx, idx] %*% B
      v <- diag(cov)
      sqrt(v)
    }))

  }else{
    f_sd <- avar_g <- NULL
  }


  structure(list(
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
    f_sd = f_sd,
    CI = CI
  ), class = "spfda.fos.gbridge")
}





