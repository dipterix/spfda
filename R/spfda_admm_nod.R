#' @import Matrix
#' @import splines
#' @export
spfda_admm = function(Y, X, W = NULL, D = NULL,
                      nknots = 100, lambda = 0.1, alpha = 0.1,max_iter = 50, use_admm = T){
  TT = ncol(Y) # T, total timepoints
  time = seq_len(TT)
  N = nrow(Y)
  P = ncol(X)

  if(is.null(W)){
    W = Diagonal(TT, 1)
  }

  # Generate knots and Bsplines
  # time = seq(1,1000, length.out = 10000); nknots = 1000
  knots = seq(min(time), max(time), length.out = nknots)
  nknots = length(knots)
  B = Matrix(t((splineDesign(knots, x = time, outer.ok = T))))
  # range(eigen(B%*%t(B))$value)
  # rowSums(B)[5]
  # sum(B[50,]^2) * nknots

  # SVD
  svd_x = svd(X)
  svd_b = svd(B%*%W)
  V = svd_x$v
  D1S = (svd_x$d)^2
  U = svd_b$u
  D2S = (svd_b$d)^2
  dxb = matrix(D1S, ncol=1) %*% matrix(D2S, nrow=1)


  # Calculate V_0
  IB = B != 0
  if(is.null(D) && alpha == 1){
    D = matrix(1, nrow = P, ncol = 1) %*% matrix(rowSums(IB), nrow = 1)
  }

  L = nrow(B)


  xtx = t(X) %*% X
  xty = t(X) %*% Y
  wwt = W %*% t(W)
  wwtbt = wwt%*%t(B)
  bwbwt = B%*%wwt%*%t(B)
  bwbwt_inv = solve(bwbwt)
  i_xtx = Diagonal(nrow(xtx))



  # Initialize ridge
  gamma = solve(xtx + i_xtx) %*% xty %*% wwtbt %*% bwbwt_inv
  eta = gamma
  theta = V %*% (gamma - eta) %*% U # I'm lazy

  # debug, recall
  # matplot(t(gamma %*% B), type = 'l')

  rhos = seq(0.1, 60, length.out = max_iter) * N

  if(alpha == 1){
    # Case 1: alpha = 1
    for(rho in rhos){
      # update gamma
      tmp = t(V) %*% (xty %*% wwtbt + use_admm * rho * eta - t(V) %*% theta %*% t(U)) %*% U
      gamma = V %*% (tmp / (rho + dxb)) %*% t(U)

      # update eta

      A = (gamma + 1/rho * (t(V) %*% theta %*% t(U)))
      eta_plus = A - lambda / rho * D
      eta_minus = A + lambda / rho * D
      eta = (eta_plus * (eta_plus > 0) + eta_minus * (eta_minus < 0)) * (eta != 0)

      theta = theta + V %*% (gamma - eta) %*% U


      # matplot(t(eta %*% B), type = 'l', lty=1)
      # matpoints(cbind(colMeans(f1),colMeans(f2),colMeans(f3)), type='l', lty=2)
    }

    re = list(
      gamma = gamma,
      theta = theta,
      eta = eta,
      B = B,
      rho = rho,
      Y = Y,
      X = X
    )
    return(re)
  }else{
    for(ii in seq_len(max_iter)){
      # zeta = alpha * apply(IB, 2, function(b){
      #   apply(gamma, 1, function(g){
      #     sum(abs(b*g))
      #   })
      # })^(alpha-1)
      # D = zeta %*% t(IB)
      D = ((abs(gamma) %*% IB)^(alpha-1)) %*% t(IB)

      re = spfda_admm(Y = Y, X = X, D = D, alpha = 1, lambda = lambda, nknots = nknots, max_iter = 5)
      gamma = re$gamma
    }

    return(re)
  }


}
