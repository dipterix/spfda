spfda_admm <- function(Y, X, N, nknots, lambda=0.1, alpha=1, max_iter = 10, D = NULL, gamma = NULL){
  TT = nrow(Y) # T, total timepoints
  time = seq_len(TT)
  N = as.vector(N)
  NN = sum(N)
  precision = Diagonal(TT)

  # Generate knots and Bsplines
  knots = seq(min(time), max(time), length.out = nknots)
  nknots = length(knots)
  B = Matrix(splineDesign(knots, x = time, outer.ok = T))

  # Initialize indices
  L = ncol(B)
  I = length(N)
  J = N
  K = nrow(X)

  # Group Y by condition
  j_ind = lapply(seq_len(I), function(i){
    J = N[i]
    (sum(N[1:i]) - J) + seq_len(J)
  })

  y_mean = sapply(j_ind, function(j){
    rowMeans(Y[,j,drop=F])
  })


  # Avoid rep calc
  sum_xi_xiT = (J * t(X)) %*% X
  sum_B_SigmaInv_yij_xi = t(B) %*% precision %*% y_mean %*% solve(sum_xi_xiT) %*% (J * t(X))
  BT_sigma_B = t(B) %*% precision %*% B
  IB = B != 0

  # Initialize gamma
  gamma = solve(Diagonal(L, NN * 0.01) + BT_sigma_B) %*% sum_B_SigmaInv_yij_xi
  eta = gamma
  theta = 0

  # matplot(y_mean, type='l', ylim = c(-5,5), lty=1)
  # matpoints(B %*% eta %*% X, type='l', lty=2)



  # ADMM for alpha = 1
  if(alpha == 1){
    if(is.null(D)){
      D = replicate(K, colSums(IB))
    }
    # rho
    rhos = seq(0.01, 20, length.out = max_iter)

    for(rho in rhos){
      # Update gamma
      gamma = solve(Diagonal(L, rho) + BT_sigma_B) %*% (sum_B_SigmaInv_yij_xi + rho * eta - theta)

      # Update eta
      eta_plus = gamma + theta / rho - lambda * D / rho
      eta_minus = gamma + theta / rho + lambda * D / rho

      eta = eta_plus * (eta_plus > 0) + eta_minus * (eta_minus < 0)

      # Update dual theta
      theta = theta + gamma - eta

      # Calculate Loss
      # print(norm(y_mean - B %*% gamma %*% X, type = 'F'))
      par(mfrow=c(1,3))
      matplot(y_mean, type='l', lty=1)
      # matpoints(B %*% eta %*% X, type='l', lty=1)
      matplot(B %*% eta %*% X, type='l', lty=1)
      matplot(B %*% eta, type='l', lty=1)
    }

    return(list(
      gamma = gamma,
      eta = eta,
      theta = theta
    ))

  }else{
    # Calculate zeta_{tk}
    err = Inf

    for(ii in seq_len(max_iter)){
      zeta = apply(gamma, 2, function(r_k){
        apply(t(IB) * r_k, 2, function(x){
          sum(abs(x))
        })
      })
      zeta = alpha * zeta^(alpha - 1)
      # Handle NAN and Inf since spline(outer=T)
      D = t(B) %*% zeta

      # Solve L1 problem
      re = spfda_admm(Y=Y, X=X, N=N, nknots=nknots, lambda=lambda, alpha=1, D = D, max_iter=10)
      gamma = re$gamma

      new_err = sum(N * apply(y_mean - B %*% re$eta %*% X, 2, norm, type='2')) / NN
      if(err < new_err){
        break
      }else{
        err = new_err
      }
    }

    return(re)
  }
}


