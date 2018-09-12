#' @importFrom glmnet glmnet
#' @import splines
#' @export
sparse_bs = function(X, Y, time, nknots = 150, niter = 4, alpha = 0.05, lambda = 0.1, plot = T, quiet = F, main = 'Data v.s. Sparse BS Recall'){
  # Multiple Y should have multiple penalties

  # Generate knots and Bsplines
  knots = seq(min(time), max(time), length.out = nknots)
  B = splineDesign(knots, x = time, outer.ok = T);


  # step 0 optimization setups
  IB = B; #cbind(1, B)
  BB = t(B) %*% B
  BY = t(B) %*% t(X)
  BI = as.matrix(colSums(B))
  A = t(B != 0)

  # Step 1 ridge
  w = solve(t(IB) %*% IB + diag(1, ncol(IB))) %*% (t(IB) %*% rowMeans(t(X)))
  # plot(
  #   time,
  #   apply(X, 2, mean),
  #   type = 'l',
  #   main = 'Data v.s. Ridge Recall',
  #   xlab = 'Time',
  #   ylab = 'Value'
  # )
  # points(time, as.vector(IB%*%w), type = 'l', col= 'red')
  # legend('topright', c('Y', 'Yhat'), col = c(1,2), lty=1)

  # step 2.1 forward step

  tau = lambda^(1/(1-alpha)) * alpha^(alpha/(1-alpha)) * (1-alpha); tau

  for(i in seq_len(niter)){
    theta = ((1-alpha) / (alpha * tau))^alpha * colSums(abs(A * as.vector(w))) ^ alpha
    ctheta = theta ^ (1 - 1/alpha)
    ctheta[!is.finite(ctheta)] = 100
    ctheta = colSums(t(A) * ctheta)
    ctheta[ctheta == 0] = runif(sum(ctheta == 0), min = -0.01, max = 0.01)
    # ctheta = ctheta / max(abs(ctheta))

    # step 2.2 lasso step
    # Calculate weighted LASSO problem
    X_w = t(t(B) / ctheta)
    glm = glmnet(X_w, as.matrix(colMeans(X)), alpha = 1, standardize = F, intercept = FALSE, lambda = 1 / nrow(X))
    w = as.vector(coef(glm))[-1] / ctheta

    # step 2.3 error
      recall = as.vector(IB %*% w)


      err = sqrt(mean(apply(t(X) - as.vector(recall), 2, norm, type = '2')^2) / nrow(X))
    if(!quiet){
      cat('Iter: ', i, ' - RMSE: ', err, '\n', sep = '')
    }
  }
  if(plot){

    x_mean = colMeans(X)
    qtiles = apply(X, 2, quantile, c(0.25,0.75))
    mid = apply(X, 2, median)
    matplot(
      time,
      t(X),
      col = 'grey80',
      type = 'l',
      main = main,
      xlab = 'Time',
      ylab = 'Value',
      ylim = range(qtiles)
    )
    matpoints(time, t(qtiles), type='l', lty = 2, col = 'yellow')
    points(time, mid, type = 'l', col = 'blue')

    points(time, x_mean, type = 'l')

    points(time, recall, type = 'l', col= 'red')
    legend(
      'topright',
      c('Y', 'YMean', 'Yhat', 'Q25-Q75', 'Median'),
      col = c('grey20', 'black', 'red', 'yellow', 'blue'),
      lty = c(1, 1, 1, 2, 1),
      cex = 0.5
    )
  }
  return(list(
    B = IB,
    w = w,
    RMSE = err,
    recall = IB %*% w
  ))
}
