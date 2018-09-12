#' Single signal estimation
#' @importFrom splines splineDesign
#' @import Matrix
#' @export
spfda_single <- function(
  X, nknots = ncol(X) / 4, alpha = 0.01, lambda = 2, niter = 10
){
  # Generate time
  time = seq_len(ncol(X))

  # Generate knots and Bsplines
  knots = seq(min(time), max(time), length.out = nknots)
  nknots = length(knots)
  B = Matrix(splineDesign(knots, x = time, outer.ok = T))


  # step 0 optimization setups
  BB = t(B) %*% B
  BY = t(B) %*% t(X)
  BI = Matrix(colSums(B), ncol = 1)
  A = t(B != 0)

  nsplines = ncol(B)
  nsignals = nrow(X)
  lambda = lambda / nknots

  x_mean = Matrix(colMeans(X), nrow = 1)

  # Step 1 warm start - ridge
  w = solve(BB + diag(1, nsplines)) %*% t(x_mean %*% B)

  # step 2.1 forward step
  tau = lambda^(1/(1-alpha)) * alpha^(alpha/(1-alpha)) * (1-alpha); tau
  eta = ((1-alpha) / (alpha * tau))^alpha
  if(is.infinite(eta)){
    eta = 1000
  }


  for(i in seq_len(niter)){
    theta = eta * colSums(abs(A * as.vector(w))) ^ alpha
    ctheta = theta ^ (1 - 1/alpha)
    ctheta[!is.finite(ctheta)] = 1000
    ctheta = colSums(t(A) * ctheta)
    #ctheta[ctheta == 0] = runif(sum(ctheta == 0), min = -0.01, max = 0.01)
    # ctheta = ctheta / max(abs(ctheta))

    # step 2.2 lasso step
    # Calculate weighted LASSO problem
    X_w = t(t(B) / ctheta)
    glm = glmnet(X_w, x_mean, alpha = 1, standardize = F, intercept = FALSE, lambda = 1 / nsignals)
    w = as.vector(coef(glm))[-1] / ctheta

    # step 2.3 error
    # recall = as.vector(B %*% w)
    #
    #
    # err = sqrt(mean(apply(t(X) - as.vector(recall), 2, norm, type = '2')^2) / nrow(X))
    # print(err)
  }

  recall = as.vector(B %*% w)
  re = list(
    W = w,
    recall = recall,
    X = X,
    knots = knots,
    B = B,
    alpha = alpha,
    lambda = lambda * nknots,
    niter = niter,
    time = time,
    err = x_mean - recall,
    .plot = T
  )

  class(re) = c('spfda_single', 'spfda_base')

  return(re)

  # if(T){
  #
  #   qtiles = apply(X, 2, quantile, c(0.25,0.75))
  #   mid = apply(X, 2, median)
  #   matplot(
  #     time,
  #     t(X),
  #     col = 'grey80',
  #     type = 'l',
  #     main = main,
  #     xlab = 'Time',
  #     ylab = 'Value',
  #     ylim = range(qtiles)
  #   )
  #   matpoints(time, t(qtiles), type='l', lty = 2, col = 'yellow')
  #   points(time, mid, type = 'l', col = 'blue')
  #
  #   points(time, x_mean, type = 'l')
  #
  #   points(time, recall, type = 'l', col= 'red')
  #   legend(
  #     'topright',
  #     c('Y', 'YMean', 'Yhat', 'Q25-Q75', 'Median'),
  #     col = c('grey20', 'black', 'red', 'yellow', 'blue'),
  #     lty = c(1, 1, 1, 2, 1),
  #     cex = 0.5
  #   )
  # }
  #

}

#' @export
plot.spfda_single <- function(
  re, quantile = c(0.25, 0.5, 0.75), xlab = 'Time', ylab = '', plot_legend = T,
  pal = c("#8DA0CB", "#66C2A5", "#FFD92F", "#FC8D62", "#E78AC3", "#A6D854", "#E5C494", "#B3B3B3"), ...){
  if(length(quantile)){
    quant = apply(re$X, 2, stats::quantile, quantile)
  }else{
    quant = NA
  }
  x_mean = colMeans(re$X)
  y_lim = range(range(quant), range(re$recall), range(x_mean), na.rm = T)

  plot(range(re$time), y_lim, type = 'n', xlab = xlab, ylab = ylab, ...)

  if(plot_legend){
    qt = paste(sprintf('%.0f%%', quantile * 100), collapse = ',')
    if(length(qt)){
      qt = paste('Quantiles -', qt)
    }else{
      qt = 'Quantiles (Hidden)'
    }

    legend('topleft', c('Recall', 'Mean', qt,
                        'Knots'), col = c(pal[1:3], 'grey80'), lty = c(1,1,2,3), lwd = c(2,1,1,0.3),
           xpd = T,
           bty = 'n')
  }

  if(length(quant) > 1){
    matpoints(re$time, t(quant), type = 'l', lty = 2, col = pal[3])
  }

  points(re$time, x_mean, type = 'l', lty = 1, col = pal[2])

  points(re$time, re$recall, type = 'l', col = pal[1], lwd = 2)

  abline(v = re$knots, lty = 3, lwd = 0.3, col = 'grey80')


}

#' @export
cat.spfda_base <- function(re, ...){
  print(re)
}

#' @export
print.spfda_base <- function(re, plot = T, ...){
  cat(
    sprintf(
      paste(
        'X: %d signals x %d time points',
        'Splines: %d knots, %d basis',
        'Penalty: \n\tlambda: %s\n\talpha: %s',
        'Results:\n\tNumber of iterations: %d\n\tPoint-wise mean distance from mean: %s',
        sep = '\n'
      ),
      nrow(re$X), ncol(re$X), length(re$knots), ncol(re$B),
      paste(sprintf('%.4f', re$lambda), collapse = ', '),
      paste(sprintf('%.4f', re$alpha), collapse = ', '),
      re$niter,
      paste(sprintf('%.4f', apply(re$err, 1, function(x){mean(abs(x))})), collapse = ', ')
    )
  )

  if(plot && re$.plot){
    plot(re)
    re$.plot = F
  }

  return(invisible(re))
}

