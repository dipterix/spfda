#' Calculates weight matrices
#' @param X design matrix
#' @param Y response matrix
#' @param bandwidth numeric band-width
#' @param part list of time point boundaries
#' @return the weight matrix
#' @export
spfda_weight <- function(X, Y, bandwidth, part){
  beta_ols <- solve(crossprod(X)) %*% crossprod(X, Y)
  err <- Y - X %*% beta_ols
  kern <- kernel("fejer", bandwidth, r = 2)
  err_smoothed <- kernapply(t(err), kern, circular = TRUE)
  t <- seq_len(ncol(Y))
  thetas <- lapply(part, function(idx){
    tha <- sapply(t[idx], function(t0){
      tmp <- (t-t0) / bandwidth
      x <- cbind(1, tmp[idx])
      xtx <- crossprod(x)
      thetas <- solve(xtx) %*% crossprod(x, err_smoothed[idx,])
      thetas[1,]
    })
    list(
      theta = tha,
      eps = rep(mean((err[,idx] - tha)^2), length(idx))
    )
  })
  theta <- do.call('cbind', lapply(thetas, '[[', 'theta'))
  eps <- mean((err - theta)^2)
  cov <- cov(theta)
  e <- eigen(cov)
  eig_val <- e$values + eps
  A <- diag(1/sqrt(eig_val)) %*% t(e$vectors)
  A
}
