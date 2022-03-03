# calculate log-likelihood
logLik_old.spfda.model <- function (object, ...) {
  res <- object$Y - object$predict(object$X)
  W <- object$W
  p <- ncol(object$X)
  N <- nrow(res)
  nT <- ncol(res)
  K <- object$K
  coef <- object$get_coef()

  if(!is.matrix(W)){
    det_w <- 1
  }else{
    det_w <- det(W)
    res <- res %*% W
  }

  # assume each row of res is uncor
  val <- sum(apply(res, 1, function(x){
    - sum(x ^ 2) / 2 + nT * (log(det_w) - log(2 * pi))
  }))

  # calculate df

  df <- sum(coef != 0)
  df_alt <- sum(abs(coef) > 0.001 * max(abs(coef)))

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
  loglik <- logLik.spfda.model(object, ...)
  df <- attr(loglik, "df")
  nobs <- attr(loglik, "nobs")
  pK <- attr(loglik, "pK")
  as.numeric((-2 * loglik + df * log(nobs) + nu * df * log(pK)) / nobs)
}



# calculate log-likelihood
#' @export
logLik.spfda.model <- function (object, eps = "auto", ...) {
  res <- object$Y - object$predict(object$X)
  W <- object$W
  p <- ncol(object$X)
  N <- nrow(res)
  nT <- ncol(res)
  K <- object$K
  coef <- object$get_coef()

  if(!is.matrix(W)){
    W <- diag(1, nrow = nT, ncol = nT)
  }else{
    res <- res %*% W
  }

  # estimate sigma^2
  wtw <- tcrossprod(W, W)
  gamma_ols <- solve(crossprod(object$X, object$X)) %*% crossprod(object$X, object$Y) %*%
    wtw %*% t(object$B) %*% solve(object$B %*% wtw %*% t(object$B))
  res_ols <- (object$Y - object$X %*% gamma_ols %*% object$B) %*% W

  sigma <- sqrt(mean(res_ols^2)); #sigma <- 1
  res <- res / (sigma)

  # assume each row of res is uncor
  val <- sum(apply(res, 1, function(x){
    - sum(x ^ 2) / 2
  }))

  # calculate df

  df <- sum(coef != 0)
  if(identical(eps, "auto")){
    eps <- 0.001 * max(abs(coef))
  }
  df_alt <- sum(abs(coef) > eps)

  attr(val, "nall") <- N
  attr(val, "nobs") <- N
  attr(val, "df") <- df
  attr(val, "df_alt") <- df_alt
  attr(val, "pK") <- p * K
  class(val) <- "logLik"
  val
}


#' @export
print.spfda.model <- function(x, digits = getOption("digits"), ...){
  if(!exists("loglik", envir = x)){
    x$loglik <- logLik(x)
  }
  if(!exists("BIC", envir = x)){
    x$BIC <- BIC(x)
  }

  cat("Model: function-on-scalar with group-bridge penalty\n",
      "Log-lik: ",
      format(x$loglik, digits = digits),
      " (df=", format(attr(x$loglik, "df")), ")\n",
      "E-BIC:   ", format(x$BIC, digits = digits), "\n",
      "RMSE :   ", x$error, "\n",
      "Parameters:\n",
      " K:      ", x$K, "\n",
      " alpha:  ", x$alpha, "\n",
      " lambda: ", x$lambda, "\n",
      sep = '')
  invisible(x)
}

#' @export
coef.spfda.model <- function(object, time = NULL, ...){
  if(!length(time) || !is.numeric(time)){
    time <- object$time
  }
  B <- object$generate_splines(time)
  coef <- object$gamma %*% B
  coef
}

#' @export
plot.spfda.model <- function(
  x, y = NULL, CI = TRUE, time = y, type = 'l', lty = 1, ylim = NULL,
  col = seq_len(ncol(x$X)), atx = NULL, labelx = atx, aty = NULL, labely = aty,
  xlab = "Time", ylab = "Coefficients",
  main = "Functional Group-Bridge Coefficients", ...){

  if(!length(time) || !is.numeric(time)){
    time <- x$time
  }
  coef <- coef(x, time = time)

  if(CI && x$CI){
    coef_orig <- x$get_coef()

    # Need to plot CI
    if(isTRUE(CI)){
      fct <- 2
      cistr <- "95"
    } else {
      fct <- abs(as.numeric(CI))
      cistr <- sprintf("%.0f", 100 - 200 * pnorm(fct, lower.tail = FALSE))
    }
    sds <- x$raw$f_sd * fct

    coef_orig_ub <- sds + coef_orig
    coef_orig_lb <- coef_orig - sds

    if(length(ylim) != 2){
      ylim <- c(min(coef_orig_lb), max(coef_orig_ub))
    }

    CI <- TRUE

  } else {
    if(length(ylim) != 2){
      ylim <- range(coef)
    }
    CI <- FALSE
  }

  axes <- get_dots('axes', ..default = FALSE, ...)
  args <- list(x = range(time), y = ylim, type = 'n', xlab = xlab, ylab = ylab, ylim = ylim, main = main, ...)
  args$axes <- axes
  do.call('plot', args)

  if( CI ){
    p <- ncol(x$X)
    for(ii in seq_len(p)){
      polygon(
        c(x$time, rev(x$time)),
        c(coef_orig_lb[ii, ], rev(coef_orig_ub[ii, ])),
        border = NA,
        col = getAlphaRGB(col[[ii]], 50)
      )
    }
  }

  matpoints(x = time, y = t(coef), type = type, lty = lty, col = col, ...)

  # add axis
  if(!length(atx)){
    atx <- pretty(time)
  }
  if(length(labelx) != length(atx)){
    labelx <- atx
  }
  axis(1, at = atx, labels = labelx, tcl = -0.2, las = 1, ...)

  if(!length(aty)){
    aty <- pretty(ylim)
  }
  if(length(labely) != length(aty)){
    labely <- aty
  }
  axis(2, at = aty, labels = labely, tcl = -0.2, las = 1, ...)


}
