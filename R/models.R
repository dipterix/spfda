# Wrapper for the model

#' @title Sparse Function-on-scalar Regression with Group Bridge Penalty
#'
#' @description \loadmathjax Function-on-scalar regression model, denote
#' \mjseqn{n} as total number of observations, \mjseqn{p} the number of
#' coefficients, \mjseqn{K} as the number of B-splines, \mjseqn{T} as total
#' time points.
#'
#' @param Y Numeric \mjseqn{n \times T} matrix, response function.
#' @param X Numeric \mjseqn{n \times p} matrix, design matrix
#' @param time Time domain, numerical length of \mjseqn{T}
#' @param nsp Integer or 'auto', number of B-splines \mjseqn{K};
#' default is 'auto'
#' @param ord B-spline order, default is \code{4}; must be \mjseqn{\geq 3}
#' @param lambda Regularization parameter \mjseqn{\gamma}
#' @param alpha Bridge parameter \mjseqn{\alpha}, default is \code{0.5}
#' @param W A \mjseqn{T \times T} weight matrix or \code{NULL}
#' (identity matrix); default is \code{NULL}
#' @param init Initial \mjseqn{\gamma}; default is \code{NULL}
#' @param max_iter Number of outer iterations
#' @param inner_iter Number of \mjseqn{ADMM} iterations (inner steps)
#' @param CI Logical, whether to calculate theoretical confidence intervals
#' @param ... Ignored
#'
#' @return A \code{spfda.model} object (environment) with following elements:
#' \describe{
#' \item{B}{B-spline basis functions used}
#' \item{error}{Root Mean Square Error ('RMSE')}
#' \item{CI}{Whether confidence intervals are calculated}
#' \item{gamma}{B-spline coefficient \mjseqn{\gamma_{p \times K}}}
#' \item{generate_splines}{Function to generate B-splines given time points}
#' \item{K}{Number of B-spline basis functions}
#' \item{knots}{B-spline knots used to fit the model}
#' \item{predict}{Function to predict responses \mjseqn{\beta(t)} given new
#'   \code{X} and/or time points}
#' \item{raw}{A list of raw variables}
#' }
#'
#'
#' @details This function implements "Functional Group Bridge for Simultaneous
#' Regression and Support Estimation" (\url{https://arxiv.org/abs/2006.10163}).
#' The model estimates functional coefficients \mjseqn{\beta(t)} under model
#' \mjsdeqn{y(t) = X\beta(t) + \epsilon(t)} with B-spline basis expansion
#' \mjsdeqn{\beta(t) = \gamma B(t) + R(t), } where \mjseqn{ R(t) } is B-spline
#' approximation error. The objective function
#' \mjsdeqn{
#' \left\| (Y-X\gamma B)W \right\|_{2}^{2} + \sum_{j,m}
#' \left\| \gamma_{j}^{T}\mathbf{1}(B^{t} > 0) \right\|_{1}^{\alpha}.
#' }
#' The input response variable is a matrix. If \mjseqn{y_{i}(t)} are observed
#' at different time points, please interpolate (e.g.
#' \code{\link[stats]{kernel}}) before feeding in.
#'
#' @examples
#'
#' dat <- spfda_simulate()
#' x <- dat$X
#' y <- dat$Y
#'
#' fit <- spfda(y, x, lambda = 5, CI = TRUE)
#'
#' BIC(fit)
#'
#' plot(fit, col = c("orange", "dodgerblue3", "darkgreen"),
#'      main = "Fitted with 95% CI", aty = c(0, 0.5, 1), atx = c(0,0.2,0.8,1))
#' matpoints(fit$time, t(dat$env$beta), type = 'l', col = 'black', lty = 2)
#' legend('topleft', c("Fitted", "Underlying"), lty = c(1,2))
#'
#' print(fit)
#' coefficients(fit)
#'
#' @export
spfda <- function(
  Y, X, lambda, time = seq(0, 1, length.out = ncol(Y)), nsp = 'auto', ord = 4,
  alpha = 0.5, W = NULL, init = NULL, max_iter = 50,
  inner_iter = 50, CI = FALSE, ...
){
  if(alpha <= 0 || alpha > 1){
    stop("alpha must be in (0,1]")
  }

  if(ord <= 2){
    stop("ord must be at least 3")
  }

  if(isTRUE(nsp == 'auto')){
    r <- ord - 2
    n_tp <- length(time)
    n_obs <- nrow(X)
    tau <- log(n_tp) / log(n_obs) * 2 * r
    if(tau > 1){
      nsp <- n_obs / 5
      nsp <- min(nsp, n_tp / 2)
    } else {
      nsp <- n_obs^(tau / (0.99 - alpha) / (2*r - 1) / 2 / r)
      nsp <- min(nsp, n_tp / 2, n_obs / 2)
    }
  }
  nsp <- ceiling(nsp)

  # if(alpha == 1){
  #   warning("The model is essentially Lasso when ", sQuote("alpha=1"), ". The implementation might be buggy under this situation.")
  #   res <- fda_glasso(Y = Y, X = X, time = time, nknots = nsp, lambda = lambda, W = W, init = init, max_iter = max_iter, ord = ord, ...)
  # }else{
  res <- fos_gp_bridge(Y = Y, X = X, time = time,
                      nknots = nsp, lambda = lambda, alpha = alpha,
                      W = W, init = init, ord = ord,
                      max_iter = max_iter, inner_iter = inner_iter, CI = CI, ...)
  # }


  # Result need to be edited
  # list(gamma = gamma, eta = eta, B = B, mse = (sqrt(mean((Y - X %*% (eta %*% B)) ^ 2))))

  env <- new.env(parent = baseenv())
  env$gamma <- res$eta
  env$knots <- c(rep(time[1], ord-1), seq(time[1], time[length(time)], length.out = nsp - ord), rep(time[length(time)], ord-1))
  env$generate_splines <- function(.time){
    if(missing(.time) || !length(.time)){
      .time <- time
    }
    return(t(splineDesign(env$knots, .time, ord = ord)))
  }
  env$B <- res$B
  env$error <- res$mse
  env$predict <- function(new_data, .time = NULL){
    B <- env$generate_splines(.time)
    new_data %*% env$gamma %*% B
  }
  env$get_coef <- function(.time = NULL){
    env$gamma %*% env$generate_splines(.time)
  }

  env$get_se <- function() {
    res$f_sd
  }

  env$raw <- res
  env$X <- X
  env$Y <- Y
  env$time <- time
  env$W <- W

  env$K <- nsp
  env$lambda <- lambda
  env$alpha <- alpha
  env$initial_gamma <- init
  env$max_iter <- max_iter
  env$CI <- res$CI

  class(env) <- c('spfda.model.gbridge', 'spfda.model')


  # calculate BIC and log-lik
  # env$BIC <- BIC(env)
  # env$loglik <- logLik(env)


  env
}


