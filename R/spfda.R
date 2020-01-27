#' @export
spfda <- function(
  # Data
  Y, X, time = seq(0, 1, length.out = ncol(Y)),
  # B-splines
  nsp = floor(length(time)/2), ord = 4,
  # params
  lambda = 0.1, alpha = 0.1, W = NULL, D = NULL,
  # other params
  init = NULL, max_iter = 50,
  # TODO: add cv and CI into the function
  ...
){
  stopifnot(alpha > 0 && alpha <= 1)

  if(alpha == 1){
    res = fda_glasso(Y = Y, X = X, time = time, nknots = nsp, lambda = lambda, W = W, D = D, init = init, max_iter = max_iter, ord = ord, ...)
  }else{
    res = fda_bridge(Y = Y, X = X, time = time,
               nknots = nsp, lambda = lambda, alpha = alpha,
               W = W, D = D, init = init, ord = ord,
               max_iter = max_iter, ...)
  }


  # Result need to be edited
  # list(gamma = gamma, eta = eta, B = B, mse = (sqrt(mean((Y - X %*% (eta %*% B)) ^ 2))))

  env = new.env(parent = baseenv())
  env$gamma = res$eta
  env$knots = c(rep(time[1], ord-1), seq(time[1], time[length(time)], length.out = nsp - ord), rep(time[length(time)], ord-1))
  env$generate_splines = function(.time){
    .time %?<-% time
    return(t(splineDesign(env$knots, .time, ord = ord)))
  }
  env$B = res$B
  env$error = res$mse
  env$predict = function(new_data, .time = NULL){
    B = env$generate_splines(.time)
    new_data %*% env$gamma %*% B
  }
  env$get_coef = function(.time = NULL){
    env$gamma %*% env$generate_splines(.time)
  }

  env$raw = res
  env$X = X
  env$Y = Y
  env$time = time
  env$W = W

  env$K = nsp
  env$lambda = lambda
  env$alpha = alpha
  env$initial_gamma = init
  env$max_iter = max_iter

  class(env) = c('spfda.model', 'environment')
  env
}
