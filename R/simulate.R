#' @title Generate toy example data
#' @param n Total number of observations
#' @param n_timepoints Total number of time points
#' @param err error level
#' @param seed random seed
#' @return A list of data generated
#' @export
spfda_simulate <- function(n = 1000, n_timepoints = 100, err = 1, seed = 1){

  n_coef = 3
  beta1 = function(t){0}
  beta2 = function(t){sin(pi * t)}
  beta3 = function(t){
    if(t < 0.2 || t >= 0.8){
      return(0)
    }else if (t >= 0.4 && t < 0.6){
      return(1)
    }else if (t < 0.4){
      return(sin((5 * t - 1) * pi / 2))
    }else{
      return(sin((5 * t - 2) * pi / 2))
    }
  }

  # Generate discrete coefficients
  time = seq(0, 1, length.out = n_timepoints)
  beta = t(sapply(1:n_coef, function(ii){
    f = get(paste0('beta', ii))
    sapply(time, f)
  }))

  set.seed(seed)
  on.exit(set.seed(-1), add = TRUE, after = TRUE)
  # X: n x 3 rnorm data
  X = rnorm(n * n_coef)
  dim(X) = c(n, n_coef)

  # have some correlation
  junk = sample(n * n_coef, n*2)
  X[junk] = X[junk] + 1
  base_curve = ((time > 0.8)*3 + (time > 0.4)*2 + (time >= 0)) / 20 * err

  err = t(replicate(n, {
    arima.sim(model=list(ar=c(.9)),n=n_timepoints) * base_curve
  }))

  Y = X %*% beta + err
  Y = Y + rnorm(length(Y))

  return(list(
    X = X,
    Y = Y,
    env = environment()
  ))

}
