#' @title Generates toy example data
#' @description Synthesized functional signals with heterogeneous error. The
#' underlying three coefficients correspond to 'dense', 'global sparse', and
#' 'local sparse' functions. See \url{https://doi.org/10.1111/biom.13684} for
#' detailed configurations.
#' @param n Total number of observations
#' @param n_timepoints Total number of time points
#' @param err Error magnitude
#' @param scale the scale of coefficients length of 1 or 3.
#' @return A list of data generated: \code{X} is scalar predictor, \code{Y} is
#' functional response.
#' @export
spfda_simulate <- function(n = 1000, n_timepoints = 100, err = 1, scale = c(1,1,1)){

  n_coef <- 3
  beta1 <- function(t){0 * scale[1]}
  beta2 <- function(t){sin(pi * t) * scale[2]}
  beta3 <- function(t){
    tmp <- function(){
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
    tmp() * scale[3]
  }

  # Generate discrete coefficients
  time <- seq(0, 1, length.out = n_timepoints)
  beta <- t(sapply(1:n_coef, function(ii){
    f <- get(paste0('beta', ii))
    sapply(time, f)
  }))

  # X: n x 3 rnorm data
  X <- rnorm(n * n_coef)
  dim(X) <- c(n, n_coef)

  # have some correlation
  junk <- sample(n * n_coef, n*2)
  X[junk] <- X[junk] + 1
  base_curve <- ((time > 0.8)*3 + (time > 0.4)*2 + (time >= 0)) / 20 * err

  err <- t(replicate(n, {
    arima.sim(model=list(ar=c(.9)),n=n_timepoints) * base_curve
  }))

  Y <- X %*% beta + err
  Y <- Y + rnorm(length(Y))

  return(list(
    X = X,
    Y = Y,
    env = environment()
  ))

}
