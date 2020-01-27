# functions to compare to other methods
#' @export
run_grmcp <- function(X, Y, K){
  n = nrow(X)
  p = ncol(X)
  fosr_result = refund::fosr.vs(Y~X-1, method = "grMCP", nbasis = K)
  fosr_mcp = fosr_result$coefficients
  se_mcp = sqrt(diag(fosr_result$vcov) / n)
  se_mcp = t(replicate(p, se_mcp))

  list(
    coefficients = fosr_mcp,
    se = se_mcp
  )
}

#' @export
run_grlasso <- function(X, Y, K){
  n = nrow(X)
  p = ncol(X)

  fosr_result = refund::fosr.vs(Y~X-1, method = "grLasso", nbasis = K)
  fosr_lasso = fosr_result$coefficients
  se_lasso = sqrt(diag(fosr_result$vcov) / n);
  se_lasso = t(replicate(p, se_lasso))


  list(
    coefficients = fosr_lasso,
    se = se_lasso
  )
}

#' @export
run_for2s <- function(X, Y, K){
  fosr_result = refund::fosr2s(Y = Y, X = X, nbasis = K)
  fosr_coef = t(fosr_result$est.func)
  se_coef = t(fosr_result$se.func);
  list(
    coefficients = fosr_coef,
    se = se_coef
  )
}


#' @export
run_ols <- function(X, Y, K){
  ols_res = spfda(Y = Y, X = X, nsp = K, ord = 4, lambda = 0, alpha = 0.5,
                  W = NULL, CI=TRUE, max_iter = 100)
  # flm = lm(Y ~ X-1)
  # flm$coefficients - coefficients
  list(
    coefficients = ols_res$get_coef(),
    se = ols_res$raw$f_sd
  )
}






