#' Ported function from \code{'refund'} package
#' @description A modified version of \code{\link[refund]{fosr.vs}}, but
#' with groups parameter to allow grouping time points rather than the
#' whole coefficient when the underlying functions are locally supported.
#' @param formula,data,nbasis,method,epsilon,max.iter_num see
#' \code{\link[refund]{fosr.vs}}
#' @param groups integer vector with length of number of time-points
#' of how time-points should be grouped; default
#' is \code{NULL}, indicating there is no local sparsity.
#' @export
fosr_vs <- function (formula, data, nbasis = 10,
          method = c("ls", "grLasso", "grMCP", "grSCAD"),
          epsilon = 1e-05, max.iter_num = 100, groups = NULL) {

  if(system.file(package = "refund") == ""){
    stop("Please install `refund` package first by running\n  install.packages('refund')")
  }

  cl <- match.call()
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data"), names(mf), 0)
  mf <- mf[c(1, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1]] <- as.name("model.frame")
  mf <- eval.parent(mf)
  Y <- stats::model.response(mf)
  W.des <- stats::model.matrix(formula, data = mf)
  N <- dim(Y)[1]
  D <- dim(Y)[2]
  K <- dim(W.des)[2]
  Theta <- splines::bs(1:D, df = nbasis, intercept = TRUE, degree = 3)
  Z <- kronecker(W.des, Theta)
  Y.vec <- as.vector(t(Y))
  ls <- stats::lm(Y.vec ~ Z + 0)
  a <- matrix(ls$residuals, nrow = N, ncol = D, byrow = TRUE)
  if (D < 15) {
    b <- cov(a)
  } else {
    w <- refund::fpca.sc(a, var = TRUE)
    b <- Reduce("+", lapply(seq_len(length(w$evalues)), function(x) {
      w$evalues[x] * w$efunctions[, x] %*% t(w$efunctions[, x])
    })) + w$sigma2 * diag(1, D, D)
  }

  if(length(groups)){
    if(length(groups) != D) {
      stop("length of groups must equal to the time length dim(Y)[2]")
    }
    groups <- as.integer(factor(groups, levels = unique(groups)))
    theta_groups <- apply(Theta, 2, function(s){
      names(which.max(table(groups[s!=0])))
    })
    theta_groups <- as.integer(factor(theta_groups))
  } else {
    theta_groups <- rep(1, nbasis)
  }
  mgn <- max(theta_groups)

  grpreg_ngroups <- lapply((seq_len(K) - 1) * mgn, function(i){
    i + theta_groups
  })
  if(isTRUE(colnames(W.des)[1] == "(Intercept)")){
    grpreg_ngroups[[1]] <- grpreg_ngroups[[1]] * 0
  }
  grpreg_ngroups <- unlist(grpreg_ngroups)


  group <- switch(
    colnames(W.des)[1],
    `(Intercept)` = {
      c(rep(0, nbasis), c(rep(1:(K - 1), each = nbasis)))
    }, {
      c(rep(1:K, each = nbasis))
    }
  )

  f <- stats::coef(ls)
  f <- replace(f, which(is.na(f)), 0)
  d <- rep(0, nbasis * K)
  num.iter <- 0
  cat("Beginning iterative algorithm \n")
  while (sum((f - d)^2) > epsilon & num.iter <= max.iter_num) {
    num.iter <- num.iter + 1
    d <- f
    c <- chol(solve(b))
    Y_new <- Y %*% t(c)
    Theta_new <- c %*% Theta
    Z_new <- kronecker(W.des, Theta_new)
    Y.vecnew <- as.vector(t(Y_new))
    if (method == "ls") {
      ls <- stats::lm(Y.vecnew ~ Z_new)
      f <- as.vector(stats::coef(ls))[-1]
    } else {
      cvlam <- grpreg::cv.grpreg(
        Z_new,
        Y.vecnew,
        group = grpreg_ngroups,
        penalty = method,
        max.iter = 5000
      )

      f <- as.vector(stats::coef(cvlam))[-1]
    }
    a <- matrix(Y.vec - Z %*% f, nrow = N, ncol = D, byrow = TRUE)
    if (D < 15) {
      b <- cov(a)
    } else {
      w <- refund::fpca.sc(a, var = TRUE)
      b <- Reduce("+", lapply(seq_len(length(w$evalues)), function(x) {
        w$evalues[x] * w$efunctions[, x] %*% t(w$efunctions[, x])
      })) + w$sigma2 * diag(1, D, D)
    }
    if (num.iter%%10 == 1)
      cat(".")
  }
  B <- matrix(f, nrow = K, byrow = TRUE)
  est <- B %*% t(Theta)
  rownames(est) <- colnames(W.des)
  fitted <- W.des %*% est
  res <- Y - fitted
  ret <- list(call = cl, formula = formula, coefficients = est,
              fitted.values = fitted, residuals = res, vcov = b, method = method)
  class(ret) <- "fosr.vs"
  return(ret)
}
