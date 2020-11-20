#' @import splines
#' @import graphics
#' @import stats
#' @importFrom grDevices col2rgb
#' @importFrom grDevices rgb
NULL

str2lang_alt <- function (s) {
  s <- sprintf("quote(%s)", gsub(pattern = "(^[ ]+)|([ ]+$)", x = s, replacement = ""))
  eval(parse(text = s))
}

str2lang <- function (s) {
  get0("str2lang", envir = baseenv(), ifnotfound = str2lang_alt)(s)
}

`%?<-%` <- function (lhs, value) {
  env <- parent.frame()
  lhs <- substitute(lhs)
  isnull <- tryCatch({
    is.null(eval(lhs, envir = env))
  }, error = function(e) {
    return(TRUE)
  })
  if (isnull) {
    eval(as.call(list(str2lang("`<-`"), lhs, value)), envir = env)
  }
}

getAlphaRGB <- function (colname, alpha) {
  c = col2rgb(colname)
  rgb(t(c), alpha = alpha, maxColorValue = 255)
}

get_dots <- function (..name, ..default = NULL, ...)
{
  call <- as.list(match.call(expand.dots = TRUE))[-1]
  call <- call[!names(call) %in% c("..name", "..default")]
  if (..name %in% names(call)) {
    idx <- which(names(call) == ..name)[1]
    return(...elt(idx))
  }
  else {
    return(..default)
  }
}

eval_dirty <- function(expr, env = parent.frame(), data = NULL){
  # if(rlang::is_quosure(expr)){
  #   expr = rlang::quo_squash(expr)
  # }
  if(!is.null(data)){
    return(base::eval(expr, enclos = env, envir = data))
  }else{
    return(base::eval(expr, envir = env))
  }
}

