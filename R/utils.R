# util functions

`%?<-%` <- function(lhs, value){
  env = parent.frame()
  lhs = substitute(lhs)

  tryCatch({
    is.null(eval(lhs, envir = env))
  }, error = function(e){
    return(TRUE)
  }) ->
    isnull

  if(isnull){
    # quo <- quo(!!lhs <- !!value)
    quo <- rlang::quo(do.call('=', list(quote(!!lhs), !!value)))
    eval_dirty(quo, env = env)   # Need to assign values, no eval_tidy
  }
}


eval_dirty <- function(expr, env = parent.frame(), data = NULL){

  if(rlang::is_quosure(expr)){
    expr = rlang::quo_squash(expr)
  }

  if(!is.null(data)){
    return(base::eval(expr, enclos = env, envir = data))
  }else{
    return(base::eval(expr, envir = env))
  }
}
