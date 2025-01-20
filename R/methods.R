##' @importFrom stats update
##' @S3method update evolaMod
update.evolaMod <- function(object, formula., ..., evaluate = TRUE) {
  if (is.null(call <- getCall(object))){stop("object should contain a 'call' component")}
  # call <- getCall(object)
  
  extras <- match.call(expand.dots = FALSE)$...
  extras[["initPop"]] <- object$pop
  extras[["simParam"]] <- object$simParam # it should be in the global environment
  
  if (!missing(formula.)){
    call$formula <- update.formula(formula(object), formula.)
  }
  if (length(extras) > 0) {
    existing <- !is.na(match(names(extras), names(call)))
    for (a in names(extras)[existing]) call[[a]] <- extras[[a]]
    if (any(!existing)) {
      call <- c(as.list(call), extras[!existing])
      call <- as.call(call)
    }
  }
  if (!evaluate) return(call)
 
  ## should be able to find model components somewhere in (1) formula env; (2) calling env;
  ##  (3) parent frame [plus its parent frames]
  ## see discusion at https://stackoverflow.com/questions/64268994/evaluate-call-when-components-may-be-scattered-among-environments
  ## FIXME: suppressWarnings(update(model)) will give
  ## Error in as.list.environment(X[[i]], ...) :
  ## promise already under evaluation: recursive default argument reference or earlier problems?
  
  ff <- environment(formula(object$call))
  pf <- parent.frame()
  sf <- sys.frames()[[1]]
  
  tryCatch(eval(call,  envir = ff),  ## try formula environment
           error = function(e) {
             tryCatch(eval(call, envir = sf),  ## try stack frame
                      error = function(e) {
                        eval(call, envir=pf) ## try parent frame
                      })
           })
  
  ##
  ## combf <- tryCatch(
  ##     do.call("c", lapply(list(ff, sf), as.list)),
  ##     error=function(e) as.list(ff)
  ## )
  ## eval(call,combf, enclos=pf)
}