## generator function yielding glue code, calling raschmodel, plmodel, gpcmodel, or rsmodel
## (analogous to raschfit, rsmfit, pcmfit, and mptfit)
generate_irtfit <- function(start = NULL, weights = NULL, estfun = FALSE, object = FALSE,
  itemtype = c("Rasch", "2PL", "3PL", "3PLu", "4PL", "1PL", "RM", "GPCM", "PCM", "RSM"),
  grouppars = FALSE, vcov = TRUE, method = "BFGS", maxit = 500L, reltol = 1e-10,
  deriv = "sum", hessian = TRUE, full = TRUE)
{
  # match itemtype argument
  type <- match.arg(itemtype)
  # copy arguments to pass to modelfit_function
  curr_call <- environment()
  # generate function
  modelfit_function <- function(y, x, start = curr_call$start, weights = curr_call$weights,
    offset = NULL, ..., estfun = curr_call$estfun, object = curr_call$object,
    itemtype = type,
    grouppars = curr_call$grouppars, vcov = curr_call$vcov,
    method = curr_call$method, maxit = curr_call$maxit, reltol = curr_call$reltol,
    deriv = curr_call$deriv, hessian = curr_call$hessian, full = curr_call$full)
  {
    # copy environment to pass arguments to plmodel
    curr_env <- environment()
    # calls for if we have a Rasch model to be estimated with CML
    if(curr_env$itemtype %in% c("Rasch")){
      # if we try to use impact, issue warning and proceed to CML
      if(!(is.null(x) || NCOL(x) == 0L || length(unique(interaction(x))) == 1)) {warning("x (impact) not used")}
      rval <- psychotools::raschmodel(y, weights = curr_env$weights, start = curr_env$start, reltol = curr_env$reltol,
                         deriv = curr_env$deriv, hessian = curr_env$hessian, maxit = curr_env$maxit,
                         full = curr_env$full, gradtol = curr_env$reltol, iterlim = curr_env$maxit)
    # calls for if we have a pl model
    } else if (curr_env$itemtype %in% c("2PL", "3PL", "3PLu", "4PL", "1PL", "RM")) {
      # if we dont have any impact set impact as NULL
      if((is.null(x) || NCOL(x) == 0L || length(unique(interaction(x))) == 1)) {
        rval <- psychotools::plmodel(y, weights = curr_env$weights, impact = NULL,
                        type = curr_env$itemtype, grouppars = curr_env$grouppars, vcov = curr_env$vcov,
                        start = curr_env$start, method = curr_env$method, maxit = curr_env$maxit, reltol = curr_env$reltol)
      # else, set impact at interactions(x)
      } else {
        rval <- psychotools::plmodel(y, weights = curr_env$weights, impact = interaction(x),
                        type = curr_env$itemtype, grouppars = curr_env$grouppars, vcov = curr_env$vcov,
                        start = curr_env$start, method = curr_env$method, maxit = curr_env$maxit, reltol = curr_env$reltol)
      }
    # calls for if we have a gpc model
    }else if(curr_env$itemtype %in% c("GPCM", "PCM")){
      # if we dont have any impact set impact as NULL
      if((is.null(x) || NCOL(x) == 0L || length(unique(interaction(x))) == 1)) {
        rval <- psychotools::gpcmodel(y, weights = curr_env$weights, impact = NULL,
                         type = curr_env$itemtype, grouppars = curr_env$grouppars, vcov = curr_env$vcov,
                         start = curr_env$start, method = curr_env$method, maxit = curr_env$maxit, reltol = curr_env$reltol,
                         nullcats = "downcode")
      # else, set impact at interactions(x)
      } else {
        rval <- psychotools::gpcmodel(y, weights = curr_env$weights, impact = interaction(x),
                         type = curr_env$itemtype, grouppars = curr_env$grouppars, vcov = curr_env$vcov,
                         start = curr_env$start, method = curr_env$method, maxit = curr_env$maxit, reltol = curr_env$reltol,
                         nullcats = "downcode")
      }
    # calls for if we have a rsm model
    }else if(curr_env$itemtype %in% c("RSM")){
      # if we try to use impact, issue warning and proceed to CML
      if(!(is.null(x) || NCOL(x) == 0L || length(unique(interaction(x))) == 1)) {warning("x (impact) not used")}
      rval <- psychotools::rsmodel(y, weights = curr_env$weights, start = curr_env$start, reltol = curr_env$reltol,
                      deriv = curr_env$deriv, hessian = curr_env$hessian,
                      maxit = curr_env$maxit, full = curr_env$full)

    }
	## Issue Warning when model did not converge in mirt
	if (!is.null(rval$mirt) && mirt::extract.mirt(rval$mirt, what="converged") == FALSE){
		warning("Parameter estimation did not converge, results should not be interpreted.")
	}
    ## Issue Warning when model did not converge in psychotools
  if (!is.null(rval$code) && rval$code != 0){
    warning("Parameter estimation did not converge, results should not be interpreted.")
  }  
	
    if(!is.null(offset)) {
      warning("offset not used")
    }
    rval <- list(
      coefficients = rval$coefficients,
      objfun = -rval$loglik,
      estfun = if (estfun) {
        if (curr_env$itemtype %in% c("Rasch")) {
          psychotools::estfun.raschmodel(rval)
        } else if (curr_env$itemtype %in% c("2PL", "3PL", "3PLu", "4PL", "1PL", "RM")) {
          psychotools::estfun.plmodel(rval)
        } else if(curr_env$itemtype %in% c("GPCM", "PCM")){
          psychotools::estfun.gpcmodel(rval)
        } else if(curr_env$itemtype %in% c("RSM")){
          psychotools::estfun.rsmodel(rval)
        }
      } else {NULL},
      object = if (object) rval else NULL)
    return(rval)

  }
  return(modelfit_function)
}
