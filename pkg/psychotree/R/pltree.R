#' Main function to run tree based DIF global testing for PL models.
#'
#' @param formula A symbolic description of the model to be fit. This should be of type \code{y ~ x1 + x2} where
#'  \code{y} should be an item response matrix and \code{x1} and \code{x2} are used as partitioning variables.
#'  For the models estimated using MML, in the case where we wish to allow for different ability parameter distributions
#'  to be estimated independently according to a group variable, the formula is modified using the operator \code{|}.
#'  For example, \code{y ~ g | x1 + x2} will use \code{x1} and \code{x2} as partitioning variables and a different ability
#'  parameter distribution will be estimated for each level in \code{g}.
#' @param data a data frame containing the variables in the model.
#' @param itemtype character string, specifying the type of IRT model to be estimated (see details below)
#' @param start an optional vector or list of starting values (see \code{\link[psychotools]{raschmodel}} or \code{\link[psychotools]{plmodel}}).
#' @param weights an optional vector of weights (interpreted as case weights).
#' @param estfun a (potentially user-defined) estfun function. (remove?)
#' @param object a (potentially user-defined) model object. (remove?)
#' @param grouppars logical. Should the estimated distributional group parameters of a multiple group model be included in the model parameters? see plmodel. I think I can remove this parameter?
#' @param vcov logical or character specifying the type of variance-covariance matrix (if any) computed for the final models when fitted using MML (see \code{\link[psychotools]{plmodel}}).
#' @param method control parameter for the optimizer employed by \code{\link[mirt]{mirt}} for the EM algorithm when models are fitted using MML (see \code{\link[psychotools]{plmodel}}).
#' @param maxit control parameter for the optimizer employed by (see \code{\link[psychotools]{raschmodel}}, \code{\link[psychotools]{plmodel}})
#' @param reltol control parameter for the optimizer employed by (see \code{\link[psychotools]{raschmodel}}, \code{\link[psychotools]{plmodel}})
#' @param deriv character. Which type of derivatives should be used for computing gradient and Hessian matrix when fitting using \code{\link[psychotools]{raschmodel}}?
#' @param hessian logical. Should the Hessian be computed for models fitted by \code{\link[psychotools]{raschmodel}}?
#' @param full logical. Should a full model object be returned for models fitted by \code{\link[psychotools]{raschmodel}}?
#' @param minsize The minimum number of observations in each node, which is passed to \code{\link[partykit]{mob_control}}. If not set, it is 300 for 2PL models and 500 for 3PL, 3PLu and 4PL models.
#' @param ... arguments passed to \code{\link[partykit]{mob_control}} for \code{pltree}, and to the underlying \code{plot} method.
#'
#' @return An object of S3 class \code{"pltree"} inheriting from class \code{"modelparty"}.
#' @importFrom partykit mob
#' @export mob
#' @export pltree
#'
#' @details Pltree implements an application of model-based recursive partitioning (implemented in \code{\link[partykit]{mob}})
#' to IRT models implemented in \code{\link[psychotools]{raschmodel}}, \code{\link[psychotools]{plmodel}}.
#'
#' Various methods are provided for \code{"pltree"} objects, most of them inherit their behavior from \code{"modelparty"} objects (e.g., \code{print}, \code{summary}).
#'
#' @examples
#' #fit a rasch tree on spisa dataset
#' library(psychotree)
#' data("SPISA")
#' fit_rasch <- pltree(spisa ~ age + gender + semester + elite + spon, 
#'                     data = SPISA, itemtype = "Rasch")
#'
#' #print the rasch tree
#' fit_rasch
#'
#' #plot the rasch tree
#' plot(fit_rasch)
#'
#' #compute summaries of the models fitted in nodes 1 and 2
#' summary(fit_rasch, 1:2)
pltree <- function(formula, data, itemtype=c("Rasch", "1PL", "2PL", "3PL", "3PLu", "4PL"),
                   start = NULL, weights = NULL, estfun = FALSE, object = FALSE,
                   grouppars = FALSE, vcov = TRUE, method = "BFGS", maxit = 500L, reltol = 1e-10,
                   deriv = "sum", hessian = TRUE, full = TRUE, minsize = NULL, ...){

  ## keep call
  cl <- match.call(expand.dots = TRUE)
  
  if (is.null(minsize) & itemtype == "2PL") {
    minsize <- 300
  }
  if (is.null(minsize) & itemtype %in% c("3PL", "3PLu", "4PL")) {
    minsize <- 500
  }
  ## use dots for setting up mob_control
  control <- partykit::mob_control(minsize = minsize,...)
  control$xtype <- "data.frame"
  control$ytype <- "matrix"
  plcontrol <- list()

  ## call mob
  m <- match.call(expand.dots = FALSE)
  ## tentative code to match itemtype with glue codes
  type <- match.arg(itemtype)
  main_call <- environment()
  m$fit <- generate_modelfit(start = main_call$start, weights = main_call$weights,
                             estfun = main_call$estfun, object = main_call$object,
                             itemtype = type,
                             grouppars = main_call$grouppars, vcov = main_call$vcov, method = main_call$method,
                             maxit = main_call$maxit, reltol = main_call$reltol,
                             deriv = main_call$deriv, hessian = main_call$hessian, full = main_call$full)
  m$control <- control
  for(n in names(plcontrol)) {
    if(!is.null(plcontrol[[n]])) m[[n]] <- plcontrol[[n]]
  }
  if("..." %in% names(m)) m[["..."]] <- NULL
  m[[1L]] <- as.name("mob")
  rval <- eval(m, parent.frame())

  ## extend class and keep original call
  rval$info$call <- cl
  class(rval) <- c("pltree", class(rval))
  return(rval)
}

## methods

#' @method print pltree
#' @rdname pltree
#' @param x an object of class \code{pltree}
#' @export
print.pltree <- function(x, ...) {
  partykit::print.modelparty(x, title = "PL Tree", objfun = "negative log-likelihood", ...)
}

#' @method plot pltree
#' @rdname pltree
#' @param x an object of class \code{pltree}
#' @param type character specifying the type of plot
#' @param terminal_panel,tp_args,tnex,drop_terminal : arguments passed to \code{\link[partykit]{mob}}.
#' @export
plot.pltree <- function(x, type = c("profile", "regions"), terminal_panel = NULL,
                        tp_args = list(...), tnex = 2L, drop_terminal = TRUE, ...)
{
  if(!is.null(terminal_panel)) {
    if(!missing(type)) warning("only one of 'type' and 'terminal_panel' should be specified")
  } else {
    terminal_panel <- switch(match.arg(type),
                             "regions" = node_regionplot,
                             "profile" = node_profileplot)
  }
  partykit::plot.modelparty(x, terminal_panel = terminal_panel,
                            tp_args = tp_args, tnex = tnex, drop_terminal = drop_terminal, ...)
}

#' @method summary pltree
#' @rdname pltree
#' @param object an object of class \code{pltree}
#' @param node numeric vector specifying for which node(s) to compute a summary (defaults to all leaves of the tree).
#' @export
summary.pltree <- function(object, node = NULL, ...) {
  summary.modelparty(object, node = node, ...)
}

#' @method coef pltree
#' @rdname pltree
#' @param object an object of class \code{pltree}
#' @param node numeric vector specifying for which node(s) to compute a summary (defaults to all leaves of the tree).
#' @export
coef.pltree <- function(object, node = NULL, ...){
  coef.modelparty(object, node = node, ...)
}

#' @method itempar pltree
#' @rdname pltree
#' @param object an object of class \code{pltree}
#' @param node numeric vector specifying for which node(s) to compute a summary (defaults to all leaves of the tree).
#' @export
itempar.pltree <- function(object, node = NULL, ...) {
  if (is.null(node))
    node <- partykit::nodeids(object, terminal = TRUE)
  cf <- apply_to_models(object, node = node, FUN = function(n) psychotools::itempar(n, ...))
  cf
}

#' @method threshpar pltree
#' @rdname pltree
#' @param object an object of class \code{pltree}
#' @param node numeric vector specifying for which node(s) to compute a summary (defaults to all leaves of the tree).
#' @export
threshpar.pltree <- function(object, node = NULL, ...) {
  if (is.null(node))
    node <- partykit::nodeids(object, terminal = TRUE)
  cf <- apply_to_models(object, node = node, FUN = function(n) psychotools::threshpar(n, ...))
  cf
}

#' @method guesspar pltree
#' @rdname pltree
#' @param object an object of class \code{pltree}
#' @param node numeric vector specifying for which node(s) to compute a summary (defaults to all leaves of the tree).
#' @export
guesspar.pltree <- function(object, node = NULL, ...) {
  if (is.null(node))
    node <- partykit::nodeids(object, terminal = TRUE)
  cf <- apply_to_models(object, node = node, FUN = function(n) psychotools::guesspar(n, ...))
  cf
}

#' @method upperpar pltree
#' @rdname pltree
#' @param object an object of class \code{pltree}
#' @param node numeric vector specifying for which node(s) to compute a summary (defaults to all leaves of the tree).
#' @export
upperpar.pltree <- function(object, node = NULL, ...) {
  if (is.null(node))
    node <- partykit::nodeids(object, terminal = TRUE)
  cf <- apply_to_models(object, node = node, FUN = function(n) psychotools::upperpar(n, ...))
  cf
}


## Imported from psychotree / partykit
apply_to_models <- function(object, node = NULL, FUN = NULL, drop = FALSE, ...) {
  if(is.null(node)) node <- nodeids(object, terminal = FALSE)
  if(is.null(FUN)) FUN <- function(object, ...) object  
  rval <- if("object" %in% object$info$control$terminal) {
    nodeapply(object, node, function(n) FUN(info_node(n)$object))
  } else {
    lapply(refit.modelparty(object, node, drop = FALSE), FUN)
  }
  names(rval) <- node
  if(drop & length(node) == 1L) rval <- rval[[1L]]
  return(rval)
}

## Imported from partykit
coef.modelparty <- function(object, node = NULL, drop = TRUE, ...) {
  if(is.null(node)) node <- nodeids(object, terminal = TRUE)
  cf <- do.call("rbind", nodeapply(object, ids = node, FUN = function(n) info_node(n)$coefficients))
  if(drop) drop(cf) else cf
}

summary.modelparty <- function(object, node = NULL, ...)
{
  ids <- if(is.null(node)) nodeids(object, terminal = TRUE) else node
  rval <- apply_to_models(object, node = ids, FUN = summary)
  if(length(ids) == 1L) rval[[1L]] else rval
}