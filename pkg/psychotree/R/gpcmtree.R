#' Main function to run tree based DIF global testing for GPCM models
#'
#' @param formula A symbolic description of the model to be fit. This should be of type \code{y ~ x1 + x2} where
#'  \code{y} should be an item response matrix and \code{x1} and \code{x2} are used as partitioning variables.
#'  In the case where we wish to allow for different ability parameter distributions
#'  to be estimated independently according to a group variable, the formula is modified using the operator \code{|}.
#'  For example, \code{y ~ g | x1 + x2} will use \code{x1} and \code{x2} as partitioning variables and a different ability
#'  parameter distribution will be estimated for each level in \code{g}.
#' @param data a data frame containing the variables in the model.
#' @param start an optional vector or list of starting values (see \code{\link[psychotools]{gpcmodel}}).
#' @param weights an optional vector of weights (interpreted as case weights).
#' @param grouppars logical. Should the estimated distributional group parameters of a multiple group model be included in the model parameters?
#' @param vcov logical or character specifying the type of variance-covariance matrix (if any) computed for the final models (see and \code{\link[psychotools]{gpcmodel}}).
#' @param nullcats see \code{\link[psychotools]{gpcmodel}}, currently only "downcode" is available.
#' @param method control parameter for the optimizer employed by \code{\link[mirt]{mirt}} for the EM algorithm (see \code{\link[psychotools]{gpcmodel}}).
#' @param maxit control parameter for the optimizer employed by \code{\link[psychotools]{gpcmodel}}.
#' @param reltol control parameter for the optimizer employed by \code{\link[psychotools]{gpcmodel}}.
#' @param minsize The minimum number of observations in each node, which is passed to \code{\link[partykit]{mob_control}}. If not set, it is 500.
#' @param ... arguments passed to \code{\link[partykit]{mob_control}} for \code{gpcmtree}, and to the underlying \code{plot} method.
#'
#' @return An object of S3 class \code{"gpcmtree"} inheriting from class \code{"modelparty"}.
#' @importFrom partykit mob
#' @export mob
#' @export gpcmtree
#'
#' @details gpcmtree implements an application of model-based recursive partitioning (implemented in \code{\link[partykit]{mob}})
#' to GPCM models implemented in \code{\link[psychotools]{gpcmodel}}.
#'
#' Various methods are provided for \code{"gpcmtree"} objects, most of them inherit their behavior from \code{"modelparty"} objects (e.g., \code{print}, \code{summary}).
#'
#'
gpcmtree <- function(formula, data,
                     weights = NULL, grouppars = FALSE, vcov = TRUE, nullcats = "downcode",
                     start = NULL, method = "BFGS", maxit = 500L, reltol = 1e-10, minsize = 500, ...){

  ## keep call
  cl <- match.call(expand.dots = TRUE)

  ## use dots for setting up mob_control
  control <- partykit::mob_control(...)
  control$xtype <- "data.frame"
  control$ytype <- "matrix"
  plcontrol <- list()

  ## call mob
  m <- match.call(expand.dots = FALSE)
  ## tentative code to match itemtype with glue codes
  main_call <- environment()
  m$fit <- generate_modelfit(start = main_call$start, weights = main_call$weights, itemtype = "GPCM",
                             grouppars = main_call$grouppars, vcov = main_call$vcov, method = main_call$method,
                             maxit = main_call$maxit, reltol = main_call$reltol)
  m$control <- control
  for(n in names(plcontrol)) {
    if(!is.null(plcontrol[[n]])) m[[n]] <- plcontrol[[n]]
  }
  if("..." %in% names(m)) m[["..."]] <- NULL
  m[[1L]] <- as.name("mob")
  rval <- eval(m, parent.frame())

  ## extend class and keep original call
  rval$info$call <- cl
  class(rval) <- c("gpcmtree", class(rval))
  return(rval)
}

## methods


#' @method print gpcmtree
#' @rdname gpcmtree
#' @param x an object of class \code{gpcmtree}
#' @param title,objfun : arguments passed to \code{\link[partykit]{mob}}.
#' @export
print.gpcmtree <- function(x,
                         title = "PL Tree",
                         objfun = "negative log-likelihood", ...) {
  partykit::print.modelparty(x, title = title, objfun = objfun, ...)
}

#' @method plot gpcmtree
#' @rdname gpcmtree
#' @param type character specifying the type of plot
#' @param terminal_panel,tp_args,tnex,drop_terminal : arguments passed to \code{\link[partykit]{mob}}.
#' @export
plot.gpcmtree <- function(x, type = c("regions", "profile"), terminal_panel = NULL,
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

#' @method summary gpcmtree
#' @rdname gpcmtree
#' @param node numeric vector specifying for which node(s) to compute a summary (defaults to all leaves of the tree).
#' @export
summary.gpcmtree <- function(object, node = NULL, ...) {
  summary.modelparty(object, node = node, ...)
}

#' @method coef gpcmtree
#' @rdname gpcmtree
#' @param object an object of class \code{gpcmtree}
#' @param node numeric vector specifying for which node(s) to compute a summary (defaults to all leaves of the tree).
#' @export
coef.gpcmtree <- function(object, node = NULL, ...){
  coef.modelparty(object, node = node, ...)
}


#' @method itempar gpcmtree
#' @rdname gpcmtree
#' @param object an object of class \code{gpcmtree}
#' @param node numeric vector specifying for which node(s) to compute a summary (defaults to all leaves of the tree).
#' @export
itempar.gpcmtree <- function(object, node = NULL, ...) {
  if (is.null(node))
    node <- partykit::nodeids(object, terminal = TRUE)
  cf <- apply_to_models(object, node = node, FUN = function(n) psychotools::itempar(n, ...))
  cf
}

#' @method threshpar gpcmtree
#' @rdname gpcmtree
#' @param object an object of class \code{gpcmtree}
#' @param node numeric vector specifying for which node(s) to compute a summary (defaults to all leaves of the tree).
#' @export
threshpar.gpcmtree <- function(object, node = NULL, ...) {
  if (is.null(node))
    node <- partykit::nodeids(object, terminal = TRUE)
  cf <- apply_to_models(object, node = node, FUN = function(n) psychotools::threshpar(n, ...))
  cf
}

#' @method guesspar gpcmtree
#' @rdname gpcmtree
#' @param object an object of class \code{gpcmtree}
#' @param node numeric vector specifying for which node(s) to compute a summary (defaults to all leaves of the tree).
#' @export
guesspar.gpcmtree <- function(object, node = NULL, ...) {
  if (is.null(node))
    node <- partykit::nodeids(object, terminal = TRUE)
  cf <- apply_to_models(object, node = node, FUN = function(n) psychotools::guesspar(n, ...))
  cf
}

#' @method upperpar gpcmtree
#' @rdname gpcmtree
#' @param object an object of class \code{gpcmtree}
#' @param node numeric vector specifying for which node(s) to compute a summary (defaults to all leaves of the tree).
#' @export
upperpar.gpcmtree <- function(object, node = NULL, ...) {
  if (is.null(node))
    node <- partykit::nodeids(object, terminal = TRUE)
  cf <- apply_to_models(object, node = node, FUN = function(n) psychotools::upperpar(n, ...))
  cf
}
