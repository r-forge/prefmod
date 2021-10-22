# Main function to run tree based DIF global testing for GPCM models
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
  m$fit <- generate_irtfit(start = main_call$start, weights = main_call$weights, itemtype = "GPCM",
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
print.gpcmtree <- function(x,
                         title = "PL Tree",
                         objfun = "negative log-likelihood", ...) {
  partykit::print.modelparty(x, title = title, objfun = objfun, ...)
}

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

itempar.gpcmtree <- function(object, node = NULL, ...) {
  if (is.null(node))
    node <- partykit::nodeids(object, terminal = TRUE)
  cf <- apply_to_models(object, node = node, FUN = function(n) psychotools::itempar(n, ...))
  cf
}

threshpar.gpcmtree <- function(object, node = NULL, ...) {
  if (is.null(node))
    node <- partykit::nodeids(object, terminal = TRUE)
  cf <- apply_to_models(object, node = node, FUN = function(n) psychotools::threshpar(n, ...))
  cf
}

guesspar.gpcmtree <- function(object, node = NULL, ...) {
  if (is.null(node))
    node <- partykit::nodeids(object, terminal = TRUE)
  cf <- apply_to_models(object, node = node, FUN = function(n) psychotools::guesspar(n, ...))
  cf
}

upperpar.gpcmtree <- function(object, node = NULL, ...) {
  if (is.null(node))
    node <- partykit::nodeids(object, terminal = TRUE)
  cf <- apply_to_models(object, node = node, FUN = function(n) psychotools::upperpar(n, ...))
  cf
}
