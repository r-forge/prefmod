## tree-based global DIF testing for parametric logistic (PL) models
npltree <- function(formula, data, type = c("Rasch", "1PL", "2PL", "3PL", "3PLu", "4PL"),
  start = NULL, weights = NULL,
  grouppars = FALSE, vcov = TRUE, method = "BFGS", maxit = 500L, reltol = 1e-10,
  deriv = "sum", hessian = TRUE, full = TRUE, minsize = NULL, ...)
{
  ## keep call
  cl <- match.call(expand.dots = TRUE)
  
  ## process defaults
  type <- match.arg(type, c("Rasch", "1PL", "2PL", "3PL", "3PLu", "4PL"))
  if (is.null(minsize) & type == "2PL") {
    minsize <- 300
  }
  if (is.null(minsize) & type %in% c("3PL", "3PLu", "4PL")) {
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
  ma_type <- match.arg(type)
  main_call <- environment()
  m$fit <- generate_irtfit(start = main_call$start, weights = main_call$weights,
                             estfun = FALSE, object = FALSE,
                             itemtype = ma_type,
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
  class(rval) <- c("npltree", class(rval))
  return(rval)
}

## methods

print.npltree <- function(x, ...) {
  partykit::print.modelparty(x, title = "PL Tree", objfun = "negative log-likelihood", ...)
}

plot.npltree <- function(x, type = c("profile", "regions"), terminal_panel = NULL,
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

itempar.npltree <- function(object, node = NULL, ...) {
  if (is.null(node))
    node <- partykit::nodeids(object, terminal = TRUE)
  cf <- apply_to_models(object, node = node, FUN = function(n) psychotools::itempar(n, ...))
  cf
}

threshpar.npltree <- function(object, node = NULL, ...) {
  if (is.null(node))
    node <- partykit::nodeids(object, terminal = TRUE)
  cf <- apply_to_models(object, node = node, FUN = function(n) psychotools::threshpar(n, ...))
  cf
}

guesspar.npltree <- function(object, node = NULL, ...) {
  if (is.null(node))
    node <- partykit::nodeids(object, terminal = TRUE)
  cf <- apply_to_models(object, node = node, FUN = function(n) psychotools::guesspar(n, ...))
  cf
}

upperpar.npltree <- function(object, node = NULL, ...) {
  if (is.null(node))
    node <- partykit::nodeids(object, terminal = TRUE)
  cf <- apply_to_models(object, node = node, FUN = function(n) psychotools::upperpar(n, ...))
  cf
}
