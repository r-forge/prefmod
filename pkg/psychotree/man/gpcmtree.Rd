\name{gpcmtree}
\alias{gpcmtree}

\alias{print.gpcmtree}
\alias{plot.gpcmtree}
\alias{itempar.gpcmtree}
\alias{threshpar.gpcmtree}
\alias{guesspar.gpcmtree}
\alias{upperpar.gpcmtree}

\title{Generalized Partial Credit Model Trees}

\description{
  Recursive partitioning (also known as trees) based on generalized partial
  credit models (GPCMs) for global testing of differential item functioning (DIF).
}

\usage{
gpcmtree(formula, data, weights = NULL,
  grouppars = FALSE, vcov = TRUE, nullcats = "downcode",
  start = NULL, method = "BFGS", maxit = 500L,
  reltol = 1e-10, minsize = 500, \dots)

\method{plot}{gpcmtree}(x, type = c("regions", "profile"), terminal_panel = NULL,
  tp_args = list(...), tnex = 2L, drop_terminal = TRUE, \dots)
}

\arguments{
\item{formula}{A symbolic description of the model to be fit. This should be of type \code{y ~ x1 + x2} where
  \code{y} should be an item response matrix and \code{x1} and \code{x2} are used as partitioning variables.
  Additionally, it is poosible to allow for impact of a group variable so that different ability distributions
  are estimated in each group. This can be specified by extending the previous \code{formula} by a group factor
  \code{g} as \code{y ~ g | x1 + x2}.}

\item{data}{a data frame containing the variables in the model.}

\item{weights}{an optional vector of weights (interpreted as case weights).}

\item{grouppars}{logical. Should the estimated distributional group parameters of a multiple group model be
  included in the model parameters?}

\item{vcov}{logical or character specifying the type of variance-covariance matrix (if any) computed for the
  final models (see \code{\link[psychotools]{gpcmodel}}).}

\item{nullcats}{character string, specifying how items with null categories (i.e., categories not observed)
  should be treated. See \code{\link[psychotools]{gpcmodel}}, currently only \code{"downcode"} is available.}

\item{start}{an optional vector or list of starting values (see \code{\link[psychotools]{gpcmodel}}).}

\item{method}{control parameter for the optimizer employed by \code{\link[mirt]{mirt}} for the EM algorithm
  (see \code{\link[psychotools]{gpcmodel}}).}

\item{maxit}{control parameter for the optimizer employed by \code{\link[psychotools]{gpcmodel}}.}

\item{reltol}{control parameter for the optimizer employed by \code{\link[psychotools]{gpcmodel}}.}

\item{minsize}{integer specification of minimum number of observations in each node, which is passed to
  \code{\link[partykit]{mob_control}}.}

\item{...}{arguments passed to \code{\link[partykit]{mob_control}} for \code{gpcmtree}, or to the
  underlying \code{plot} method, respectively.}

\item{x}{an object of class \code{gpcmtree}.}

\item{type}{character specifying the type of plot.}

\item{terminal_panel, tp_args, tnex, drop_terminal}{arguments passed to \code{\link[partykit]{mob}}.}
}


\details{
Generalized partial credit model (GPCM) trees are an application of model-based recursive partitioning
(implemented in \code{\link[partykit]{mob}}) to GPCM models (implemented in \code{\link[psychotools]{gpcmodel}}).
For technical and algorithmic details, see the documentation of the two core functions
linked above as well as \code{vignette("mob", package = "partykit")}.

Various methods are provided for \code{"gpcmtree"} objects, most of them inherit their behavior from
\code{"modelparty"} objects (e.g., \code{print}, \code{summary}). Additionally, dedicated extractor
functions or provided for the different groups of model parameters in each node of the tree: 
\code{\link[psychotools]{itempar}} (item parameters),
\code{\link[psychotools]{threshpar}} (threshold parameters),
\code{\link[psychotools]{guesspar}} (guessing parameters),
\code{\link[psychotools]{upperpar}} (upper asymptote parameters).
}

\value{
An object of S3 class \code{"gpcmtree"} inheriting from class \code{"modelparty"}.
}

\seealso{\code{\link[partykit]{mob}}, \code{\link[psychotools]{plmodel}},
  \code{\link{rstree}}, \code{\link{pctree}}, \code{\link{raschtree}}, \code{\link{npltree}}}
