\name{gpcmtree}

\alias{gpcmtree}
\alias{print.gpcmtree}
\alias{plot.gpcmtree}
\alias{itempar.gpcmtree}
\alias{threshpar.gpcmtree}
\alias{guesspar.gpcmtree}
\alias{upperpar.gpcmtree}

\title{Main function to run tree based DIF global testing for generalized partial credit models}

\description{
   Recursive partitioning (also known as trees) based on generalized partial credit models.
}

\usage{
gpcmtree(formula, data, weights = NULL,
  grouppars = FALSE, vcov = TRUE, nullcats = "downcode",
  start = NULL, method = "BFGS", maxit = 500L,
  reltol = 1e-10, minsize = 500, \dots)

\method{print}{gpcmtree}(x, title = "PL Tree", objfun = "negative log-likelihood", ...)

\method{plot}{gpcmtree}(x, type = c("regions", "profile"), terminal_panel = NULL,
  tp_args = list(...), tnex = 2L, drop_terminal = TRUE,
  \dots
)

\method{itempar}{gpcmtree}(object, node = NULL, \dots)

\method{threshpar}{gpcmtree}(object, node = NULL, \dots)

\method{guesspar}{gpcmtree}(object, node = NULL, \dots)

\method{upperpar}{gpcmtree}(object, node = NULL, \dots)
}

\arguments{
\item{formula}{A symbolic description of the model to be fit. This should be of type \code{y ~ x1 + x2} where
\code{y} should be an item response matrix and \code{x1} and \code{x2} are used as partitioning variables.
In the case where we wish to allow for different ability parameter distributions
to be estimated independently according to a group variable, the formula is modified using the operator \code{|}.
For example, \code{y ~ g | x1 + x2} will use \code{x1} and \code{x2} as partitioning variables and a different ability
parameter distribution will be estimated for each level in \code{g}.}

\item{data}{a data frame containing the variables in the model.}

\item{weights}{an optional vector of weights (interpreted as case weights).}

\item{grouppars}{logical. Should the estimated distributional group parameters of a multiple group model be included in the model parameters?}

\item{vcov}{logical or character specifying the type of variance-covariance matrix (if any) computed for the final models (see \code{\link[psychotools]{gpcmodel}}).}

\item{nullcats}{character string, specifying how items with null categories (i.e., categories not observed) should be treated. See \code{\link[psychotools]{gpcmodel}}, currently only "downcode" is available.}

\item{start}{an optional vector or list of starting values (see \code{\link[psychotools]{gpcmodel}}).}

\item{method}{control parameter for the optimizer employed by \code{\link[mirt]{mirt}} for the EM algorithm (see \code{\link[psychotools]{gpcmodel}}).}

\item{maxit}{control parameter for the optimizer employed by \code{\link[psychotools]{gpcmodel}}.}

\item{reltol}{control parameter for the optimizer employed by \code{\link[psychotools]{gpcmodel}}.}

\item{minsize}{The minimum number of observations in each node, which is passed to \code{\link[partykit]{mob_control}}. If not set, it is 500.}

\item{...}{arguments passed to \code{\link[partykit]{mob_control}} for \code{gpcmtree}, and to the underlying \code{plot} method.}

\item{x}{an object of class \code{gpcmtree}}

\item{title, objfun}{: arguments passed to \code{\link[partykit]{mob}}.}

\item{type}{character specifying the type of plot}

\item{terminal_panel, tp_args, tnex, drop_terminal}{: arguments passed to \code{\link[partykit]{mob}}.}

\item{object}{an object of class \code{gpcmtree}}

\item{node}{numeric vector specifying for which node(s) to compute a summary (defaults to all leaves of the tree).}
}


\details{
\code{gpcmtree} implements an application of model-based recursive partitioning (implemented in \code{\link[partykit]{mob}})
to GPCM models implemented in \code{\link[psychotools]{gpcmodel}}.

Various methods are provided for \code{"gpcmtree"} objects, most of them inherit their behavior from \code{"modelparty"} objects (e.g., \code{print}, \code{summary}).
}

\value{
An object of S3 class \code{"gpcmtree"} inheriting from class \code{"modelparty"}.
}

\seealso{\code{\link[partykit]{mob}}, \code{\link[psychotools]{plmodel}},
  \code{\link{rstree}}, \code{\link{pctree}}, \code{\link{raschtree}}, \code{\link{pltree}}}
