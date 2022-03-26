\name{pltree}
\alias{pltree}

\alias{print.pltree}
\alias{plot.pltree}
\alias{itempar.pltree}
\alias{threshpar.pltree}
\alias{guesspar.pltree}
\alias{upperpar.pltree}

\title{Main function to run tree based DIF global testing for PL models.}

\description{
  Recursive partitioning (also known as trees) based on PL models.
}

\usage{
pltree(formula, data, type = c("Rasch", "1PL", "2PL", "3PL", "3PLu", "4PL"),
  start = NULL, weights = NULL, grouppars = FALSE,
  vcov = TRUE, method = "BFGS", maxit = 500L,
  reltol = 1e-10, deriv = "sum", hessian = TRUE,
  full = TRUE, minsize = NULL, \dots)

\method{print}{pltree}(x, \dots)

\method{plot}{pltree}(
  x, type = c("profile", "regions"), terminal_panel = NULL,
  tp_args = list(...), tnex = 2L, drop_terminal = TRUE,
  \dots
)

\method{itempar}{pltree}(object, node = NULL, \dots)

\method{threshpar}{pltree}(object, node = NULL, \dots)

\method{guesspar}{pltree}(object, node = NULL, \dots)

\method{upperpar}{pltree}(object, node = NULL, \dots)
}

\arguments{
\item{formula}{A symbolic description of the model to be fit. This should be of type \code{y ~ x1 + x2} where
\code{y} should be an item response matrix and \code{x1} and \code{x2} are used as partitioning variables.
For the models estimated using MML, in the case where we wish to allow for different ability parameter distributions
to be estimated independently according to a group variable, the formula is modified using the operator \code{|}.
For example, \code{y ~ g | x1 + x2} will use \code{x1} and \code{x2} as partitioning variables and a different ability
parameter distribution will be estimated for each level in \code{g}.}

\item{data}{a data frame containing the variables in the model.}

\item{type}{character, specifying either the type of IRT model in \code{pltree} (see also \code{\link[psychotools]{plmodel}}) or the type of visualization to be used in the \code{plot} method.}

\item{start}{an optional vector or list of starting values (see \code{\link[psychotools]{raschmodel}} or \code{\link[psychotools]{plmodel}}).}

\item{weights}{an optional vector of weights (interpreted as case weights).}

\item{grouppars}{logical. Should the estimated distributional group parameters of a multiple group model be included in the model parameters? (see \code{\link[psychotools]{plmodel}}).}

\item{vcov}{logical or character specifying the type of variance-covariance matrix (if any) computed for the final models when fitted using MML (see \code{\link[psychotools]{plmodel}}).}

\item{method}{control parameter for the optimizer used by \code{\link[mirt]{mirt}} for the EM algorithm when models are fitted using MML (see \code{\link[psychotools]{plmodel}}).}

\item{maxit}{control parameter for the optimizer used by \code{\link[psychotools]{raschmodel}} or \code{\link[psychotools]{plmodel}} (see \code{\link[psychotools]{raschmodel}}, \code{\link[psychotools]{plmodel}})}

\item{reltol}{control parameter for the optimizer used by \code{\link[psychotools]{raschmodel}} or \code{\link[psychotools]{plmodel}} (see \code{\link[psychotools]{raschmodel}}, \code{\link[psychotools]{plmodel}})}

\item{deriv}{character. Which type of derivatives should be used for computing gradient and Hessian matrix when fitting Rasch models  with the CML method (see \code{\link[psychotools]{raschmodel}})?}

\item{hessian}{logical. Should the Hessian be computed for Rasch models fitted with the CML method (see \code{\link[psychotools]{raschmodel}})?}

\item{full}{logical. Should a full model object be returned for Rasch models fitted with the CML method (see \code{\link[psychotools]{raschmodel}})?}

\item{minsize}{The minimum number of observations in each node, which is passed to \code{\link[partykit]{mob_control}}. If not set, it is 300 for 2PL models and 500 for 3PL, 3PLu and 4PL models.}

\item{...}{arguments passed to \code{\link[partykit]{mob_control}} for \code{pltree}, and to the underlying \code{plot} method.}

\item{x}{an object of class \code{pltree}}

\item{terminal_panel, tp_args, tnex, drop_terminal}{: arguments passed to \code{\link[partykit]{mob}}.}

\item{object}{an object of class \code{pltree}}

\item{node}{numeric vector specifying for which node(s) to compute a summary (defaults to all leaves of the tree).}
}

\details{
\code{pltree} implements an application of model-based recursive partitioning (implemented in \code{\link[partykit]{mob}})
to IRT models implemented in \code{\link[psychotools]{raschmodel}}, \code{\link[psychotools]{plmodel}}.

Various methods are provided for \code{"pltree"} objects, most of them inherit their behavior from \code{"modelparty"} objects (e.g., \code{print}, \code{summary}).
}

\value{
An object of S3 class \code{"pltree"} inheriting from class \code{"modelparty"}.
}

\seealso{\code{\link[partykit]{mob}}, \code{\link[psychotools]{plmodel}},
  \code{\link{rstree}}, \code{\link{pctree}}, \code{\link{raschtree}}, \code{\link{gpcmtree}}}

\examples{
o <- options(digits = 4)

# fit a Rasch (1PL) tree on the SPISA data set
library("psychotree")
data("SPISA", package = "psychotree")
plt <- pltree(spisa[, 1:9] ~ age + gender + semester + elite + spon, 
  data = SPISA, type = "Rasch")
plt

# visualize
plot(plt)

# compute summaries of the models fitted in nodes 1 and 2
summary(plt, 1:2)

options(digits = o$digits)
}
