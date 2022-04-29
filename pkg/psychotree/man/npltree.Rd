\name{npltree}
\alias{npltree}

\alias{print.npltree}
\alias{plot.npltree}
\alias{itempar.npltree}
\alias{threshpar.npltree}
\alias{guesspar.npltree}
\alias{upperpar.npltree}

\title{Parametric Logisitic (n-PL) IRT Model Trees}

\description{
  Recursive partitioning (also known as trees) based on parametric logistic
  (n-PL) item response theory (IRT) models for global testing of differential
  item functioning (DIF).
}

\usage{
npltree(formula, data, type = c("Rasch", "1PL", "2PL", "3PL", "3PLu", "4PL"),
  start = NULL, weights = NULL, grouppars = FALSE,
  vcov = TRUE, method = "BFGS", maxit = 500L,
  reltol = 1e-10, deriv = "sum", hessian = TRUE,
  full = TRUE, minsize = NULL, \dots)

\method{plot}{npltree}(x, type = c("profile", "regions"), terminal_panel = NULL,
  tp_args = list(...), tnex = 2L, drop_terminal = TRUE, \dots)
}

\arguments{
\item{formula}{A symbolic description of the model to be fit. This should be of type \code{y ~ x1 + x2} where
  \code{y} should be an item response matrix and \code{x1} and \code{x2} are used as partitioning variables.
  For the models estimated using marginal maximum likelihood (MML), it is additionally poosible to allow for
  impact of a group variable so that different ability distributions are estimated in each group.
  This can be specified by extending the previous \code{formula} by a group factor \code{g} as \code{y ~ g | x1 + x2}.}

\item{data}{a data frame containing the variables in the model.}

\item{type}{character, specifying either the type of IRT model in \code{npltree} (see also \code{\link[psychotools]{nplmodel}})
  or the type of visualization to be used in the \code{plot} method, respectively.}

\item{start}{an optional vector or list of starting values (see \code{\link[psychotools]{raschmodel}}
  or \code{\link[psychotools]{nplmodel}}).}

\item{weights}{an optional vector of weights (interpreted as case weights).}

\item{grouppars}{logical. Should the estimated distributional group parameters of a multiple-group model
  be included in the model parameters? (See \code{\link[psychotools]{nplmodel}}.)}

\item{vcov}{logical or character specifying the type of variance-covariance matrix (if any) computed for
  the final models when fitted using MML (see \code{\link[psychotools]{nplmodel}}).}

\item{method}{control parameter for the optimizer used by \code{\link[mirt]{mirt}} for the EM algorithm
  when models are fitted using MML (see \code{\link[psychotools]{nplmodel}}).}

\item{maxit}{control parameter for the optimizer used by \code{\link[psychotools]{raschmodel}} or
  \code{\link[psychotools]{nplmodel}} (see \code{\link[psychotools]{raschmodel}}, \code{\link[psychotools]{nplmodel}}).}

\item{reltol}{control parameter for the optimizer used by \code{\link[psychotools]{raschmodel}} or
  \code{\link[psychotools]{nplmodel}} (see \code{\link[psychotools]{raschmodel}}, \code{\link[psychotools]{nplmodel}}).}

\item{deriv}{character. Which type of derivatives should be used for computing gradient and Hessian matrix
  when fitting Rasch models  with the conditional maximum likelihood (CML) method (see \code{\link[psychotools]{raschmodel}})?}

\item{hessian}{logical. Should the Hessian be computed for Rasch models fitted with the CML method
  (see \code{\link[psychotools]{raschmodel}})?}

\item{full}{logical. Should a full model object be returned for Rasch models fitted with the CML method
  (see \code{\link[psychotools]{raschmodel}})?}

\item{minsize}{The minimum number of observations in each node, which is passed to \code{\link[partykit]{mob_control}}.
  If not set, it is 300 for 2PL models and 500 for 3PL, 3PLu, and 4PL models.}

\item{...}{arguments passed to \code{\link[partykit]{mob_control}} for \code{npltree}, and to the underlying
  \code{plot} method.}

\item{x}{an object of class \code{npltree}.}

\item{terminal_panel, tp_args, tnex, drop_terminal}{arguments passed to \code{\link[partykit]{mob}}.}

}

\details{
Parametric logistic (n-PL) model trees are an application of model-based recursive partitioning
(implemented in \code{\link[partykit]{mob}}) to item response theory (IRT) models (implemented in
\code{\link[psychotools]{raschmodel}} and \code{\link[psychotools]{nplmodel}}). While the \code{"Rasch"}
model is estimated by conditional maximum likelihood (CML) all other n-PL models are estimated by
marginal maximum likelihood (MML) via the standard EM algorithm. The latter allow the specification
of multiple-group model to capture group impact on the ability distributions.

Various methods are provided for \code{"npltree"} objects, most of them inherit their behavior from
\code{"modelparty"} objects (e.g., \code{print}, \code{summary}). Additionally, dedicated extractor
functions or provided for the different groups of model parameters in each node of the tree: 
\code{\link[psychotools]{itempar}} (item parameters),
\code{\link[psychotools]{threshpar}} (threshold parameters),
\code{\link[psychotools]{guesspar}} (guessing parameters),
\code{\link[psychotools]{upperpar}} (upper asymptote parameters).
}

\value{
An object of S3 class \code{"npltree"} inheriting from class \code{"modelparty"}.
}

\seealso{\code{\link[partykit]{mob}}, \code{\link[psychotools]{nplmodel}},
  \code{\link{rstree}}, \code{\link{pctree}}, \code{\link{raschtree}}, \code{\link{gpcmtree}}}

\examples{
o <- options(digits = 4)

# fit a Rasch (1PL) tree on the SPISA data set
library("psychotree")
data("SPISA", package = "psychotree")
nplt <- npltree(spisa[, 1:9] ~ age + gender + semester + elite + spon, 
  data = SPISA, type = "Rasch")
nplt

# visualize
plot(nplt)

# compute summaries of the models fitted in nodes 1 and 2
summary(nplt, 1:2)

options(digits = o$digits)
}
