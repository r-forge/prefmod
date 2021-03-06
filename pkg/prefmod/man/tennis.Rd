\encoding{UTF-8}
\name{tennis}
\alias{tennis}
\docType{data}
\title{Data (paired comparisons): Preferred Interview Partner}
\description{%
  The data describes results from a paired comparison study where 68 male and 96 female students were asked whom they would prefer to interview.
  The potential interview partners were Bonnie Blair, Jackie Joyner, and Jennifer Capriati.
}
\usage{tennis}
\format{%
  A data frame with 16 observations on the following 5 variables.
  \describe{
    \item{\code{n}}{counts of response pattern (\code{C1, C2, C3})}
    \item{\code{C1}}{Blair vs. Joyner: (1) Blair preferred, (\eqn{-1}) Joyner preferred}
    \item{\code{C2}}{Blair vs. Capriati: (1) Blair preferred, (\eqn{-1}) Capriati preferred)}
    \item{\code{C3}}{Joyner vs. Capriati: (1) Joyner preferred, (\eqn{-1}) Capriati preferred)}
    \item{\code{SEX}}{a numeric vector: (1) male, (2) female}
  }
}
\references{%
  Böckenholt, U., & Dillon, W. R., (1997). Modeling within-subject dependencies in ordinal paired comparison data. \emph{Psychometrika, 62}(3), 411--434.
}
\examples{
tdat <- expand.mat(tennis[, -1], tennis[, 1])
head(tdat)
}
\keyword{datasets}
