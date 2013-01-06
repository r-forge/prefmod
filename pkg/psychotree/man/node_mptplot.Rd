\name{node_mptblot}
\alias{node_mptblot}

\title{Panel-Generating Function for Visualizing MPT Tree Models}

\description{
  Panel-generating function for visualizing the model parameters
  from the nodes in MPT tree models.
}

\usage{
node_mptblot(mobobj, id = TRUE,
  names = TRUE, abbreviate = TRUE, index = TRUE, ref = TRUE,
  col = "black", linecol = "lightgray", cex = 0.5, pch = 19,
  xscale = NULL, yscale = NULL, ylines = 1.5)

}

\arguments{
  \item{mobobj}{an object of class \code{"mob"} based on Bradley-Terry models
    fitted by \code{\link{btReg}}.}
  \item{id}{logical. Should the node ID be displayed?}
  \item{names}{logical. Should the names for the objects be displayed?}
  \item{abbreviate}{logical or numeric. Should object names be abbreviated?
    If numeric this controls the length of the abbreviation.}
  \item{index}{logical. Should different indexes for different stimuli be used?}
  \item{ref}{logical. Should a horizontal line for the reference level be drawn?}
  \item{col, cex, pch}{graphical appearance of plotting symbols.}
  \item{linecol}{line color for reference line (if \code{ref}).}
  \item{xscale, yscale}{x and y axis limits.}
  \item{ylines}{numeric. Number of lines used for y-axis labels.}
}

\details{
  The panel-generating function \code{node_mptblot} is called by the \code{plot}
  method for \code{"mpttree"} objects and does not have to be called by the
  user directly.
}

\value{
  A panel function which can be supplied to the \code{plot} method for
  \code{"mob"} objects.
}

\seealso{\code{\link{mptree}}}

\keyword{hplot}