\name{node_profileplot}
\alias{node_profileplot}

\title{Panel-Generating Function for Visualizing IRT Tree Models}

\description{
  Panel-generating function for visualizing profiles (points and lines)
  of the parameters from the nodes in IRT tree models.
}

\usage{
node_profileplot(
  mobobj,
  what = c("items", "thresholds", "discriminations"),
  parg = list(type = NULL, ref = NULL, alias = TRUE),
  id = TRUE,
  names = FALSE,
  abbreviate = TRUE,
  index = TRUE,
  ref = TRUE,
  col = "black",
  border = col,
  linecol = "black",
  refcol = "lightgray",
  bg = "white",
  cex = 0.5,
  pch = 21,
  xscale = NULL,
  yscale = NULL,
  ylines = 2,
  \dots
)
}

\arguments{
\item{mobobj}{an object of class \code{"pltree"} or class \code{"mob"}
fitted by \code{\link[psychotree]{pltree}}}

\item{what}{specifying the type of parameters to be plotted}

\item{parg}{supplementary arguments for \code{"what"}}

\item{id}{logical. Should the node ID be displayed?}

\item{names}{logical or character. If \code{TRUE}, the names of
the items are displayed on the x-axis. If \code{FALSE}, numbers of
items are shown. Alternatively a character vector of the same
length as the number of items can be supplied.}

\item{abbreviate}{logical. Should item names be abbreviated?
If numeric this controls the length of the abbreviation.}

\item{index}{logical. Should different indexes for different items be used?}

\item{ref}{logical. Should a horizontal line for the reference level be drawn?}

\item{col, border, pch, cex}{graphical appearance of plotting symbols.}

\item{linecol, refcol}{character, specifying the line color to use for the
profile lines and reference line, respectively.}

\item{bg}{color for background filling.}

\item{xscale, yscale}{x and y axis limits.}

\item{ylines}{numeric. Number of lines used for y-axis labels.}

\item{...}{further arguments currently not used.}
}

\details{
The panel-generating function \code{node_regionplot} is called by the \code{plot}
method of \code{"gpcmtree"} object by default and does not have to be called by the user directly.
See \code{\link[psychotools]{regionplot}} for details and references
of the drawn region plots and possible values and their meaning for
the argument \code{type} (taken by \code{node_regionplot}).
}

\value{
A panel function which can be supplied to the \code{plot} method for \code{"pltree"} objects
or \code{"mob"} objects fitted by \code{\link[psychotree]{pltree}} or \code{\link[psychotree]{gpcmtree}}.
}

\keyword{hplot}
