\name{node_regionplot}

\alias{node_regionplot}

\title{Panel-Generating Function for Visualizing IRT Tree Models}

\usage{
node_regionplot(
  mobobj,
  names = FALSE,
  abbreviate = TRUE,
  type = c("mode", "median", "mean"),
  ref = NULL,
  ylim = NULL,
  off = 0.1,
  col_fun = gray.colors,
  bg = "white",
  uo_show = TRUE,
  uo_col = "red",
  uo_lty = 2,
  uo_lwd = 1.25,
  ylines = 2
)
}

\arguments{
\item{mobobj}{an object of class \code{"pltree"} or class \code{"mob"}
fitted by \code{\link[psychotree]{pltree}}}

\item{names}{logical or character. If \code{TRUE}, the names of the items are displayed
on the x-axis. If \code{FALSE}, numbers of items are shown.
Alternatively a character vector of the same length as the number of items can be supplied.}

\item{abbreviate}{logical. Should item names be abbreviated?
If numeric this controls the length of the abbreviation.}

\item{type}{character, specifying which type of threshold parameters
are to be used to mark the category regions per item in the plot
(see \code{\link[psychotools]{regionplot}} for details).}

\item{ref}{a vector of labels or position indices of item parameters
which should be used as restriction/for normalization.
If \code{NULL} (the default), all items are used (sum zero restriction).
 See \code{\link[psychotools]{threshpar}} for more details.}

\item{ylim}{y axis limits}

\item{off}{numeric, the distance (in scale units) between two item rectangles.}

\item{col_fun}{function. Function to use for creating the color
palettes for the rectangles. Per default \code{gray.colors} is
used. Be aware that \code{col_fun} should accept as first argument an
integer specifying the number of colors to create.}

\item{bg}{color for background filling.}

\item{uo_show}{logical. If set to \code{TRUE} (the default),
disordered absolute item threshold parameters are indicated by a
horizontal line (only if \code{type} is set to \code{"mode"}).}

\item{uo_col}{character, color of indication lines (if \code{uo_show}).}

\item{uo_lty}{numeric. Line typ of indication lines (if \code{uo_show}).}

\item{uo_lwd}{numeric. Line width of indication lines (if \code{uo_show}).}

\item{ylines}{numeric. Number of lines used for y-axis labels.}
}

\value{
A panel function which can be supplied to the \code{plot} method for \code{"pltree"} objects
 or \code{"mob"} objects fitted by \code{\link[psychotree]{pltree}}.
}

\description{
Panel-Generating Function for Visualizing IRT Tree Models
}
