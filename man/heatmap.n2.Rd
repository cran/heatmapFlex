% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/heatmap.n2.R
\name{heatmap.n2}
\alias{heatmap.n2}
\title{Main heatmap function}
\usage{
heatmap.n2(
  x,
  main = "Heatmap",
  ...,
  na.color = "gray80",
  sidebars = NULL,
  picketdata = NULL,
  r.cex = 1,
  c.cex = 1,
  titleheight = lcm(1.2),
  dendroheight = lcm(3),
  dendrowidth = lcm(3),
  labelheight = lcm(1),
  labelwidth = lcm(1),
  picketheight = 0.4,
  sidebarwidth = lcm(1),
  sidebar.cex = 1,
  colorkeylabels = FALSE,
  legendcorner = "bottomleft",
  plot = TRUE,
  factorpalettefn = "Pastel1",
  add.sig = FALSE,
  pv = NULL,
  order_list = TRUE,
  genes2highl = NULL
)
}
\arguments{
\item{x}{(\code{numeric}). Numeric matrix.}

\item{main}{(\code{character}). Main plot title.}

\item{...}{Additional arguments passed to \code{prepare_heatmap_data2}.}

\item{na.color}{(\code{character}). Colour code or name used for missing values. Defaults to \code{gray80}.}

\item{sidebars}{(\code{list}). List of data frames for each of the four sides of the heatmap ("bottom", "left", "top", "right").
Each column in the data frame represents some statistics or other numerical variable describing the rows or
columns of the data matrix.}

\item{picketdata}{(\code{data.frame} or \code{numeric}). Usually a data frame with columns consisting of 0 and 1.
The columns describe categories that the samples (rows in the picketdata) belong to, such as sex or a health stage.
The categories must be binary, so the data frame would for example have two rows, 'female' and 'male',
and the rows (=samples) are coded with 1 (=TRUE) and 0 (=FALSE). See \code{\link[Heatplus]{picketPlot}} for a
description of the general purpose of a picket plot.}

\item{r.cex}{(\code{numeric}). Font size for row labels.}

\item{c.cex}{(\code{numeric}). Font size for column labels.}

\item{titleheight}{(\code{numeric} or \code{lcm(numeric)}). Relative or absolute height (using \code{\link[graphics]{lcm}}, see \code{\link[graphics]{layout}}) of the plot title.}

\item{dendroheight}{(\code{numeric} or \code{lcm(numeric)}). Relative or absolute height (using \code{\link[graphics]{lcm}}, see \code{\link[graphics]{layout}}) of the dendrogram.}

\item{dendrowidth}{(\code{numeric} or \code{lcm(numeric)}). Relative or absolute width (using \code{\link[graphics]{lcm}}, see \code{\link[graphics]{layout}}) of the dendrogram.}

\item{labelheight}{(\code{numeric} or \code{lcm(numeric)}). Relative or absolute height (using \code{\link[graphics]{lcm}}, see \code{\link[graphics]{layout}}) of the labels.}

\item{labelwidth}{(\code{numeric} or \code{lcm(numeric)}). Relative or absolute width (using \code{\link[graphics]{lcm}}, see \code{\link[graphics]{layout}}) of the labels.}

\item{picketheight}{(\code{numeric} or \code{lcm(numeric)}). Relative or absolute height (using \code{\link[graphics]{lcm}}, see \code{\link[graphics]{layout}}) of the picket plots.}

\item{sidebarwidth}{(\code{numeric} or \code{lcm(numeric)}). Relative or absolute width (using \code{\link[graphics]{lcm}}, see \code{\link[graphics]{layout}}) of the sidebars.}

\item{sidebar.cex}{(\code{numeric}). Font size for sidebar labels.}

\item{colorkeylabels}{(\code{logical}). Should colour key labels be added?}

\item{legendcorner}{(\code{character}). Position of the legend.}

\item{plot}{(\code{logical}). Draw the plot? Defaults to \code{TRUE}.}

\item{factorpalettefn}{(\code{character}). Name of the colour palette.}

\item{add.sig}{(\code{logical}). Should significance asterixes be drawn?}

\item{pv}{(\code{numeric}). Matrix of the same dimensions as \code{x} with P-Values for correlation significance.}

\item{order_list}{(\code{logical}). Should the order of the correlation matrix be reversed? Meaningful if the order
of input variables should be preserved because \code{\link[graphics]{image}} turns the input matrix.}

\item{genes2highl}{(\code{character}). Vector of gene symbols (or whatever labels are used) to be highlighted.
If not \code{NULL} will draw a semi-transparent rectangle around the labels and rows or columns in the heatmap
labels.}
}
\value{
In addition to its side-effect of plotting the heatmap (if \code{plot} is \code{TRUE}, which is the default),
    the function will invisibly returns the prepared "display list", a \code{list} of all parts of the heatmap used by
    \code{hlayout()} to actually draw the final heatmap.
}
\description{
This is the main function to be called be end users. It accepts a numeric matrix and draws a heatmap.
}
\details{
The result can be used for zooming. For simple basic cases, picketvar can be given directly. Factor
    sidebars are supported, but legends are only shown for the first two.
    Arguments that will be passed on the \code{prepare_heatmap_data2}:
\tabular{ll}{
\tab labRow (\code{character}). Custom row labels\cr
\tab labCol (\code{character}). Custom column labels\cr
\tab Rowv (\code{dendrogram} or \code{integer}). Custom dendrogram object or integer vector giving the ordering index for rows\cr
\tab Colv (\code{dendrogram} or \code{integer}). Custom dendrogram object or integer vector giving the ordering index for columns\cr
\tab reorder (\code{logical}). Boolean of length 2 for rows and columns. Should rows and/or columns be reordered according to \code{Rowv}/\code{Colv}\cr
\tab distfun (\code{function}). Distance function. Defaults to \code{\link[stats]{dist}}\cr
\tab hclustfun (\code{function}). Hierarchical clustering function. Defaults to \code{\link[stats]{hclust}}\cr
\tab rowMembers (\code{character}). Group vector to split rows by\cr
\tab colMembers (\code{character}). Group vector to split columns by\cr
\tab spacer (\code{integer}). Length 1 integer used as spacer to separate groups (see \code{rowMembers}/\code{colMembers}). Defaults to 1\cr
\tab scale (\code{character}). One of "row", "column" or "none": By which dimension should data be scaled? Defaults to "none"\cr
\tab trim (\code{numeric}). Value to "cut off" data distribution. Values and both ends of the distribution, larger or smaller, respectively, will be made equal to \code{+/-trim}. Defaults to \code{NULL}, no trimming\cr
\tab zlim (\code{numeric}). Value to set the ranges for plotting \code{x} and \code{y} values. Supports both symmetrical (with default \code{NULL}) and asymmetrical (\code{NA}) automatic \code{zlim.} Used for zooming\cr
\tab col (\code{character}) Name (acronym) of the colour palette to use. Can be one of "RdBkGn" (c("green", "black", "red")), "BuYl" (c("blue", "yellow")), "BuWtRd" (c("blue", "white", "red")) or a valid name used by \code{\link[RColorBrewer]{brewer.pal}}. Defaults to NULL which will set "RdBkGn"\cr
\tab filter (\code{integer} or \code{logical}.) For filtering rows/columns: TRUE (=1.0) = remove rows/columns with only NAs, 0.5 = remove if >= 50 per cent NAs, etc. Defaults to \code{c(TRUE, TRUE)}\cr
}
}
\examples{
# Generate a random 10x10 matrix and plot it using default values
mat <- matrix(rnorm(100), nrow = 10)
heatmap.n2(mat)

# Generate a random 10x10 matrix with two distinct sets, order it using
# default clustering methods, split it into each two groups along both rows
# and columns and adjust colour palette and dendrogram dimensions:
mat <- matrix(c(rnorm(50, mean = 1), rnorm(50, mean = -1)), nrow = 10)
heatmap.n2(mat, col = "BuWtRd", rowMembers=rep(1:2, each=5),
    colMembers=rep(1:2, each=5), labRow=paste0("gene-", 1:10),
    labCol=paste0(c("A", "B"), rep(1:5, 2)), r.cex=0.8,
    dendroheight = lcm(2.2), dendrowidth = lcm(2.4))
}
\seealso{
\code{\link[graphics]{layout}}

\code{\link[RColorBrewer]{brewer.pal}}

\code{\link[grDevices]{colorRampPalette}}

\code{\link[stats]{dist}}

\code{\link[stats]{hclust}}

\code{\link[Heatplus]{picketPlot}}
}
