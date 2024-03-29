% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/scHOT_plot.R
\name{plotOrderedExpression}
\alias{plotOrderedExpression}
\title{plotOrderedExpression}
\usage{
plotOrderedExpression(
  scHOT,
  genes,
  positionType = NULL,
  branches = NULL,
  ranked_by = NULL,
  xvals = NULL,
  subsetBranch = NULL,
  facet = FALSE,
  positionColData = NULL,
  assayName = NULL,
  point_size = 3
)
}
\arguments{
\item{scHOT}{A scHOT object, where the expression
data is stored in the assay slot, with assay name "expression".}

\item{genes}{is a character vector for gene names}

\item{positionType}{A string indicates the position type,
either trajectory or spatial}

\item{branches}{A character indicates that the colnames stored
the branch information in colData}

\item{ranked_by}{A character indicates that the colnames stored
the ranking information of the cells in colData, such as trajectory time}

\item{xvals}{A character indicates that the colnames stored in
colData of the x-values associated with the samples in branchData}

\item{subsetBranch}{subsetBranch is a character vector containing
the names of the branches to be plotted. If NULL it will plot all branches}

\item{facet}{can either be FALSE, "branch", "gene", or "both"}

\item{positionColData}{A vector indicates column names of colData that
stored the postion informaton (for spatial type of data)}

\item{assayName}{the name of the assay that are used to plot.}
}
\value{
\code{ggplot} a ggplot object for ribbon plot with points
}
\description{
the plotOrderedExpression function plots expression vectors
along branches and genes as ribbon plots
}
