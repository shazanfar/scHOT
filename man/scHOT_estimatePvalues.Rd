% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/scHOT_functions.R
\name{scHOT_estimatePvalues}
\alias{scHOT_estimatePvalues}
\title{scHOT_estimatePvalues}
\usage{
scHOT_estimatePvalues(
  scHOT,
  usenperm_estimate = FALSE,
  nperm_estimate = 10000,
  maxDist = 0.1,
  plot = FALSE,
  verbose = FALSE
)
}
\arguments{
\item{scHOT}{A scHOT object}

\item{usenperm_estimate}{Logical (default FALSE) if
number of neighbouring permutations should be used, or
if difference of global higher order statistic should be used}

\item{nperm_estimate}{Number of neighbouring permutations to
use for p-value estimation}

\item{maxDist}{max difference of global higher order statistic
to use for p-value estimation (default 0.1)}

\item{plot}{A logical input indicating whether the results are plotted}

\item{verbose}{A logical input indicating whether the intermediate
steps will be printed}
}
\value{
scHOT A scHOT object with results stored in scHOT_output slot
}
\description{
Estimate p-values based on already run permutation tests
}
\examples{
 data(MOB_subset)
 sce_MOB_subset <- MOB_subset$sce_MOB_subset
 scHOT_spatial <- scHOT_buildFromSCE(sce_MOB_subset,
                                     assayName = "logcounts",
                                    positionType = "spatial",
                                     positionColData = c("x", "y"))
pairs <- matrix(c("Arrb1", "Mtor", "Dnm1l", "Gucy1b3"), ncol = 2, byrow = TRUE)
rownames(pairs) <- apply(pairs,1,paste0,collapse = "_")

scHOT_spatial <- scHOT_addTestingScaffold(scHOT_spatial, pairs)

scHOT_spatial <- scHOT_setWeightMatrix(scHOT_spatial,
                                       positionColData = c("x","y"),
                                        positionType = "spatial",
                                        nrow.out = NULL,
                                        span = 0.05)
scHOT_spatial <- scHOT_calculateGlobalHigherOrderFunction(
  scHOT_spatial,
  higherOrderFunction = weightedSpearman,
  higherOrderFunctionType = "weighted")
scHOT_spatial <- scHOT_setPermutationScaffold(scHOT_spatial,
                                              numberPermutations = 100)
scHOT_spatial <- scHOT_calculateHigherOrderTestStatistics(
  scHOT_spatial,
  higherOrderSummaryFunction = sd)

scHOT_spatial <- scHOT_performPermutationTest(
  scHOT_spatial,
  verbose = TRUE,
  parallel = FALSE)

scHOT_spatial <- scHOT_estimatePvalues(scHOT_spatial)

}
