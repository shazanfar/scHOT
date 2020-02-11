#' A wrap up function to perform scHOT
#'
#' @title scHOT
#'
#' @param scHOT A scHOT object
#' @param testingScaffold A matrix with rows for each testing combination
#' @param weightMatrix A matrix indicates the weight matrix for scHOT analysis
#' @param positionType A string indicates the position type, either trajectory or spatial
#' @param positionColData If positionType is "trajectory" then positionColData should be a sortable vector
#' if positionType is "spatial" then positionColData should be a matrix type object.
#' It should be stored in colData (SHILA to check).
#' @param nrow.out SHILA to input
#' @param averageAcrossTrajectoryTies SHILA to input
#' @param higherOrderFunction A function object indicates the higher order function
#' @param higherOrderFunctionType is weighted or unweighted, determines if there
#' is a weighting argument in the higher order function
#' @param numberPermutations The number of permutatuion, set as 1000 by default
#' @param numberScaffold The number of scaffold, set as 100 by default
#' @param storePermutations  a logical flag on whether
#' @param higherOrderSummaryFunction A functon indicates the higher order summary function (SHILA to check)
#' @param parallel A logical input indicates whether run the permutation test using multiple cores in parallel.
#' @param BPPARAM  A \code{BiocParallelParam} class object from the \code{BiocParallel} package is used. Default is SerialParam().
#' @param usenperm SHILA to write
#' @param nperm number of permutation
#' @param maxDist SHILA to write
#' @param plot A logical input indicates whether the results are plotted
#' @param verbose A logical input indicates whether the intermediate steps will be printed
#'
#' @param ... parameters for function trajectoryWeightMatrix or spatialWeightMatrix
#'
#' @importFrom stats sd
#'
#' @export


scHOT <- function(scHOT,
                  testingScaffold = NULL,
                  weightMatrix = NULL,
                  positionType = NULL,
                  positionColData = NULL,
                  nrow.out = NULL,
                  averageAcrossTrajectoryTies = FALSE,
                  higherOrderFunction = NULL,
                  higherOrderFunctionType = NULL,
                  numberPermutations = 1000,
                  numberScaffold = 100,
                  storePermutations = TRUE,
                  higherOrderSummaryFunction = sd,
                  parallel = FALSE,
                  BPPARAM = BiocParallel::SerialParam(),
                  usenperm = FALSE,
                  nperm = 10000,
                  maxDist = 0.1,
                  plot = FALSE,
                  verbose = TRUE,
                  ...
) {


  # add testing scaffold


  if (ncol(testingScaffold) != 2) {
    stop("testingScaffold must be a matrix with two columns \n")
  }

  rownames(testingScaffold) <- paste(testingScaffold[, 1], testingScaffold[, 2])

  if (verbose) {
    cat("Adding testing scaffold \n")
  }

  scHOT <- scHOT_addTestingScaffold(scHOT,
                                    testingScaffold = testingScaffold)


  # set weight matrix

  if (verbose) {
    cat("Set weight matrix \n")
  }

  scHOT <- scHOT_setWeightMatrix(scHOT,
                                 weightMatrix = weightMatrix,
                                 positionColData = positionColData,
                                 positionType = positionType,
                                 nrow.out = nrow.out,
                                 averageAcrossTrajectoryTies = averageAcrossTrajectoryTies,
                                 ...)

  if (verbose) {
    cat("Calculate gobal higher order function \n")
  }

  scHOT <- scHOT_calculateGlobalHigherOrderFunction(scHOT,
                                                    higherOrderFunction = higherOrderFunction,
                                                    higherOrderFunctionType = higherOrderFunctionType)


  if (verbose) {
    cat("Calculate Higher Order Test Statistics \n")
  }

  scHOT <- scHOT_setPermutationScaffold(scHOT,
                                        numberPermutations = numberPermutations,
                                        numberScaffold = numberScaffold,
                                        storePermutations = storePermutations)




  scHOT <- scHOT_calculateHigherOrderTestStatistics(scHOT,
                                                    higherOrderSummaryFunction = higherOrderSummaryFunction)

  if (verbose) {
    cat("Perform Permutation Test \n")
  }




  scHOT <- scHOT_performPermutationTest(scHOT,
                                        verbose = verbose,
                                        parallel = parallel,
                                        BPPARAM = BPPARAM)



  if (verbose) {
    cat("Estimating p-values \n")
  }


  scHOT <- scHOT_estimatePvalues(scHOT,
                                 usenperm = usenperm,
                                 nperm = nperm,
                                 maxDist = maxDist,
                                 plot = plot,
                                 verbose = verbose)

  return(scHOT)
}
