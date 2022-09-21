#' A wrapper function to perform scHOT
#'
#' @title scHOT
#'
#' @param scHOT A scHOT object
#' @param testingScaffold A matrix with rows for each testing combination
#' @param weightMatrix A matrix indicates the weight matrix for scHOT analysis
#' @param positionType A string indicating the position type,
#' either "trajectory" or "spatial"
#' @param positionColData Either trajectory or spatial information
#' for each sample. If positionType is "trajectory"
#' then positionColData should be a character or numeric indicating
#' the subset of colData of the scHOT object.
#' If positionType is "spatial" then positionColData should be
#' a character or numeric vector indicating the subset of colData that
#' give the full spatial coordinates.
#' @param nrow.out The number of weightings to include for testing,
#' a smaller value is faster for computation
#' @param averageAcrossTrajectoryTies Logical indicating whether ties
#' in the trajectory should be given the same local weights
#' @param higherOrderFunction A function object indicates the
#' higher order function
#' @param higherOrderFunctionType is "weighted" or "unweighted",
#' determines if there
#' is a weighting argument in the higher order function
#' @param numberPermutations The number of permutations,
#' set as 1000 by default
#' @param numberScaffold The number of testing scaffolds to
#'  perform permutations, set as 100 by default
#' @param storePermutations  a logical flag on whether
#' permutation values should be saved
#' @param higherOrderSummaryFunction A functon indicating the higher order
#' summary function (default is standard deviation `sd`)
#' @param parallel A logical input indicating whether to run
#' the permutation test using multiple cores in parallel.
#' @param BPPARAM  A \code{BiocParallelParam} class object from
#' the \code{BiocParallel} package is used. Default is SerialParam().
#' @param usenperm_estimate Logical (default FALSE) if number of neighbouring
#' permutations should be used to estimate P-values, or
#' if difference of global higher order statistic should be used
#' @param nperm_estimate Number of neighbouring permutations to
#' use for p-value estimation
#' @param maxDist max difference of global higher order statistic to
#' use for p-value estimation (default 0.1)
#' @param plot A logical input indicating whether the results are plotted
#' @param verbose A logical input indicating whether the intermediate
#' steps will be printed
#'
#' @param ... parameters for function trajectoryWeightMatrix or spatialWeightMatrix
#'
#' @return A scHOT object
#'
#' @examples
#' data(MOB_subset)
#' sce_MOB_subset <- MOB_subset$sce_MOB_subset
#' scHOT_spatial <- scHOT_buildFromSCE(sce_MOB_subset,
#'                                     assayName = "logcounts",
#'                                     positionType = "spatial",
#'                                     positionColData = c("x", "y"))
#' pairs <- matrix(c("Arrb1", "Mtor", "Dnm1l", "Gucy1b3"), ncol = 2, byrow = TRUE)
#' rownames(pairs) <- apply(pairs,1,paste0,collapse = "_")
#'
#' scHOT_spatial <- scHOT(scHOT_spatial,
#'                        testingScaffold = pairs,
#'                        positionType = "spatial",
#'                        positionColData = c("x", "y"),
#'                        nrow.out = NULL,
#'                        higherOrderFunction = weightedSpearman,
#'                        higherOrderFunctionType = "weighted",
#'                        numberPermutations = 100,
#'                        higherOrderSummaryFunction = sd,
#'                        parallel = TRUE,
#'                        BPPARAM = MulticoreParam(2),
#'                        verbose = TRUE,
#'                        span = 0.05)
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
                  usenperm_estimate = FALSE,
                  nperm_estimate = 10000,
                  maxDist = 0.1,
                  plot = FALSE,
                  verbose = TRUE,
                  ...
) {
  
  
  # add testing scaffold
  
  
  # if (ncol(testingScaffold) != 2) {
  #   stop("testingScaffold must be a matrix with two columns \n")
  # }
  
  # rownames(testingScaffold) <- paste(testingScaffold[, 1], testingScaffold[, 2])
  # if (nrow(testingScaffold) == 1) {
  # rownames(testingScaffold) <- apply(testingScaffold,1,paste0, collapse = "_")
  # } else {
  rownames(testingScaffold) <- apply(testingScaffold, 1, paste, collapse = "_")
  # }
  
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
                                 averageAcrossTrajectoryTies =
                                   averageAcrossTrajectoryTies,
                                 ...)
  
  if (verbose) {
    cat("Calculating global higher order function \n")
  }
  
  scHOT <- scHOT_calculateGlobalHigherOrderFunction(scHOT,
                                                    higherOrderFunction =
                                                      higherOrderFunction,
                                                    higherOrderFunctionType =
                                                      higherOrderFunctionType)
  
  
  if (verbose) {
    cat("Calculating higher order test statistics \n")
  }
  
  scHOT <- scHOT_setPermutationScaffold(scHOT,
                                        numberPermutations = numberPermutations,
                                        numberScaffold = numberScaffold,
                                        storePermutations = storePermutations)
  
  
  
  
  scHOT <- scHOT_calculateHigherOrderTestStatistics(scHOT,
                                                    higherOrderSummaryFunction =
                                                      higherOrderSummaryFunction)
  
  if (verbose) {
    cat("Perform Permutation Test \n")
  }
  
  
  
  
  scHOT <- scHOT_performPermutationTest(scHOT,
                                        verbose = verbose,
                                        parallel = parallel,
                                        BPPARAM = BPPARAM)
  
  if (storePermutations) {
    
    if (verbose) {
      cat("Estimating p-values \n")
    }
    
    scHOT <- scHOT_estimatePvalues(scHOT,
                                   usenperm_estimate = usenperm_estimate,
                                   nperm_estimate = nperm_estimate,
                                   maxDist = maxDist,
                                   plot = plot,
                                   verbose = verbose)
  } else {
    if (verbose) {
      cat("No permutations stored, not estimating p-values \n")
    }
  }
  
  return(scHOT)
}
