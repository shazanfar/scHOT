
###################################################################


#' Calculate and store the higherOrderSequence and higherOrderTestStatistic
#'
#' @title scHOT_calculateHigherOrderTestStatistics
#'
#' @param scHOT A scHOT object
#' @param higherOrderSummaryFunction A functon indicates the higher order summary function (SHILA to check)
#' @param ... parameters for higherOrderSummaryFunction
#'
#'
#' @importFrom SummarizedExperiment assay
#' @importFrom IRanges NumericList
#'
#' @export





scHOT_calculateHigherOrderTestStatistics <- function(
  scHOT,
  higherOrderSummaryFunction = NULL,
  ...) {
  # calculate and store the higherOrderSequence and higherOrderTestStatistic

  higherOrderFunctionType <- scHOT@params$higherOrderFunctionType
  if (is.null(higherOrderFunctionType)) {
    stop("No higherOrderFunctionType found in scHOT object's params slot, please provide one
         or run scHOT_calculateGlobalHigherOrderFunction!")
  }

  higherOrderFunction <- scHOT@params$higherOrderFunction
  if (is.null(higherOrderFunction)) {
    stop("No higherOrderFunction found in scHOT object's params slot, please provide one
         or run scHOT_calculateGlobalHigherOrderFunction!")
  }

  if (is.null(higherOrderSummaryFunction)) {
    higherOrderSummaryFunction = scHOT@params$higherOrderSummaryFunction
    if (is.null(higherOrderSummaryFunction)) {
      stop(paste0("No higherOrderSummaryFunction provided or stored in scHOT object's params slot,
                  please provide one!"))
    }
  } else {
    message("higherOrderSummaryFunction will replace any stored param")
    scHOT@params$higherOrderSummaryFunction = higherOrderSummaryFunction
  }

  testingScaffold = scHOT@testingScaffold

  if (is.null(testingScaffold)) {
    stop("No testingScaffold found in scHOT,
         please provide one or run scHOT_addTestingScaffold function!")
  }

  expressionData = SummarizedExperiment::assay(scHOT, "expression")
  if (is.null(expressionData)) {
    stop(paste0("No expressionData found, please provide an",
                " assay(scHOT, \"expression\") slot"))
  }

  weightMatrix = scHOT@weightMatrix
  if (is.null(weightMatrix)) {
    stop("No weightMatrix found in scHOT,
         please provide one or run scHOT_setWeightMatrix!")
  }

  if (higherOrderFunctionType == "unweighted") {
    higherOrderSequence = functionOverScaffoldWithWeights(
      testingScaffold,
      expressionData,
      higherOrderFunction,
      weightMatrix)
  } else {
    # for weighted
    higherOrderSequence = weightedFunctionOverScaffold(
      testingScaffold,
      expressionData,
      higherOrderFunction,
      weightMatrix)
  }

  if (nrow(scHOT@scHOT_output) == 0) {
    scHOT@scHOT_output = DataFrame(
      testingScaffold
    )
  }

  scHOT@scHOT_output$higherOrderSequence = split(
    higherOrderSequence, 1:nrow(higherOrderSequence))

  scHOT@scHOT_output$higherOrderSequence <- IRanges::NumericList(
    scHOT@scHOT_output$higherOrderSequence)

  scHOT@scHOT_output$higherOrderStatistic = unlist(
    lapply(scHOT@scHOT_output$higherOrderSequence,
           higherOrderSummaryFunction, ...))

  return(scHOT)
}
