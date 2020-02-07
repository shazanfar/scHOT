
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
#' @return scHOT A scHOT object with results stored in scHOT_output slot
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

  testingScaffold <- scHOT@testingScaffold

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



###################################################################


#' Perform permutation test
#'
#' @title scHOT_performPermutationTest
#'
#' @param scHOT A scHOT object
#' @param verbose A logical input indicates whether the intermediate steps will be printed
#' @param parallel A logical input indicates whether run the permutation test using multiple cores in parallel.
#' @param BPPARAM  A \code{BiocParallelParam} class object from the \code{BiocParallel} package is used. Default is SerialParam().
#'
#' @importFrom SummarizedExperiment assay
#' @importFrom IRanges NumericList
#' @importFrom stats p.adjust na.omit
#' @importFrom BiocParallel bplapply SerialParam
#'
#' @return scHOT A scHOT object with results stored in scHOT_output slot
#'
#' @export



scHOT_performPermutationTest <- function(scHOT,
                                         verbose = FALSE,
                                         parallel = FALSE,
                                         BPPARAM = BiocParallel::SerialParam()) {

  # do this in a for loop

  higherOrderFunctionType <- scHOT@params$higherOrderFunctionType
  if (is.null(higherOrderFunctionType)) {
    stop("No higherOrderFunctionType found in scHOT object's params slot, please provide one
         or run scHOT_calculateGlobalHigherOrderFunction!")
  }


  higherOrderFunction <- scHOT@params$higherOrderFunction
  if (is.null(higherOrderFunction)) {
    stop("No higherOrderFunction found in scHOT object's params slot, please provide one
         or run scHOT_calculateGlobalHigherOrderFunction!")  }


  higherOrderSummaryFunction <- scHOT@params$higherOrderSummaryFunction
  if (is.null(higherOrderSummaryFunction)) {
    stop(paste0("No higherOrderSummaryFunction provided or stored in scHOT object's params slot,
                  please provide one!"))
  }

  expressionData <- SummarizedExperiment::assay(scHOT, "expression")
  if (is.null(expressionData)) {
    stop(paste0("No expressionData found, please provide",
                " an assay(scHOT, \"expression\") slot"))
  }

  weightMatrix = scHOT@weightMatrix
  if (is.null(weightMatrix)) {
    stop("No weightMatrix found in scHOT object,
         please provide one or run scHOT_setWeightMatrix function!")
  }



  DF <- scHOT@scHOT_output

  if (nrow(DF) == 0) {
    stop("Need scHOT@scHOT_output object")
  }

  if (any(DF$storePermutations)) {
    DF$permutations = split(rep(NA, nrow(DF)), 1:nrow(DF))
  }

  DF$pvalPermutations = rep(NA, nrow(DF))

  for (i in 1:nrow(DF)) {

    niter = DF[i,"numberPermutations"]

    if (niter == 0) next

    if (verbose) {
      if (nrow(DF) > 100) {
        if (i == 1 | i %% 10 == 0 | i == nrow(DF)) {
          message(paste0("Permutation testing combination ", i, "...\n"))
        }
      } else {
        message(paste0("Permutation testing combination ", i, "...\n"))
      }
    }

    store = DF[i,"storePermutations"]
    obs = DF[i,"higherOrderStatistic"]

    scaffold = scHOT@testingScaffold[i,,drop = FALSE]

    exprs = expressionData[scaffold[1,],,drop = FALSE]
    nc = ncol(exprs)


    if (parallel) {

      permutations <- BiocParallel::bplapply(seq_len(niter), function(i) {

        if (higherOrderFunctionType == "unweighted") {
          higherOrderSequence <- functionOverScaffoldWithWeights(
            scaffold,
            exprs[,sample(nc),drop = FALSE],
            higherOrderFunction,
            weightMatrix)
        } else {
          # for weighted
          higherOrderSequence <- weightedFunctionOverScaffold(
            scaffold,
            exprs[,sample(nc),drop = FALSE],
            higherOrderFunction,
            weightMatrix)
        }
        higherOrderSummaryFunction(higherOrderSequence)

      }, BPPARAM = BPPARAM)

      permutations <- unlist(permutations)

    } else {

      permutations <- replicate(niter, {

        if (higherOrderFunctionType == "unweighted") {
          higherOrderSequence <- functionOverScaffoldWithWeights(
            scaffold,
            exprs[,sample(nc),drop = FALSE],
            higherOrderFunction,
            weightMatrix)
        } else {
          # for weighted
          higherOrderSequence <- weightedFunctionOverScaffold(
            scaffold,
            exprs[,sample(nc),drop = FALSE],
            higherOrderFunction,
            weightMatrix)
        }
        higherOrderSummaryFunction(higherOrderSequence)
      })

    }


    if (store) {
      DF$permutations[[i]] <- permutations
    }

    pval = mean(permutations >= obs, na.rm = TRUE)
    if (pval == 0) {
      pval <- 1/(length(stats::na.omit(permutations)) + 1)
    }

    DF[i,"pvalPermutations"] <- pval

  }

  DF$FDRPermutations = stats::p.adjust(DF$pvalPermutations, method = "BH")

  scHOT@scHOT_output <- DF

  if (store) {
    scHOT@scHOT_output$permutations <- IRanges::NumericList(
      scHOT@scHOT_output$permutations)
  }

  return(scHOT)

}



###################################################################
#'
#' @importFrom graphics plot
#' @importFrom reshape sort_df
#'

estimatePvalues <- function(stats,
                            globalCors,
                            permstats,
                            usenperm = FALSE,
                            nperm = 10000,
                            maxDist = 0.1,
                            plot = FALSE,
                            verbose = FALSE) {

  if (plot) {
    graphics::plot(globalCors[names(permstats)],
                   unlist(lapply(permstats,
                                 function(x) mean(unlist(x)))),
                   xlab = "Global Spearman correlation",
                   ylab = "Mean of permuted statistics")
  }


  # melt permstats list
  permstatsDF = data.frame(
    genepair = rep(names(permstats),
                   times = unlist(
                     lapply(permstats, function(x) length(unlist(x))))),
    stat = unlist(permstats))

  # add globalCor in the permstats data frame
  permstatsDF$globalCor = globalCors[as.character(permstatsDF$genepair)]


  # this can be do it in parallel
  results = sapply(names(stats), function(genepair) {

    if (verbose) {
      print(genepair)
    }

    stat = stats[genepair]
    globalCor = globalCors[genepair]

    permstatsDF$dist = abs(permstatsDF$globalCor - globalCor)

    if (usenperm) {
      if (nperm >= nrow(permstatsDF)) {
        message(paste0("nperm given is larger than or equal to total",
                       " number of permutations... using all",
                       " permutations"))
        # permstatsDF_sorted_sub = permstatsDF
        permstatsDF_sorted_sub = na.omit(permstatsDF)
      }
      else {
        permstatsDF_sorted = reshape::sort_df(na.omit(permstatsDF), "dist")
        permstatsDF_sorted_sub = permstatsDF_sorted[1:nperm, ]
      }
    } else {
      if (is.null(maxDist)) {
        message(paste0("usenperm set as FALSE, but maxDist not",
                       " given... defaulting to maxDist = 0.1"))
        maxDist = 0.1
      }

      permstatsDF_sorted_sub = subset(na.omit(permstatsDF), dist <
                                        maxDist)
      if (nrow(permstatsDF_sorted_sub) == 0) {
        message(paste0("no permutations found within given",
                       " maxDist... using all permutations instead"))
        permstatsDF_sorted_sub = na.omit(permstatsDF)
      }
    }

    permstats_sub = permstatsDF_sorted_sub$stat


    if (length(permstats_sub) == 0)
      stop("please set a larger maxDist, or set a larger usenperm")

    pval = mean(permstats_sub >= stat)
    if (pval == 0) {
      pval = 1/(1 + length(permstats_sub))
    }
    nperm = nrow(permstatsDF_sorted_sub)
    globalCorMin = min(permstatsDF_sorted_sub$globalCor)
    globalCorMax = max(permstatsDF_sorted_sub$globalCor)
    return(data.frame(genepair = genepair,
                      stat = stat,
                      globalCor = globalCor,
                      pval = pval,
                      nperm = nperm,
                      globalCorMin = globalCorMin,
                      globalCorMax = globalCorMax))
  }, simplify = FALSE)

  resultsDF = do.call(rbind, results)
  return(resultsDF)
}



###################################################################

#' Estimate p-values based on already run permutation tests
#'
#' @title scHOT_estimatePvalues
#'
#' @param scHOT A scHOT object
#' @param usenperm SHILA to write
#' @param nperm number of permutation
#' @param maxDist SHILA to write
#' @param plot A logical input indicates whether the results are plotted
#' @param verbose A logical input indicates whether the intermediate steps will be printed
#'
#' @return scHOT A scHOT object with results stored in scHOT_output slot
#'
#'
#' @importFrom stats p.adjust
#'
#' @export
#'

scHOT_estimatePvalues <- function(scHOT,
                                  usenperm = FALSE,
                                  nperm = 10000,
                                  maxDist = 0.1,
                                  plot = FALSE,
                                  verbose = FALSE) {

  stats = scHOT@scHOT_output$higherOrderStatistic
  globalCors = scHOT@scHOT_output$globalHigherOrderFunction
  permstats = as.list(scHOT@scHOT_output$permutations)


  # check the input

  if (is.null(stats)) {
    stop("No higherOrderStatistic found in scHOT object's scHOT_output slot,
         please run function scHOT_calculateHigherOrderTestStatistics!")
  }

  if (is.null(globalCors)) {
    stop("No globalHigherOrderFunction found in scHOT object's scHOT_output slot,
         please run function scHOT_calculateGlobalHigherOrderFunction!")
  }

  if (is.null(permstats)) {
    stop("No permutations found in scHOT object's scHOT_output slot,
         please run function scHOT_performPermutationTest!")
  }




  if (is.null(names(stats))) {
    names(stats) <- paste0("stat_", 1:length(stats))
  }

  names(globalCors) <- names(stats)
  names(permstats) <- names(stats)

  out = estimatePvalues(stats,
                        globalCors,
                        permstats,
                        usenperm,
                        nperm,
                        maxDist,
                        plot,
                        verbose)


  scHOT@scHOT_output$numberPermutationsEstimated = out$nperm
  scHOT@scHOT_output$globalLowerRangeEstimated = out$globalCorMin
  scHOT@scHOT_output$globalUpperRangeEstimated = out$globalCorMax
  scHOT@scHOT_output$pvalEstimated = out$pval
  scHOT@scHOT_output$FDREstimated = stats::p.adjust(
    scHOT@scHOT_output$pvalEstimated, method = "BH")

  return(scHOT)
}

