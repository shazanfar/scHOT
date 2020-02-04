
#' Create scHOT object from a matrix
#'
#' @title scHOT_buildFromMatrix
#'
#' @param mat A matrix with rows for genes and columns for cells
#' @param cellData A dataframe or DataFrame object with rows for cells
#'
#' @importFrom S4Vectors DataFrame SimpleList
#' @importFrom SingleCellExperiment SingleCellExperiment
#' @importFrom methods as
#'
#' @examples
#'
#' dat <- rbind(rnorm(50), rnorm(50), rnorm(50))
#' colnames(dat) <- paste0("cell_", 1:ncol(dat))
#' rownames(dat) <- c("T","Gata1", "Tal1")
#'
#' scHOT <-scHOT_buildFromMatrix(dat, cellData = data.frame(1:ncol(dat)))
#'
#' @export



scHOT_buildFromMatrix <- function(mat, cellData = NULL) {

  if (is.null(cellData)) {
    cellData <- S4Vectors::DataFrame(row.names = colnames(mat))
  }

  if (!class(cellData) == "DataFrame") {
    cellData <- S4Vectors::DataFrame(cellData)
  }

  SCE <- SingleCellExperiment::SingleCellExperiment(assays = S4Vectors::SimpleList(expression = mat),
                                                    colData = cellData)

  methods::as(SCE, "scHOT")
}

###################################################################



#' Create scHOT object from a SingleCellExperiment object
#'
#' @title scHOT_buildFromSCE
#'
#' @param sce A SingleCellExperiment object
#' @param assayName is a single assay to pull out from sce
#'
#' @importFrom S4Vectors SimpleList
#' @importFrom SingleCellExperiment SingleCellExperiment
#' @importFrom SummarizedExperiment colData assay
#' @importFrom methods as
#' @examples
#'
#' library(SingleCellExperiment)
#' dat <- rbind(rnorm(50), rnorm(50), rnorm(50))
#' colnames(dat) <- paste0("cell_", 1:ncol(dat))
#' rownames(dat) <- c("T","Gata1", "Tal1")
#'
#' sce <- SingleCellExperiment::SingleCellExperiment(assays = S4Vectors::SimpleList(counts = dat))
#' scHOT <- scHOT_buildFromSCE(sce)
#'
#' @export



scHOT_buildFromSCE <- function(sce, assayName = "counts") {

  # will keep the colData but nothing else

  SCE <- SingleCellExperiment::SingleCellExperiment(assays = S4Vectors::SimpleList(expression = SummarizedExperiment::assay(sce, assayName)),
                                                    colData = SummarizedExperiment::colData(sce))

  methods::as(SCE, "scHOT")
}



#' Add a testing scaffold to a scHOT object
#'
#' @title scHOT_addTestingScaffold
#'
#' @param scHOT A scHOT object
#' @param testingScaffold A matrix with rows for each testing combination
#'
#' @export


scHOT_addTestingScaffold <- function(scHOT, testingScaffold) {
  # scaffold is a matrix with rows for each testing combination
  # and columns for level of dimensionality (1 for single gene etc.)

  # check if any are not found in data
  genesFound <- rowSums(apply(testingScaffold, 1:2, function(x) x %in% rownames(scHOT)))

  if (any(genesFound != ncol(testingScaffold))) {
    warning("Some genes not found in data, removing these..")
  }

  # Remove the genes that are not in data
  testingScaffold <- testingScaffold[genesFound == ncol(testingScaffold), , drop = FALSE]

  colnames(testingScaffold) <- paste0("gene_", 1:ncol(testingScaffold))

  scHOT@testingScaffold <- testingScaffold

  return(scHOT)
}


###################################################################


#' Create weight matrix for trajectory data
#'
#' @title trajectoryWeightMatrix
#'
#' @param n indicates the number of cels
#' @param type Type of weight matrix, one of "harmonic", "triangular", "block"
#' @param span proportion of samples to include on either side, default is 0.5
#'
#' @examples
#'
#' trajectoryWeightMatrix(100)
#' trajectoryWeightMatrix(100, type = "triangular")
#' trajectoryWeightMatrix(100, type = "block")
#' trajectoryWeightMatrix(100, type = "harmonic")
#'
#' @export



trajectoryWeightMatrix <- function(n, type = NULL, span = NULL) {
  # n indicates the number of cels


  if (is.null(type)) {
    type = "triangular"
    message("type not specified, defaulting to triangular")
  }

  if (!type %in% (c("harmonic", "triangular", "block"))) {
    type = "triangular"
    message("type given is not one of \"harmonic\",\"triangular\", or \"block\", defaulting to \"triangular\"")
  }

  if (type != "harmonic") {
    if (is.null(span)) {
      span = 0.5
      message("span not specified, defaulting to 0.5")
    }
    if (span < 0 | span > 1) {
      span = 0.5
      message("span specified not between 0 and 1, defaulting to 0.5")
    }
  }

  W = matrix(0, nrow = n, ncol = n)
  spanSamples = ceiling(span * n)

  if (type == "harmonic") {
    for (i in 1:n) {
      for (j in 1:n) {
        W[i, j] <- 1/(1 + abs(i - j))
      }
    }
  }

  if (type == "triangular") {
    for (i in 1:n) {
      for (j in 1:n) {
        W[i, j] <- (abs(j - i) <= spanSamples) * ((spanSamples +
                                                     1) - abs(j - i))/(spanSamples + 1)
      }
    }
  }

  if (type == "block") {
    for (i in 1:n) {
      for (j in 1:n) {
        W[i, j] <- (abs(j - i) <= spanSamples) * 1
      }
    }
  }

  return(W)
}

###################################################################


#' Create weight matrix for spatial data
#'
#' @title spatialWeightMatrix
#'
#' @param x a matrix with rows corresponding to cells and columns corresponding to dimensions to calculate distance
#' @param span proportion of samples to include on either side, default is 0.5
#'
#' @importFrom stats dist
#'
#' @export

spatialWeightMatrix <- function(x, span = NULL) {

  # x is a matrix with rows corresponding to cells and columns corresponding to dimensions to calculate distance

  n = nrow(x)

  if (is.null(span)) {
    span = 0.5
    message("span not specified, defaulting to 0.5")
  }

  # calculate euclidean distance of points on 2D
  # coords = cbind(x,y)
  coords = as.matrix(x)
  d = as.matrix(stats::dist(coords))

  # extract a weights vector per cell

  W_raw = sapply(seq_len(n), function(cell) {
    dvec = d[cell,]
    vals = rep(0, n)
    vals[order(dvec)[1:ceiling(span*n)]] = seq(1, 0, length.out = ceiling(span*n))
    return(vals)
  }, simplify = FALSE)

  W = do.call(rbind, W_raw)

  return(W)

}

###################################################################


#' The thin function extracts the rows of a matrix evenly so that roughly n number of rows remain. Used for thinning down the weight matrix to speed up overall computation.
#'
#' @title thin
#' @param W matrix
#' @param n rough number of rows to keep
#' @return \code{matrix} of thinned matrix keeping only roughly n rows.

#' @examples
#'
#' W = trajectoryWeightMatrix(500)
#' W_small = thin(W, n = 100)
#'
#' @export

thin = function(W, n = 100) {
  # given a matrix W, extract the rows so that
  # the total number of rows is roughly n

  N = nrow(W)

  if (n == N) return(W)

  if (n > N) stop("need to set n to less than nrow(W)")

  new_n = floor(N/n)

  index = seq(from = 1, to = N, by = new_n)

  return(W[index,])
}


###################################################################


#' Create scHOT object from a SingleCellExperiment object
#'
#' @title scHOT_setWeightMatrix
#'
#' @param scHOT A scHOT object
#' @param weightMatrix A matrix indicates the weight matrix for scHOT analysis
#' @param positionType A string indicates the position type, either trajectory or spatial
#' @param cellPosition If positionType is "trajectory" then cellPosition should be a sortable vector
#' if positionType is "spatial" then cellPosition should be a matrix type object.
#' It should be stored in colData (SHILA to check).
#' @param nrow.out SHILA to input
#' @param averageAcrossTrajectoryTies SHILA to input
#' @param ... parameters for function trajectoryWeightMatrix or spatialWeightMatrix
#'
#' @importFrom methods as
#'
#' @export

scHOT_setWeightMatrix <- function(scHOT,
                                  weightMatrix = NULL,
                                  positionType = "trajectory",
                                  cellPosition = NULL,
                                  nrow.out = NULL,
                                  averageAcrossTrajectoryTies = FALSE,
                                  # type = NULL,
                                  # span = NULL,
                                  ...) {

  # either set own weight matrix in object or generate a new one using parameters
  # if positionType is "trajectory" then cellPosition should be a sortable vector
  # if positionType is "spatial" then cellPosition should be a matrix type object


  # input check

  if (class(scHOT) != "scHOT") {
    stop("scHOT needs to be scHOT class")
  }

  if (!is.null(weightMatrix)) {

    if (ncol(weightMatrix) != ncol(scHOT)) {
      stop("Provided weightMatrix must have same number of columns as scHOT object!")
    }

  } else {


    message("weightMatrix not provided, generating one using parameter settings...")


    n = ncol(scHOT)

    # if (!positionType %in% c("trajectory","spatial")) {
    #   stop("select either \"trajectory\", or \"spatial\" as positionType")
    # }

    positionType <- match.arg(positionType, c("trajectory","spatial"), several.ok = FALSE)

    if (!all(cellPosition %in% colnames(colData(scHOT)))) {
      stop("at least one cellPosition column names not found in colData(scHOT)")
    }

    if (positionType == "trajectory") {

      if (length(cellPosition) > 1) {
        stop("only one cellPosition column name to be given for trajectory")
      }

      weightMatrix = trajectoryWeightMatrix(n, ...)
      colnames(weightMatrix) <- colnames(scHOT)
      weightMatrix = weightMatrix[ , order(colData(scHOT)[, cellPosition], na.last = NA), drop = FALSE]

      if (averageAcrossTrajectoryTies) {
        # average across ties
        weightMatrixTied = t(apply(weightMatrix,1,function(x, y = colData(scHOT)[, cellPosition]){
          unsplit(tapply(x,y,mean),y)
        }))
        colnames(weightMatrixTied) <- colnames(weightMatrix)
        weightMatrix = weightMatrixTied
      }
    }

    if (positionType == "spatial") {
      weightMatrix = spatialWeightMatrix(as.matrix(colData(scHOT)[, cellPosition]), ...)
      colnames(weightMatrix) <- colnames(scHOT)
    }

    if (!is.null(nrow.out)) {
      weightMatrix = thin(weightMatrix, n = nrow.out)
    }

  }

  weightMatrix <- methods::as(weightMatrix, "dgCMatrix")
  scHOT@weightMatrix <- weightMatrix

  return(scHOT)
}

###################################################################


functionOverScaffoldWithWeights = function(testingScaffold,
                                           expressionData,
                                           higherOrderFunction,
                                           weightMatrix) {

  # functionOverScaffoldWithWeights takes in a set of weights and
  # applies a function to it selectively in such a way that it's
  # equivalent to a block weight matrix set up
  # higherOrderFunction must have same number of arguments as testingScaffold

  sapply(1:nrow(weightMatrix), function(i) {
    apply(testingScaffold, 1, function(x) {
      argvals = split(expressionData[x, weightMatrix[i, , drop = FALSE] > 0, drop = FALSE],
                      1:length(x))
      names(argvals) <- names(formals(higherOrderFunction))[1:length(argvals)]
      return(do.call(higherOrderFunction, args = argvals))
    })
  })
}

###################################################################

weightedFunctionOverScaffold = function(testingScaffold,
                                        expressionData,
                                        weightedHigherOrderFunction,
                                        weightMatrix) {

  # weightedFunctionOverScaffold takes in a set of weights
  # and applies a weighted function to it
  #  weightedFunction must have sample weighting as its first argument

  sapply(1:nrow(weightMatrix), function(i) {
    apply(testingScaffold, 1, function(x) {
      argvals = split(rbind(w = weightMatrix[i, , drop = FALSE],
                            expressionData[x, , drop = FALSE]),
                      1:(length(x) + 1))
      names(argvals) <- names(formals(weightedHigherOrderFunction))[1:length(argvals)]
      return(do.call(weightedHigherOrderFunction, args = argvals))
    })
  })
}


###################################################################


#' Calculates the global higher order function
#'
#' @title scHOT_calculateGlobalHigherOrderFunction
#'
#' @description this calculates the global higher order function and stores it in the output
#' if these aren't found in the params slot then they need to be specified here
#'
#' @param scHOT A scHOT object
#' @param higherOrderFunction A function object indicates the higher order function
#' @param higherOrderFunctionType is weighted or unweighted, determines if there
#' is a weighting argument in the higher order function
#'
#' @importFrom SummarizedExperiment assay
#'
#'
#' @export




scHOT_calculateGlobalHigherOrderFunction <- function(scHOT,
                                                    higherOrderFunction = NULL,
                                                    higherOrderFunctionType = NULL) {

  # this calculates the global higher order function and stores it in the output
  # if these aren't found in the params slot then they need to be specified here
  # higherOrderFunctionType is weighted or unweighted, determines if there
  # is a weighting argument in the higher order function

  if (!is.null(higherOrderFunctionType)) {
    message("higherOrderFunctionType given will replace any stored param")
    scHOT@params$higherOrderFunctionType <- higherOrderFunctionType
  }

  if (!higherOrderFunctionType %in% c("weighted","unweighted")) {
    stop("higherOrderFunctionType must be either \"weighted\" or \"unweighted\"")
  }

  if (is.null(higherOrderFunctionType)) {
    if (is.null(scHOT@params$higherOrderFunctionType)) {
      stop("higherOrderFunctionType not given and not saved in params, please provide as either \"weighted\" or \"unweighted\"")
    } else {
      higherOrderFunctionType <- scHOT@params$higherOrderFunctionType
    }
  }

  if (!is.null(higherOrderFunction)) {
    message("higherOrderFunction given will replace any stored param")
    scHOT@params$higherOrderFunction <- higherOrderFunction
  }
  if (class(higherOrderFunction) != "function") {
    stop("higherOrderFunction must be a function object")
  }
  if (is.null(higherOrderFunction)) {
    if (is.null(scHOT@params$higherOrderFunction)) {
      stop("higherOrderFunction not given and not saved in params, please provide as function object")
    } else {
      higherOrderFunction <- scHOT@params$higherOrderFunction
    }
  }

  testingScaffold = scHOT@testingScaffold

  if (is.null(testingScaffold)) {
    stop("No testingScaffold found in the scHOT object, please provide one!")
  }

  expressionData = SummarizedExperiment::assay(scHOT, "expression")

  if (is.null(expressionData)) {
    stop("No expressionData found, please provide an assay(scHOT, \"expression\") slot")
  }

  weightMatrix_trivial = matrix(1,nrow = 1, ncol = ncol(expressionData))

  if (higherOrderFunctionType == "unweighted") {
    global = functionOverScaffoldWithWeights(testingScaffold,
                                             expressionData,
                                             higherOrderFunction,
                                             weightMatrix_trivial)
  } else {
    # for weighted
    global = weightedFunctionOverScaffold(testingScaffold,
                                          expressionData,
                                          higherOrderFunction,
                                          weightMatrix_trivial)
  }

  if (nrow(scHOT@scHOT_output) == 0) {
    scHOT@scHOT_output = DataFrame(
      testingScaffold
    )
  }
  scHOT@scHOT_output$globalHigherOrderFunction <- global

  return(scHOT)
}

###################################################################

stratifiedSample = function(stats, length = 100) {
  nsamps = 5
  ranges = range(stats[is.finite(stats)])
  fac = cut(stats, seq(from = ranges[1] - 1e-04, to = ranges[2],
                       length.out = ceiling(length/nsamps)))
  sampleindices = unlist(tapply(1:length(stats), fac, function(x) {
    if (length(x) < nsamps)
      return(x)
    sample(x, nsamps, replace = FALSE)
  }))
  sampleindices = unique(sampleindices[!is.na(sampleindices)])
  return(sampleindices)
}




###################################################################


#' Set permutation scaffold
#'
#' @title scHOT_setPermutationScaffold
#'
#' @param scHOT A scHOT object
#' @param numberPermutations The number of permutatuion, set as 1000 by default
#' @param numberScaffold The number of scaffold, set as 100 by default
#' @param storePermutations  a logical flag on whether
#' Permutations should be stored, or discarded
#'
#' @export




scHOT_setPermutationScaffold = function(scHOT,
                                        numberPermutations = 1000,
                                        numberScaffold = 100,
                                        storePermutations = TRUE) {

  # if you want all combinations to do permutations then set
  # numberScaffold much higher than the testingScaffold
  # storePermutations is a logical flag on whether
  # Permutations should be stored,
  # or discarded

  testingScaffold = scHOT@testingScaffold

  if (is.null(testingScaffold)) {
    stop("No testingScaffold found in the scHOT object, please provide one!")
  }

  if (nrow(scHOT@scHOT_output) == 0) {
    scHOT@scHOT_output = DataFrame(
      testingScaffold
    )
  }

  if (numberScaffold > nrow(scHOT@testingScaffold)) {
    message(paste0("numberScaffold set higher than the scaffold,",
                   " performing permutations for all tests"))
    scHOT@scHOT_output$numberPermutations = numberPermutations
  } else {

    if (is.null(scHOT@scHOT_output$globalHigherOrderFunction)) {
      stop(paste("need scHOT@scHOT_output$globalHigherOrderFunction to",
                 " take random stratified sample"))
    }

    scHOT@scHOT_output$numberPermutations = 0
    scHOT@scHOT_output$numberPermutations[
      stratifiedSample(
        scHOT@scHOT_output$globalHigherOrderFunction,
        length = numberScaffold)] <- numberPermutations
  }
  scHOT@scHOT_output$storePermutations = storePermutations

  return(scHOT)
}
