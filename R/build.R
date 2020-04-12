
#' Create scHOT object from a matrix
#'
#' @title scHOT_buildFromMatrix
#'
#' @param mat A matrix with rows for genes and columns for cells
#' @param cellData A dataframe or DataFrame object with rows for cells
#' @param positionType A string indicating the position type,
#' either "trajectory" or "spatial"
#' @param positionColData Strings indicate the position information
#' stored in colData. If positionType is "trajectory" then
#' positionColData should be a sortable vector
#' if positionType is "spatial" then positionColData should
#' be a matrix type object.
#'
#' @return A scHOT object
#'
#' @importFrom S4Vectors DataFrame SimpleList
#' @importFrom SingleCellExperiment SingleCellExperiment
#' @importFrom methods as validObject is
#'
#' @examples
#'
#' dat <- rbind(rnorm(50), rnorm(50), rnorm(50))
#' colnames(dat) <- paste0("cell_", 1:ncol(dat))
#' rownames(dat) <- c("gene_1","gene_2", "gene_2")
#'
#' scHOT <-scHOT_buildFromMatrix(dat, cellData = data.frame(1:ncol(dat)))
#'
#' @export



scHOT_buildFromMatrix <- function(mat, cellData = NULL, positionType = NULL,
                                  positionColData = NULL) {
  if (is.null(cellData)) {
    cellData <- S4Vectors::DataFrame(row.names = colnames(mat))
  }

  if (!is(cellData, "DataFrame")) {
    cellData <- S4Vectors::DataFrame(cellData)
  }

  SCE <- SingleCellExperiment(assays = S4Vectors::SimpleList(expression = mat),
                              colData = cellData)

  scHOT <- methods::as(SCE, "scHOT")

  scHOT@positionType <- positionType

  scHOT@positionColData <- positionColData

  methods::validObject(scHOT)

  return(scHOT)
}

###################################################################



#' Create scHOT object from a SingleCellExperiment object
#'
#' @title scHOT_buildFromSCE
#'
#' @param sce A SingleCellExperiment object
#' @param assayName is a single assay to pull out from sce as
#' the expression matrix input of scHOT
#' @param positionType A string indicates the position type,
#' either trajectory or spatial
#' @param positionColData Strings indicate the position information
#' stored in colData.
#' If positionType is "trajectory" then positionColData
#' should be a sortable vector
#' if positionType is "spatial" then positionColData
#' should be a matrix type object.
#'
#'
#' @return A scHOT object
#'
#'
#' @importFrom S4Vectors SimpleList
#' @importFrom SingleCellExperiment SingleCellExperiment
#' @importFrom SummarizedExperiment colData assay
#' @importFrom methods as validObject
#' @examples
#'
#' library(SingleCellExperiment)
#' dat <- rbind(rnorm(50), rnorm(50), rnorm(50))
#' colnames(dat) <- paste0("cell_", 1:ncol(dat))
#' rownames(dat) <- c("gene_1","gene_2", "gene_2")
#'
#' sce <- SingleCellExperiment::SingleCellExperiment(assays =
#' S4Vectors::SimpleList(counts = dat))
#' scHOT <- scHOT_buildFromSCE(sce)
#'
#' @export



scHOT_buildFromSCE <- function(sce,
                               assayName = "counts",
                               positionType = NULL,
                               positionColData = NULL) {

  # will keep the colData but nothing else

  SCE <- SingleCellExperiment(assays = SimpleList(expression =
                                                    assay(sce, assayName)),
                              colData = colData(sce))


  scHOT <- methods::as(SCE, "scHOT")

  scHOT@positionType <- positionType

  scHOT@positionColData <- positionColData

  methods::validObject(scHOT)

  return(scHOT)
}



#' Add a testing scaffold to a scHOT object
#'
#' @title scHOT_addTestingScaffold
#'
#' @param scHOT A scHOT object
#' @param testingScaffold A matrix with rows for each testing combination,
#' and columns for level of dimensionality (1 for single gene etc.)
#'
#' @return A scHOT object with slot testingScaffold saved
#'
#' @examples
#'
#' data(MOB_subset)
#' sce_MOB_subset <- MOB_subset$sce_MOB_subset
#' scHOT_spatial <- scHOT_buildFromSCE(sce_MOB_subset,
#' assayName = "logcounts",
#' positionType = "spatial",
#' positionColData = c("x", "y"))
#' pairs <- matrix(c("Arrb1", "Mtor", "Dnm1l", "Gucy1b3"), ncol = 2, byrow = TRUE)
#' scHOT_spatial <- scHOT_addTestingScaffold(scHOT_spatial, pairs)
#'
#'
#' @export


scHOT_addTestingScaffold <- function(scHOT, testingScaffold) {
  # scaffold is a matrix with rows for each testing combination
  # and columns for level of dimensionality (1 for single gene etc.)

  # check if any are not found in data
  genesFound <- rowSums(apply(testingScaffold, seq_len(2),
                              function(x) x %in% rownames(scHOT)))

  if (any(genesFound != ncol(testingScaffold))) {
    warning("Some genes not found in data, removing these..")
  }

  # Remove the genes that are not in data
  testingScaffold <- testingScaffold[genesFound == ncol(testingScaffold), , drop = FALSE]

  colnames(testingScaffold) <- paste0("gene_", seq_len(ncol(testingScaffold)))

  scHOT@testingScaffold <- testingScaffold

  return(scHOT)
}


###################################################################


#' Create weight matrix for trajectory data
#'
#' @title trajectoryWeightMatrix
#'
#' @param n indicates the number of cels
#' @param type Type of weight matrix, one of "triangular" (default),
#' "block", and "harmonic"
#' @param span proportion of samples to include on either side, default is 0.25
#'
#' @return A weighted matrix
#'
#' @examples
#'
#' W <- trajectoryWeightMatrix(100)
#' W <- trajectoryWeightMatrix(100, type = "triangular")
#' W <- trajectoryWeightMatrix(100, type = "block")
#' W <- trajectoryWeightMatrix(100, type = "harmonic")
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
    message("type given is not one of \"harmonic\",\"triangular\",
            or \"block\", defaulting to \"triangular\"")
  }

  if (type != "harmonic") {
    if (is.null(span)) {
      span = 0.25
      message("span not specified, defaulting to 0.25")
    }
    if (span < 0 | span > 1) {
      span = 0.25
      message("span specified not between 0 and 1, defaulting to 0.25")
    }
  }

  W = matrix(0, nrow = n, ncol = n)
  spanSamples = ceiling(span * n)

  if (type == "harmonic") {
    for (i in seq_len(n)) {
      j = seq_len(n)
      W[i, j] <- 1/(1 + abs(i - j))
    }
  }

  if (type == "triangular") {
    for (i in seq_len(n)) {
      j = seq_len(n)
      W[i, j] <- (abs(j - i) <= spanSamples) *
        ((spanSamples + 1) - abs(j - i))/(spanSamples + 1)
    }
  }

  if (type == "block") {
    for (i in seq_len(n)) {
      j = seq_len(n)
      W[i, j] <- (abs(j - i) <= spanSamples) * 1
    }
  }

  return(W)
}

###################################################################


#' Create weight matrix for spatial data
#'
#' @title spatialWeightMatrix
#'
#' @param x a matrix with rows corresponding to cells and columns
#' corresponding to dimensions to calculate Euclidean distance
#' @param span proportion of samples to include on either side,
#' default is 13/(number of rows in `x`), corresponding roughly
#' to points within a diamond shape distance away
#'
#' @return A weighted matrix
#'
#' @examples
#'
#' spat_x <- rnorm(50)
#' spat_y <- rnorm(50)
#' spat_coord <- cbind(spat_x, spat_y)
#' W <- spatialWeightMatrix(spat_coord)
#'
#' @importFrom stats dist
#'
#' @export

spatialWeightMatrix <- function(x, span = NULL) {

  # x is a matrix with rows corresponding to cells and columns
  # corresponding to dimensions to calculate distance

  n = nrow(x)

  if (is.null(span)) {
    span = 13/n
    if (span > 1) span = 0.5
    message("span not specified, defaulting to ", round(span, 2))
  }

  # calculate euclidean distance of points on 2D
  # coords = cbind(x,y)
  coords = as.matrix(x)
  d = as.matrix(stats::dist(coords))

  # extract a weights vector per cell

  W_raw = sapply(seq_len(n), function(cell) {
    dvec = d[cell,]
    vals = rep(0, n)
    #vals[order(dvec)[1:ceiling(span*n)]] = seq(1, 0, length.out = ceiling(span*n))
    vals[order(dvec)[seq_len(ceiling(span*n))]] =
      seq(1, 0, length.out = ceiling(span*n))
    return(vals)
  }, simplify = FALSE)

  W = do.call(rbind, W_raw)

  return(W)

}

###################################################################


#' The thin function extracts the rows of a matrix evenly so that
#' roughly n number of rows remain. Used for thinning down
#' the weight matrix to speed up overall computation.
#'
#' @title thin
#' @param W matrix
#' @param n rough number of rows to keep
#'
#' @return \code{matrix} of thinned matrix keeping only roughly n rows.
#'
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

  W <- W[index,]
  rownames(W) <- index
  return(W)
}


###################################################################


#' Create scHOT object from a SingleCellExperiment object
#'
#' @title scHOT_setWeightMatrix
#'
#' @param scHOT A scHOT object
#' @param weightMatrix A matrix indicating the weight matrix for
#' scHOT analysis, such as the output from `trajectoryWeightMatrix`
#' or `spatialWeightMatrix`. If this is not NULL then other parameters
#' are ignored.
#' @param positionType A string indicating the position type,
#' either "trajectory" or "spatial"
#' @param positionColData Either trajectory or spatial information
#' for each sample. If positionType is "trajectory"
#' then positionColData should be a character or
#' numeric indicating the subset of colData of the scHOT object.
#' If positionType is "spatial" then positionColData
#' should be a character or numeric vector indicating
#' the subset of colData that give the full spatial coordinates.
#' @param nrow.out The number of weightings to include for testing,
#' a smaller value is faster for computation
#' @param averageAcrossTrajectoryTies Logical indicating whether ties
#' in the trajectory should be given the same local weights
#' @param ... parameters for function trajectoryWeightMatrix
#'  or spatialWeightMatrix
#'
#' @return A scHOT object with slot weightMatrix saved
#'
#' @examples
#'
#'  data(MOB_subset)
#'  sce_MOB_subset <- MOB_subset$sce_MOB_subset
#'  scHOT_spatial <- scHOT_buildFromSCE(sce_MOB_subset,
#'                                      assayName = "logcounts",
#'                                     positionType = "spatial",
#'                                      positionColData = c("x", "y"))
#'  pairs <- matrix(c("Arrb1", "Mtor", "Dnm1l", "Gucy1b3"), ncol = 2, byrow = TRUE)
#'  rownames(pairs) <- apply(pairs,1,paste0,collapse = "_")
#'  scHOT_spatial <- scHOT_addTestingScaffold(scHOT_spatial, pairs)
#'  scHOT_spatial <- scHOT_setWeightMatrix(scHOT_spatial,
#'                                        positionColData = c("x","y"),
#'                                         positionType = "spatial",
#'                                         nrow.out = NULL,
#'                                         span = 0.05)
#'
#' @importFrom methods as is
#'
#' @export

scHOT_setWeightMatrix <- function(scHOT,
                                  weightMatrix = NULL,
                                  positionType = NULL,
                                  positionColData = NULL,
                                  nrow.out = NULL,
                                  averageAcrossTrajectoryTies = FALSE,
                                  ...) {

  # either set own weight matrix in object or generate a new one using parameters
  # if positionType is "trajectory" then positionColData should be a sortable vector
  # if positionType is "spatial" then positionColData should be a matrix type object


  # input check

  if (!is(scHOT, "scHOT")) {
    stop("scHOT needs to be scHOT class")
  }

  if (is.null(positionType)) {

    if (is.null(scHOT@positionType)) {
      stop("Both positionType and scHOT@positionType are NULL.")
    } else {
      positionType <- scHOT@positionType
    }

  }

  positionType <- match.arg(positionType, c("trajectory","spatial"), several.ok = FALSE)



  if (is.null(positionColData)) {

    if (is.null(scHOT@positionColData)) {
      stop("Both positionColData and scHOT@positionColData are NULL.")
    } else {
      positionColData <- scHOT@positionColData
    }

  }

  if (!all(positionColData %in% colnames(colData(scHOT)))) {
    stop("at least one positionColData column names not found in colData(scHOT)")
  }


  if (!is.null(weightMatrix)) {

    if (ncol(weightMatrix) != ncol(scHOT)) {
      stop("Provided weightMatrix must have same number of columns as scHOT object!")
    }

    rownames(weightMatrix) <- seq_len(nrow(weightMatrix))

  } else {


    message("weightMatrix not provided, generating one using parameter settings...")


    n = ncol(scHOT)

    # if (!positionType %in% c("trajectory","spatial")) {
    #   stop("select either \"trajectory\", or \"spatial\" as positionType")
    # }


    if (positionType == "trajectory") {

      if (length(positionColData) > 1) {
        stop("only one positionColData column name to be given for trajectory")
      }

      weightMatrix = trajectoryWeightMatrix(n, ...)
      colnames(weightMatrix) <- colnames(scHOT)
      weightMatrix = weightMatrix[ , order(colData(scHOT)[, positionColData],
                                           na.last = NA), drop = FALSE]

      if (averageAcrossTrajectoryTies) {
        # average across ties
        weightMatrixTied = t(apply(weightMatrix,
                                   1, function(x,
                                               y = colData(scHOT)[, positionColData]){
                                     unsplit(tapply(x,y,mean),y)
                                   }))
        colnames(weightMatrixTied) <- colnames(weightMatrix)
        weightMatrix = weightMatrixTied
      }
    }

    if (positionType == "spatial") {
      weightMatrix = spatialWeightMatrix(as.matrix(colData(scHOT)[, positionColData]), ...)
      colnames(weightMatrix) <- colnames(scHOT)
      # rownames(weightMatrix) <- colnames(scHOT)
    }

    if (!is.null(nrow.out)) {
      weightMatrix <- thin(weightMatrix, n = nrow.out)
    } else {
      rownames(weightMatrix) <- seq_len(nrow(weightMatrix))
    }

  }

  weightMatrix <- methods::as(weightMatrix, "dgCMatrix")

  scHOT@weightMatrix <- weightMatrix


  scHOT@positionType <- positionType
  scHOT@positionColData <- positionColData

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

  sapply(seq_len(nrow(weightMatrix)), function(i) {
    apply(testingScaffold, 1, function(x) {
      argvals = split(expressionData[x,
                                     weightMatrix[i, ,
                                                  drop = FALSE] > 0,
                                     drop = FALSE],
                      seq_len(length(x)))
      #names(argvals) <- names(formals(higherOrderFunction))[1:length(argvals)]
      names(argvals) <- names(formals(higherOrderFunction))[seq_len(length(argvals))]
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
  #  weightedFunction must have sample weighting as its last argument

  sapply(seq_len(nrow(weightMatrix)), function(i) {
    apply(testingScaffold, 1, function(x) {
      # argvals = split(rbind(expressionData[x, , drop = FALSE],
      #                       w = weightMatrix[i, , drop = FALSE]),
      #                 1:(length(x) + 1))
      argvals = split(rbind(expressionData[x, , drop = FALSE],
                            w = weightMatrix[i, , drop = FALSE]),
                      seq_len((length(x) + 1)))
      #names(argvals) <- names(formals(weightedHigherOrderFunction))[1:length(argvals)]
      names(argvals) <- names(formals(weightedHigherOrderFunction))[seq_len(length(argvals))]
      return(do.call(weightedHigherOrderFunction, args = argvals))
    })
  })
}


###################################################################


#' Calculates the global higher order function
#'
#' @title scHOT_calculateGlobalHigherOrderFunction
#'
#' @description this calculates the global higher order function
#' and stores it in the output
#' if these aren't found in the params slot
#' then they need to be specified here
#'
#' @param scHOT A scHOT object
#' @param higherOrderFunction A function object indicating the higher order function
#' @param higherOrderFunctionType is "weighted" or "unweighted", determines if there
#' is a weighting argument in the higher order function
#'
#' @return A scHOT object with scHOT_output$globalHigherOrderFunction in
#' slot scHOT_output saved
#' @examples
#'
#'  data(MOB_subset)
#'  sce_MOB_subset <- MOB_subset$sce_MOB_subset
#'  scHOT_spatial <- scHOT_buildFromSCE(sce_MOB_subset,
#'                                      assayName = "logcounts",
#'                                     positionType = "spatial",
#'                                      positionColData = c("x", "y"))
#' pairs <- matrix(c("Arrb1", "Mtor", "Dnm1l", "Gucy1b3"), ncol = 2, byrow = TRUE)
#' rownames(pairs) <- apply(pairs,1,paste0,collapse = "_")
#' scHOT_spatial <- scHOT_addTestingScaffold(scHOT_spatial, pairs)
#'  scHOT_spatial <- scHOT_setWeightMatrix(scHOT_spatial,
#'                                        positionColData = c("x","y"),
#'                                         positionType = "spatial",
#'                                         nrow.out = NULL,
#'                                         span = 0.05)
#' scHOT_spatial <- scHOT_calculateGlobalHigherOrderFunction(
#'   scHOT_spatial,
#'   higherOrderFunction = weightedSpearman,
#'   higherOrderFunctionType = "weighted")
#'
#' @importFrom SummarizedExperiment assay
#' @importFrom methods is
#'
#'
#' @export




scHOT_calculateGlobalHigherOrderFunction <- function(
  scHOT,
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
    stop("higherOrderFunctionType must be either \"weighted\"
         or \"unweighted\"")
  }

  if (is.null(higherOrderFunctionType)) {
    if (is.null(scHOT@params$higherOrderFunctionType)) {
      stop("higherOrderFunctionType not given and not saved in params,
           please provide as either \"weighted\" or \"unweighted\"")
    } else {
      higherOrderFunctionType <- scHOT@params$higherOrderFunctionType
    }
  }

  if (!is.null(higherOrderFunction)) {
    message("higherOrderFunction given will replace any stored param")
    scHOT@params$higherOrderFunction <- higherOrderFunction
  }
  if (!is(higherOrderFunction, "function")) {
    stop("higherOrderFunction must be a function object")
  }
  if (is.null(higherOrderFunction)) {
    if (is.null(scHOT@params$higherOrderFunction)) {
      stop("higherOrderFunction not given and not saved in params,
           please provide as function object")
    } else {
      higherOrderFunction <- scHOT@params$higherOrderFunction
    }
  }

  testingScaffold = scHOT@testingScaffold

  if (is.null(testingScaffold)) {
    stop("No testingScaffold found in the scHOT object,
         please provide one!")
  }

  expressionData = SummarizedExperiment::assay(scHOT, "expression")

  if (is.null(expressionData)) {
    stop("No expressionData found, please provide an assay(scHOT,
         \"expression\") slot")
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

  scHOT <- scHOT_stripOutput(scHOT, force = FALSE)

  scHOT@scHOT_output$globalHigherOrderFunction <- global

  return(scHOT)
}

###################################################################

stratifiedSample = function(stats, length = 100) {
  nsamps = 5
  ranges = range(stats[is.finite(stats)])
  fac = cut(stats, seq(from = ranges[1] - 1e-04, to = ranges[2],
                       length.out = ceiling(length/nsamps)))
  # sampleindices = unlist(tapply(1:length(stats), fac, function(x) {
  #   if (length(x) < nsamps)
  #     return(x)
  #   sample(x, nsamps, replace = FALSE)
  # }))
  sampleindices = unlist(tapply(seq_len(length(stats)), fac, function(x) {
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
#' @param numberPermutations The number of permutations, set as 1000 by default
#' @param numberScaffold The number of scaffold, set as 100 by default, minimum 6.
#' if you want all combinations to do permutations then set,
#' numberScaffold much higher than the testingScaffold
#' @param storePermutations  a logical flag on whether
#' Permutations should be stored, or discarded once used
#'
#' @return A scHOT object with storePermutations in slot scHOT_output saved
#'
#' @examples
#'  data(MOB_subset)
#'  sce_MOB_subset <- MOB_subset$sce_MOB_subset
#'  scHOT_spatial <- scHOT_buildFromSCE(sce_MOB_subset,
#'                                      assayName = "logcounts",
#'                                     positionType = "spatial",
#'                                      positionColData = c("x", "y"))
#' pairs <- matrix(c("Arrb1", "Mtor", "Dnm1l", "Gucy1b3"), ncol = 2, byrow = TRUE)
#' rownames(pairs) <- apply(pairs,1,paste0,collapse = "_")
#' scHOT_spatial <- scHOT_addTestingScaffold(scHOT_spatial, pairs)
#'  scHOT_spatial <- scHOT_setWeightMatrix(scHOT_spatial,
#'                                        positionColData = c("x","y"),
#'                                         positionType = "spatial",
#'                                         nrow.out = NULL,
#'                                         span = 0.05)
#' scHOT_spatial <- scHOT_calculateGlobalHigherOrderFunction(
#'   scHOT_spatial,
#'   higherOrderFunction = weightedSpearman,
#'   higherOrderFunctionType = "weighted")
#' scHOT_spatial <- scHOT_setPermutationScaffold(scHOT_spatial,
#'                                               numberPermutations = 100)
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

  scHOT <- scHOT_stripOutput(scHOT, force = FALSE)

  if (numberScaffold > nrow(scHOT@testingScaffold)) {
    message("numberScaffold set higher than the scaffold, setting permutation number for all tests")
    scHOT@scHOT_output$numberPermutations = numberPermutations
  } else {

    if (numberScaffold < 6) {
      message("numberScaffold set lower than 6, resetting to 6")
      numberScaffold <- 6
    }

    if (is.null(scHOT@scHOT_output$globalHigherOrderFunction)) {
      stop("need scHOT@scHOT_output$globalHigherOrderFunction to take random stratified sample")
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


###################################################################


#' Strip the scHOT output
#'
#' @title scHOT_stripOutput
#'
#' @param scHOT A scHOT object
#' @param force A logical indicates whther forcing stripping the scHOT output
#' @param store A logical flag on whether the scHOT should be stored as .rds file
#' @param file_name A string indicates the file name of the scHOT will be stored
#' @examples
#' data(MOB_subset)
#' sce_MOB_subset <- MOB_subset$sce_MOB_subset
#' scHOT_spatial <- scHOT_buildFromSCE(sce_MOB_subset,
#'                                      assayName = "logcounts",
#'                                     positionType = "spatial",
#'                                      positionColData = c("x", "y"))
#'
#' scHOT_spatial <- scHOT_stripOutput(scHOT_spatial)
#'
#' @return A scHOT object with scHOT_output striped
#'
#' @examples
#'
#'
#' @export



scHOT_stripOutput <- function(scHOT, force = TRUE,
                              store = FALSE,
                              file_name = NULL) {

  scHOT_output <- scHOT@scHOT_output

  if (!is.null(scHOT_output)) {
    # current_testing <- paste(scHOT@testingScaffold[, 1],
    # scHOT@testingScaffold[, 2], sep = "_")

    if (nrow(scHOT@testingScaffold) == 1) {
      current_testing <- apply(scHOT@testingScaffold, 1, paste, collapse = "_")
      current_output_testing <- apply(scHOT_output[, grepl("gene_[0-9]",
                                                           colnames(scHOT_output)),
                                                   drop = FALSE], 1, paste,
                                      collapse = "_")
    } else {
      current_testing <- apply(scHOT@testingScaffold, 1, paste, sep = "_")
      current_output_testing <- apply(scHOT_output[, grepl("gene_[0-9]",
                                                           colnames(scHOT_output)),
                                                   drop = FALSE], 1, paste,
                                      sep = "_")
    }

    # current_testing <- apply(scHOT@testingScaffold, 1, paste, sep = "_")
    # current_output_testing <- paste(scHOT_output$gene_1,
    # scHOT_output$gene_2, sep = "_")
    # current_output_testing <- apply(scHOT_output[, grepl("gene_[0-9]",
    # colnames(scHOT_output)), drop = FALSE], 1, paste, sep = "_")

    if (force | !all(current_testing %in% current_output_testing) | length(current_testing) != length(current_output_testing)) {

      warnings("The current testing scaffold does not match with the
               current scHOT_output. The scHOT_output will be stripped out.")
      scHOT_output <- DataFrame(NULL)

      if (store) {

        if (is.null(file_name)) {
          file_name <- "scHOT.Rds"
        } else {
          file_name <- paste0(file_name, ".Rds")
        }

        message("The current scHOT wil be stored as ", file_name)
        saveRDS(scHOT, file = file_name)

      }

    }

  } else {
    scHOT_output <- DataFrame(NULL)
  }



  scHOT@scHOT_output <- scHOT_output

  return(scHOT)
}

