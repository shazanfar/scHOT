



#' @importClassesFrom Matrix dgCMatrix
#'

setClassUnion("matrixORdgCMatrix", c("matrix", "dgCMatrix"))

#' @importClassesFrom S4Vectors DataFrame DFrame
#'

setClassUnion("data.frameORDataFrame", c("data.frame", "DataFrame"))


setClassUnion("characterORNULL", c("character", "NULL"))



check_validity <- function(object) {

  if (!is.null(object@positionType)) {
    if (!object@positionType %in% c("trajectory", "spatial") |
        length(object@positionType) > 1) {
      errors <- c("positionType should be either trajectory or spatial")
      return(errors)
    }
  }

  if (!is.null(object@positionColData)) {
    if (!all(object@positionColData %in% colnames(colData(object)))) {
      errors <- c("positionColData should be one of the colnames of
                  colData of the scHOT object")
      return(errors)
    }
  }

  return(TRUE)

}


#' scHOT class
#'
#' @slot testingScaffold A matrix with rows for each testing combination
#' @slot weightMatrix A matrix or dgCMatrix indicates the weight matrix
#' @slot scHOT_output A data.frame or DtatFrame to store output from scHOT
#' @slot params A list of parameters
#' @slot positionType A character indicates the type of the position,
#' either trajectory or spatial
#' @slot positionColData A vector indicates column names of colData
#' that stored the postion informaton
#'
#' @importClassesFrom SingleCellExperiment SingleCellExperiment
#'
#' @export



setClass("scHOT",
         slots = c(testingScaffold = "matrix",
                   weightMatrix = "matrixORdgCMatrix",
                   scHOT_output = "data.frameORDataFrame",
                   params = "list",
                   positionType = "characterORNULL",
                   positionColData = "characterORNULL"),
         contains = "SingleCellExperiment")


#' @importFrom S4Vectors coolcat
#' @importFrom methods callNextMethod
#'
#'

setMethod("show", "scHOT", function(object) {
  methods::callNextMethod()
  cat("testingScaffold dim:", dim(object@testingScaffold), "\n")
  cat("weightMatrix dim:", dim(object@weightMatrix), "\n")
  S4Vectors::coolcat("scHOT_output colnames (%d): %s\n",
                     colnames(object@scHOT_output))
  S4Vectors::coolcat("param names (%d): %s\n", names(object@params))
  cat("position type:", as.character(object@positionType),"\n")
})


#' @importFrom S4Vectors setValidity2
#'
setValidity2("scHOT", check_validity)


