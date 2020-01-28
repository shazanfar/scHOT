



#' @importClassesFrom Matrix dgCMatrix
#'

setClassUnion("matrixORdgCMatrix", c("matrix", "dgCMatrix"))

#' @importClassesFrom S4Vectors DataFrame DFrame
#'

setClassUnion("data.frameORDataFrame", c("data.frame", "DataFrame"))




#' scHOT class
#'
#' @slot testingScaffold A matrix with rows for each testing combination
#' @slot weightMatrix A matrix or dgCMatrix indicates the weight matrix
#' @slot scHOT_output A data.frame or DtatFrame to store output from scHOT
#' @slot params A list of parameters
#'
#' @importClassesFrom SingleCellExperiment SingleCellExperiment
#'
#' @export



setClass("scHOT",
         slots = c(testingScaffold = "matrix",
                   weightMatrix = "matrixORdgCMatrix",
                   scHOT_output = "data.frameORDataFrame",
                   params = "list"),
         contains = "SingleCellExperiment")


#' @importFrom S4Vectors coolcat
#' @importFrom methods callNextMethod
#'

setMethod("show", "scHOT", function(object) {
  methods::callNextMethod()
  cat("testingScaffold dim:", dim(object@testingScaffold), "\n")
  cat("weightMatrix dim:", dim(object@weightMatrix), "\n")
  S4Vectors::coolcat("scHOT_output colnames (%d): %s\n", colnames(object@scHOT_output))
  S4Vectors::coolcat("param names (%d): %s\n", names(object@params))
})
