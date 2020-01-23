



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

