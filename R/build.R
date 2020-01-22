#' @title Create scHOT object from a matrix
#'
#' @param mat A matrix with rows for genes and columns for cells
#' @param cellData A dataframe or DataFrame object with rows for cells
#'
#' @importFrom S4Vectors DataFrame SimpleList
#' @importFrom SingleCellExperiment SingleCellExperiment
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

  as(SCE, "scHOT")
}

#' @title Create scHOT object from a SingleCellExperiment object
#'
#' @param sce A SingleCellExperiment object
#' @param assayName is a single assay to pull out from sce
#'
#' @importFrom S4Vectors SimpleList
#' @importFrom SingleCellExperiment SingleCellExperiment
#' @importFrom SummarizedExperiment colData assay
#' @export



scHOT_buildFromSCE <- function(sce, assayName = "counts") {

  # will keep the colData but nothing else

  SCE <- SingleCellExperiment::SingleCellExperiment(assays = S4Vectors::SimpleList(expression = SummarizedExperiment::assay(sce, assayName)),
                                                    colData = SummarizedExperiment::colData(sce))

  as(SCE, "scHOT")
}
