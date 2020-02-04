# yes 	 plotColouredExpression
# yes 	 plotEgoNetwork
# yes 	 plotNetworkPathway
# yes 	 plotOrderedExpression
# yes 	 plotWCorLine
# yes 	 plotWcorsClusterPathway





##############################################

#' the plotColouredExpression function plots an n-panel scatterplot of the gene pairs split by early, mid, and late in the sample ordering.
#'
#' @title plotColouredExpression
#'
#' @param scHOT A scHOT object, where the expression data is stored in the assay slot, with assay name "expression".
#' @param genepair is either a single character string with an underscore, or a two length character vector
#' @param branches A character indicates that the colnames stored the branch information in colData
#' @param ranked_by A character indicates that the colnames stored the ranking information of the cells in colData, such as trajectory time,
#' If it is NULL, it will be ranked based on the branch information.
#' @param subsetBranch subsetBranch is a character vector containing the names of the branches to be plotted. If NULL it will plot all branches
#' @param n number of panels to split ranked samples into, default 3.
#' @param fittedline logical default TRUE, add a lm straight line to the plot
#' @return \code{ggplot} a ggplot object of scatterplots of expression split by sample ordering
#'
#'
#'
#' @importFrom SummarizedExperiment assay colData
#' @import ggplot2
#'
#'
#' @export
#'

plotColouredExpression = function(scHOT,
                                  genepair,
                                  branches = NULL,
                                  ranked_by = NULL,
                                  subsetBranch = NULL,
                                  n = 3,
                                  fittedline = TRUE) {


  if (length(genepair) == 1) {
    genepair = unlist(strsplit(genepair, "_"))
  }
  else {
    genepair = genepair[1:2]
  }

  if (length(genepair) != 2) {
    stop("genepair is either a single character string with an underscore, or a two length character vector")
  }

  branchData <- SummarizedExperiment::assay(scHOT, "expression")




  if (is.null(ranked_by)) {
    message("ranked_by information is not provided, the expression data is ranked by the branches")
  } else {

    if (!(ranked_by %in% colnames(SummarizedExperiment::colData(scHOT)))) {
      stop("branches provided is not found in colData(scHOT)")
    }

    branchData <- branchData[, order(SummarizedExperiment::colData(scHOT)[, ranked_by])]
  }



  if (is.null(branches)) {
    message("branches information is not provided")
    branchData = list(Branch = branchData)
  } else {
    if (!(branches %in% colnames(colData(scHOT)))) {
      stop("branches provided is not found in colData(scHOT)")
    }

    branch_info <- SummarizedExperiment::colData(scHOT)[, branches]

    branchData <- lapply(split(seq_along(branch_info), branch_info),
                         function(idx) branchData[, idx])[order(unique(branch_info))]
  }



  if (!all(unlist(lapply(branchData, function(x) genepair %in% rownames(x))))) {
    stop("At least one gene name not exists in the expression data")
  }



  gdf = do.call(rbind, lapply(branchData, function(branch) {
    gdf_list_1 = data.frame(Sample = colnames(branch), order = 1:ncol(branch),
                            ExpressionGene1 = branch[genepair[1], ], ExpressionGene2 = branch[genepair[2], ])
    if (n > 1) {
      # gdf_list_1$ordercut = cut(gdf_list_1$order, n, labels = unlist(ifelse(n == 3,
      #                                                                       list(c("Early", "Middle", "Late")),
      #                                                                       list(paste0("Group ", 1:n)))))


      gdf_list_1$ordercut = cut(gdf_list_1$order, n, labels = unlist(ifelse(n == 3,
                                                                     list(c("Early", "Middle", "Late")),
                                                                     ifelse(n == 2, list(c("Early", "Late")),
                                                                            list(paste0("Group ", 1:n))))))
    } else {
      gdf_list_1$ordercut = rep("All cells", times = nrow(gdf_list_1))
    }
    return(gdf_list_1)
  }))

  gdf$branch = rep(names(branchData), times = unlist(lapply(branchData,
                                                            ncol)))
  if (!is.null(subsetBranch)) {
    gdf_sub = subset(gdf, gdf$branch %in% subsetBranch)
    if (nrow(gdf_sub) == 0)
      stop("no branches with names in subsetBranch, please re-run with correct names (should match names of branchData)")
  }
  else {
    gdf_sub = gdf
  }
  g = ggplot2::ggplot(gdf_sub, aes(x = gdf_sub$ExpressionGene1, y = gdf_sub$ExpressionGene2,
                          color = order)) +
    geom_point(show.legend = FALSE, alpha = 0.7) +
    facet_grid(branch ~ ordercut, scales = "free_y") +
    scale_color_gradientn(colours = c("orange", "blue")) +
    xlab(genepair[1]) + ylab(genepair[2]) + theme_minimal() +
    NULL
  if (fittedline) {
    g = g +
      geom_smooth(colour = "black", fill = NA, linetype = "dashed",
                  method = "lm") +
      NULL
  }
  return(g)
}


