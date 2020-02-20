

##############################################

#' the plotColouredExpression function plots an n-panel scatterplot of the gene pairs split by early, mid, and late in the sample ordering.
#'
#' @title plotColouredExpression
#'
#' @param scHOT A scHOT object.
#' @param genepair is either a single character string with an underscore, or a two length character vector
#' @param branches A character indicates that the colnames stored the branch information in colData
#' @param ranked_by A character indicates that the colnames stored the ranking information of the cells in colData, such as trajectory time,
#' If it is NULL, it will be ranked based on the branch information.
#' @param subsetBranch subsetBranch is a character vector containing the names of the branches to be plotted. If NULL it will plot all branches
#' @param n number of panels to split ranked samples into, default 3.
#' @param fittedline logical default TRUE, add a lm straight line to the plot
#' @param assayName the name of the assay that are used to plot.
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
                                  fittedline = TRUE,
                                  assayName = NULL) {


  if (length(genepair) == 1) {
    genepair = unlist(strsplit(genepair, "_"))
  }
  else {
    genepair = genepair[1:2]
  }

  if (length(genepair) != 2) {
    stop("genepair is either a single character string with an underscore, or a two length character vector")
  }

  if (is.null(assayName)) {
    assayName <- "expression"
  }


  branchData <- SummarizedExperiment::assay(scHOT, assayName)




  if (is.null(ranked_by)) {
    message("ranked_by information is not provided, the expression data is ranked by the branches")
  } else {

    if (!(ranked_by %in% colnames(SummarizedExperiment::colData(scHOT)))) {
      stop("ranked_by provided is not found in colData(scHOT)")
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



##############################################

#' the plotOrderedExpression function plots expression vectors along branches and genes as ribbon plots
#'
#' @title plotOrderedExpression
#' @param scHOT A scHOT object, where the expression data is stored in the assay slot, with assay name "expression".
#' @param genepair is either a single character string with an underscore, or a two length character vector
#' @param positionType A string indicates the position type, either trajectory or spatial
#' @param branches A character indicates that the colnames stored the branch information in colData
#' @param ranked_by A character indicates that the colnames stored the ranking information of the cells in colData, such as trajectory time
#' @param xvals A character indicates that the colnames stored in colData of the x-values associated with the samples in branchData
#' @param subsetBranch subsetBranch is a character vector containing the names of the branches to be plotted. If NULL it will plot all branches
#' @param facet can either be FALSE, "branch", "gene", or "both"
#' @param positionColData A vector indicates column names of colData that stored the postion informaton (for spatial type of data)
#' @param assayName the name of the assay that are used to plot.
#'
#' @return \code{ggplot} a ggplot object for ribbon plot with points
#'
#'
#' @importFrom SummarizedExperiment assay colData
#' @import ggplot2
#' @importFrom reshape melt
#'
#' @export
#'
plotOrderedExpression = function(scHOT,
                                 genepair,
                                 positionType = NULL,
                                 branches = NULL,
                                 ranked_by = NULL,
                                 xvals = NULL,
                                 subsetBranch = NULL,
                                 facet = FALSE,
                                 positionColData = NULL,
                                 assayName = NULL) {

  # branchData is a list containing matrices of the cell expression per branch
  # assumed that the columns of each matrix in branchData is ordered by pseudotime
  # if branchData is not given as a list, it will be converted into a list containing branchData
  # gene is either a single character string with an underscore, or a two length character vector
  # xvals is a list containing the x-values associated with the samples in branchData (if NULL, samples will just be plotted against their rank)
  # subsetBranch is a character vector containing the names of the branches to be plotted. If NULL it will plot all branches
  # facet can either be FALSE, "branch", "gene", or "both"


  if (is.null(positionType)) {

    if (is.null(scHOT@positionType)) {
      stop("Both positionType and scHOT@positionType are NULL.")
    } else {
      positionType <- scHOT@positionType
    }

  }

  if (length(genepair) == 1) {
    genepair = unlist(strsplit(genepair, "_"))
  }
  else {
    genepair = genepair[1:2]
  }

  if (length(genepair) != 2) {
    stop("genepair is either a single character string with an underscore, or a two length character vector")
  }

  if (is.null(assayName)) {
    assayName <- "expression"
  }

  branchData <- SummarizedExperiment::assay(scHOT, assayName)

  if (scHOT@positionType == "trajectory") {


    if (is.null(ranked_by)) {
      message("ranked_by information is not provided, the expression data is ranked by the branches")
    } else {

      if (!(ranked_by %in% colnames(SummarizedExperiment::colData(scHOT)))) {
        stop("ranked_by provided is not found in colData(scHOT)")
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




    gdf_list = sapply(genepair, function(g) {

      gdf = do.call(rbind,lapply(branchData, function(branch){
        gdf_list_1 = data.frame(
          Sample = colnames(branch),
          order = 1:ncol(branch),
          ExpressionGene = branch[g,],
          gene = g
        )
        return(gdf_list_1)
      }))
      gdf$branch = rep(names(branchData), times = unlist(lapply(branchData, ncol)))
      return(gdf)
    }, simplify = FALSE)

    gdf = do.call(rbind, gdf_list)

    if (!is.null(subsetBranch)) {
      gdf_sub = subset(gdf, gdf$branch %in% subsetBranch)
      if (nrow(gdf_sub) == 0) stop("no branches with names in subsetBranch,
                                 please re-run with correct names (should match names of branchData)")
    } else {
      gdf_sub = gdf
    }


    if (!is.null(xvals)) {

      if (!(xval %in% colnames(SummarizedExperiment::colData(scHOT)))) {
        stop("xval provided is not found in colData(scHOT)")
      }
      xval <- lapply(split(seq_along(branch_info), branch_info),
                     function(idx) xval[, idx])[order(unique(branch_info))]
      xval <- apply(as.matrix(gdf_sub), 1, function(x)xvals[[x["branch"]]][as.numeric(x["order"])])
      gdf_sub$order <- xval

    }

    g = ggplot(gdf_sub, aes(x = order, y = gdf_sub$ExpressionGene, colour = gdf_sub$gene,
                            fill = gdf_sub$gene, linetype = gdf_sub$branch, shape = gdf_sub$branch)) +
      geom_point() +
      labs(fill = "Gene", col = "Gene", linetype = "Branch", shape = "Branch") +
      theme_minimal() + geom_smooth() +
      ylab("Expression") +
      ggtitle(paste0(genepair, collapse = ", ")) +
      NULL

    if (facet == "branch") {
      g = g + facet_grid(~branch)
    }

    if (facet == "gene") {
      g = g + facet_grid(gene~.)
    }

    if (facet == "both") {
      g = g + facet_grid(gene~branch)
    }
  }

  if (scHOT@positionType == "spatial") {

    if (is.null(positionColData)) {

      if (is.null(scHOT@positionColData)) {
        stop("Both positionColData and scHOT@positionColData are NULL.")
      } else {
        positionColData <- scHOT@positionColData
      }

    }

    coords_info <- data.frame(SummarizedExperiment::colData(scHOT)[, positionColData])
    colnames(coords_info) <- positionColData

    branch_long <- reshape::melt(cbind(coords_info, t(branchData[genepair, , drop = FALSE])),
                                 id.vars = positionColData)


    colnames(branch_long) <- c("x", "y", "genepair", "value")

    g <- ggplot(branch_long,  aes(x = branch_long$x, y = branch_long$y, color = branch_long$value)) +
      geom_point(size = 3) +
      # geom_point(size = 0.5, colour = "black") +
      theme_minimal() +
      facet_wrap(~genepair) +
      scale_alpha_continuous(range = c(0,0.5)) +
      scale_color_viridis_c(breaks = c(0,max(branch_long$value)),
                            limits = c(0,max(branch_long$value)),
                            labels = c("Low","High")) +
      theme(panel.grid = element_blank()) +
      theme(axis.ticks = element_blank()) +
      theme(axis.text = element_blank()) +
      xlab("") + ylab("") +
      coord_fixed() +
      labs(col = "Expression value") +
      NULL
  }




  return(g)
}


##############################################

#' the plotEgoNetwork function plots network graphs with edges coloured by weights in the network
#'
#' @title plotEgoNetwork
#' @param scHOT a scHOT object
#' @param hubnode is a character vector of node(s) to include as hub nodes
#' @param network is an igraph network
#' @param weight A string indicates the column name stored in scHOT_output slot that are used as the weights of the network
#' @param subset is a logical asking if you should subset based on the weight (default FALSE)
#' @param thresh is the subset weight threshold
#' @return \code{igraph} object containing the network graphed. Produces an igraph plot
#'
#'
#' @import igraph
#' @importFrom graphics plot
#'
#' @export


plotEgoNetwork = function(scHOT, hubnode, network,
                          weight = "higherOrderStatistic",
                          subset = FALSE,
                          thresh = NULL) {


  # hubnode is a character vector of node(s) to include as hub nodes
  # g is an igraph network, with E(g)[[weight]] given as DCARS test statistics
  # weight is a character vector containing edge weights, associated with different branches
  # subset is a logical asking if you should subset based on the weight
  # thresh is the subset weight threshold

  if (is.null(weight)) {
    igraph::E(network)$weight <- rep(1, ncol(scHOT))
  } else {
    if (!(weight %in% colnames(scHOT@scHOT_output))) {
      stop("weight provided is not found in scHOT_output")
    }

    igraph::E(network)$weight <- scHOT@scHOT_output[, weight]
  }

  if (!hubnode %in% names(V(network))) {
    stop("The provided hubnode is not a node in the network")
  }


  nodes = unique(names(unlist(igraph::neighborhood(network, nodes = hubnode))))
  subego = igraph::induced.subgraph(network, vids = nodes)
  subego = igraph::simplify(subego, edge.attr.comb = "mean")

  if (subset) {
    if (!all(weight %in% names(edge_attr(subego)))) {
      stop("at least one weight missing from edge attributes, either re-specify weights or rerun with subset = FALSE")
    }
    if (is.null(thresh)) {
      message("no threshold given, using 0 as default")
      thresh = 0
    }
    keepedges = apply(sapply(weight, function(w) igraph::edge_attr(subego)[[w]] > thresh,
                             simplify = TRUE),1,any)
    subego = igraph::subgraph.edges(subego, which(keepedges))
  }

  igraph::V(subego)$color = "beige"
  igraph::V(subego)$label.color = "black"
  igraph::V(subego)$label.family = "sans"
  igraph::V(subego)$label.cex = 0.7
  igraph::V(subego)$size = 20
  igraph::V(subego)$frame.color = "black"
  igraph::V(subego)$frame.size = 5
  lyout = igraph::layout.davidson.harel(subego)
  width = 5

  maxval = ceiling(max(50*unlist(edge_attr(subego)[["weight"]])))
  colvals = grDevices::colorRampPalette(c("grey","red"))(maxval)

  graphics::plot(subego, layout = lyout,
                 edge.color = colvals[ceiling(50*edge_attr(subego)[["weight"]])],
                 edge.width = width,
                 # xlab = "weight",
                 main = hubnode)

  # par(mfrow = c(1,length(weight)))
  #
  # for (i in weight) {
  #   plot(subego, layout = lyout,
  #        edge.color = colvals[ceiling(50*edge_attr(subego)[["weight"]])],
  #        edge.width = width,
  #        xlab = i,
  #        main = hubnode)
  # }
  return(subego)
}


##############################################


#' the plotHigherOrderSequence function plots weighted correlation vectors (stored in higherOrderSequence) as line plots
#'
#' @title plotHigherOrderSequence
#' @param scHOT A scHOT object with higherOrderSequence in scHOT_output slot
#' @param gene is either a logical vector matching rows of entries in wcorsList, or a character of a gene
#' @param positionType A string indicates the position type, either trajectory or spatial
#' @param branches A character indicates that the colnames stored the branch information in colData (for trajectory type of data)
#' @param positionColData A vector indicates column names of colData that stored the postion informaton (for spatial type of data)
#' @return \code{ggplot} object with line plots
#'
#'
#'
#'
#' @importFrom reshape melt
#' @import ggplot2
#' @importFrom ggforce geom_voronoi_tile
#'
#' @export
#'
plotHigherOrderSequence <- function(scHOT,
                                    gene,
                                    positionType = NULL,
                                    branches = NULL,
                                    positionColData = NULL) {

  # wcorsList is a list of matrices, with each matrix gene pair x samples weighted correlation vectors,
  # assumed that they have same number of rows
  # gene is either a logical vector matching rows of entries in wcorsList, or a character of a gene
  # matchExact matches gene names by splitting instead of using grep, but is slower

  if (!("higherOrderSequence" %in% colnames(scHOT@scHOT_output))) {
    stop("higherOrderSequence is not found in scHOT_output")
  }

  wcor <- as.matrix(scHOT@scHOT_output$higherOrderSequence)
  rownames(wcor) <- paste(scHOT@scHOT_output$gene_1, scHOT@scHOT_output$gene_2, sep = "_")
  colnames(wcor) <- NULL

  if (ncol(wcor) != ncol(scHOT)) {
    warning("Not all the cell position has higherOrderSequence statistics, set nrow.out = NULL in scHOT_setWeightMatrix to calculate higherOrderSequence for all positions!")
  }

  if (is.null(positionType)) {

    if (is.null(scHOT@positionType)) {
      stop("Both positionType and scHOT@positionType are NULL.")
    } else {
      positionType <- scHOT@positionType
    }

  }

  positionType <- match.arg(positionType, c("trajectory","spatial"), several.ok = FALSE)




  if (positionType == "trajectory") {

    if (is.null(branches)) {
      message("branches information is not provided")
      wcorsList = list(Branch = wcor)
    } else {
      if (!(branches %in% colnames(colData(scHOT)))) {
        stop("branches provided is not found in colData(scHOT)")
      }

      branch_info <- SummarizedExperiment::colData(scHOT)[, branches]

      wcorsList <- lapply(split(seq_along(branch_info), branch_info),
                          function(idx) wcor[, idx])[order(unique(branch_info))]
    }



    if (is.null(names(wcorsList))) {
      names(wcorsList) <- paste0("Branch_",1:length(wcorsList))
    }

    if (is.logical(gene[1])) {

      if (length(unique(unlist(lapply(wcorsList,nrow)))) > 1) {
        stop("cannot use logical subset when weighted correlation matrices have differing rows")
      }

      if (length(gene) != nrow(wcorsList[[1]])) {
        stop("cannot use logical subset when length of gene doesn't match nrow of wcorsList matrices")
      }

      wcors_longList = lapply(wcorsList,function(branch){
        reshape::melt(t(branch[gene,]))
      })

      gene = ""
    } else {

      gene = paste0(sort(gene), collapse = "|")

      wcors_longList = lapply(wcorsList,function(branch){
        reshape::melt(t(branch[grepl(gene, rownames(branch)), , drop = FALSE]))
      })

    }

    branch_long = do.call(rbind, wcors_longList)
    branch_long = cbind(
      rep(names(wcors_longList), unlist(lapply(wcors_longList, nrow))),
      branch_long)
    # colnames(branch_long) = c("branch","SampleOrder", "GenePair","WeightedCorrelation")
    colnames(branch_long) = c("branch","SampleOrder", "GenePair", "WeightedCorrelation")



    if (max(abs(branch_long$WeightedCorrelation)) < 1) {
      ylimit <- ylim(c(-1, 1))
    } else {
      ylimit <- NULL
    }

    g <- ggplot(branch_long,
                aes(x = branch_long$SampleOrder,
                    y = branch_long$WeightedCorrelation,
                    group = branch_long$GenePair,
                    col = branch_long$GenePair)) +
      geom_line(size = 2, alpha = 0.6) +
      facet_grid(~branch, scales = "free_x") +
      theme_minimal() +
      ylimit +
      geom_hline(yintercept = 0, size = 1, colour = "grey") +
      ggtitle(gene) +
      xlab("Sample Order") +
      ylab("Weighted Correlation") +
      labs(col = "Gene Pair") +
      NULL


  }


  # spatial case

  if (positionType == "spatial") {
    if (is.logical(gene[1])) {

      if (length(unique(unlist(lapply(wcorsList, nrow)))) > 1) {
        stop("cannot use logical subset when weighted correlation matrices have differing rows")
      }

      if (length(gene) != nrow(wcor)) {
        stop("cannot use logical subset when length of gene doesn't match nrow of wcorsList matrices")
      }


      gene = ""
    } else {


      if (!all(gene %in% rownames(wcor))) {

        if (length(gene) == 2) {
          if (!paste0(sort(gene), collapse = "_") %in% rownames(wcor)) {
            stop("gene pairs has no higherOrderSequence ")
          } else {
            gene <- paste0(sort(gene), collapse = "_")
          }
        } else {
          stop("gene pairs has no higherOrderSequence ")
        }
      }
    }

    if (is.null(positionColData)) {

      if (is.null(scHOT@positionColData)) {
        stop("Both positionColData and scHOT@positionColData are NULL.")
      } else {
        positionColData <- scHOT@positionColData
      }

    }

    coords_info <- data.frame(SummarizedExperiment::colData(scHOT)[, positionColData])
    colnames(coords_info) <- positionColData

    wcor_all <- matrix(NA, nrow = length(gene), ncol = ncol(scHOT))
    rownames(wcor_all) <- gene
    colnames(wcor_all) <- seq_len(ncol(scHOT))
    wcor_all[gene, rownames(scHOT@weightMatrix)] <- wcor[gene, , drop = FALSE]


    branch_long <- reshape::melt(cbind(coords_info, t(wcor_all)),
                                 id.vars = positionColData)

    colnames(branch_long) <- c("x", "y", "genepair", "value")

    g <- ggplot(branch_long,  aes(x = branch_long$x, y = branch_long$y, fill = branch_long$value)) +
      ggforce::geom_voronoi_tile(max.radius = 1) +
      geom_point(size = 0.5, colour = "black") +
      theme_minimal() +
      facet_wrap(~genepair) +
      scale_alpha_continuous(range = c(0,0.5)) +
      scale_fill_gradient2(low = "blue", mid = "white", high = "red", limits = c(-1,1)) +
      theme(panel.grid = element_blank()) +
      theme(axis.ticks = element_blank()) +
      theme(axis.text = element_blank()) +
      xlab("") + ylab("") +
      coord_fixed() +
      labs(fill = "Weighted Correlation") +
      NULL
  }

  return(g)
}



##############################################

#' the scHOT_plotPermutationDistributions function plots the permutation test statistics
#' as a diagnostic plot for estimating p-values
#'
#' @title scHOT_plotPermutationDistributions
#' @param scHOT a scHOT object
#' @return \code{ggplot} graph of global higher order function and the
#' permutation scHOT test statistics. This should have a continuous pattern
#' to be reliably used for p-value estimation
#'
#'
#' @import ggplot2
#' @importFrom ggplot2 ggplot
#'
#' @export

scHOT_plotPermutationDistributions = function(scHOT) {

  scHOT_output = scHOT@scHOT_output

  permstatsDF = data.frame(
    test = rep(seq_len(nrow(scHOT_output)),
               times = unlist(
                 lapply(scHOT_output$permutations, function(x) length(unlist(x))))),
    stat = unlist(scHOT_output$permutations),
    globalHigherOrderFunction = rep(scHOT_output$globalHigherOrderFunction,
                                    times = unlist(
                                      lapply(scHOT_output$permutations, function(x) length(unlist(x)))))
  )

  quantileDF = data.frame(
    test = rep(seq_len(nrow(scHOT_output)),
               times = lapply(scHOT_output$permutations, function(x) length(unlist(x)) > 0)),
    quantile_0.9 = unlist(lapply(scHOT_output$permutations, function(x) quantile(x, 0.9, na.rm = TRUE))),
    globalHigherOrderFunction = rep(scHOT_output$globalHigherOrderFunction,
                                    times = lapply(scHOT_output$permutations, function(x) length(unlist(x)) > 0))
  )

  quantileDF$quantile_0.9_fitted <- NA
  quantileDF$quantile_0.9_fitted[!is.na(quantileDF$quantile_0.9)] = loess(quantile_0.9 ~ globalHigherOrderFunction, data = quantileDF)$fitted

  gBase = ggplot(permstatsDF, aes(x = globalHigherOrderFunction, y = stat))

  if (require(scattermore) & require(scales)) {
    gBase = gBase + geom_scattermore()
  }

  g = gBase +
    geom_line(aes(y = quantile_0.9_fitted),
              data = reshape::sort_df(quantileDF, "globalHigherOrderFunction")) +
    theme_classic() +
    ylab("Permuted scHOT test statistics") +
    NULL

  return(g)
}
