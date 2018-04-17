#' Plot The Pathway-CEG Network for A Drug
#'
#' This function allows user to plot the pathway-CEG network for a specified pathway in a specified drug.
#' @importFrom igraph graph.adjacency plot.igraph
#' @param path.analysis.res The pathway analysis results. Output of function pathway.analysis().
#' @param pathway.FDR.cutoff The FDR threshold to select significant drug specific pathways and the default is 0.1.
#' @param drug.name The name of the drug user would like to analyze.The specified drug name should among those that are listed in the output table of function network.graph().
#' @param pathway.name The name of the pathway user would like to analyze. The specified pathway name should among those that are listed in the output table of function network.graph().
#' @param show.plot True or False, specifying if you want to show the pathway-CEG plot for a specific drug and pathway.
#' @keywords Drug pathway CEG network
#' @export
#' @examples
#'
#' ## get the Dr. Insight drug identification results
#' drug.ident.res = drug.ident(query.data = example.disease, cmap.ref.profiles = example.drug.profiles,
#'                  repurposing.unit = "treatment", connectivity = "negative")
#'
#' ## load in example pathway data
#' data("example.pathway")
#'
#' ## Performe pathway analysis (for the drugs that are identified by ident.drug())
#' path.analysis.res = pathway.analysis(drug.ident.res = drug.ident.res,
#'                     pathway.list = example.pathway,drug.FDR.cutoff = 0.5)
#'
#'path.CEG.network = CEG.network.graph(path.analysis.res = path.analysis.res,
#'                   pathway.FDR.cutoff = 0.5,drug.name = "drug1",
#'                   pathway.name = "pathway5", show.plot = TRUE)



##*******************************
CEG.network.graph = function(path.analysis.res = path.analysis.res,pathway.FDR.cutoff = 0.1,
                             drug.name = NULL, pathway.name = NULL, show.plot = TRUE){
  pvals = path.analysis.res$drug.pathway.pvalues
  drug.CEGs = path.analysis.res$drug.CEGs
  pathways = path.analysis.res$processed.pathway.list

  if(!(drug.name %in% colnames(pvals))){
    stop(simpleError("Drug name did not found in the pathway analysis result. Please specify a drug that is being analyzed in pathway analysis."))
  } else{
    if(!(pathway.name %in% rownames(pvals))){
      stop(simpleError("Pathway name did not found in the pathway analysis result. Please specify a pathway that exists in the pathway list."))
    } else{
      CEG = drug.CEGs$CEGs[[drug.name]]
      path.gene = pathways[[pathway.name]]
      up.common = intersect(names(CEG$zscore.up.CEGs),path.gene)
      down.common = intersect(names(CEG$zscore.down.CEGs), path.gene)
      N = 1+length(up.common) + length(down.common)

      dat = matrix(0, nrow = N, ncol = N)
      rownames(dat) = colnames(dat) = c(pathway.name,up.common,down.common)
      dat[1,2:ncol(dat)] = dat[2:ncol(dat),1] = 1
    }
  }

  label_path = pathway.name
  label_gene = c(up.common,down.common)


  g = graph.adjacency(dat, mode = 'undirected',diag = F)

  if(show.plot){
    plot.igraph(g,
                vertex.shape = c(rep("square",length(label_path)),rep("circle",length(label_gene))),
                vertex.size = c(20,rep(15,length(label_gene))),
                vertex.color = c("orange",rep("lightcoral",length(up.common)),rep("lightgreen",length(label_gene))),
                vertex.label = c(label_path,label_gene),
                vertex.label.color = "black",
                vertex.label.font = 2,
                vertex.label.cex= c(rep(1.2,length(label_path)),rep(1.2,length(label_gene))),
                # layout = layout_nicely(g, dim = 2),
                edge.width = 2,
                vertex.label.dist = 0.2
    )
    # legend('topleft',legend=c("pathway","up-regulated CEGs","down-regualted CEGs"),
    #        col=c("orange","lightcoral","lightgreen"),
    #        pch=c(22,19,19), pt.bg='orange',pt.cex = 1.5)
  }


  res = list()
  res[[1]] = list()
  names(res) = drug.name
  res[[1]][[1]] = list()
  names(res[[1]]) = pathway.name
  res[[1]][[1]]$up.CEGs = up.common
  res[[1]][[1]]$down.CEGs = down.common
  return(res)
}
