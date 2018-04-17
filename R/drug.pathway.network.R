#' Plot Drug and Pathway Interaction Network
#'
#' This function allows user to plot the drug-pathway interaction network and returns the drug-pathway interaction contigency table.
#' @importFrom igraph graph.adjacency plot.igraph
#' @param path.analysis.res The pathway analysis results. Output of function pathway.analysis().
#' @param pathway.FDR.cutoff The FDR threshold to select significant drug specific pathways and the default is 0.1.
#' @param pathway.label True or False specofying if show pathway labels in the graph.
#' @param drug.label.size The number indicating the size of drug labels in the graph.
#' @param path.label.size The number indicating the size of pathway labels in the graph.
#' @param return.adj.table True or False specifying if return the resulted drug-pathway contigency table.
#' In the returned table, 1's represent that the pathway is up-regulated by the drug; -1's reprensent down-regulation and 0's no-regulation.
#' @keywords drug pathway network
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
#' drug.pathway.network = network.graph(path.analysis.res, pathway.FDR.cutoff = 0.5,
#'                        return.adj.table = TRUE, pathway.label = TRUE)


##*******************************
network.graph = function(path.analysis.res = path.analysis.res,pathway.FDR.cutoff = 0.1,
                             pathway.label = FALSE, drug.label.size = 1.2, path.label.size = 0.8,
                             return.adj.table = TRUE){
  pvals = pval.table = path.analysis.res$drug.pathway.pvalues
  idx1 = which(pvals <= pathway.FDR.cutoff)
  idx2 = which(pvals > pathway.FDR.cutoff)
  pvals[idx1] = 1
  pvals[idx2] = 0

  pvals = pvals[apply(pvals,1,function(x){sum(x == 1) >= 1}),]

  N = nrow(pvals) + ncol(pvals)
  drug.pathway = matrix(NaN, nrow = N, ncol = N)
  rownames(drug.pathway) = colnames(drug.pathway) = c(rownames(pvals),colnames(pvals))

  drug.pathway[1:nrow(pvals),1:nrow(pvals)] = 0
  drug.pathway[1:nrow(pvals),(nrow(pvals)+1):nrow(drug.pathway)] = pvals
  drug.pathway[(nrow(pvals)+1):ncol(drug.pathway),1:nrow(pvals)] = t(pvals)
  drug.pathway[(nrow(pvals)+1):ncol(drug.pathway),(nrow(pvals)+1):nrow(drug.pathway)] = 0


  label_drug = colnames(pvals)
  label_path = rownames(pvals)


  g = graph.adjacency(drug.pathway, mode = 'undirected',diag = F)
  if(pathway.label == T){
    plot.igraph(g,
                vertex.shape = c(rep("square",length(label_path)),rep("circle",length(label_drug))),
                vertex.size = c(rep(8,length(label_path)),rep(12,length(label_drug))),
                vertex.color = c(rep("orange",length(label_path)),rep("lightblue2",length(label_drug))),
                vertex.label = c(label_path,label_drug),
                vertex.label.color = "black",
                vertex.label.font = 2,
                vertex.label.cex= c(rep(path.label.size,length(label_path)),rep(drug.label.size,length(label_drug))),
                # layout = layout_nicely(g, dim = 2),
                edge.width = 2,
                vertex.label.dist = 0.2
    )
  } else if(pathway.label == F){
    plot.igraph(g,
                vertex.shape = c(rep("square",length(label_path)),rep("circle",length(label_drug))),
                vertex.size = c(rep(8,length(label_path)),rep(12,length(label_drug))),
                vertex.color = c(rep("orange",length(label_path)),rep("lightblue2",length(label_drug))),
                vertex.label = c(rep(NA,length(label_path)),label_drug),
                vertex.label.color = "black",
                vertex.label.font = 2,
                vertex.label.cex= c(rep(path.label.size,length(label_path)),rep(drug.label.size,length(label_drug))),
                # layout = layout_nicely(g, dim = 2),
                edge.width = 2,
                vertex.label.dist = 0.2
    )
  }
  # legend('topleft',legend=c("drug","pathway"),
  #        col=c("lightblue","orange"),
  #        pch=c(19,22), pt.bg='orange',pt.cex = 1.5)
  if(return.adj.table){
    p.signs = path.analysis.res$pathway.pvalue.signs
    p.signs = p.signs[rownames(pvals),]
    pval.table = pvals
    pval.table[p.signs == "-"] = -(pval.table[p.signs == "-"])
    return(pval.table)
  }
}
