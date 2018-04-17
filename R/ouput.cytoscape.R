#' Drug and Pathway connection output files for Cytoscape visulization
#'
#' This fucniton allows user to get the two files needed for Cytoscape to visulize the drug-pathway network
#' @param path.analysis.res The pathway analysis results. Output of pathway.analysis().
#' @param pathway.FDR.cutoff The FDR threshold to select significant drug specific pathways and the default is 0.1.
#' @keywords pathway analysis output
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
#'                     pathway.list = example.pathway, drug.FDR.cutoff = 0.5)
#'
#' ## get the pathway analysis ouput that can be loaded into Cytoscape for visulization
#' network.cytoscape = make.cytoscape.network(path.analysis.res = path.analysis.res,
#' pathway.FDR.cutoff = 0.5)

##*******************************
## make data frame of Dr.Inight identified drug-pathway interactions (both min and max sides, drug level)
make.cytoscape.network = function(path.analysis.res = path.analysis.res,
                                  pathway.FDR.cutoff = 0.1){

  pvals = path.analysis.res$drug.pathway.pvalues
  p.signs = path.analysis.res$pathway.pvalue.signs
  idx1 = which(pvals <= pathway.FDR.cutoff)
  idx2 = which(pvals > pathway.FDR.cutoff)
  pvals[idx1] = 1
  pvals[idx2] = 0
  pvals = pvals[apply(pvals,1,function(x){sum(x == 1) >= 1}),]
  p.signs = p.signs[rownames(pvals),]

  ## make drug pathway connection file
  tmp = colnames(pvals)
  Dr_pathway = data.frame(drug = NA,pathway = NA, directions = NA)
  for(i in 1:length(tmp)){
    tmp2 = which(pvals[,i] == 1)
    tmp3 = p.signs[tmp2,i]
    Dr_pathway_tmp = data.frame(drug = rep(tmp[i],length(tmp2)), pathway = rownames(pvals)[tmp2],
                                directions = tmp3)
    Dr_pathway = rbind(Dr_pathway,Dr_pathway_tmp)
  }
  Dr_pathway = Dr_pathway[-1,]
  Dr_pathway$directions[Dr_pathway$directions == "+"] = "up"
  Dr_pathway$directions[Dr_pathway$directions == "-"] = "down"


  ## make node attributes file
  tmp2 = data.frame(unique(Dr_pathway$pathway))
  colnames(tmp2) = "node"
  tmp2$is.drug = FALSE
  tmp2$is.pathway = TRUE
  tmp2$degree = sapply(tmp2$node,function(x){length(Dr_pathway[Dr_pathway$pathway == x,]$pathway)})

  tmp = data.frame(unique(Dr_pathway$drug))
  colnames(tmp) = "node"
  tmp$is.drug = TRUE
  tmp$is.pathway = FALSE
  tmp$degree = NA

  tmp = rbind(tmp,tmp2)

  network.cytoscape = list()
  network.cytoscape$drug.pathway.network = Dr_pathway
  network.cytoscape$node.attributes = tmp

  # dir.create(paste0(file.path,"/cytoscape/"))
  # write.csv(network.cytoscape$drug.pathway.network,
  #           file = paste0(file.path,"/cytoscape/drug.pathway.network.csv"),
  #           row.names = F,
  #           quote = F)
  # write.csv(network.cytoscape$node.attributes,
  #           file = paste0(file.path,"/cytoscape/drug.pathway.network.node.attributes.csv"),row.names = F,
  #           quote = F)
  return(network.cytoscape)
}

