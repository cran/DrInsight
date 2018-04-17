#' Read In and Process CMap Reference Profiles
#'
#' This function allows user to load in the CMap drug rank matrix user downloaded from CMap website: https://portals.broadinstitute.org/cmap/ (data matrix in the "downloads" section).
#' @param cmap.data.path The local path and the name of the downloaded data matrix file. The data matrix file should be in txt format.
#' @param probe.to.genes The ID converter between Affymetrix probe IDs (the IDs used in CMap data matrix) and official gene symbol.
#' The packge comes with an embeded probe.to.genes file that can be directly used. User can use their own converter file with two columns named "ID", and "Gene.Symbol'.
#' @param drug.info The drug instance information from CMap data. The deafault drug.info file comes with DrInsight package.
#' @keywords CMap Reference Profiles
#' @export


get.cmap.ref = function(cmap.data.path = NULL, probe.to.genes = probe.to.genes, drug.info = drug.info){
  cat("\n")
  cat("\n")
  message("Loading CMap drug matrix. This may take some time ... \n")
  cmap.drug.rank = read.table(cmap.data.path,row.names = 1, header = T)

  cmap.drug.rank = cmap.drug.rank[probe.to.genes$ID,]

  rownames(cmap.drug.rank) = probe.to.genes$Gene.Symbol

  cmap.ref.profiles = list(drug.info = drug.info, drug.rank.matrix = cmap.drug.rank)
  return(cmap.ref.profiles)
}
