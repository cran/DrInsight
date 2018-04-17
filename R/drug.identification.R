#' Dr. Insight Drug Identification
#'
#' This function allows user to use the differential expression data of their interested disease phenotype or drug profile to query against CMap reference drug expression profiles.
#' This function returns an obejct of three elements: a table of drug p-values reflecting the connectivity between query data and CMap drugs; a table contains all the meta information of CMap drugs; and an object of p-values of identified CEGs for each drug instance.
#' @importFrom stats ks.test median pbeta phyper qnorm p.adjust
#' @importFrom utils read.table write.csv
#' @param query.data User provided differential expression analysis results(e.g. t-test statistic scores) of querying data (either disease phenotype data, or drug expression data). The column names for gene symboles and statistic scores must be "geneSymbol" and "score".
#' @param cmap.ref.profiles Dr. Insight provided CMap drug rank matrix containing all 6100 instances of CMap data set.The example of reference data can be loaded with data("example.profiles"). The real cmap data can be loaded with the function 'get.cmap.ref' (See instructions in vignette).
#' @param repurposing.unit The parameter of either "treatment" or "drug", which indicates if user want the algorithm to test drug repurposing p value at treatment level or drug level. The default is "treatment", which treats the drug data from different cell lines separately.
#' @param CEG.threshold The p value threshold to select the consistently differential expressed genes (CEGs).
#' @param connectivity The type of connectivity, either "negative" or "positive". Negative connectivity is used when the query data is the differential scores from disease data, and Dr. Insight will repurpose drugs that can potentially reverse the query disease phenotype. Positive connectivity is used when the query data is from a drug profile, and Dr. Insight will return the drugs that are similar to the query drug.
#' @keywords Drug Indetification
#' @export
#' @examples
#'
#' ## load in the example query disease data
#' data("example.disease")
#' data("example.drug.profiles")
#'
#' ## get the Dr. Insight drug identification results
#' drug.ident.res = drug.ident(query.data = example.disease, cmap.ref.profiles = example.drug.profiles,
#'                  repurposing.unit = "treatment", connectivity = "negative")
#'
#' drug.pvals = drug.ident.res$drug.pvals


drug.ident = function(query.data = NULL,cmap.ref.profiles = NULL,
                      repurposing.unit = "treatment",CEG.threshold = 0.05, connectivity = "negative"){
  cmap.drug.rank = cmap.ref.profiles$drug.rank.matrix
  e1 = simpleError("Did not find the column named 'geneSymbol' in query data that contains the gene symbols in it.")
  e2 = simpleError("Did not find the column named 'score' in query data that contains the test statistics or any values that you would like to rank the genes.")

  cat("\n")
  cat("\n")
  message("Data preprocessing ...\n")
  cat("\n")
  if("score" %in% colnames(query.data)){
    if("geneSymbol" %in% colnames(query.data)){
      tmp = data_preprocess(query.data, cmap.drug.rank,connectivity = connectivity)
      query.data = tmp[[1]]
      cmap.drug.rank = tmp[[2]]
      rm(tmp)
    } else{
      stop(e1)
    }
  } else{
    stop(e2)
  }

  message("Identifying drug instance CEGs...\n")
  cat("\n")
  p_min = get_gene_pval('min',cmap.drug.rank,query.data)
  p_max = get_gene_pval('max',cmap.drug.rank,query.data)

  ## select the smallest p value (between 2 p values) as the p value of the gene
  p_score = pmin(p_min,p_max)
  z_score = qnorm(p_score,lower.tail = F)
  CEG.pvals = get_CEGs(p_min, p_max, z_score,threshold = CEG.threshold)

  message("Calculating drug connectivity p values ...\n")
  cat("\n")
  drug.info = cmap.ref.profiles$drug.info
  if(repurposing.unit == "drug"){
    drug.info$drug = drug.info$cmap_name
  } else if(repurposing.unit == "treatment"){
    drug.info$drug = drug.info$treatment
  } else{
    stop(simpleError("Please set the repurposing unit to either 'drug' or 'treatment'."))
  }

  drug.pvals = get_drug_pval(CEGsum = CEG.pvals$CEG.sumz.scores,drug.info = drug.info)

  drug.pvals = drug.pvals[order(drug.pvals$pval),]
  drug.pvals$drug = rownames(drug.pvals)
  rownames(drug.pvals) = NULL
  drug.pvals = drug.pvals[,c(2,1)]
  drug.pvals$FDR = p.adjust(drug.pvals$pval,method = "fdr")

  res = list(drug.pvals,drug.info,CEG.pvals$CEG.pvals)
  names(res) = c("drug.pvals","drug.info","CEG.pvals")
  return(res)
}

