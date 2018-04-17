#' Drug Mechanism of Action Pathway Analysis
#'
#' This function allows user to run pathway analysis on the identified significat drugs, therefore to detect drug mechanism related pathways, and the CEGs in the pathways.
#' @importFrom qusage read.gmt
#' @param drug.ident.res Dr. Insight drug identification analysis results. Output of function drug.ident().
#' @param pathway.list The pathways used to analyze drug mechanism. Should be a list of pathways where the names of lists are pathway names.
#' Dr.Insight provides NIH pathway interaction database (PID) pathways, which can be called by data("pathway.PID").
#' If you want to analyze with other pathways (from GSEA MsigDB), please use the parameter: pathway.list.path.
#' @param pathway.list.path The local path and file name of the pathways downloaded from GSEA MsigDB http://software.broadinstitute.org/gsea/msigdb/collections.jsp, if pathway.list is NULL.
#' Should be in full directory and file name (See instructions in vignette).
#' @param drug.FDR.cutoff The FDR threshold to select significant repurposed drugs by Dr. Insight for pathway analysis.
#' @param CEG.threshold The p value threhold for selecting significant consistently differential expressed genes (CEGs).
#' @keywords pathway analysis
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



pathway.analysis = function(drug.ident.res = drug.ident.res,
                            pathway.list = NULL,
                            pathway.list.path = NULL,
                            drug.FDR.cutoff = 0.1, CEG.threshold = 0.05){
  if(is.null(pathway.list) & is.null(pathway.list.path)){
    stop(simpleError("Please provide a pathway list for pathway analysis."))
  } else{
    if(is.null(pathway.list)){
      cat("\n")
      cat("\n")
      message("Using pathways provided in 'pathway.list.path' for pathway analysis.")
      pathway.list = qusage::read.gmt(pathway.list.path)
    } else if(!is.null(pathway.list) & (!is.null(pathway.list.path))){
      cat("\n")
      cat("\n")
      message("Warning: 'pathway.list' is provided. Ignore the 'pathway.list.path' given.")
    }
    options(warn=-1)
    cat("\n")
    drug.pvals = drug.ident.res$drug.pvals
    drug.info = drug.ident.res$drug.info
    ## pathway preprocess
    pathway_list_processed = lapply(pathway.list,function(x){x[x %in% rownames(drug.ident.res$CEG.pvals$up)]})
    ## remove empty pathways
    pathway_list_processed = pathway_list_processed[setdiff(seq(length(pathway_list_processed)),(which(sapply(pathway_list_processed,length) <= 1)))]

    ident_drugs = drug.pvals[drug.pvals$FDR < drug.FDR.cutoff,]$drug
    if(length(ident_drugs) <= 1){
      simpleError("Not enough drugs pass the drug FDR cutoff. Please select a less stringent cutoff threshold.")
    } else{
      ## aggregte results to drug level
      ident_drugs = unique(sapply(ident_drugs,function(x){strsplit(x,split = "_")[[1]][1]}))

      message("Generating CEGs for identified significant drugs ...\n")
      cat("\n")
      p_min = drug.ident.res$CEG.pvals$down
      p_max = drug.ident.res$CEG.pvals$up

      z_min = qnorm(p_min,lower.tail = F)
      z_max = qnorm(p_max,lower.tail = F)

      drug.CEGs = get_target_drug_CEG(z_min, z_max, ident_drugs, drug.info, CEG.threshold)

      message("Pathway significance tests ...\n")
      cat("\n")
      ## pathway tests
      path_hyper_p_min = getPathHyperP(drug.CEGs$down_CEG_zscore,pathway_list_processed,CEG.threshold)
      path_hyper_p_max = getPathHyperP(drug.CEGs$up_CEG_zscore,pathway_list_processed,CEG.threshold)

      path_ks_p_min = getPathwayP(z_min,pathway_list_processed, drug.info, ident_drugs)
      path_ks_p_max = getPathwayP(z_max,pathway_list_processed, drug.info, ident_drugs)

      ## FDR correction for both pathway tests
      path_hyper_p_min = apply(path_hyper_p_min,2,function(x){p.adjust(x,method = "fdr")})
      path_hyper_p_max = apply(path_hyper_p_max,2,function(x){p.adjust(x,method = "fdr")})

      path_ks_p_min = t(apply(path_ks_p_min,1,function(x){p.adjust(x,method = "fdr")}))
      path_ks_p_max = t(apply(path_ks_p_max,1,function(x){p.adjust(x,method = "fdr")}))


      ## select maximnum p value of two tests and cut with path.threshold
      path_p_down = pmax(path_hyper_p_min,path_ks_p_min)
      path_p_up = pmax(path_hyper_p_max,path_ks_p_max)

      ## some pathways are significant at both sides in one drug, select the more significant one to present
      path_p = path_p_signs =path_p_up
      idx = which(path_p_down < path_p_up)
      path_p[idx] = path_p_down[idx]
      path_p_signs[idx] = "-"
      path_p_signs[path_p_signs != "-"] = "+"

      path.analysis.res = list(path_p,path_p_signs, drug.CEGs, pathway_list_processed)
      names(path.analysis.res) = c("drug.pathway.pvalues","pathway.pvalue.signs","drug.CEGs","processed.pathway.list")

      return(path.analysis.res)
    }
  }
}

