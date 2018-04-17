## pathway significance test with both hypergeometric test and K-S test

## drug level CEG aggregation by instance median z score
get_target_drug_CEG = function(z_min, z_max, ident_drugs, drug.info, CEG.threshold){
  ## produce CEG z score matrix using median of muti-instance drug z scores
  z_min_drug = z_max_drug = matrix(NaN,nrow = nrow(z_min),ncol=length(ident_drugs))
  rownames(z_min_drug) = rownames(z_max_drug) = rownames(z_min)
  colnames(z_min_drug) = colnames(z_max_drug) = ident_drugs
  for(i in 1:ncol(z_min_drug)){
    ins = drug.info[drug.info$cmap_name == ident_drugs[i],]$instance_id
    z_min_drug[,i] = 'if'((length(ins) == 1),(z_min[,ins]),(apply(z_min[,ins],1,median)))
    z_max_drug[,i] = 'if'((length(ins) == 1),(z_max[,ins]),(apply(z_max[,ins],1,median)))
  }

  ## produce CEG list of target drugs
  CEG_min = apply(z_min_drug,2,function(x){rownames(z_min_drug)[x > qnorm(CEG.threshold,lower.tail = F)]})
  CEG_max = apply(z_max_drug,2,function(x){rownames(z_max_drug)[x > qnorm(CEG.threshold,lower.tail = F)]})
  CEGs = list()
  for(i in 1:length(CEG_min)){
    CEGs[[i]] = list()
    CEGs[[i]]$zscore.up.CEGs = z_max_drug[CEG_max[[i]],i]
    CEGs[[i]]$zscore.down.CEGs = z_min_drug[CEG_min[[i]],i]
  }
  names(CEGs) = names(CEG_min)

  res = list(z_max_drug,z_min_drug,CEGs)
  names(res) = c("up_CEG_zscore","down_CEG_zscore","CEGs")
  return(res)
}

## use hypergeometric test for pathway significance within each drug
getPathHyperP = function(z_order_drug, pathway_list_processed, CEG.threshold){
  ## calculate number of CEGs (z_score_drug != 0) in each drug
  treat_hit = apply(z_order_drug,2,function(x){sum(x >= qnorm(CEG.threshold,lower.tail = F))})

  path_hyper_p = matrix(NaN,nrow = length(pathway_list_processed),ncol=ncol(z_order_drug))
  for (i in 1:length(pathway_list_processed)){
    z_set = z_order_drug[rownames(z_order_drug) %in% pathway_list_processed[[i]],]
    ## calculate the number of CEG in each pathway
    set_hit = colSums(ifelse(z_set >= qnorm(CEG.threshold,lower.tail = F), 1,0))
    ## hyper test #hit-1 result samilar with fisher test #hit
    path_hyper_p[i,] = phyper(set_hit-1,nrow(z_set),nrow(z_order_drug)-nrow(z_set),treat_hit,lower.tail = F)
  }
  colnames(path_hyper_p) = colnames(z_order_drug)
  rownames(path_hyper_p) = names(pathway_list_processed)

  return(path_hyper_p)
}



## use ks test to test against null scores, across cmap drugs
getPathwayP = function(z_order, pathway_list_processed, drug.info, ident_drugs){
  pathway_score= matrix(0,nrow = length(pathway_list_processed),ncol=ncol(z_order))
  for (i in 1:length(pathway_list_processed)){
    z_set = z_order[rownames(z_order) %in% pathway_list_processed[[i]],]
    if(is.matrix(z_set)){ ## deal with pathways that only have one gene (can't be matrix)
      pathway_score[i,] = apply(z_set,2,function(x){sum(x)})
    } else {
      pathway_score[i,] = z_set
    }
  }
  colnames(pathway_score) = colnames(z_order)
  rownames(pathway_score) = names(pathway_list_processed)

  pathway_score_data = pathway_score[,drug.info[drug.info$cmap_name %in% ident_drugs,]$instance_id]
  pathway_score_rest = pathway_score[,drug.info[!(drug.info$cmap_name %in% ident_drugs),]$instance_id]

  pathway_p = matrix(0,nrow = nrow(pathway_score),ncol = length(ident_drugs))
  for (j in 1:length(ident_drugs)){
    ins = drug.info[drug.info$cmap_name == ident_drugs[j],]$instance_id
    ## use K-S test, one drug v.s. others (per pathway)
    for(i in 1:nrow(pathway_p)){
      pathway_p[i,j] = ks.test(pathway_score_data[i,ins],pathway_score_rest[i,],alternative = "less")$p.value
    }
  }
  colnames(pathway_p) = ident_drugs
  rownames(pathway_p) = rownames(pathway_score_data)
  return(pathway_p)
}

