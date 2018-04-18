## ----eval=F--------------------------------------------------------------
#  library("DrInsight")
#  # load in and process the downloaded CMap drug rank matrix
#  cmap.ref.profiles = get.cmap.ref(cmap.data.path = '/path_to/rankMatrix.txt',
#                                   probe.to.genes = probe.to.genes, drug.info = drug.info)

## ----eval=F--------------------------------------------------------------
#  # select the drug p value table
#  drug.pvals = drug.ident.res$drug.pvals

