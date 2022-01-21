#' Run DRIM-Seq
#'
#' Runs DRIM-seq to see to see differential transcript useage (DTU).
#'
#' @param cts A dataframe of counts as produced by the load DRIM functions for RSEM and Salmon. First two columns should be feature_id and gene_id, the transcript and gene names.
#' @param meta A dataframe with meta data. Needs a column sample_id whose names match the names of samples in cts.
#' @param form The formula to use for testing, using columns in meta.
#' @param coef The coefficient to return
#' @param min_samps_feature_expr The number of samples that a given trascript must be expressed in.
#' @return A dmTest object from DRIM-seq
#' @export
RunDRIM<-function(cts,meta,form,coef=2,min_samps_feature_expr=3)
{
print("Set up")
cts=cts[rowSums(cts[,3:dim(cts)[2]])>0,]
d <- dmDSdata(counts=cts, samples=meta)
print("Filter")
d<- dmFilter(d,min_samps_feature_expr=min_samps_feature_expr, min_feature_expr=10,min_gene_expr = 10,min_samps_gene_expr=(dim(cts)[2]-2))
design_full <- model.matrix(form, data=DRIMSeq::samples(d))
print(colnames(design_full))
print("Calc precision!")
d <- dmPrecision(d, design=design_full)
print("Fit!")
d <- dmFit(d, design=design_full)
print("Test!")
d <- dmTest(d, coef=coef)
print("Return Results!")
return(d)
}
