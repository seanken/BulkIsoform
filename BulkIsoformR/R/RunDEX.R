#' Run DEX-Seq on Isoforms
#'
#' A method to run DEX-Seq for ISOform level DTU testing. Implements the approach (including copying some code) from Love et al 2018 (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6178912/)
#'
#' @param cts A dataframe of counts as produced by the load DRIM functions for RSEM and Salmon. First two columns should be feature_id and gene_id, the transcript and gene names.
#' @param meta A dataframe with meta data.
#' @param form The formula to use for testing, using columns in meta.
#' @param form_red the reduced form that is tested against
#' @param min_samps_feature_expr The number of samples that a given trascript must be expressed in.
#' @return A list with 2 entries: the first the DEX results at Isoform leve, the second the results at Gene level
#' @export
RunDexIso<-function(cts,meta,form=~sample + exon + condition:exon,form_red=~sample + exon,min_samps_feature_expr=3)
{
print("Set up")
cts=cts[rowSums(cts[,3:dim(cts)[2]])>0,]
d <- dmDSdata(counts=cts, samples=meta)
d<- dmFilter(d,min_samps_feature_expr=min_samps_feature_expr, min_feature_expr=10,min_gene_expr = 10,min_samps_gene_expr=(dim(cts)[2]-2))
print("Make DEXSeq object")
sample.data <- DRIMSeq::samples(d)
count.data <- round(as.matrix(counts(d)[,-c(1:2)]))
dxd <- DEXSeqDataSet(countData=count.data,sampleData=sample.data,design=form,featureID=counts(d)$feature_id,groupID=counts(d)$gene_id)
print("Fit Model")
dxd <- estimateSizeFactors(dxd)
dxd <- estimateDispersions(dxd, quiet=TRUE)
dxd <- testForDEU(dxd, reducedModel=form_red)
dxr <- DEXSeqResults(dxd, independentFiltering=FALSE)
print("get gene level and return")
qval <- perGeneQValue(dxr)
lst=list()
lst[["Results"]]=dxr
lst[["PerGene"]]=qval
return(lst)
}
