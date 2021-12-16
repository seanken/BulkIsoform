
##
##lst is list of bams
##
#' Run DiffUTR
#'
#' Command to run diffUTR
#'
#' @param lst A list of bams
#' @param meta A dataframe of meta data
#' @param gtf The GTF to use
#' @param useDEX If true use DEX
#' @param usediffSpliceDGE If true use diffSpliceDGE
#' @return The results of diffUTR
#' @export
RunDiffUTR<-function(lst=c(),meta=c(),gtf="",useDEX=F,usediffSpliceDGE=F)
{
print("Make Bins")
bins=prepareBins(gtf,stranded=T)
print("Count Reads")
print(lst)
print("hi")
counts=countFeatures(lst,bins,2,isPairedEnd=TRUE)
print(str(counts))
colData(counts)[,"condition"]=meta[,"condition"]
print("Run DE Test!")
if(!useDEX)
{
if(usediffSpliceDGE)
{
counts=diffSpliceDGEWrapper(counts,~condition)
return(counts)
}
print("Run DiffSplice")
counts <- diffSpliceWrapper( counts, design = ~condition )
}
else{
print("Run DEXSeq")
counts <- DEXSeqWrapper(counts)
}
return(counts)
}


