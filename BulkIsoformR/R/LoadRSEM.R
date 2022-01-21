#' Load RSEM output
#'
#' A function to load a list of data from RSEM in tximport format
#'
#' @param files A name array of isoform or gene level files from RSEM
#' @param gene If true load like genes, else load like isoforms
#' @return 
#' @export
loadRSEM<-function(files,gene=T)
{
res=c()
if(gene)
{
res=tximport(files, type = "rsem", txIn = FALSE, txOut = FALSE)
}
else{
res=tximport(files, type = "rsem", txIn = TRUE, txOut = TRUE)
res[["Gene2Tx"]]=read.table(files[1],header=T,stringsAsFactors=F,sep="\t")[,1:2]
}
return(res)
}

#' Make RSEM DRIM Ready
#'
#' Makes RSEM data ready for DTU testing
#'
#' @param counts A matrix of count data with rownames being transcript names
#' @param gene2trans A dataframe with one entry (the first) being transcript name, the second being the gene name
#' @return Ready from DRIM or other testing
#' @export
prepRSEMforDrim<-function(counts,gene2trans)
{
dat=data.frame(gene_id=gene2trans[,2],feature_id=gene2trans[,1],counts)
colnames(dat)=c("gene_id","feature_id",colnames(counts))
return(dat)
}

#' Load RSEM in testing format
#'
#' A function to load a list of data from RSEM as count matrix ready for DTU testing
#'
#' @param files A name array of isoform files from RSEM
#' @param nams An array of names, if not given just uses files as names. Used as sample names.
#' @return Ready from DRIM or other testing
#' @export
loadRSEMForDrim<-function(files,nams=c())
{
txi=loadRSEM(files,gene=F)
cts=txi[["counts"]]
if(length(nams)!=dim(cts)[2]){nams=files}
colnames(cts)=nams
dat=prepRSEMforDrim(cts,txi[["Gene2Tx"]])
if(sum(rownames(dat)!=dat$feature_id)!=0)
{print("Mismatch issue!")}
rownames(dat)=NULL
return(dat)
}

