#' Load Salmon output
#'
#' A function to load a list of data from Salmon in tximport format
#'
#' @param files A name array of sf files from Salmon
#' @return 
#' @export
loadSalmon<-function(files)
{
res=tximport(files,type="salmon", txOut=TRUE,countsFromAbundance="scaledTPM")
return(res)
}

#' Make Salmon DRIM Ready
#'
#' Makes Salmon data ready for DTU testing
#'
#' @param counts A matrix of count data with rownames being transcript names
#' @param gene2trans A dataframe with one entry (the first) being transcript name, the second being the gene name
#' @return Ready from DRIM or other testing
#' @export
prepSalmonforDrim<-function(counts,gene2trans)
{
rownames(gene2trans)=gene2trans[,"feature_id"]
gene2trans=gene2trans[rownames(counts),]
dat=data.frame(gene_id=gene2trans[,2],feature_id=gene2trans[,1],counts)
colnames(dat)=c("gene_id","feature_id",colnames(counts))
return(dat)
}

#' Load Salmon in testing format
#'
#' A function to load a list of data from Salmon as count matrix ready for DTU testing
#'
#' @param files A name array of sf files from Salmont
#' @param gene2trans A dataframe with one entry (the first) being transcript name, the second being the gene name. If NULL will try to load assuming data from our pipeline.
#' @param nams An array of names, if not given just uses files as names. Used as sample names.
#' @return Ready from DRIM or other testing
#' @export
loadSalmonForDrim<-function(files,gene2trans=NULL,nams=c())
{
if(is.null(gene2trans))
{
fil=sub("quant.sf$","trans2gene.txt",lst[1])
gene2trans=read.table(fil,sep="\t",stringsAsFactors=F)
colnames(gene2trans)=c("feature_id","gene_id")
}

txi=loadSalmon(files)
cts=txi[["counts"]]
if(length(nams)!=dim(cts)[2]){nams=files}
colnames(cts)=nams
print(dim(cts))
print(dim(gene2trans))
dat=prepSalmonforDrim(cts,gene2trans)
if(sum(rownames(dat)!=dat$feature_id)!=0)
{print("Mismatch issue!")}
rownames(dat)=NULL
return(dat)
}


