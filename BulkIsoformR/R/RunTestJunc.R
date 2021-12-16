



##
##
#' Test Junction usage changes
#'
#' Tests for changes between conditions of the use of junctions in each gene.
#'
#' @param dat_peak Peak data from regtools
#' @param meta Meta data matching peak data, rownames must match that col names of dat_peak
#' @param form Formula to fit model for
#' @param dat_gene A gene by sample matrix, if not given produced by collapsing dat_peak by gene
#' @param numReads The min number of reads before filtering
#' @param toTest The coefficient to test
#' @return A per peak p-value and effect size estimate
#' @export
TestPeakOutput<-function(dat_peak,meta,form=~condition,dat_gene=NULL,numReads=10,toTest="conditionKO")
{
collapse_genes=F
if(is.null(dat_gene))
{
dat_gene=dat_peak
collapse_genes=T
}

if(collapse_genes)
{
print("Collapse genes")
dat=data.frame(dat_gene)
print(head(dat))
genes=as.character(lapply(rownames(dat_peak),function(x){strsplit(x,":")[[1]][1]}))
dat["gene"]=genes
dat<-dat %>% gather(Sample,Count,-gene) %>% group_by(gene,Sample) %>% summarise(Count=sum(Count)) %>% spread(Sample,Count,fill=0) %>% as.data.frame()
rownames(dat)=dat[,"gene"]
dat=dat[,2:dim(dat)[2]]
dat_gene=as.matrix(dat)
}





print("Make Peak gene matrix")
genes=as.character(lapply(rownames(dat_peak),function(x){strsplit(x,":")[[1]][1]}))
dat_peak=dat_peak[genes %in% rownames(dat_gene),]
dat_gene=dat_gene[genes[genes %in% rownames(dat_gene)],]
rownames(dat_gene)=rownames(dat_peak)
dat_gene=dat_gene-dat_peak
print("Filter!")
minGene=apply(dat_gene,1,min)
num_peak=apply(dat_peak,1,sum)
dat_gene=dat_gene[minGene>5 & num_peak>numReads,]
dat_peak=dat_peak[minGene>5 & num_peak>numReads,]
print(dim(dat_gene))
print(dim(dat_peak))
print(head(dat_gene))
print(head(dat_peak))
#print("Get Percent Intronic Per Sample")
#ratio=apply(dat_int,2,sum)/apply(dat_ex,2,sum)
#dat_perc=matrix(rep(ratio,digenedat_ex)[1]),ncol=length(ratio),byrow=T)
print("Format")
dat_peak=as.matrix(dat_peak)
dat_gene=as.matrix(dat_gene)
dat_gene=log(dat_gene)
#dat_perc=log(dat_perc)
#dat_ex=dat_ex+dat_perc
print("Fit!")
print(sample)
#meta=meta[!duplicated(meta[,sample]),]
#print(meta)
#rownames(meta)=meta[,sample]
meta=meta[colnames(dat_peak),]
print(dim(dat_peak))
fit=glm_gp(dat_peak,form,meta,size_factors = 1,offset=dat_gene,overdispersion_shrinkage = T)
#return(fit)
print("Test!")
out=test_de(fit,toTest)
out["cond"]="HT"
out=out[order(out$pval),]
out["gene"]=as.character(lapply(out[,1],function(x){strsplit(x,":")[[1]][1]}))
return(out)
}
