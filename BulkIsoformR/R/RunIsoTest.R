#' Run alternate DTU test
#'
#' Runs an alternate test for DTU
#'
#' @param cts A dataframe of counts as produced by the load DRIM functions for RSEM and Salmon. First two columns should be feature_id and gene_id, the transcript and gene names.
#' @param meta A dataframe with meta data used for testing. Rownames should be the sample names.
#' @param form The formula to use for testing, using columns in meta.
#' @param contrasts List of contrasts to test. If not supplied will test all coeficients except the intercept.
#' @param numSamp The number of samples a transcript needs to be expressed in to be tested
#' @return Returns a list with the DTU results for each contrast tested.
#' @export 
RunISO<-function(cts,meta,form,contrasts=c(),numSamp=4)
{
print("Set up")
cts=cts[rowSums(cts[,3:dim(cts)[2]]>0)>numSamp,]
rownames(cts)=cts[,"feature_id"]
tx2gene=cts[,1:2]
tab<-cts %>% gather(Sample,Count,-gene_id,-feature_id) %>% group_by(gene_id,Sample) %>% summarise(Tot=sum(Count)) %>% as.data.frame()
tab<-inner_join(tx2gene,tab)%>% spread(Sample,Tot,fill=0)
rownames(tab)=tab[,"feature_id"]
cts=cts[,3:dim(cts)[2]]
tab=tab[rownames(cts),colnames(cts)]
tab2=tab-cts
readsPerGene=rowMeans(tab>10)
readsOut=rowMeans(tab2>0)
tab2=tab2[readsPerGene==1 & readsOut==1,]
cts=cts[readsPerGene==1 & readsOut==1,]
tab2=log(tab2)
cts=as.matrix(floor(cts))
tab2=as.matrix(tab2)
print("Run GlmGamPoi")
meta=meta[colnames(cts),]
fit=glm_gp(cts,design=form,col_data=meta,offset=tab2,size_factors=1)
print(colnames(fit$Beta))
print("Get DE")
if(length(contrasts)==0)
{
contrasts=colnames(fit$Beta)
contrasts=contrasts[contrasts!="Intercept"]
}
res=lapply(contrasts,function(x){print(x);DE=test_de(fit,contrast=x);colnames(DE)[1]="feature_id";DE=inner_join(tx2gene,DE);DE=DE[order(DE$pval),];DE["Contrast"]=x;return(DE)})
names(res)=contrasts
print("Return results!")
return(res)
}

#' Gene level stats
#'
#' Takes results of RunISO and aggregates to gene level using FWER correction.
#'
#' @param res The list of results from RunISO. Alternately, a data frame with just one such result.
#' @return A list with one entry for each entry in res, with gene level results
#' @export
GeneLevel<-function(res)
{
if(class(res)=="data.frame")
{
lst=list()
lst[["Test"]]=res
res=lst
}

res=lapply(res,function(x){tab=x %>% group_by(gene_id,Contrast) %>% summarise(Num=length(feature_id),pval=min(1,pval*Num)) %>% as.data.frame();tab["padj"]=p.adjust(tab[,"pval"],"fdr");tab=tab[order(tab$pval),]})
return(res)

}
