##
##Same as TestIntrons_nb.R except has an offset term to deal with changes in intron levels between samples
##

##
##Test to see if differences in intron levl between conditions
##dat_all is count data with introns+exons
##dat_exon is count data with just exonic reads
##meta idata meta data
##form is the formula to use
##
#' Test for increase/decrease in intron levels
#'
#' Tests each gene to see if that gene has more/less introns than expected.
#'
#' @param dat_all A gene by sample matrix of counts, including both introns and exons
#' @param dat_ex A gene by sample matrix of counts, including onlyexons
#' @param meta The meta data to use
#' @param offset A boolea variable that determines if there should be an offset based on the percent intronic reads in each sample
#' @param form The formula to be tested (~condition by default)
#' @param test The coefficient to be tested
#' @return The results of testing each gene for changes in the ratio between the intronic and exonic reads
#' @export 
TestIntrons<-function(dat_all,dat_exon,meta,form=~condition,test="conditionKO",offset=T)
{
print("Process data")
dat_intron=dat_all-dat_exon
dat_intron=data.frame(dat_intron)
dat_exon=data.frame(dat_exon)
mn1=apply(dat_exon,1,sum)
mn2=apply(dat_intron,1,sum)

mn3=apply(dat_exon>0,1,sum)
mn4=apply(dat_intron>0,1,sum)

pos=apply(dat_intron,1,min)
dat_all=dat_all[mn1>10 & mn2>10 & pos>0 & mn3>1 & mn4>1,]
dat_intron=dat_intron[mn1>10 & mn2>10 & pos>0 & mn3>1 & mn4>1,]

dat1=dat_intron
dat2=log(dat_all)

dat1=as.matrix(dat1)
dat2=as.matrix(dat2)

if(offset)
{
tot_introns=colSums(dat_intron)
tot_reads=colSums(dat_all)
ratio=tot_introns/tot_reads
print(ratio)
dat3=matrix(rep(ratio,dim(dat1)[1]),ncol=length(ratio),byrow=T)
dat3=log(dat3)
dat2=dat2+dat3
rownames(dat2)=rownames(dat1)
}



fit=glm_gp(dat1,form,meta,size_factors = 1,offset=dat2,overdispersion_shrinkage = T)
#return(fit)
dat2=test_de(fit,test)
dat2["test"]=test
dat2=dat2[order(dat2$pval),]
print(head(dat2))
return(dat2)

}



