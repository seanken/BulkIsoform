#' Load Feature Counts
#'
#' A function to load a list of count data from featureCount
#'
#' @param lst A name array of count files produced by featureCount
#' @return 
#' @export
loadFeature=function(lst)
{
dat=do.call(cbind,lapply(lst,function(x){tab=read.table(x,header=T,sep="\t",stringsAsFactors=F);col=dim(tab)[2];tab=tab[colnames(tab)[col]]}))
#col=7;dat=do.call(cbind,lapply(lst,function(x){tab=read.table(x,header=T,sep="\t",stringsAsFactors=F);tab=tab[colnames(tab)[col]]}))
colnames(dat)=names(lst)
dat=data.frame(dat)
x=lst[1]
rownames(dat)=read.table(x,header=T,sep="\t",stringsAsFactors=F)[,1]
return(dat)
}


##lst is a named list
#' Load PICARD QC data
#'
#' A function to load QC from PICARD
#'
#' @param lst A name array of files produced by RNASeqMetrics
#' @return 
#' @export
loadQC<-function(lst)
{
out=lapply(names(lst),function(x){fil=lst[x];dat=read.table(fil,stringsAsFactors=F,header=T,nrows=1,sep="\t");dat=data.frame(dat);dat["Name"]=x;return(dat)})
dat=do.call(rbind,out)
for(i in colnames(dat)){if(i!="Name"){dat[i]=as.numeric(dat[,i])}}
return(dat)
}

##lst is list of directories
#' Load Pipeline Output
#'
#' Loads the output of the pipeline into R
#'
#' @param lst A named list of directories produced by the pipeline (one directory per sample)
#' @param QCFilter True if you want to perform QC filtering, False otherwise
#' @return A list of data loaded, with entires All for feature count with introns+exons, Exons for counts with only exons, and QC for QC information
#' @export
getAll<-function(lst,QCFilter=T) #,introns=F)
{
lst1=sub("$","/Counts/counts.exons.txt",lst)
lst2=sub("$","/Counts/counts.introns.txt",lst)
lst3=sub("$","/QC/output.QC.txt",lst)
print("Get Exons")
dat_ex=loadFeature(lst1)
print("Get All Reads")
dat_all=loadFeature(lst2)
#meta["SampleType"]=as.character(lapply(lst,function(x){s=strsplit(x,"-")[[1]];paste(s[1:(length(s)-1)],collapse="_")}))
print("Get QC")
names(lst3)=colnames(dat_all)
QC=loadQC(lst3)
print("Put together!")
input=list()
input[["All"]]=dat_all
input[["Exons"]]=dat_ex
#input[["Meta"]]=meta
input[["QC"]]=QC
#if(introns)
#{
#input[["Introns"]]=loadData(lst)
#}
if(QCFilter)
{
print("Filter")
input=BasicQC(input)
}
return(input)
}

#' Removes cells
#'
#' Given the input produced by getAll filters out cells
#'
#' @param input A list with entries All, Exons, Meta, and QC as produced by getAll
#' @param to_keep A list of samples to keep
#' @param to_tem A list of camples to remove
#' @return The input list with bad samples removed
#' @export
cleanUp<-function(input,to_keep=c(),to_rem=c())
{
if(length(to_keep)<3)
{
to_keep=setdiff(colnames(input[["All"]]),to_rem)
}
input[["All"]]=input[["All"]][,to_keep]
input[["Exons"]]=input[["Exons"]][,to_keep]
#input[["Meta"]]=input[["Meta"]][to_keep,]
input[["QC"]]=input[["QC"]][input[["QC"]][,"Name"] %in% to_keep,]
return(input)
}


#' Basic QC filtering
#'
#' Does basic QC filtering on data
#'
#' @param input A list with entries All, Exons, Meta, and QC as produced by getAll
#' @param numSD Number of standard deviations to use for filtering
#' @return The input list after filtering based on QC
#' @export
BasicQC<-function(input,numSD=3)
{
QC=input[["QC"]]
QC=QC[,apply(QC,2,function(x){!is.na(x[1])})]
rem=do.call(c,apply(QC[,1:(dim(QC)[2]-1)],2,function(x){mn=median(x);sd_val=sd(x);if(sd_val==0){return(c())};QC[!(x<mn+numSD*sd_val & x>mn-numSD*sd_val),"Name"]}))
rem=unique(rem)
input=cleanUp(input,to_rem=rem)
return(input)
}

