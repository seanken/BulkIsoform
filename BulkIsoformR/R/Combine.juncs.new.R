##combine junctions with minDist of each other
#' Merge junctions
#'
#' Merges Junctions together is they are within short distance of one another
#'
#' @param dat Data frame of junction count information
#' @param minDist The distance to collapse
#' @return A new junction dataframe with collapsed junctions
#' @export
combinedJunctions=function(dat,minDist=50)
{
dat_orig=dat
#dat=dat[order(dat$End,decreasing=T),]
print(dim(dat))
dat=dat %>% group_by(Start,Sign,Chr) %>% summarise() %>% as.data.frame()
print(dim(dat))

dat=dat[order(dat$Start,decreasing=F),]
dat=dat[order(dat$Chr),]
dat=dat[order(dat$Sign),]
start=rep(0,dim(dat)[1])
i=0
chr=""
strand=""
curRange=c(0,0)
print("Cluster Junction Starts")
for(j in 1:dim(dat)[1]){
curchr=dat[j,"Chr"]
curstrand=dat[j,"Sign"]
if(curRange[2]>dat[j,"Start"] & curchr==chr & curstrand==strand)
{
curRange[2]=dat[j,"Start"]+minDist
}
else
{
chr=curchr
strand=curstrand
i=i+1
curRange=c(dat[j,"Start"]-minDist,dat[j,"Start"]+minDist)
}
start[j]=i
}
dat["Start_clust"]=start
print("Join")
dat_orig=inner_join(dat_orig,dat)

dat=dat_orig
print(dim(dat))
dat=dat %>% group_by(End,Sign,Chr) %>% summarise() %>% as.data.frame()
print(dim(dat))

dat=dat[order(dat$End,decreasing=F),]
dat=dat[order(dat$Chr),]
dat=dat[order(dat$Sign),]
end=rep(0,dim(dat)[1])
chr=""
strand=""
i=0
print("Cluster Junction End")
curRange=c(0,0)
for(j in 1:dim(dat)[1]){
curstrand=dat[j,"Sign"]
curchr=dat[j,"Chr"]
if(curRange[2]>dat[j,"End"] & curchr==chr & curstrand==strand)
{
curRange[2]=dat[j,"End"]+minDist
}
else
{
chr=curchr
strand=curstrand
i=i+1
curRange=c(dat[j,"End"]-minDist,dat[j,"End"]+minDist)
}
end[j]=i
}
dat["End_clust"]=end
print("Join")
dat_orig=inner_join(dat_orig,dat)
return(dat_orig)

}

##load Data for junction/genes
#' Load Junctions
#'
#' Loads the annotated junctions produced by the pipeline
#'
#' @param lst A list of annotated junction files
#' @return A data frame with junction count information
#' @export
loadData_Junc=function(lst)
{
print("Load data!")
out=lapply(lst,function(x){
dat=read.table(x,sep="\t",stringsAsFactors=F)
print(x)
print("Reformat and filter")
dat=dat[,c(1,2,3,4,5,6,16)]
colnames(dat)=c("Chr","Start","End","Name","Count","Sign","Gene")
tab=table(dat[,"Name"])
dat=dat[tab[dat[,"Name"]]==1,]
dat["Sample"]=x
return(dat)
})
dat=do.call(rbind,out)
print(head(dat))
dat=dat %>% group_by(Chr,Start,End,Sign,Gene,Sample) %>% summarise(Count=sum(Count)) %>% as.data.frame()
print("Merge Junctions")
dat=combinedJunctions(dat)
print("Finish!")
print(head(dat))
dat<-dat %>% unite(Name,Gene,Chr,Start_clust,End_clust,Sign,sep=":") 
tab<-dat %>% group_by(Name) %>% summarise(minStart=min(Start),maxStart=max(Start),minEnd=min(End),maxEnd=max(End)) %>% unite(FullName,Name,minStart,maxStart,minEnd,maxEnd,sep="_",remove=F)
tab=tab[,c("FullName","Name")]
dat<-inner_join(dat,tab)
dat<-dat %>% group_by(Sample,FullName) %>% summarise(Count=sum(Count)) %>% spread(Sample,Count,fill=0) %>% as.data.frame()
rownames(dat)=dat[,1]
dat=dat[,2:dim(dat)[2]]
return(dat)
}

