#!/usr/bin/env Rscript
Args <- commandArgs()
if ( length(Args) < 8  ){
fnm<-Args[4]
fnm<-substring(fnm,8,nchar(as.character(fnm)))
print(paste("Usage - RScript labels_in.csv eigenanat.csv outname.csv "))
q()
}
ARGIND<-6
labelsIn<-c(as.character(Args[ARGIND]))
ARGIND<-ARGIND+1
eiganatMatrixFN<-c(as.character(Args[ARGIND]))
ARGIND<-ARGIND+1
outname<-c(as.character(Args[ARGIND]))
ARGIND<-ARGIND+1
# to report regions that are included in each eigenanatomy component
labList<-read.csv(labelsIn) # the text names and numbers of each region 
eiganatMatrix<-read.csv(eiganatMatrixFN) # the text names and numbers of each region 
labels<-as.numeric(eiganatMatrix[1,]) # get list of labels 
neiganats<-nrow(eiganatMatrix)-1 # count eigenanatomy components 
# volumes <- apply(labels>0, 1, sum) # get volume of non-zero component in each row
# ulabs<-sort(unlist(unique(c(labels))),index.return=T)$x # get unique labels 
# get volume of each label in label map
ulabs<-labList[,1]
alllabelvols<-rep(0,length(ulabs))
for ( x in c(1:length(ulabs)) ) alllabelvols[x]<-sum(labels==ulabs[x])
# first fill the matrix of volumes : N-EigAnat vs M-anatomical regions
volmatrix<-matrix(NA,length(ulabs),neiganats)
for ( x in c(1:neiganats) ) 
  {
    evec<-( as.numeric(eiganatMatrix[x+1,]) > 0 )
    loclabs<-labels*(evec)
    for ( y in c(1:length(ulabs)) ) volmatrix[y,x]<-sum(loclabs==ulabs[y])
    if ( sum(volmatrix[,x]) > 0 )
      volmatrix[,x]<-volmatrix[,x]/sum(volmatrix[,x])
}

# now print out the results 
formattedoutput<-rep("",neiganats)
for ( x in c(1:neiganats) ) 
  {
    preout<-paste("eigenanatomy",x-1,sep='')
    anatlist<-""
    for ( y in c(1:length(ulabs)) )
      {
        if ( !is.na(volmatrix[y,x]) )
        if ( volmatrix[y,x] > 0 )
          anatlist<-c(anatlist,as.character(labList[y,2]))
      }
    anatlist<-toString(anatlist)
    anatlist<-gsub(",","",anatlist)
    anatlist<-paste(preout,anatlist,sep='')
    anatlist<-gsub(" ",",",anatlist)
    formattedoutput[x]<-as.character(anatlist)
  }

df<-data.frame(formattedoutput,t(volmatrix))
colnames(df)<-c("AllAnat",as.character(paste("Ratio_",labList[,2],sep='')))
write.csv(df,outname)
q()

adjmat<-matrix(NA,length(ulabs),length(ulabs))
 x<-5 
 for ( y in c(1:length(ulabs)) ) 
 for ( z in c(1:length(ulabs)) ) 
   if ( volmatrix[y,x] > 0 &&  volmatrix[z,x] > 0 && y != z ) adjmat[y,z]<-1  
 g1 <- graph.adjacency( adjmat,mode=c("undirected")) #  mode=c("directed", "undirected", "max","min", "upper", "lower", "plus"))
 plot(g1)
 coords <- layout.fruchterman.reingold(g1, dim=2)
 rglplot(g1, layout=coords)
 
# simulate some data 
 pvals<-rep(0.01,neiganats)
 anatnames<-gsub(" ","",paste("A",c(1:length(ulabs))))
 # now add pvals to that 
 pvalsAndVols<-t(cbind(pvals,t(volmatrix)))
 # now add totvols to that 
 pvalsAndVols<-t( cbind(volumes[2:length(volumes)],t(pvalsAndVols)))
 rownames<-rep("",nrow(pvalsAndVols))
 rownames[1]<-"TotalVolume"
 rownames[2]<-"Significance"
 rownames[3:length(rownames)]<-anatnames
# make a data frame to hold this information 
 colnames<-paste("EigenAnat",c(1:neiganats))
 dfout<-data.frame(pvalsAndVols) 
 colnames(prop)<-colnames
 rownames(prop)<-rownames


# edgeList <- matrix( c("foo", "bar", "bar", "bone","bone","foo","bone","foo"), nc=2, byrow=TRUE)
# ntwk <- graph.edgelist(edgeList,directed=FALSE)
# pdf('temp.pdf')
# plot(ntwk)
# dev.off()

# create a graph that shows the connections within an evec
# rglplot(g, layout=coords)
# g<-g13
# g <- graph.lattice( c(5,5,5) )
# plot(g)
# coords <- layout.fruchterman.reingold(g, dim=3)
# rglplot(g, layout=coords) 

# now format/print out the results 
# get anatomy names formatted
 anat<-toString(paste(anatnames,"&"))
 anat<-gsub(",","",anat)
 print(paste("EigenAnatomy & Volume & q-value & Anatomy "),quote=F)
 for ( x in c(2:neiganats) ) 
 {
   if ( pvals[x-1] < 0.05001 )
   {
   labelvols<-volmatrix[,x-1]
   if ( sum(labelvols[2:length(labelvols)]) > 0 )
   {
   labelvolsnorm<-labelvols/alllabelvols
   labelVolsStr<-toString(labelvols)
   labelVolsStr<-gsub(",","& ",labelVolsStr)
   labelVolsStr<-gsub(" 0& "," . & ",labelVolsStr)
   locanat<-c("")
   for ( y in 2:length(ulabs) )
   {
     if ( labelvols[y-1] > 0 ) locanat<-paste(toString(locanat)," + ",anatnames[y-1])
   }
#    print(paste(x-1,volumes[x],pvals[x-1],labelVolsStr),quote=F,zero.print = ".")
   print(paste(x-1,volumes[x],pvals[x-1],locanat),quote=F,zero.print = ".")
   }
   }
 }
# q()
