#!/usr/bin/env Rscript
Args <- commandArgs()
if ( length(Args) < 8  ){
fnm<-Args[4]
fnm<-substring(fnm,8,nchar(as.character(fnm)))
print(paste("Usage - RScript eigenanat_descriptor_file.csv eigenanat.csv outprefix "))
q()
}
ARGIND<-6
eigenanatDesriptorFN<-c(as.character(Args[ARGIND]))
ARGIND<-ARGIND+1
eigenanatProjFN<-c(as.character(Args[ARGIND]))
ARGIND<-ARGIND+1
outname<-c(as.character(Args[ARGIND]))
ARGIND<-ARGIND+1
 library('MASS')
 library(methods) 
 library('Cairo')
 library('igraph')
 library('misc3d')
 library('rgl')
# to report regions that are included in each eigenanatomy component
eigenanatDesriptor<-read.csv(eigenanatDesriptorFN) # the text names and numbers of each region
a<-eigenanatDesriptor
eigenanatProj<-read.csv(eigenanatProjFN) # the text names and numbers of each region 
neig<-ncol(eigenanatProj)
for ( y in c(1:neig) )
  {
  vy<-eigenanatProj[,y]
  forheatmap<-matrix(0,neig,neig)
  ewt<-0
  for ( z in c(1:neig) )
  {
    vz<-eigenanatProj[,z]
    if ( sd(vz) > 0 & sd(vy) > 0 )
    {
    corval<-cor.test(vz,vy)$est
    if ( y != z )
      {
        forheatmap[y,z]<-corval
        ewt<-c(ewt,corval)
      }
    }
  }
  ewt<-ewt[2:length(ewt)]
  thresh<-0.25 # (max(forheatmap))*0.9
  ewt<-ewt[ ewt > thresh ] 
#  print('edge weight ') ; print(ewt)
  adjmat<-as.matrix( forheatmap > thresh , nrow=neig,ncol=neig)
  g1<-graph.adjacency( adjmat,mode=c("undirected")) #  mode=c("directed", "undirected", "max","min", "upper", "lower", "plus"))
  g1 <- set.edge.attribute(g1 , "weights", value=ewt*10)
#heatmap(forheatmap,Rowv=NA,Colv=NA,scale="none") # ,labRow=a$ROIName,labCol=a$ROIName)
  if ( y == 1 ){
    coords <- layout.fruchterman.reingold(g1, dim=3)
   # from ImageMath ROIStatistics 
   coords[,3]<-as.numeric(a$comX)
   coords[,1]<-as.numeric(a$comY)
   coords[,2]<-as.numeric(a$comZ)
  }
# coords<-as.matrix(read.csv('NirepManualGraphCoordsAxial.csv'))[,2:3]
# plot(g1,layout=coords) # instead use Cairo
# a$ROINumber
 colors<-rep("red",neig)
 colors[1:neig/2]="green"
  outfn<-paste(outname,y,"cnx.png",sep='')
  print(paste("write",outfn))
  png(outfn)
  #
  plot(g1,vertex.color=colors,vertex.label=a$ROINumber,layout=coords,edge.width=E(g1)$weights)
  dev.off()
}
warnings()
q()


                                        # plot(g1, layout=coords,
#       vertex.label=names(a), vertex.size=13, edge.color="black",
#       vertex.color=colors)
# coords <- layout.fruchterman.reingold(g1, dim=2)
# rglplot(g1, layout=coords)
 
# simulate some data 
# pvals<-rep(0.01,neiganats)
# anatnames<-gsub(" ","",paste("A",c(1:length(ulabs))))
 # now add pvals to that 
# pvalsAndVols<-t(cbind(pvals,t(volmatrix)))
 # now add totvols to that 
 # pvalsAndVols<-t( cbind(volumes[2:length(volumes)],t(pvalsAndVols)))
 # rownames<-rep("",nrow(pvalsAndVols))
 # rownames[1]<-"TotalVolume"
 # rownames[2]<-"Significance"
 # rownames[3:length(rownames)]<-anatnames
# make a data frame to hold this information 
# colnames<-paste("EigenAnat",c(1:neiganats))
# dfout<-data.frame(pvalsAndVols) 
# colnames(prop)<-colnames
# rownames(prop)<-rownames


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
# anat<-toString(paste(anatnames,"&"))
# anat<-gsub(",","",anat)
# print(paste("EigenAnatomy & Volume & q-value & Anatomy "),quote=F)
# for ( x in c(2:neiganats) ) 
# {
#   if ( pvals[x-1] < 0.05001 )
#   {
#   labelvols<-volmatrix[,x-1]
#   if ( sum(labelvols[2:length(labelvols)]) > 0 )
#   {
#   labelvolsnorm<-labelvols/alllabelvols
#   labelVolsStr<-toString(labelvols)
#   labelVolsStr<-gsub(",","& ",labelVolsStr)
#   labelVolsStr<-gsub(" 0& "," . & ",labelVolsStr)
#   locanat<-c("")
#   for ( y in 2:length(ulabs) )
#   {
#     if ( labelvols[y-1] > 0 ) locanat<-paste(toString(locanat)," + ",anatnames[y-1])
#   }
##    print(paste(x-1,volumes[x],pvals[x-1],labelVolsStr),quote=F,zero.print = ".")
#   print(paste(x-1,volumes[x],pvals[x-1],locanat),quote=F,zero.print = ".")
#   }
#   }
# }
# q()
