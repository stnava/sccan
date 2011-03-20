#!/usr/bin/Rscript
#################################
# simulation studies for pscca  # 
#################################
# define the CCA progam 
CCA<-"~/code/sccan/bin/sccan "
# define the parameters for the simulation 
 sparseness<-(-0.33) # for X, Y matrices
 testspatiallocalization<-1 # tests non-overlapping signals 
 nsub<-40 ; nvoxy<-117*(2-testspatiallocalization) ; nvoxx<-225*(2-testspatiallocalization) ; # size of simulated images 
 noise<-0.02 # increase this to get extra noise
 totalSimulations<-1 # number of random repeats 
 nvoxz<-1 # ignore this 
ptot1<-rep(1,totalSimulations)
ptot2<-rep(1,totalSimulations)
for ( sim in c(1:totalSimulations) ) {
# simulate true signal Z
Z<-matrix(c(1:nsub)/nsub,nrow=nsub,ncol=1)
Z<-(Z-mean(Z))/sd(Z)
# random signal
ysig<-matrix(rnorm(nvoxy,0,1),nrow=1,ncol=nvoxy)
Y<-Z%*%ysig+matrix(rnorm(nsub*nvoxy,0,1),nrow=nsub,ncol=nvoxy)
xsig<-matrix(rnorm(nvoxx,0,1),nrow=1,ncol=nvoxx) 
X<-Z%*%xsig+matrix(rnorm(nsub*nvoxx,0,1),nrow=nsub,ncol=nvoxx)

# produce corrupted 'true' signal - covariate
zsig<-matrix(rnorm(nsub,0,1),nrow=1,ncol=nsub) 
Zs<-Z+t(zsig)*noise
Z<-Zs

# repeat for different true signal Z2
Z2<-matrix(c(1:nsub)/nsub,nrow=nsub,ncol=1)
Z2<-exp(Z2*Z2)# *Z2
zsig<-matrix(rnorm(nsub,0,1),nrow=1,ncol=nsub) 
Z2<-Z2+t(zsig)*noise
# plot(Z2)
Z2<-sample(Z2) # permute the signal 

# random signal
ysig<-matrix(rnorm(nvoxy,0,1),nrow=1,ncol=nvoxy)
Y2<-Z2%*%ysig+matrix(rnorm(nsub*nvoxy,0,1),nrow=nsub,ncol=nvoxy)
xsig<-matrix(rnorm(nvoxx,0,1),nrow=1,ncol=nvoxx) 
X2<-Z2%*%xsig+matrix(rnorm(nsub*nvoxx,0,1),nrow=nsub,ncol=nvoxx)

# concatenate the matrices made of 2 signals
if ( testspatiallocalization == 1 ) {
X<-matrix(c(X,X2),nrow=nsub,ncol=nvoxx*2) 
Y<-matrix(c(Y,Y2),nrow=nsub,ncol=nvoxy*2)
# Z<-matrix(c(Z,Z2),nrow=nsub,ncol=2)
nvoxx<-nvoxx*2
nvoxy<-nvoxy*2
}

# Y<-Y[sample(c(1:nsub)),]
# X<-X[sample(c(1:nsub)),]

# write out the simulated data
#
DIR<-"/Users/brianavants/Data/PSCCA_simulation/"
pre<-paste(DIR,"Z",sep='')
oname<-paste(pre,".mhd",sep='')
fn<-paste(pre,".raw",sep='')
mhd<-paste("ObjectType = Image
NDims = 2
BinaryData = True
BinaryDataByteOrderMSB = False
CompressedData = False
TransformMatrix = 1 0 0 1
Offset = 0 0
CenterOfRotation = 0 0
ElementSpacing = 1 1
DimSize = ",nsub,ncol(Z),"
AnatomicalOrientation = ??
ElementType = MET_FLOAT
ElementDataFile = ",fn)
# print(oname)
write.table(mhd,oname,row.names=F,col.names=F,sep=" ",quote=FALSE) 
writeBin(c(as.real(Z)),fn,size=4,endian = .Platform$endian)

pre<-paste(DIR,"X",sep='')
oname<-paste(pre,".mhd",sep='')
fn<-paste(pre,".raw",sep='')
mhd<-paste("ObjectType = Image
NDims = 2
BinaryData = True
BinaryDataByteOrderMSB = False
CompressedData = False
TransformMatrix = 1 0 0 1
Offset = 0 0
CenterOfRotation = 0 0
ElementSpacing = 1 1
DimSize = ",nsub,nvoxx,"
AnatomicalOrientation = ??
ElementType = MET_FLOAT
ElementDataFile = ",fn)
# print(oname)
write.table(mhd,oname,row.names=F,col.names=F,sep=" ",quote=FALSE) 
writeBin(c(as.real(X)),fn,size=4,endian = .Platform$endian)


pre<-paste(DIR,"Y",sep='')
oname<-paste(pre,".mhd",sep='')
fn<-paste(pre,".raw",sep='')
mhd<-paste("ObjectType = Image
NDims = 2
BinaryData = True
BinaryDataByteOrderMSB = False
CompressedData = False
TransformMatrix = 1 0 0 1
Offset = 0 0
CenterOfRotation = 0 0
ElementSpacing = 1 1
DimSize = ",nsub,nvoxy,"
AnatomicalOrientation = ??
ElementType = MET_FLOAT
ElementDataFile = ",fn)
# print(oname)
write.table(mhd,oname,row.names=F,col.names=F,sep=" ",quote=FALSE) 
writeBin(c(as.real(Y)),fn,size=4,endian = .Platform$endian)


exe<-"for N in X Y Z ; do StackSlices ${N}mask.nii.gz 0 -1 -1 $N.mhd ; ThresholdImage 2 ${N}mask.nii.gz  ${N}mask.nii.gz  -9.e9 9.e9 ; SmoothImage 2 ${N}.mhd 0 ${N}.mhd ;  done "
bb<-try(system(exe, intern = TRUE, ignore.stderr = TRUE))

exe1<-paste(CCA," --scca two-view[X.mhd,Y.mhd,Xmask.nii.gz,Ymask.nii.gz,",sparseness,",",sparseness,"]   -o TEST1.nii.gz -p 100   " )
exe2<-paste(CCA," --scca partial[X.mhd,Y.mhd,Z.mhd,Xmask.nii.gz,Ymask.nii.gz,Zmask.nii.gz,",sparseness,",",sparseness,", -1]   -o TEST2.nii.gz -p 0  --partial-scca-option PminusRQminusR ")

# print(exe1)
 print(exe2,quote=F)
# bb1<-try(system(exe1, intern = TRUE, ignore.stderr = TRUE))
bb2<-try(system(exe2, intern = TRUE, ignore.stderr = TRUE))
print(bb2)
q()

pv1<-unlist(strsplit(bb1[110], " ", fixed = TRUE))
pv2<-unlist(strsplit(bb2[132], " ", fixed = TRUE))
# print(paste(" scca: p ",pv1[3]," corr ",pv1[7] ," pscca : ",pv2[3]," corr ",pv2[7]))
ptot1[sim]<-as.real(pv1[3])
ptot2[sim]<-as.real(pv2[3])
# print(ptot2)
print(paste("AVG results at sim",sim," with sparseness ",sparseness," - scca: p ",mean(ptot1[1:sim])," pscca : p ",mean(ptot2[1:sim])))
}
