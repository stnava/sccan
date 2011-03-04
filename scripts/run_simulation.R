#################################
# simulation studies for pscca  # 
#################################
sparseness<-0.5 # for X, Y matrices
testspatiallocalization<-1 # tests non-overlapping signals 
 nsub<-40 ; nvoxy<-300 ; nvoxx<-180 ; 
 nvoxz<-1
 noise<-0.05
# simulate true signal Z
Z<-matrix(c(1:nsub)/nsub,nrow=nsub,ncol=1)
Z<-(Z-mean(Z))/sd(Z)
# random signal
ysig<-matrix(rnorm(nvoxy,0,1),nrow=1,ncol=nvoxy)
Y<-Z%*%ysig+matrix(rnorm(nsub*nvoxy,0,1),nrow=nsub,ncol=nvoxy)
xsig<-matrix(rnorm(nvoxx,0,1),nrow=1,ncol=nvoxx) 
X<-Z%*%xsig+matrix(rnorm(nsub*nvoxx,0,1),nrow=nsub,ncol=nvoxx)

# repeat for different true signal Z2
Z2<-matrix(c(1:nsub)/nsub,nrow=nsub,ncol=1)
Z2<-exp(Z2)# *Z2
zsig<-matrix(rnorm(nsub,0,1),nrow=1,ncol=nsub) 
Z2<-Z2+t(zsig)*noise
plot(Z2)
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
nvoxx<-nvoxx*2
nvoxy<-nvoxy*2
}
# produce corrupted 'true' signal 
zsig<-matrix(rnorm(nsub,0,1),nrow=1,ncol=nsub) 
Zs<-Z+t(zsig)*noise
Z<-Zs

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
DimSize = ",nsub,1,"
AnatomicalOrientation = ??
ElementType = MET_FLOAT
ElementDataFile = ",fn)
print(oname)
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
print(oname)
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
print(oname)
write.table(mhd,oname,row.names=F,col.names=F,sep=" ",quote=FALSE) 
writeBin(c(as.real(Y)),fn,size=4,endian = .Platform$endian)


# define the CCA progam 
CCA<-"~/code/sccan/bin/sccan "

exe<-"for N in X Y Z ; do StackSlices ${N}mask.nii.gz 0 -1 -1 $N.mhd ; ThresholdImage 2 ${N}mask.nii.gz  ${N}mask.nii.gz  -9.e9 9.e9 ; SmoothImage 2 ${N}.mhd 0.5 ${N}.mhd ;  done "
bb<-try(system(exe, intern = TRUE, ignore.stderr = TRUE))

exe1<-paste(CCA," --scca two-view[X.mhd,Y.mhd,Xmask.nii.gz,Ymask.nii.gz,",sparseness,",",sparseness,"]   -o TEST1.nii.gz -p 100   " )
exe2<-paste(CCA," --scca partial[X.mhd,Y.mhd,Z.mhd,Xmask.nii.gz,Ymask.nii.gz,Zmask.nii.gz,",sparseness,",",sparseness,", -1]   -o TEST2.nii.gz -p 100  --partial-scca-option PminusRQminusR ")

print(exe1)
print(exe2)
bb1<-try(system(exe1, intern = TRUE, ignore.stderr = TRUE))
print(bb1[111])
bb2<-try(system(exe2, intern = TRUE, ignore.stderr = TRUE))
print(bb2[133])
