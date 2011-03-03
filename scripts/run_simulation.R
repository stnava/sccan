
 nsub<-50 ; nvoxy<-200 ; nvoxx<-300
 nvoxz<-1
# simulate true signal Z
Z<-matrix(c(1:nsub)/nsub-0.5,nrow=nsub,ncol=1)
Z<-(Z-mean(Z))/sd(Z)
# random signal
ysig<-matrix(rnorm(nvoxy,0,1),nrow=1,ncol=nvoxy)
Y<-Z%*%ysig+matrix(rnorm(nsub*nvoxy,0,1),nrow=nsub,ncol=nvoxy)
xsig<-matrix(rnorm(nvoxx,0,1),nrow=1,ncol=nvoxx) 
X<-Z%*%xsig+matrix(rnorm(nsub*nvoxx,0,1),nrow=nsub,ncol=nvoxx)

# produce corrupted 'true' signal 
zsig<-matrix(rnorm(nsub,0,1),nrow=1,ncol=nsub) 
Zs<-Z+t(zsig)*0.05
plot(Z,Zs)
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

spar<-0.5
exe1<-paste(CCA," --scca two-view[X.mhd,Y.mhd,X.mhd,Y.mhd,",spar,",",spar,"]   -o TEST1.nii.gz -p 100   " )
exe2<-paste(CCA," --scca partial[X.mhd,Y.mhd,Z.mhd,X.mhd,Y.mhd,Z.mhd,",spar,",",spar,", -1]   -o TEST2.nii.gz -p 100  --partial-scca-option PminusRQminusR ")

bb1<-try(system(exe1, intern = TRUE, ignore.stderr = TRUE))
print(bb1[111])
bb2<-try(system(exe2, intern = TRUE, ignore.stderr = TRUE))
print(bb2[133])
