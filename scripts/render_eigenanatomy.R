#!/usr/bin/env Rscript
Args <- commandArgs()
if ( length(Args) < 8  ){
fnm<-Args[4]
fnm<-substring(fnm,8,nchar(as.character(fnm)))
print(paste("Usage - RScript glassbrain.nii eigenanat_description_file.csv output_prefix "))
q()
}
ARGIND<-6
glassBrainFN<-c(as.character(Args[ARGIND]))
ARGIND<-ARGIND+1
eiganatViewsFN<-c(as.character(Args[ARGIND]))
ARGIND<-ARGIND+1
outname<-c(as.character(Args[ARGIND]))
ARGIND<-ARGIND+1
# libaries for rendering
library('AnalyzeFMRI') ; library('misc3d') ; library('rgl')
perms<-c(1,2,3)
img<-f.read.nifti.volume(glassBrainFN)
img <- aperm(img,c(perms,4))[,,1:dim(img)[3],1]
bimg<-(img > 0 )
# bimg<-GaussSmoothArray(bimg, voxdim=dim(bimg), ksize=11, sigma=diag(6, 3), mask=NULL, var.norm=FALSE)
# bimg<-(bimg > 0.25 )
brain <- contour3d(bimg, level = 1, alpha = 0.2, draw = FALSE,smooth=1,material="metal")
mylist<-list(brain)
mycolorlist<-colors()[ round(randu*length(colors()))[1:max(img),2] ]
mycolorlist<-c("red","blue","green","violet","yellow","tomato1","turquoise1","steelblue2","rosybrown2","skyblue4","orange4","maroon","gold","darkseagreen","burlywood","azure")
# list of files
projfiles<-read.csv(eiganatViewsFN,h=F)
for ( x in  1:nrow(projfiles) )
  {
    fn<-(as.character(projfiles[x,1]))
    print(paste("fn",fn))
    img<-f.read.nifti.volume(fn)
    img <- aperm(img,c(perms,4))[,,1:dim(img)[3],1]
    if ( max(img) > 0 ) {
      bimg<-(img > 0 )
      loccolor<-mycolorlist[x %% length(mycolorlist) ]
      lab <- contour3d(bimg, level = 1, alpha = 1, draw = F,smooth=1,color = loccolor )
      loclist<-list(lab)
      mylist <- c(mylist,loclist)
    }
  }
drawScene.rgl(c(mylist))
# write out some frames --- note that movie3d can be used instead to make a full movie of the results 
play3d(spin3d(axis=c(0,0,1), rpm=20), duration=0.9)
play3d(spin3d(axis=c(0,1,0), rpm=20), duration=0.2)
play3d(spin3d(axis=c(1,0,0), rpm=20), duration=0.05)
rgl.snapshot( paste(outname,'renderLeft.png',sep=''), fmt="png", top=TRUE )
# play3d(spin3d(axis=c(0,0,1), rpm=20), duration=1.5)
play3d(spin3d(axis=c(0,1,0), rpm=20), duration=0.8)
play3d(spin3d(axis=c(0,0,1), rpm=20), duration=0.85)
rgl.snapshot(  paste(outname,'renderAx.png',sep=''), fmt="png", top=TRUE )
# play3d(spin3d(axis=c(0,1,0), rpm=20), duration=1.5)
play3d(spin3d(axis=c(0,0,1), rpm=20), duration=0.8)
play3d(spin3d(axis=c(0,1,0), rpm=20), duration=0.8)
rgl.snapshot( paste(outname,'renderRight.png',sep=''), fmt="png", top=TRUE )
play3d(spin3d(axis=c(0,0,1), rpm=20), duration=(2.2))
rgl.snapshot( paste(outname,'renderCor.png',sep=''), fmt="png", top=TRUE )
