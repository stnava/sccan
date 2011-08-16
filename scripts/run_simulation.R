#!/usr/bin/Rscript
#################################
# simulation studies for pscca  # 
#################################
# define the CCA progam 
CCA<-"~/code/sccan/bin/sccan "
# define the parameters for the simulation 
 sparseness<-(-0.33) # for X, Y matrices
 testspatiallocalization<-1 # tests non-overlapping signals 
 nsub<-75 ; nvoxy<-100*(2-testspatiallocalization) ; nvoxx<-80*(2-testspatiallocalization) ; # size of simulated images 
# nsub<-20 ; nvoxy<-1 ; nvoxx<-5 ; # size of simulated images 
 noise<-0.05 # increase this to get extra noise
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
Z2<-exp(Z2*Z2*Z2)
# Z2<-(Z2)^(1/3)*(-1.)
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
Z<-matrix(c(Z),nrow=nsub,ncol=1)
# Z<-matrix(c(Z,Z2),nrow=nsub,ncol=2)
nvoxx<-nvoxx*2
nvoxy<-nvoxy*2
}

# Y<-Y[sample(c(1:nsub)),]
# X<-X[sample(c(1:nsub)),]

# write out the simulated data
#
write.csv(X,"X.csv",row.names =F)
write.csv(Y,"Y.csv",row.names =F)
write.csv(Z,"Z.csv",row.names =F)

exe1<-paste(CCA," --scca two-view[X.csv,Y.csv,na,na,",sparseness,",",sparseness,"]   -o TEST1.nii.gz -p 100 -i 50   " )
exe2<-paste(CCA," --scca partial[X.csv,Y.csv,Z.csv,na,na,na,",sparseness,",",sparseness,", -1]   -o TEST2.nii.gz  -e 0 -n 5 -i 50 -r 0 -p 20  --partial-scca-option PQminusR")

# print(exe1)
 print(exe2,quote=F)
q()
# bb1<-try(system(exe1, intern = TRUE, ignore.stderr = TRUE))
bb2<-try(system(exe2, intern = TRUE, ignore.stderr = TRUE))
print(bb2)

pv1<-unlist(strsplit(bb1[110], " ", fixed = TRUE))
pv2<-unlist(strsplit(bb2[132], " ", fixed = TRUE))
# print(paste(" scca: p ",pv1[3]," corr ",pv1[7] ," pscca : ",pv2[3]," corr ",pv2[7]))
ptot1[sim]<-as.real(pv1[3])
ptot2[sim]<-as.real(pv2[3])
# print(ptot2)
print(paste("AVG results at sim",sim," with sparseness ",sparseness," - scca: p ",mean(ptot1[1:sim])," pscca : p ",mean(ptot2[1:sim])))
}
