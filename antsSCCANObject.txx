/*=========================================================================


  Program:   Advanced Normalization Tools
  Module:    $RCSfile: antsSCCANObject.txx,v $
  Language:  C++
  Date:      $Date: $
  Version:   $Revision: $

  Copyright (c) ConsortiumOfANTS. All rights reserved.
  See accompanying COPYING.txt or
  http://sourceforge.net/projects/advants/files/ANTS/ANTSCopyright.txt
  for details.

  This software is distributed WITHOUT ANY WARRANTY; without even
  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
  PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#include <vnl/vnl_random.h>
#include <vnl/algo/vnl_matrix_inverse.h>
#include "antsSCCANObject.h"

namespace itk {
namespace ants {

template <class TInputImage, class TRealType>
antsSCCANObject<TInputImage, TRealType>::antsSCCANObject( ) 
{
  this->m_CorrelationForSignificanceTest=0;
  this->m_SpecializationForHBM2011=false;
  this->m_AlreadyWhitened=false;
  this->m_PinvTolerance=1.e-5;
  this->m_PercentVarianceForPseudoInverse=0.95;
  this->m_MaximumNumberOfIterations=100;
  this->m_MaskImageP=NULL;
  this->m_MaskImageQ=NULL;
  this->m_MaskImageR=NULL;
  this->m_KeepPositiveP=true;
  this->m_KeepPositiveQ=true;
  this->m_KeepPositiveR=true;
  this->m_FractionNonZeroP=0.5;
  this->m_FractionNonZeroQ=0.5;
  this->m_FractionNonZeroR=0.5;
  this->m_ConvergenceThreshold=1.e-6;
} 

template <class TInputImage, class TRealType>
typename antsSCCANObject<TInputImage, TRealType>::VectorType
antsSCCANObject<TInputImage, TRealType>
::InitializeV( typename antsSCCANObject<TInputImage, TRealType>::MatrixType p ) 
{
  VectorType w_p( p.columns() );
  vnl_random randgen(time(0));
  for (unsigned long i=0; i < p.columns(); i++)
    { 
//      w_p(i)=randgen.drand32();//1.0/p.rows();//
      w_p(i)=1.0+fabs(randgen.drand32())*1.e-3; 
    }
  w_p=w_p/w_p.two_norm();
  return w_p;
}

template <class TInputImage, class TRealType>
typename antsSCCANObject<TInputImage, TRealType>::MatrixType
antsSCCANObject<TInputImage, TRealType>
::NormalizeMatrix( typename antsSCCANObject<TInputImage, TRealType>::MatrixType p ) 
{
  vnl_random randgen(time(0));
  MatrixType np( p.rows() , p.columns() );
  for (unsigned long i=0; i < p.columns(); i++)
  { 
    VectorType wpcol=p.get_column(i);
    VectorType wpcol2=wpcol-wpcol.mean();
    double sd=wpcol2.squared_magnitude();
    sd=sqrt( sd/(p.rows()-1) );
    if ( sd <= 0 ) {
      std::cout << " bad-row " << i <<  wpcol << std::endl;
      for (unsigned long j=0; j < wpcol.size(); j++)
	wpcol2(j)=randgen.drand32();
      wpcol2=wpcol2-wpcol2.mean();
      sd=wpcol2.squared_magnitude();
      sd=sqrt( sd/(p.rows()-1) );
    }
    wpcol=wpcol2/sd;
    np.set_column(i,wpcol);
  }
  return np;
}


template <class TInputImage, class TRealType>
typename antsSCCANObject<TInputImage, TRealType>::MatrixType
antsSCCANObject<TInputImage, TRealType>
::WhitenMatrixOrGetInverseCovarianceMatrix( typename antsSCCANObject<TInputImage, TRealType>::MatrixType rin , bool whiten_else_invcovmatrix ) 
{
  double regularization=1.e-0;
  if (  rin.columns() > rin.rows() ) 
    {
      MatrixType dd=rin*rin.transpose();
      dd.set_identity(); 
      dd=dd*regularization+rin*rin.transpose();
      vnl_svd_economy<RealType> eig(dd);
//      vnl_svd_economy<RealType> eig(rin*rin.transpose());
      VectorType eigvals=eig.lambdas();
      RealType eigsum=eigvals.sum();
      RealType total=0; 
      unsigned int eigct=0;
      while ( total/eigsum < this->m_PercentVarianceForPseudoInverse ) 
	{
	  total+=eigvals(eigct);
	  eigct++;
	}
      DiagonalMatrixType r_diag_inv(eigct);
      MatrixType r_eigvecs(eig.V().get_column(0).size(),eigct);
      r_eigvecs.fill(0);
      for (unsigned int j=0; j<eigct; j++) 
	{
	  RealType eval=eigvals(j);
	  if ( eval > this->m_PinvTolerance )  {// FIXME -- check tolerances against matlab pinv
	   if ( whiten_else_invcovmatrix )
	    r_diag_inv(j)=1/sqrt(eval);// need sqrt for whitening 
	   else 
	    r_diag_inv(j)=1/(eval);// need eval for inv cov 
	    r_eigvecs.set_column(j,eig.V().get_column(j));
	  }
	  else r_diag_inv(j)=0;// need sqrt for whitening 
	}
      MatrixType evecs=rin.transpose()*r_eigvecs;
     if ( whiten_else_invcovmatrix ) 
       return (rin*evecs)*(r_diag_inv*evecs.transpose());
     else // return inv cov matrix 
       return (evecs)*(r_diag_inv*evecs.transpose());
    }
  else 
    {     
      MatrixType dd=rin.transpose()*rin;
      dd.set_identity(); 
      dd=dd*regularization+rin.transpose()*rin;
      vnl_svd_economy<RealType> eig(dd);
      VectorType eigvals=eig.lambdas();
      RealType eigsum=eigvals.sum();
      RealType total=0; 
      unsigned int eigct=0;
      while ( total/eigsum < this->m_PercentVarianceForPseudoInverse ) 
	{
	  total+=eigvals(eigct);
	  eigct++;
	}
      DiagonalMatrixType r_diag_inv(eigct);
      MatrixType r_eigvecs(eig.V().get_column(0).size(),eigct);
      r_eigvecs.fill(0);
      for (unsigned int j=0; j<eigct; j++) 
	{
	  RealType eval=eigvals(j);
	  if ( eval > this->m_PinvTolerance )  {// FIXME -- check tolerances against matlab pinv
	   if ( whiten_else_invcovmatrix )
	    r_diag_inv(j)=1/sqrt(eval);// need sqrt for whitening 
	   else 
	    r_diag_inv(j)=1/(eval);// need eval for inv cov
	    r_eigvecs.set_column(j,eig.V().get_column(j));
	  }
	  else r_diag_inv(j)=0;// need sqrt for whitening 
	}
      MatrixType wmatrix=r_eigvecs*r_diag_inv*r_eigvecs.transpose();
      
     if ( whiten_else_invcovmatrix ) 
       return (rin*wmatrix);
      else // return inv cov matrix 
       return wmatrix;
    }

}



template <class TInputImage, class TRealType>
typename antsSCCANObject<TInputImage, TRealType>::VectorType
antsSCCANObject<TInputImage, TRealType>
::SoftThreshold( typename antsSCCANObject<TInputImage, TRealType>::VectorType
 v_in, TRealType fractional_goal , bool allow_negative_weights )
{
//  std::cout <<" allow neg weights? " << allow_negative_weights << std::endl;
  VectorType v_out(v_in);
  if ( fractional_goal > 1 ) return v_out;
  RealType minv=v_in.min_value();
  RealType maxv=v_in.max_value();
  if ( fabs(v_in.min_value()) > maxv ) maxv=fabs(v_in.min_value());
  minv=0;
  RealType lambg=1.e-3;
  RealType frac=0;
  unsigned int its=0,ct=0;
  RealType soft_thresh=lambg;
 
  RealType minthresh=0,minfdiff=1;
  unsigned int maxits=1000;
  for ( its=0; its<maxits; its++) 
  {
    soft_thresh=(its/(float)maxits)*maxv;
    ct=0; 
    for ( unsigned int i=0; i<v_in.size(); i++) {
      RealType val=v_in(i);
      if ( allow_negative_weights ) val=fabs(val);
      if ( val < soft_thresh ) 
      {
	v_out(i)=0;
	ct++;
      }
      else v_out(i)=v_in(i);
    }
    frac=(float)(v_in.size()-ct)/(float)v_in.size();
    //    std::cout << " cur " << frac << " goal "  << fractional_goal << " st " << soft_thresh << " th " << minthresh << std::endl;
    if ( fabs(frac - fractional_goal) < minfdiff ) {
      minthresh=soft_thresh;
      minfdiff= fabs(frac - fractional_goal) ;
    }
  }
  ct=0;
  for ( unsigned int i=0; i<v_in.size(); i++) {
    RealType val=v_in(i);
    if ( allow_negative_weights ) val=fabs(val);
    if ( val < minthresh ) 
      {
	v_out(i)=0;
	ct++;
      }
    else v_out(i)=v_in(i);
  }
  frac=(float)(v_in.size()-ct)/(float)v_in.size();
  // std::cout << " frac non-zero " << frac << " wanted " << fractional_goal << " allow-neg " << allow_negative_weights << std::endl;
  return v_out;
}


template <class TInputImage, class TRealType>
typename antsSCCANObject<TInputImage, TRealType>::VectorType
antsSCCANObject<TInputImage, TRealType>
::TrueCCAPowerUpdate( TRealType penalty1,  typename antsSCCANObject<TInputImage, TRealType>::MatrixType p , typename antsSCCANObject<TInputImage, TRealType>::VectorType  w_q ,  typename antsSCCANObject<TInputImage, TRealType>::MatrixType q, bool keep_pos ,   typename antsSCCANObject<TInputImage, TRealType>::VectorType covariate, bool factorOutR )
{
  RealType norm=0;
  // recall that the matrices below have already be whitened ....
  // we bracket the computation and use associativity to make sure its done efficiently 
  //vVector wpnew=( (CppInv.transpose()*p.transpose())*(CqqInv*q) )*w_q;
  VectorType wpnew;
  if ( factorOutR ){
    VectorType temp=q*w_q;
    wpnew=p.transpose()*( temp - this->m_MatrixRRt*temp ); 
  }
  else {
    VectorType temp=q*w_q;
// we now orthogonalize temp wrt the covariate 
    if ( covariate.two_norm() > 0 ) { //gram-schmidt 
      double a=inner_product(temp,covariate);
      double b=inner_product(covariate,covariate);
      double ratio=a/b;
      temp=temp-ratio*covariate;
    }
    wpnew=p.transpose()*temp;
  }
  wpnew=this->SoftThreshold( wpnew , penalty1 , !keep_pos );
  norm=wpnew.two_norm();
  if ( norm > 0 ) wpnew=wpnew/(norm);
  return wpnew;
}

template <class TInputImage, class TRealType>
void
antsSCCANObject<TInputImage, TRealType>
::UpdatePandQbyR( ) 
{
// R is already whitened 
    switch( this->m_SCCANFormulation )
    {
    case PQ:
      {
// do nothing
      break;
      }
    case PminusRQ:
      {
        this->m_MatrixP=(this->m_MatrixP-this->m_MatrixRRt*this->m_MatrixP);
      break;
      }
    case PQminusR:
      {
        this->m_MatrixQ=(this->m_MatrixQ-this->m_MatrixRRt*this->m_MatrixQ);
      break;
      }
    case PminusRQminusR :
      {
/** P_R =   P - R_w R_w^T P */
/** Q_R =   Q - R_w R_w^T Q */
        this->m_MatrixP=(this->m_MatrixP-this->m_MatrixRRt*this->m_MatrixP);
        this->m_MatrixQ=(this->m_MatrixQ-this->m_MatrixRRt*this->m_MatrixQ);
      break;
    case PQR :
      std::cout <<" You should call mscca not pscca " << std::endl;
      break;
      }
    }

}

template <class TInputImage, class TRealType>
TRealType
antsSCCANObject<TInputImage, TRealType>
::RunSCCAN2multiple( unsigned int n_vecs ) 
{
  unsigned long projection_iteration_limit=500;
  RealType truecorr=0;
  unsigned int nr1=this->m_MatrixP.rows();
  unsigned int nr2=this->m_MatrixQ.rows();
  if ( nr1 != nr2 ) 
  {
    std::cout<< " P rows " << this->m_MatrixP.rows() << " cols " << this->m_MatrixP.cols() << std::endl;
    std::cout<< " Q rows " << this->m_MatrixQ.rows() << " cols " << this->m_MatrixQ.cols() << std::endl;
    std::cout<< " R rows " << this->m_MatrixR.rows() << " cols " << this->m_MatrixR.cols() << std::endl;
    std::cout<<" N-rows for MatrixP does not equal N-rows for MatrixQ " << nr1 << " vs " << nr2 << std::endl;
    exit(1);
  }
  if (  this->m_CovariatesP.size() == 0 ){
    this->m_CovariatesP.set_size( nr1 );
    this->m_CovariatesQ.set_size( nr2 );
    this->m_CovariatesP.fill(0);
    this->m_CovariatesQ.fill(0);
  }
  this->m_VariatesP.set_size(this->m_MatrixP.cols(),n_vecs);
  this->m_VariatesQ.set_size(this->m_MatrixQ.cols(),n_vecs);
  this->m_VariatesP.fill(0);
  this->m_VariatesQ.fill(0);
  TRealType corr=0;
  unsigned int j=0; 
  bool notdone=true;
  while ( notdone ) {
    this->m_MatrixP=this->NormalizeMatrix(this->m_OriginalMatrixP);  
    this->m_MatrixQ=this->NormalizeMatrix(this->m_OriginalMatrixQ);  
    this->m_MatrixP=this->WhitenMatrix(this->m_MatrixP);  
    this->m_MatrixQ=this->WhitenMatrix(this->m_MatrixQ);
    if ( this->m_MatrixR.size() > 0 ) {
      this->m_MatrixR=this->NormalizeMatrix(this->m_MatrixR);  
      this->m_MatrixR=this->WhitenMatrix(this->m_MatrixR);  
      this->m_MatrixRRt=this->m_MatrixR*this->m_MatrixR.transpose(); 
      this->UpdatePandQbyR( );
    }
    std::cout << " get canonical variate number " << j+1 << std::endl;
    double initcorr=1.e-5;
    truecorr=initcorr;
    double deltacorr=1,lastcorr=initcorr*0.5;
    this->m_WeightsP=this->InitializeV(this->m_MatrixP);
    this->m_WeightsQ=this->InitializeV(this->m_MatrixQ);
    unsigned long its=0;
    while ( its < this->m_MaximumNumberOfIterations && deltacorr > this->m_ConvergenceThreshold  )
    {
      unsigned long proj_lim=projection_iteration_limit;
      if ( its == 0 && j > 0) { 
        proj_lim=1;
	MatrixType temp;
	MatrixType pp;
	pp.set_size(this->m_MatrixP.rows(),j);
	for (unsigned int kk=0; kk<j; kk++) pp.set_column(kk,this->m_MatrixP*this->m_VariatesP.get_column(kk));
        temp=this->NormalizeMatrix(pp);  
        temp=this->WhitenMatrix(temp);  
        temp=temp*temp.transpose(); 
        this->m_MatrixP=(this->m_MatrixP-temp*this->m_MatrixP);
	MatrixType qq;
	qq.set_size(this->m_MatrixQ.rows(),j);
	for (unsigned int kk=0; kk<j; kk++) qq.set_column(kk,this->m_MatrixQ*this->m_VariatesQ.get_column(kk));
        temp=this->NormalizeMatrix(qq);  
        temp=this->WhitenMatrix(temp);  
        temp=temp*temp.transpose(); 
        this->m_MatrixQ=(this->m_MatrixQ-temp*this->m_MatrixQ);
      }
      {
      RealType norm=0;
      double ip=1; unsigned long ct=0;
      VectorType pj=this->m_MatrixP.transpose()*(this->m_MatrixQ*this->m_WeightsQ );
      pj=pj/pj.two_norm();
      while ( fabs(ip) > 1.e-2 && ct < proj_lim ) {
        ip=0;
        VectorType bb;
        for (unsigned int kk=0; kk<j; kk++) {
          bb=this->m_VariatesP.get_column(kk);
          pj=pj - inner_product(this->m_MatrixP*pj,this->m_MatrixP*bb)/inner_product(this->m_MatrixP*bb,this->m_MatrixP*bb)*bb;
        }
        pj=this->SoftThreshold( pj , this->m_FractionNonZeroP , !this->m_KeepPositiveP );
        norm=pj.two_norm();
        if ( norm > 0 ) this->m_WeightsP=pj/(norm);
        else this->m_WeightsP.fill(0);
        for (unsigned int kk=0; kk<j; kk++) {
          bb=this->m_VariatesP.get_column(kk);
          ip+=inner_product(this->m_MatrixP*pj,this->m_MatrixP*bb);
        }
	ip/=(RealType)j;
	ct++;
 	std::cout << " pip " << ip << " ct " << ct <<  std::endl;
      }
      this->m_VariatesP.set_column(j,this->m_WeightsP);
      }

      {
      RealType norm=0;
      double ip=1; unsigned long ct=0;
      VectorType qj=this->m_MatrixQ.transpose()*(this->m_MatrixP*this->m_WeightsP);
      qj=qj/qj.two_norm();
      while ( fabs(ip) > 1.e-2 && ct < proj_lim ) {
        ip=0;
        VectorType bb;
        for (unsigned int kk=0; kk<j; kk++) {
          bb=this->m_VariatesQ.get_column(kk);
          qj=qj - inner_product(this->m_MatrixQ*qj,this->m_MatrixQ*bb)/inner_product(this->m_MatrixQ*bb,this->m_MatrixQ*bb)*bb;
        }
        qj=this->SoftThreshold( qj , this->m_FractionNonZeroQ , !this->m_KeepPositiveQ );
        norm=qj.two_norm();
        if ( norm > 0 ) this->m_WeightsQ=qj/(norm);
        else this->m_WeightsQ.fill(0);
        for (unsigned int kk=0; kk<j; kk++) {
          bb=this->m_VariatesQ.get_column(kk);
          ip+=inner_product(this->m_MatrixQ*qj,this->m_MatrixQ*bb);
        }
	ip/=(RealType)j;
	ct++;
 	std::cout << " qip " << ip << " ct " << ct << std::endl;
      }
      this->m_VariatesQ.set_column(j,this->m_WeightsQ);
      }
      truecorr=this->PearsonCorr( this->m_MatrixP*this->m_WeightsP , this->m_MatrixQ*this->m_WeightsQ );
      deltacorr=(truecorr-lastcorr);
      lastcorr=truecorr;
      ++its;
      std::cout << " p-ip-a " << inner_product( this->m_MatrixP*this->m_VariatesP.get_column(1) , this->m_MatrixP*this->m_VariatesP.get_column(0) )  << " q-ip " << inner_product( this->m_MatrixQ*this->m_VariatesQ.get_column(1) , this->m_MatrixQ*this->m_VariatesQ.get_column(0) ) << " corr " << truecorr << std::endl;
      std::cout << " p-ip-b " << inner_product( this->m_MatrixP*this->m_VariatesP.get_column(2) , this->m_MatrixP*this->m_VariatesP.get_column(0) )  << " q-ip " << inner_product( this->m_MatrixQ*this->m_VariatesQ.get_column(2) , this->m_MatrixQ*this->m_VariatesQ.get_column(0) ) << " corr " << truecorr << std::endl;
      std::cout << " p-ip-c " << inner_product( this->m_MatrixP*this->m_VariatesP.get_column(2) , this->m_MatrixP*this->m_VariatesP.get_column(1) )  << " q-ip " << inner_product( this->m_MatrixQ*this->m_VariatesQ.get_column(2) , this->m_MatrixQ*this->m_VariatesQ.get_column(1) ) << " corr " << truecorr << std::endl;
    }// inner_it
    std::cout << "  canonical variate number " << j+1 << " corr " << truecorr << std::endl;
    if ( fabs(truecorr) < 1.e-2 || (j+1) == n_vecs ) notdone=false;
    else {
     if ( j == 0 )
       corr+=fabs(truecorr);
     j++;
    }
  }   
/*
  std::cout <<"  this->m_VariatesP " << this->m_VariatesP  << std::endl;
  MatrixType projp=this->m_MatrixP*this->m_VariatesP;
  MatrixType projq=this->m_MatrixQ*this->m_VariatesQ;
  MatrixType cr=projp.transpose()*projq;
  MatrixType denom=this->InverseCovarianceMatrix(projp)*this->InverseCovarianceMatrix(projq);
  cr=cr*denom;
  std::cout << " cr " << cr << std::endl;
*/
  this->m_WeightsP=this->m_VariatesP.get_column(0);
  this->m_WeightsQ=this->m_VariatesQ.get_column(0);

  return corr;
}

template <class TInputImage, class TRealType>
TRealType
antsSCCANObject<TInputImage, TRealType>
::RunSCCAN2( ) 
{
  RealType truecorr=0;
  unsigned int nr1=this->m_MatrixP.rows();
  unsigned int nr2=this->m_MatrixQ.rows();
  if ( nr1 != nr2 ) 
  {
    std::cout<< " P rows " << this->m_MatrixP.rows() << " cols " << this->m_MatrixP.cols() << std::endl;
    std::cout<< " Q rows " << this->m_MatrixQ.rows() << " cols " << this->m_MatrixQ.cols() << std::endl;
    std::cout<< " R rows " << this->m_MatrixR.rows() << " cols " << this->m_MatrixR.cols() << std::endl;
    std::cout<<" N-rows for MatrixP does not equal N-rows for MatrixQ " << nr1 << " vs " << nr2 << std::endl;
    exit(1);
  }
  else {
//  std::cout << " P-positivity constraints? " <<  this->m_KeepPositiveP << " frac " << this->m_FractionNonZeroP << " Q-positivity constraints?  " << m_KeepPositiveQ << " frac " << this->m_FractionNonZeroQ << std::endl;
  }
  if (  this->m_CovariatesP.size() == 0 ){
    this->m_CovariatesP.set_size( nr1 );
    this->m_CovariatesQ.set_size( nr2 );
    this->m_CovariatesP.fill(0);
    this->m_CovariatesQ.fill(0);
  }
  else {
//    std::cout <" P Covar " << this->m_CovariatesP << std::endl;
//    std::cout <" Q Covar " << this->m_CovariatesQ << std::endl;
  }
  this->m_WeightsP=this->InitializeV(this->m_MatrixP);
  this->m_WeightsQ=this->InitializeV(this->m_MatrixQ);

  if ( !this->m_AlreadyWhitened ) {
  std::cout <<" norm P " << std::endl;
    this->m_MatrixP=this->NormalizeMatrix(this->m_MatrixP);  
  std::cout <<" norm Q " << std::endl;
    this->m_MatrixQ=this->NormalizeMatrix(this->m_MatrixQ);  
    if ( this->m_MatrixR.size() > 0 ) {
      this->m_MatrixR=this->NormalizeMatrix(this->m_MatrixR);  
      this->m_MatrixR=this->WhitenMatrix(this->m_MatrixR);  
      this->m_MatrixRRt=this->m_MatrixR*this->m_MatrixR.transpose(); 
      this->UpdatePandQbyR( );
    }
    this->m_MatrixP=this->WhitenMatrix(this->m_MatrixP);  
    this->m_MatrixQ=this->WhitenMatrix(this->m_MatrixQ);
    this->m_AlreadyWhitened=true;
  }  
  for (unsigned int outer_it=0; outer_it<2; outer_it++) {
  truecorr=0;
  double deltacorr=1,lastcorr=1;
  unsigned long its=0;
  while ( its < this->m_MaximumNumberOfIterations && deltacorr > this->m_ConvergenceThreshold  )
  {
    this->m_WeightsP=this->TrueCCAPowerUpdate(this->m_FractionNonZeroP,this->m_MatrixP,this->m_WeightsQ,this->m_MatrixQ,this->m_KeepPositiveP,this->m_CovariatesQ,false);
    this->m_WeightsQ=this->TrueCCAPowerUpdate(this->m_FractionNonZeroQ,this->m_MatrixQ,this->m_WeightsP,this->m_MatrixP,this->m_KeepPositiveQ,this->m_CovariatesP,false);    
    truecorr=this->PearsonCorr( this->m_MatrixP*this->m_WeightsP , this->m_MatrixQ*this->m_WeightsQ );
    deltacorr=fabs(truecorr-lastcorr);
    lastcorr=truecorr;
    ++its;
  }// inner_it 
  }//outer_it 
  this->m_CorrelationForSignificanceTest=truecorr;
  if ( this->m_SpecializationForHBM2011 ){
    truecorr=this->SpecializedCorrelation();
  }
  this->m_CorrelationForSignificanceTest=truecorr;
  return truecorr;
}

template <class TInputImage, class TRealType>
TRealType
antsSCCANObject<TInputImage, TRealType>
::RunSCCAN3( ) 
{
  unsigned int nc1=this->m_MatrixP.rows();
  unsigned int nc2=this->m_MatrixQ.rows();
  unsigned int nc3=this->m_MatrixR.rows();
  if ( nc1 != nc2 || nc1 != nc3 || nc3 != nc2) 
  {
    std::cout<< " P rows " << this->m_MatrixP.rows() << " cols " << this->m_MatrixP.cols() << std::endl;
    std::cout<< " Q rows " << this->m_MatrixQ.rows() << " cols " << this->m_MatrixQ.cols() << std::endl;
    std::cout<< " R rows " << this->m_MatrixR.rows() << " cols " << this->m_MatrixR.cols() << std::endl;
    std::cout<<" N-rows do not match "  << std::endl;
    exit(1);
  }

  this->m_WeightsP=this->InitializeV(this->m_MatrixP);
  this->m_WeightsQ=this->InitializeV(this->m_MatrixQ);
  this->m_WeightsR=this->InitializeV(this->m_MatrixR);
  if ( !this->m_AlreadyWhitened ) {
  this->m_MatrixP=this->NormalizeMatrix(this->m_MatrixP);  
  this->m_MatrixP=this->WhitenMatrixOrGetInverseCovarianceMatrix(this->m_MatrixP);  
  this->m_MatrixQ=this->NormalizeMatrix(this->m_MatrixQ);  
  this->m_MatrixQ=this->WhitenMatrixOrGetInverseCovarianceMatrix(this->m_MatrixQ);  
  this->m_MatrixR=this->NormalizeMatrix(this->m_MatrixR);  
  this->m_MatrixR=this->WhitenMatrixOrGetInverseCovarianceMatrix(this->m_MatrixR);  
  if ( this->m_SpecializationForHBM2011 ) this->SoftThresholdByRMask();
  this->m_AlreadyWhitened=true;
  }
  RealType truecorr=0;
  RealType norm=0,deltacorr=1,lastcorr=1;
  unsigned long its=0;
  while ( its < this->m_MaximumNumberOfIterations && deltacorr > this->m_ConvergenceThreshold  )
  {
  /** for sparse mcca 
   *     w_i \leftarrow \frac{ S( X_i^T ( \sum_{j \ne i} X_j w_j  ) }{norm of above } 
   */
    this->m_WeightsP=this->m_MatrixP.transpose()*(this->m_MatrixQ*this->m_WeightsQ+this->m_MatrixR*this->m_WeightsR);
    this->m_WeightsP=this->SoftThreshold( this->m_WeightsP , this->m_FractionNonZeroP,!this->m_KeepPositiveP);
    norm=this->m_WeightsP.two_norm();
    this->m_WeightsP=this->m_WeightsP/(norm);

    this->m_WeightsQ=this->m_MatrixQ.transpose()*(this->m_MatrixP*this->m_WeightsP+this->m_MatrixR*this->m_WeightsR);
    this->m_WeightsQ=this->SoftThreshold( this->m_WeightsQ , this->m_FractionNonZeroQ,!this->m_KeepPositiveQ);
    norm=this->m_WeightsQ.two_norm();
    this->m_WeightsQ=this->m_WeightsQ/(norm);

    this->m_WeightsR=this->m_MatrixR.transpose()*(this->m_MatrixP*this->m_WeightsP+this->m_MatrixQ*this->m_WeightsQ);
    this->m_WeightsR=this->SoftThreshold( this->m_WeightsR , this->m_FractionNonZeroR,!this->m_KeepPositiveR);
    norm=this->m_WeightsR.two_norm();
    this->m_WeightsR=this->m_WeightsR/(norm);

    VectorType pvec=this->m_MatrixP*this->m_WeightsP;
    VectorType qvec=this->m_MatrixQ*this->m_WeightsQ;
    VectorType rvec=this->m_MatrixR*this->m_WeightsR;

    double corrpq=this->PearsonCorr( pvec , qvec );
    double corrpr=this->PearsonCorr( pvec , rvec );
    double corrqr=this->PearsonCorr( rvec , qvec );
    truecorr=corrpq+corrpr+corrqr;
    deltacorr=fabs(truecorr-lastcorr);
    lastcorr=truecorr;
   // std::cout << " correlation of projections: pq " << corrpq << " pr " << corrpr << " qr " << corrqr << " at-it " << its << std::endl;
    its++;

  }
  //  std::cout << " PNZ-Frac " << this->CountNonZero(this->m_WeightsP) << std::endl;
  //  std::cout << " QNZ-Frac " << this->CountNonZero(this->m_WeightsQ) << std::endl;
  //  std::cout << " RNZ-Frac " << this->CountNonZero(this->m_WeightsR) << std::endl;

  if ( this->m_SpecializationForHBM2011 ) truecorr=this->SpecializedCorrelation();
  this->m_CorrelationForSignificanceTest=truecorr;
  return truecorr;

}

} // namespace ants
} // namespace itk
