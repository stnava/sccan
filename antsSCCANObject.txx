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
#include "antsSCCANObject.h"

namespace itk {
namespace ants {

template <class TInputImage, class TRealType>
typename antsSCCANObject<TInputImage, TRealType>::VectorType
antsSCCANObject<TInputImage, TRealType>
::InitializeV( typename antsSCCANObject<TInputImage, TRealType>::MatrixType p ) 
{
  VectorType w_p( p.columns() );
  vnl_random randgen(time(0));
  for (unsigned long i=0; i < p.columns(); i++)
    { 
      w_p(i)=randgen.drand32();//1.0/p.rows();//
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
::WhitenMatrix( typename antsSCCANObject<TInputImage, TRealType>::MatrixType rin ) 
{

  if (  rin.columns() > rin.rows() ) 
    {
      vnl_svd_economy<RealType> eig(rin*rin.transpose());
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
	    r_diag_inv(j)=1/sqrt(eval);// need sqrt for whitening 
	    r_eigvecs.set_column(j,eig.V().get_column(j));
	  }
	  else r_diag_inv(j)=0;// need sqrt for whitening 
	}
      MatrixType evecs=rin.transpose()*r_eigvecs;
      return (rin*evecs)*(r_diag_inv*evecs.transpose());

    }
  else 
    {     
      vnl_svd_economy<RealType> eig(rin.transpose()*rin);
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
	    r_diag_inv(j)=1/sqrt(eval);// need sqrt for whitening 
	    r_eigvecs.set_column(j,eig.V().get_column(j));
	  }
	  else r_diag_inv(j)=0;// need sqrt for whitening 
	}
      MatrixType wmatrix=r_eigvecs*r_diag_inv*r_eigvecs.transpose();
      return (rin*wmatrix);

    }

}

template <class TInputImage, class TRealType>
typename antsSCCANObject<TInputImage, TRealType>::VectorType
antsSCCANObject<TInputImage, TRealType>
::SoftThreshold( typename antsSCCANObject<TInputImage, TRealType>::VectorType
 v_in, TRealType fractional_goal , bool allow_negative_weights )
{
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
  //    std::cout << " frac non-zero " << frac << std::endl;
  return v_out;
}


template <class TInputImage, class TRealType>
typename antsSCCANObject<TInputImage, TRealType>::VectorType
antsSCCANObject<TInputImage, TRealType>
::TrueCCAPowerUpdate( TRealType penalty1,  typename antsSCCANObject<TInputImage, TRealType>::MatrixType p , typename antsSCCANObject<TInputImage, TRealType>::VectorType  w_q ,  typename antsSCCANObject<TInputImage, TRealType>::MatrixType q, bool keep_pos )
{
  RealType norm=0;
  // inverse covar is symmetric but we transpose anyway for clarity
  // we bracket the computation and use associativity to make sure its done efficiently 
  //vVector wpnew=( (CppInv.transpose()*p.transpose())*(q*CqqInv) )*w_q;
  VectorType wpnew=p.transpose()*(q*w_q);
  wpnew=this->SoftThreshold( wpnew , penalty1 , !keep_pos );
  norm=wpnew.two_norm();
  if ( norm > 0 )
    return wpnew/(norm);
  else return wpnew;
}

template <class TInputImage, class TRealType>
TRealType
antsSCCANObject<TInputImage, TRealType>
::RunSCCAN2( ) 
{
  unsigned int nc1=this->m_MatrixP.rows();
  unsigned int nc2=this->m_MatrixQ.rows();
  if ( nc1 != nc2 ) 
  {
    std::cout<< " P rows " << this->m_MatrixP.rows() << " cols " << this->m_MatrixP.cols() << std::endl;
    std::cout<< " Q rows " << this->m_MatrixQ.rows() << " cols " << this->m_MatrixQ.cols() << std::endl;
    std::cout<< " R rows " << this->m_MatrixR.rows() << " cols " << this->m_MatrixR.cols() << std::endl;
    std::cout<<" N-rows for MatrixP does not equal N-rows for MatrixQ " << nc1 << " vs " << nc2 << std::endl;
    exit(1);
  }
  else {
  std::cout << " P-positivity constraints? " <<  this->m_KeepPositiveP << " frac " << this->m_FractionNonZeroP << " Q-positivity constraints?  " << m_KeepPositiveQ << " frac " << this->m_FractionNonZeroQ << std::endl;
  }
  this->m_WeightsP=this->InitializeV(this->m_MatrixP);
  this->m_WeightsQ=this->InitializeV(this->m_MatrixQ);
  this->m_MatrixP=this->NormalizeMatrix(this->m_MatrixP);  
  this->m_MatrixQ=this->NormalizeMatrix(this->m_MatrixQ);  
  this->m_MatrixP=this->WhitenMatrix(this->m_MatrixP);  
  this->m_MatrixQ=this->WhitenMatrix(this->m_MatrixQ);  
  RealType truecorr=0;
  //  for (unsigned int loop=0; loop<maxccaits; loop++) 
  double deltacorr=1,lastcorr=1;
  unsigned long its=0;
  while ( its < this->m_MaximumNumberOfIterations && deltacorr > this->m_ConvergenceThreshold  )
  {
    this->m_WeightsP=this->TrueCCAPowerUpdate(this->m_FractionNonZeroP,this->m_MatrixP,this->m_WeightsQ,this->m_MatrixQ,this->m_KeepPositiveP);
    this->m_WeightsQ=this->TrueCCAPowerUpdate(this->m_FractionNonZeroQ,this->m_MatrixQ,this->m_WeightsP,this->m_MatrixP,this->m_KeepPositiveQ);
    truecorr=this->PearsonCorr( this->m_MatrixP*this->m_WeightsP , this->m_MatrixQ*this->m_WeightsQ );
    deltacorr=fabs(truecorr-lastcorr);
    lastcorr=truecorr;
    ++its;
  }

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
  this->m_MatrixP=this->NormalizeMatrix(this->m_MatrixP);  
  this->m_MatrixQ=this->NormalizeMatrix(this->m_MatrixQ);  
  this->m_MatrixR=this->NormalizeMatrix(this->m_MatrixR);  
  this->m_MatrixP=this->WhitenMatrix(this->m_MatrixP);  
  this->m_MatrixQ=this->WhitenMatrix(this->m_MatrixQ);  
  this->m_MatrixR=this->WhitenMatrix(this->m_MatrixR);  
  RealType truecorr=0;
  RealType norm=0,deltacorr=1,lastcorr=1;
  unsigned long its=0;
  while ( its < this->m_MaximumNumberOfIterations && deltacorr > this->m_ConvergenceThreshold  )
  {
  /** for sparse mcca 
   *     w_i \leftarrow \frac{ S( X_i^T ( \sum_{j \ne i} X_j w_j  ) }{norm of above } 
   */
    this->m_WeightsP=this->m_MatrixP.transpose()*(this->m_MatrixQ*this->m_WeightsQ+this->m_MatrixR*this->m_WeightsR);
    this->m_WeightsP=this->SoftThreshold( this->m_WeightsP , this->m_FractionNonZeroP,this->m_KeepPositiveP);
    norm=this->m_WeightsP.two_norm();
    this->m_WeightsP=this->m_WeightsP/(norm);

    this->m_WeightsQ=this->m_MatrixQ.transpose()*(this->m_MatrixP*this->m_WeightsP+this->m_MatrixR*this->m_WeightsR);
    this->m_WeightsQ=this->SoftThreshold( this->m_WeightsQ , this->m_FractionNonZeroQ,this->m_KeepPositiveQ);
    norm=this->m_WeightsQ.two_norm();
    this->m_WeightsQ=this->m_WeightsQ/(norm);

    this->m_WeightsR=this->m_MatrixR.transpose()*(this->m_MatrixP*this->m_WeightsP+this->m_MatrixQ*this->m_WeightsQ);
    this->m_WeightsR=this->SoftThreshold( this->m_WeightsR , this->m_FractionNonZeroR,this->m_KeepPositiveR);
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
    std::cout << " correlation of projections: pq " << corrpq << " pr " << corrpr << " qr " << corrqr << " at-it " << its << std::endl;
    its++;

  }

  return truecorr;

}

} // namespace ants
} // namespace itk
