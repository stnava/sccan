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
#include <vnl/algo/vnl_generalized_eigensystem.h>
#include "antsSCCANObject.h"

namespace itk {
namespace ants {

template <class TInputImage, class TRealType>
antsSCCANObject<TInputImage, TRealType>::antsSCCANObject( ) 
{
  this->m_Debug=false;
  this->m_CorrelationForSignificanceTest=0;
  this->m_SpecializationForHBM2011=false;
  this->m_AlreadyWhitened=false;
  this->m_PinvTolerance=1.e-2;
  this->m_PercentVarianceForPseudoInverse=0.9;
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
  w_p.fill(0);
  for (unsigned int its=0; its<1; its++) {
  vnl_random randgen(time(0));
  for (unsigned long i=0; i < p.columns(); i++)
    { 
//      w_p(i)+=randgen.normal();
//      w_p(i)=randgen.drand32();//1.0/p.rows();//
      w_p(i)=1.0/p.rows(); //+randgen.normal()*1.e-3;
//1.0+fabs(randgen.drand32())*1.e-3; 
    }
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
      if ( this->m_Debug ) 
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
::VNLPseudoInverse( typename antsSCCANObject<TInputImage, TRealType>::MatrixType rin , bool take_sqrt ) {
    double pinvTolerance=this->m_PinvTolerance;
      MatrixType dd=rin;
      unsigned int ss=dd.rows();
      if ( dd.rows() > dd.columns() ) ss=dd.columns();
      vnl_svd<RealType> eig(dd,pinvTolerance);
      for (unsigned int j=0; j<ss; j++) 
	{
	  RealType eval=eig.W(j,j);
	  if ( eval > pinvTolerance )  {// FIXME -- check tolerances against matlab pinv
	    eig.W(j,j)=1/(eval);// need eval for inv cov 
	    if ( take_sqrt ) eig.W(j,j)=1/sqrt(eval);
	  }
	  else eig.W(j,j)=0; 
	} 
  return ( eig.recompose() ).transpose();
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
  if ( v_out.two_norm() > 0 ) return v_out/v_out.two_norm();
  return v_out;
}


template <class TInputImage, class TRealType>
typename antsSCCANObject<TInputImage, TRealType>::VectorType
antsSCCANObject<TInputImage, TRealType>
::TrueCCAPowerUpdate( TRealType penalty1,  typename antsSCCANObject<TInputImage, TRealType>::MatrixType p , typename antsSCCANObject<TInputImage, TRealType>::VectorType  w_q ,  typename antsSCCANObject<TInputImage, TRealType>::MatrixType q, bool keep_pos , bool factorOutR )
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
void antsSCCANObject<TInputImage, TRealType>
::RunDiagnostics( unsigned int n_vecs ) 
{ 
  std::cout << "Quantitative diagnostics: "<<std::endl;
  std::cout << "Type 1: correlation from canonical variate to confounding vector "<<std::endl;
  std::cout << "Type 2: correlation from canonical variate to canonical variate "<<std::endl;
  RealType corrthresh=0.3;
  MatrixType omatP=this->NormalizeMatrix(this->m_OriginalMatrixP);
  MatrixType omatQ=this->NormalizeMatrix(this->m_OriginalMatrixQ);
  if (this->m_OriginalMatrixR.size()>0){
   for (unsigned int wv=0; wv<n_vecs; wv++)
     for (unsigned int col=0; col<this->m_MatrixR.columns(); col++) 
      {
      RealType a=this->PearsonCorr(omatP*this->m_VariatesP.get_column(wv) , this->m_OriginalMatrixR.get_column(col) ); 
      RealType b=this->PearsonCorr(omatQ*this->m_VariatesQ.get_column(wv) , this->m_OriginalMatrixR.get_column(col) );
      std::cout << "Pvec " << wv << " confound " << col << " : " << a <<std::endl; 
      std::cout << "Qvec " << wv << " confound " << col << " : " << b <<std::endl; 
      if ( fabs(a) > corrthresh && fabs(b) > corrthresh ) {
       std::cout << " correlation with confound too high for variate " << wv << " corrs " << a << " and " << b <<  std::endl;        
  //     this->m_CanonicalCorrelations[wv]=0;
      }
      }
  }
  for (unsigned int wv=0; wv<n_vecs; wv++)
    for (unsigned int yv=wv+1; yv<n_vecs; yv++)
      {
      RealType a=this->PearsonCorr(omatP*this->m_VariatesP.get_column(wv) ,omatP*this->m_VariatesP.get_column(yv)); 
      if ( fabs(a) > corrthresh ) {
      std::cout << " not orthogonal p " <<  a << std::endl; 
 //       this->m_CanonicalCorrelations[yv]=0; 
      }
      RealType b=this->PearsonCorr(omatQ*this->m_VariatesQ.get_column(wv) ,omatQ*this->m_VariatesQ.get_column(yv));
      if ( fabs(b) > corrthresh ) { 
      std::cout << " not orthogonal q " <<  a << std::endl; 
   //     this->m_CanonicalCorrelations[yv]=0; 
      }
      std::cout << "Pvec " << wv << " Pvec " << yv << " : " << a <<std::endl; 
      std::cout << "Qvec " << wv << " Qvec " << yv << " : " << b <<std::endl; 
      }
 
}


struct my_sccan_sort_class {
  bool operator() (double i, double j) { return (i>j);}
} my_sccan_sort_object;

template <class TInputImage, class TRealType>
TRealType antsSCCANObject<TInputImage, TRealType>
::SparseCCA(unsigned int nvecs)
{
  std::cout <<" ed sparse cca " << std::endl;
  unsigned int nsubj=this->m_MatrixP.rows();
  this->m_MatrixP=this->NormalizeMatrix(this->m_OriginalMatrixP);  
  this->m_MatrixQ=this->NormalizeMatrix(this->m_OriginalMatrixQ);

  RealType tau=0.01;
  if (this->m_Debug) std::cout << " inv view mats " << std::endl;
  MatrixType inviewcovmatP=( (this->m_MatrixP*this->m_MatrixP.transpose())*(this->m_MatrixP*this->m_MatrixP.transpose()) )*(1-tau)+  ( this->m_MatrixP*this->m_MatrixP.transpose() )*tau*(RealType)nsubj;
  MatrixType inviewcovmatQ=( (this->m_MatrixQ*this->m_MatrixQ.transpose())*(this->m_MatrixQ*this->m_MatrixQ.transpose()) )*(1-tau)+( this->m_MatrixQ*this->m_MatrixQ.transpose() )*tau*(RealType)nsubj;

/** standard cca */
//  MatrixType CppInv=this->PseudoInverseCovMat(this->m_MatrixP);
//  MatrixType CqqInv=this->PseudoInverseCovMat(this->m_MatrixQ);
//  MatrixType cov=(CppInv*this->m_MatrixP).transpose()*(CqqInv*this->m_MatrixQ);
//  MatrixType TT=(this->m_MatrixP.transpose()*this->m_MatrixQ)*(this->m_MatrixP.transpose()*this->m_MatrixQ).transpose();
//  MatrixType TT=( (CppInv*this->m_MatrixP).transpose()*(CqqInv*this->m_MatrixQ))*
//                         (this->m_MatrixP.transpose()*this->m_MatrixQ).transpose();

/** dual cca */
  MatrixType CppInv=this->PseudoInverse(inviewcovmatP); 
  MatrixType CqqInv=this->PseudoInverse(inviewcovmatQ);
 
  /** need the eigenvectors of this reduced matrix */ 
  //  eMatrix pccain=CppInv*Cpq*CqqInv*Cqp;
  MatrixType ccap=( (CppInv*(this->m_MatrixP*this->m_MatrixP.transpose())).transpose() *
                  (CqqInv*(this->m_MatrixQ*this->m_MatrixQ.transpose()))*
                        ( (this->m_MatrixP*this->m_MatrixP.transpose() ).transpose()*
                          (this->m_MatrixQ*this->m_MatrixQ.transpose() ) ).transpose() );
// convert to eigen3 format
  eMatrix pccain=this->mVtoE(ccap);
 
  typedef Eigen::EigenSolver<eMatrix> eigsolver;
  eigsolver pEG( pccain );
//  eigsolver qEG( qccain );
  eMatrix pccaVecs = pEG.pseudoEigenvectors();
  eMatrix pccaSquaredCorrs=pEG.pseudoEigenvalueMatrix();

// call this function to check we are doing conversions correctly , matrix-wise 
  this->mEtoV(pccaVecs);

// copy to stl vector so we can sort the results  
  std::vector<TRealType> evals(pccaSquaredCorrs.cols(),-9.e9);
  for ( long j=0; j<pccaSquaredCorrs.cols(); ++j){
    RealType val=pccaSquaredCorrs(j,j);
    evals[j]=val;
  }

// sort and reindex the eigenvectors/values 
  sort (evals.begin(), evals.end(), my_sccan_sort_object); 
  std::vector<int> sorted_indices(nvecs,-1);
  for (unsigned int i=0; i<evals.size(); i++) {
  for (unsigned int j=0; j<evals.size(); j++) {
    if ( evals[i] == pccaSquaredCorrs(j,j) &&  sorted_indices[i] == -1 ) {
      sorted_indices[i]=j;
      pccaSquaredCorrs(j,j)=0;
    }
  }}
  this->m_CanonicalCorrelations.set_size(nvecs);
  this->m_CanonicalCorrelations.fill(0); 
// map the variates back to P, Q space
  this->m_VariatesP.set_size(this->m_MatrixP.cols(),nvecs);
  this->m_VariatesQ.set_size(this->m_MatrixQ.cols(),nvecs);
  for (unsigned int i=0; i<nvecs; i++) {
    VectorType temp=this->vEtoV( pccaVecs.col(  sorted_indices[i] ) );
    VectorType tempq;
    tempq=( CqqInv*( (this->m_MatrixQ*this->m_MatrixQ.transpose() )* (this->m_MatrixP*this->m_MatrixP.transpose() ) ))*temp;
    VectorType pvar=this->SoftThreshold(  temp*this->m_MatrixP , this->m_FractionNonZeroP , !this->m_KeepPositiveP ); 
    VectorType qvar=this->SoftThreshold( tempq*this->m_MatrixQ , this->m_FractionNonZeroQ , !this->m_KeepPositiveQ );
    this->m_VariatesP.set_column( i, pvar  );
    this->m_VariatesQ.set_column( i, qvar  );
  }

for (unsigned int i=0; i<nvecs; i++) {
    this->m_CanonicalCorrelations[i]=
      this->PearsonCorr(this->m_MatrixP*this->GetVariateP(i),this->m_MatrixQ*this->GetVariateQ(i) );
    std::cout << "correlation of mapped back data " << this->m_CanonicalCorrelations[i] <<  " eval " << evals[i] << std::endl;
  }
  for (unsigned int i=0; i<nvecs-1; i++) {
    std::cout << "inner prod of projections " <<  this->PearsonCorr( this->m_MatrixP*this->GetVariateP(i) ,  this->m_MatrixP*this->GetVariateP(i+1) ) << std::endl;
  }

  return this->m_CanonicalCorrelations[0];

}


template <class TInputImage, class TRealType>
TRealType antsSCCANObject<TInputImage, TRealType>
::SparsePartialCCA(unsigned int nvecs)
{
  std::cout <<" ed sparse partial cca " << std::endl;
  unsigned int nsubj=this->m_MatrixP.rows();
  this->m_MatrixP=this->NormalizeMatrix(this->m_OriginalMatrixP);  
  this->m_MatrixQ=this->NormalizeMatrix(this->m_OriginalMatrixQ);
  this->m_MatrixR=this->NormalizeMatrix(this->m_OriginalMatrixR);
  
  RealType tau=0.01;
  if (this->m_Debug) std::cout << " inv view mats " << std::endl;
  this->m_MatrixRRt=this->ProjectionMatrix(this->m_OriginalMatrixR);
  MatrixType PslashR=this->m_MatrixP-(this->m_MatrixRRt*this->m_MatrixP);
  MatrixType QslashR=this->m_MatrixQ-this->m_MatrixRRt*this->m_MatrixQ;
  if (this->m_Debug) {
    std::cout <<" corr-pre " << this->PearsonCorr( this->m_MatrixP.get_column(0) , this->m_OriginalMatrixR.get_column(0)  ) << std::endl  ; std::cout <<" corr-post " << this->PearsonCorr( PslashR.get_column(0) , this->m_OriginalMatrixR.get_column(0)  ) << std::endl; 
    std::cout <<" corr-pre " << this->PearsonCorr( this->m_MatrixQ.get_column(0) , this->m_OriginalMatrixR.get_column(0)  ) << std::endl  ; std::cout <<" corr-post " << this->PearsonCorr( QslashR.get_column(0) , this->m_OriginalMatrixR.get_column(0)  ) << std::endl;
  }
  MatrixType inviewcovmatP=( (PslashR*PslashR.transpose())*(PslashR*PslashR.transpose()) )*(1-tau)+  ( PslashR*PslashR.transpose() )*tau*(RealType)nsubj;
  MatrixType inviewcovmatQ=( (QslashR*QslashR.transpose())*(QslashR*QslashR.transpose()) )*(1-tau)+( QslashR*QslashR.transpose() )*tau*(RealType)nsubj;

/** dual cca */
  if (this->m_Debug) std::cout << " inv view mats pinv " << std::endl;
  MatrixType CppInv=this->PseudoInverse(inviewcovmatP); 
  MatrixType CqqInv=this->PseudoInverse(inviewcovmatQ);
  if (this->m_Debug) std::cout << " inv view mats pinv done " << std::endl;
 
  /** need the eigenvectors of this reduced matrix */ 
  MatrixType ccap=( (CppInv*(PslashR*PslashR.transpose())).transpose() *
                  (CqqInv*(QslashR*QslashR.transpose()))*
                        ( (PslashR*PslashR.transpose() ).transpose()*
                          (QslashR*QslashR.transpose() ) ).transpose() );
// convert to eigen3 format
  eMatrix pccain=this->mVtoE(ccap);
 
  typedef Eigen::EigenSolver<eMatrix> eigsolver;
  eigsolver pEG( pccain );
  eMatrix pccaVecs = pEG.pseudoEigenvectors();
  eMatrix pccaSquaredCorrs=pEG.pseudoEigenvalueMatrix();
// call this function to check we are doing conversions correctly , matrix-wise 
  this->mEtoV(pccaVecs);
// map the variates back to P, Q space and sort them 
  this->m_CanonicalCorrelations.set_size(nvecs);
  this->m_CanonicalCorrelations.fill(0); 
// copy to stl vector so we can sort the results  
  std::vector<TRealType> evals(pccaSquaredCorrs.cols(),0);
  std::vector<TRealType> oevals(pccaSquaredCorrs.cols(),0);
  for ( long j=0; j<pccaSquaredCorrs.cols(); ++j){
    RealType val=pccaSquaredCorrs(j,j);
    evals[j]=val;
    oevals[j]=val;
    if ( val > 0 ){
      VectorType temp=this->vEtoV( pccaVecs.col(  j ) );
      VectorType tempq=( CqqInv*( (QslashR*QslashR.transpose() )* (PslashR*PslashR.transpose() ) ))*temp;
      VectorType pvar=this->SoftThreshold(  temp*PslashR , this->m_FractionNonZeroP , !this->m_KeepPositiveP ); 
      VectorType qvar=this->SoftThreshold( tempq*QslashR , this->m_FractionNonZeroQ , !this->m_KeepPositiveQ );
      evals[j]=this->PearsonCorr(PslashR*pvar,QslashR*qvar);
      oevals[j]=evals[j];
    }
  }

// sort and reindex the eigenvectors/values 
  sort (evals.begin(), evals.end(), my_sccan_sort_object); 
  std::vector<int> sorted_indices(nvecs,-1);
  for (unsigned int i=0; i<evals.size(); i++) {
  for (unsigned int j=0; j<evals.size(); j++) {
    if ( evals[i] == oevals[j] &&  sorted_indices[i] == -1 ) {
      sorted_indices[i]=j;
      oevals[j]=0;
    }
  }}

  this->m_VariatesP.set_size(PslashR.cols(),nvecs);
  this->m_VariatesQ.set_size(QslashR.cols(),nvecs);
  for (unsigned int i=0; i<nvecs; i++) {
    VectorType temp=this->vEtoV( pccaVecs.col(  sorted_indices[i] ) );
    VectorType tempq=( CqqInv*( (QslashR*QslashR.transpose() )* (PslashR*PslashR.transpose() ) ))*temp;
    VectorType pvar=this->SoftThreshold(  temp*PslashR , this->m_FractionNonZeroP , !this->m_KeepPositiveP ); 
    VectorType qvar=this->SoftThreshold( tempq*QslashR , this->m_FractionNonZeroQ , !this->m_KeepPositiveQ );
    this->m_VariatesP.set_column( i, pvar  );
    this->m_VariatesQ.set_column( i, qvar  );
  }

  for (unsigned int i=0; i<nvecs; i++) {
    this->m_CanonicalCorrelations[i]=
      this->PearsonCorr(PslashR*this->GetVariateP(i),QslashR*this->GetVariateQ(i) );
    std::cout << "correlation of mapped back data " << this->m_CanonicalCorrelations[i] <<  " eval " << evals[i] << std::endl;
  }
  for (unsigned int i=0; i<nvecs-1; i++) {
    std::cout << "inner prod of projections " <<  this->PearsonCorr( PslashR*this->GetVariateP(i) ,  PslashR*this->GetVariateP(i+1) ) << std::endl;
  }

  return this->m_CanonicalCorrelations[0];

}



template <class TInputImage, class TRealType>
void antsSCCANObject<TInputImage, TRealType>
::WhitenDataSetForRunSCCANMultiple(unsigned int nvecs)
{
    if ( this->m_Debug ) std::cout << " now whiten and apply R " << std::endl;
    if ( this->m_OriginalMatrixR.size() > 0 || nvecs > 0  ) {
      this->m_MatrixP=this->NormalizeMatrix(this->m_OriginalMatrixP);  
      if ( this->m_VariatesP.size() > 0 ) {
        this->m_MatrixRp.set_size(this->m_MatrixP.rows(),this->m_OriginalMatrixR.cols()+nvecs);
	this->m_MatrixRp.set_columns(0,this->m_OriginalMatrixR);
	this->m_MatrixRp.set_columns(this->m_OriginalMatrixR.cols(), 
          this->m_MatrixP*(this->m_VariatesP.get_n_columns(0,nvecs)));
      } 
      else {
        this->m_MatrixRp=this->NormalizeMatrix(this->m_OriginalMatrixR);  
      }
      MatrixType temp;
      this->m_MatrixRp=ProjectionMatrix(this->m_MatrixRp);
      temp=this->m_MatrixP-this->m_MatrixRp*this->m_MatrixP;
      this->m_MatrixP=temp;
//      this->m_MatrixP=this->WhitenMatrixByAnotherMatrix(this->m_MatrixP,temp);
  if (this->m_Debug && this->m_OriginalMatrixR.size()>0 ) {
    std::cout <<" corr-pre " << this->PearsonCorr( this->m_MatrixP.get_column(0) , this->m_OriginalMatrixR.get_column(0)  ) << std::endl  ; std::cout <<" corr-post " << this->PearsonCorr( temp.get_column(0) , this->m_OriginalMatrixR.get_column(0)  ) << std::endl; 
  }
  if (this->m_Debug && this->m_VariatesP.cols() > 1 ) {
    std::cout <<" corr-pre " << this->PearsonCorr( (this->m_OriginalMatrixP*this->m_VariatesP).get_column(0) , (this->m_OriginalMatrixP*this->m_VariatesP).get_column(1) ) << std::endl  ; std::cout <<" corr-post " << this->PearsonCorr( (this->m_MatrixP*this->m_VariatesP).get_column(0) , (this->m_MatrixP*this->m_VariatesP).get_column(1) ) << std::endl; 
  }
  

      this->m_MatrixQ=this->NormalizeMatrix(this->m_OriginalMatrixQ);  
      if ( this->m_VariatesQ.size() > 0 ) {
        this->m_MatrixRq.set_size(this->m_MatrixQ.rows(),this->m_OriginalMatrixR.cols()+nvecs);
	this->m_MatrixRq.set_columns(0,this->m_OriginalMatrixR);
	this->m_MatrixRq.set_columns(this->m_OriginalMatrixR.cols(),
          this->m_MatrixQ*(this->m_VariatesQ.get_n_columns(0,nvecs)));
      } 
      else {
        this->m_MatrixRq=this->NormalizeMatrix(this->m_OriginalMatrixR); 
      }
      this->m_MatrixRq=ProjectionMatrix(this->m_MatrixRq);
      temp=this->m_MatrixQ-this->m_MatrixRq*this->m_MatrixQ;
      this->m_MatrixQ=temp;
//      this->m_MatrixQ=this->WhitenMatrixByAnotherMatrix( this->m_MatrixQ,temp);
    }
    else {
      this->m_MatrixP=this->NormalizeMatrix(this->m_OriginalMatrixP);  
      this->m_MatrixQ=this->NormalizeMatrix(this->m_OriginalMatrixQ);
      this->m_MatrixP=this->WhitenMatrix(this->m_MatrixP);  
      this->m_MatrixQ=this->WhitenMatrix(this->m_MatrixQ);
    }


    if ( this->m_Debug ) std::cout << "  whiten and apply R done " << std::endl;
}


template <class TInputImage, class TRealType>
TRealType
antsSCCANObject<TInputImage, TRealType>
::RunSCCAN2multiple( unsigned int n_vecs ) 
{
  std::cout << " power iteration (partial) scca " <<std::endl;
  this->m_CanonicalCorrelations.set_size(n_vecs);
  this->m_CanonicalCorrelations.fill(0); 
  RealType truecorr=0;
  unsigned int nr1=this->m_MatrixP.rows();
  unsigned int nr2=this->m_MatrixQ.rows();
  this->m_VariatesP.set_size(0,0);
  this->m_VariatesQ.set_size(0,0);
  if ( nr1 != nr2 ) 
  {
    std::cout<< " P rows " << this->m_MatrixP.rows() << " cols " << this->m_MatrixP.cols() << std::endl;
    std::cout<< " Q rows " << this->m_MatrixQ.rows() << " cols " << this->m_MatrixQ.cols() << std::endl;
    std::cout<< " R rows " << this->m_MatrixR.rows() << " cols " << this->m_MatrixR.cols() << std::endl;
    std::cout<<" N-rows for MatrixP does not equal N-rows for MatrixQ " << nr1 << " vs " << nr2 << std::endl;
    exit(1);
  }
  if  (  !this->m_AlreadyWhitened  ) {
    if ( this->m_Debug ) std::cout << " whiten " << std::endl;
    this->WhitenDataSetForRunSCCANMultiple();    
    this->m_AlreadyWhitened=true; 
    if ( this->m_Debug ) std::cout << " whiten done " << std::endl;
  } 
  this->m_VariatesP.set_size(this->m_MatrixP.cols(),n_vecs);
  this->m_VariatesQ.set_size(this->m_MatrixQ.cols(),n_vecs);
  for (unsigned int kk=0;kk<n_vecs; kk++) {
    this->m_VariatesP.set_column(kk,this->InitializeV(this->m_MatrixP));
    this->m_VariatesQ.set_column(kk,this->InitializeV(this->m_MatrixQ));
  }
  MatrixType p_evecs_factor;
  MatrixType q_evecs_factor;
  unsigned int which_e_vec=0; 
  bool notdone=true;
  while ( notdone ) {
     if ( this->m_Debug ) std::cout << " get canonical variate number " << which_e_vec+1 << std::endl;
    double initcorr=1.e-5;
    truecorr=initcorr;
    double deltacorr=1,lastcorr=initcorr*0.5;
    this->m_WeightsP= this->m_VariatesP.get_column(which_e_vec);
    this->m_WeightsQ= this->m_VariatesQ.get_column(which_e_vec);
    unsigned long its=0, min_its=3;
    if ( this->m_Debug ) std::cout << " Begin " << std::endl;
    while ( its < this->m_MaximumNumberOfIterations && deltacorr > this->m_ConvergenceThreshold || its < min_its )
    {
      if ( its == 0 &&  which_e_vec > 0) this->WhitenDataSetForRunSCCANMultiple(which_e_vec);

      {
        VectorType proj=this->m_MatrixQ*this->m_WeightsQ;
        if ( this->m_MatrixRq.size() > 0 ) proj=proj-this->m_MatrixRq*proj;
        for (unsigned int kk=0; kk<which_e_vec; kk++) 
          proj=this->Orthogonalize(proj,this->m_OriginalMatrixQ*this->m_VariatesQ.get_column(kk));
        this->m_WeightsP=this->m_MatrixP.transpose()*(proj);
        this->m_WeightsP=this->SoftThreshold( this->m_WeightsP , this->m_FractionNonZeroP , !this->m_KeepPositiveP );
        this->m_VariatesP.set_column(which_e_vec,this->m_WeightsP);
      }

      {
        VectorType proj=this->m_MatrixP*this->m_WeightsP;
        if ( this->m_MatrixRp.size() > 0 ) proj=proj-this->m_MatrixRp*proj;
        for (unsigned int kk=0; kk<which_e_vec; kk++) 
          proj=this->Orthogonalize(proj,this->m_OriginalMatrixP*this->m_VariatesP.get_column(kk));
        this->m_WeightsQ=this->m_MatrixQ.transpose()*(proj);
        this->m_WeightsQ=this->SoftThreshold( this->m_WeightsQ , this->m_FractionNonZeroQ , !this->m_KeepPositiveQ );
        this->m_VariatesQ.set_column(which_e_vec,this->m_WeightsQ);
      }

      truecorr=this->PearsonCorr( this->m_MatrixP*this->m_WeightsP , this->m_MatrixQ*this->m_WeightsQ );
      if ( this->m_Debug ) std::cout << " corr " << truecorr << std::endl;
      deltacorr=fabs(truecorr-lastcorr);
      lastcorr=truecorr;
      ++its;
    }// inner_it
    this->m_CanonicalCorrelations[which_e_vec]=truecorr; 
    std::cout << "  canonical variate number " << which_e_vec+1 << " corr " << this->m_CanonicalCorrelations[which_e_vec]  << std::endl;
    if ( fabs(truecorr) < 1.e-2 || (which_e_vec+1) == n_vecs ) notdone=false;
    else which_e_vec++;
  }     
  this->RunDiagnostics(n_vecs);
  RealType corrsum=0;
  for ( unsigned int i=0; i < this->m_CanonicalCorrelations.size(); i++) 
    corrsum+=fabs(this->m_CanonicalCorrelations[i]);
  return corrsum;
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
  this->m_WeightsP=this->InitializeV(this->m_MatrixP);
  this->m_WeightsQ=this->InitializeV(this->m_MatrixQ);

//  if ( !this->m_AlreadyWhitened ) 
  {
   if ( this->m_Debug ) std::cout <<" norm P " << std::endl;
    this->m_MatrixP=this->NormalizeMatrix(this->m_MatrixP);  
   if ( this->m_Debug ) std::cout <<" norm Q " << std::endl;
    this->m_MatrixQ=this->NormalizeMatrix(this->m_MatrixQ);  
    if ( this->m_OriginalMatrixR.size() > 0 ) {
      this->m_MatrixR=this->NormalizeMatrix(this->m_OriginalMatrixR);  
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
    this->m_WeightsP=this->TrueCCAPowerUpdate(this->m_FractionNonZeroP,this->m_MatrixP,this->m_WeightsQ,this->m_MatrixQ,this->m_KeepPositiveP,false);
    this->m_WeightsQ=this->TrueCCAPowerUpdate(this->m_FractionNonZeroQ,this->m_MatrixQ,this->m_WeightsP,this->m_MatrixP,this->m_KeepPositiveQ,false);    
    truecorr=this->PearsonCorr( this->m_MatrixP*this->m_WeightsP , this->m_MatrixQ*this->m_WeightsQ );
    deltacorr=fabs(truecorr-lastcorr);
    lastcorr=truecorr;
    ++its;
  }// inner_it 
  }//outer_it 
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
  this->m_MatrixP=this->WhitenMatrix(this->m_MatrixP);  
  this->m_MatrixQ=this->NormalizeMatrix(this->m_MatrixQ);  
  this->m_MatrixQ=this->WhitenMatrix(this->m_MatrixQ);  
  this->m_MatrixR=this->NormalizeMatrix(this->m_MatrixR);  
  this->m_MatrixR=this->WhitenMatrix(this->m_MatrixR);  
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

  this->m_CorrelationForSignificanceTest=truecorr;
  return truecorr;

}

} // namespace ants
} // namespace itk




/*

      if (its == 0 )
      if ( which_e_vec > 0   && false ) 
       {
// here, factor out previous evecs globally 
	MatrixType temp;
	MatrixType pp;
	unsigned int basect=this->m_MatrixR.columns(); 
	basect=0;
	pp.set_size(this->m_MatrixP.rows(),which_e_vec+basect);
//	for (unsigned int kk=0; kk<this->m_MatrixR.columns(); kk++) 
//          pp.set_column(kk,this->m_MatrixR.get_column(kk));
	unsigned int colcount=0; //this->m_MatrixR.columns();
	for (unsigned int kk=which_e_vec-1; kk<which_e_vec; kk++) { 
          pp.set_column(colcount,this->m_MatrixP*this->m_VariatesP.get_column(kk));
	  colcount++;
        }
        temp=this->NormalizeMatrix(pp);  
        temp=this->WhitenMatrix(temp);  
        temp=temp*temp.transpose(); 
        this->m_MatrixP=(this->m_MatrixP-temp*this->m_MatrixP);

	MatrixType qq;
	qq.set_size(this->m_MatrixQ.rows(),which_e_vec+this->m_MatrixR.columns());
	for (unsigned int kk=0; kk<this->m_MatrixR.columns(); kk++) 
          qq.set_column(kk,this->m_MatrixR.get_column(kk));
	colcount=this->m_MatrixR.columns();
	for (unsigned int kk=which_e_vec-1; kk<which_e_vec; kk++) { 
          qq.set_column(colcount,this->m_MatrixQ*this->m_VariatesQ.get_column(kk));
	  colcount++;
        }
        temp=this->NormalizeMatrix(qq);  
        temp=this->WhitenMatrix(temp);  
        temp=temp*temp.transpose(); 
        this->m_MatrixQ=(this->m_MatrixQ-temp*this->m_MatrixQ);
      }




//m_Debug=true; 
      double ip=1; unsigned long ct=0,max_ip_its=50;
      double deltaip=1,lastip=0;
      while ( (deltaip) > 1.e-3 && ct < max_ip_its && which_e_vec > 0 || ct < 4 ) {
        ip=0;
        this->m_WeightsP=this->SoftThreshold( this->m_WeightsP , this->m_FractionNonZeroP , !this->m_KeepPositiveP );
        VectorType ptem=this->m_WeightsP;
	if ( which_e_vec >= 1 ) 
          ptem=this->Orthogonalize(ptem,this->m_VariatesP.get_column(0),&this->m_MatrixP);
  	if ( which_e_vec >= 2 ) 
          ptem=this->Orthogonalize(ptem,this->m_VariatesP.get_column(1),&this->m_MatrixP);
	this->m_WeightsP=ptem;
        this->m_WeightsP=this->SoftThreshold( this->m_WeightsP , this->m_FractionNonZeroP , !this->m_KeepPositiveP );
        ip+=this->PearsonCorr(this->m_MatrixP*this->m_WeightsP,this->m_MatrixP*this->m_VariatesP.get_column(0));
	if ( which_e_vec >= 2) 
          ip+=this->PearsonCorr(this->m_MatrixP*this->m_WeightsP,this->m_MatrixP*this->m_VariatesP.get_column(1));
	deltaip=fabs(lastip)-fabs(ip);
	lastip=ip;
	ct++;
         if ( this->m_Debug ) std::cout << " pip-b " << ip << " delt " << deltaip << std::endl;
        }

       ip=1; ct=0;
       deltaip=1;lastip=0;
      while ( (deltaip) > 1.e-3 && ct < max_ip_its && which_e_vec > 0  || ct < 4 ) {
        ip=0;
        this->m_WeightsQ=this->SoftThreshold( this->m_WeightsQ , this->m_FractionNonZeroQ , !this->m_KeepPositiveQ );
        VectorType ptem=this->m_WeightsQ;
	if ( which_e_vec >= 1 ) 
          ptem=this->Orthogonalize(ptem,this->m_VariatesQ.get_column(0),&this->m_MatrixQ);
  	if ( which_e_vec >= 2 ) 
          ptem=this->Orthogonalize(ptem,this->m_VariatesQ.get_column(1),&this->m_MatrixQ);
	this->m_WeightsQ=ptem;
        this->m_WeightsQ=this->SoftThreshold( this->m_WeightsQ , this->m_FractionNonZeroQ , !this->m_KeepPositiveQ );
        ip+=this->PearsonCorr(this->m_MatrixQ*this->m_WeightsQ,this->m_MatrixQ*this->m_VariatesQ.get_column(0));
	if ( which_e_vec >= 2) 
          ip+=this->PearsonCorr(this->m_MatrixQ*this->m_WeightsQ,this->m_MatrixQ*this->m_VariatesQ.get_column(1));
	deltaip=fabs(lastip)-fabs(ip);
	lastip=ip;
	ct++;
         if ( this->m_Debug ) std::cout << " qip-b " << ip << " delt " << deltaip << std::endl;
        }


// alternative tools for factoring out evecs 
//      if ( its == 0 && which_e_vec > 0  ) {
//        this->WhitenDataSetForRunSCCANMultiple();   
        q_evecs_factor=this->m_MatrixQ*this->m_VariatesQ.get_n_columns(0,which_e_vec);
        q_evecs_factor=this->NormalizeMatrix( q_evecs_factor );  
        q_evecs_factor=this->WhitenMatrix(q_evecs_factor);  
        q_evecs_factor=q_evecs_factor*q_evecs_factor.transpose(); 
 //       MatrixType temp=this->m_MatrixP-q_evecs_factor*this->m_MatrixP;
 //       temp=this->InverseCovarianceMatrix(temp,&this->m_MatrixP);  
 //       this->m_MatrixP=temp;

        p_evecs_factor=this->m_MatrixP*this->m_VariatesP.get_n_columns(0,which_e_vec);
        p_evecs_factor=this->NormalizeMatrix( p_evecs_factor );  
        p_evecs_factor=this->WhitenMatrix(p_evecs_factor);  
        p_evecs_factor=p_evecs_factor*p_evecs_factor.transpose();
//        temp=this->m_MatrixQ-p_evecs_factor*this->m_MatrixQ;
//        temp=this->InverseCovarianceMatrix(temp,&this->m_MatrixQ);  
//        this->m_MatrixQ=temp;

      }

*/