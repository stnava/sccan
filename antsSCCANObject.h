/*=========================================================================

  Program:   Advanced Normalization Tools
  Module:    $RCSfile: antsSCCANObject.h,v $
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
#ifndef __antsSCCANObject_h
#define __antsSCCANObject_h

#include <vnl/algo/vnl_matrix_inverse.h>
#include <vnl/algo/vnl_cholesky.h>
#include "itkImageToImageFilter.h"
/** Custom SCCA implemented with vnl and ITK: Flexible positivity constraints, image ops, permutation testing, etc. */
namespace itk {
namespace ants {

template<class TInputImage, class TRealType = double>
class ITK_EXPORT antsSCCANObject :
    public ImageToImageFilter<TInputImage, TInputImage>
{
public:
  /** Standard class typdedefs. */
  typedef antsSCCANObject                     Self;
  typedef ImageToImageFilter<TInputImage, TInputImage>  Superclass;
  typedef SmartPointer<Self>                                 Pointer;
  typedef SmartPointer<const Self>                           ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Run-time type information (and related methods). */
  itkTypeMacro( antsSCCANObject, ImageToImageFilter );

  /** Dimension of the images. */
  itkStaticConstMacro( ImageDimension, unsigned int,
                       TInputImage::ImageDimension );

  itkStaticConstMacro( MatrixDimension, unsigned int, 2 );

  /** Typedef support of input types. */
  typedef TInputImage                                 ImageType;
  typedef typename ImageType::Pointer                 ImagePointer;
  typedef typename ImageType::PixelType               PixelType;
  typedef typename ImageType::IndexType               IndexType;

  /** Some convenient typedefs. */
  typedef TRealType                                      RealType;
  typedef Image<RealType,
    itkGetStaticConstMacro( ImageDimension )>         RealImageType;

  /** note, eigen for pseudo-eigenvals  */
  typedef vnl_matrix<RealType>        MatrixType;
  typedef vnl_vector<RealType>        VectorType;
  typedef MatrixType                  VariateType;
  typedef vnl_diag_matrix<RealType>   DiagonalMatrixType;

  enum SCCANFormulationType{ PQ , PminusRQ ,  PQminusR ,  PminusRQminusR , PQR  };

  /** ivars Set/Get functionality */
  itkSetMacro( MaximumNumberOfIterations, unsigned int );
  itkGetConstMacro( MaximumNumberOfIterations, unsigned int );
  itkSetMacro( ConvergenceThreshold, RealType );
  itkGetConstMacro( ConvergenceThreshold, RealType );
  itkGetConstMacro( CurrentConvergenceMeasurement, RealType );
  itkGetConstMacro( ElapsedIterations, unsigned int );
  itkSetMacro( SCCANFormulation, SCCANFormulationType );
  itkGetConstMacro( SCCANFormulation, SCCANFormulationType );

  void SetPseudoInversePercentVariance( RealType p ) { this->m_PercentVarianceForPseudoInverse=p; }

  itkSetMacro( FractionNonZeroP, RealType );
  itkSetMacro( KeepPositiveP, bool );
  void SetMaskImageP( ImagePointer mask ) { this->m_MaskImageP=mask; }
  void SetMatrixP(  MatrixType matrix ) { this->m_OriginalMatrixP.set_size(matrix.rows(),matrix.cols());  this->m_MatrixP.set_size(matrix.rows(),matrix.cols()); this->m_OriginalMatrixP.update(matrix); this->m_MatrixP.update(matrix); }
  MatrixType GetMatrixP(  ) { return this->m_MatrixP; }

  MatrixType PseudoInverse( MatrixType rin ) {
    double pinvTolerance=1.e-5;
    double delta=1.e-0; 
    //    vnl_svd<double> mysvd(rin,pinvTolerance); 
    //  return mysvd.pinverse();

    //  see Limit Relations in http://en.wikipedia.org/wiki/Moore%E2%80%93Penrose_pseudoinverse
      // we use this equation 
      // $   A^+ = limit as \delta approaches 0  ( A^T A + delta Id )^{-1} A^T $
      // or if cols < rows 
      // $   A^+ = limit as \delta approaches 0  A^T ( A A^T + delta Id )^{-1}  $
    if (  rin.columns() > rin.rows() ) 
    {
      // first add the identity to 
      MatrixType dd=rin*rin.transpose();
      dd.set_identity(); 
      dd=dd*delta+rin*rin.transpose();
      vnl_svd_economy<RealType> eig(dd);
      VectorType eigvals=eig.lambdas();
      RealType eigsum=eigvals.sum();
      RealType total=0; 
      unsigned int eigct=0;
      while ( total/eigsum < 1 ) 
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
	  if ( eval > pinvTolerance )  {// FIXME -- check tolerances against matlab pinv
	    r_diag_inv(j)=1/(eval);// need eval for inv cov 
	    r_eigvecs.set_column(j,eig.V().get_column(j));
	  }
	  else r_diag_inv(j)=0; 
	}
      // $   A^+ = limit as \delta approaches 0  A^T ( A A^T + delta Id )^{-1}  $
      // $   A^+ \approx A^T ( A A^T + delta Id )^{-1}  $
      // svd approx to inverse is  = V \Sigma^{-1} U^T
      // $   A^+ \approx A^T ( V r_diag_inv V^T )   $
      MatrixType wmatrix= rin.transpose()* ( r_eigvecs*r_diag_inv*(r_eigvecs.transpose() ));
      // std::cout << "row rank-deficient wm rows " << wmatrix.rows() << " cols " << wmatrix.cols() <<std::endl;
      // std::cout << wmatrix * rin  << std::endl;
      return wmatrix;

    }
  else 
    {     
      MatrixType dd=rin.transpose()*rin;
      dd.set_identity(); 
      dd=dd*delta+rin.transpose()*rin;
      vnl_svd_economy<RealType> eig(dd);
      VectorType eigvals=eig.lambdas();
      RealType eigsum=eigvals.sum();
      RealType total=0; 
      unsigned int eigct=0;
      while ( total/eigsum < 1 ) 
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
	  if ( eval > pinvTolerance )  {// FIXME -- check tolerances against matlab pinv
	    r_diag_inv(j)=1/(eval);// need eval for inv cov 
	    r_eigvecs.set_column(j,eig.V().get_column(j));
	  }
	  else r_diag_inv(j)=0; 
	}
      // $   A^+ = limit as \delta approaches 0  ( A^T A + delta Id )^{-1} A^T $
      MatrixType wmatrix=(  r_eigvecs*r_diag_inv*(r_eigvecs.transpose() ) ) * rin.transpose();
      //      std::cout << "full rank wm rows " << wmatrix.rows() << " cols " << wmatrix.cols() <<std::endl;
      //  std::cout << wmatrix * rin  << std::endl;
      return wmatrix;
    }


  }

  VectorType OrthogonalizeVector(MatrixType P, MatrixType varP , unsigned int col1, unsigned int col2 ,  RealType fnz, bool keeppos )
  {
    MatrixType tt=(P*varP).transpose();
    double ratio=inner_product(tt.get_row(col2),tt.get_row(col1))/inner_product(tt.get_row(col1),tt.get_row(col1));
    VectorType  ortho=tt.get_row(col2)-tt.get_row(col1)*ratio;
    VectorType x = this->PseudoInverse(P)*ortho; 
    //    x=this->SoftThreshold( x , fnz , !keeppos);
    return x;
  }

  VectorType Orthogonalize(VectorType Mvec, VectorType V )
  {
      double ratio=inner_product(Mvec,V)/inner_product(V,V);
      VectorType  ortho=Mvec-V*ratio;
    return ortho;
  }

  MatrixType OrthogonalizeMatrix(MatrixType M, VectorType V )
  {
    for ( unsigned int j = 0 ; j < M.cols() ; j++ ) {
      VectorType Mvec=M.get_column(j);
      double ratio=inner_product(Mvec,V)/inner_product(V,V);
      VectorType  ortho=Mvec-V*ratio;
      M.set_column(j,ortho);
    }
    return M;
  }
  
  itkSetMacro( FractionNonZeroQ, RealType );
  itkSetMacro( KeepPositiveQ, bool );
  void SetMaskImageQ( ImagePointer mask ) { this->m_MaskImageQ=mask; }
  void SetMatrixQ(  MatrixType  matrix ) {  this->m_OriginalMatrixQ.set_size(matrix.rows(),matrix.cols());  this->m_MatrixQ.set_size(matrix.rows(),matrix.cols()); this->m_OriginalMatrixQ.update(matrix); this->m_MatrixQ.update(matrix);}
  MatrixType GetMatrixQ(  ) { return this->m_MatrixQ; }

  itkSetMacro( FractionNonZeroR, RealType );
  itkSetMacro( KeepPositiveR, bool );
  void SetMaskImageR( ImagePointer mask ) { this->m_MaskImageR=mask; }
  void SetMatrixR(  MatrixType matrix ) {  this->m_OriginalMatrixR.set_size(matrix.rows(),matrix.cols());  this->m_MatrixR.set_size(matrix.rows(),matrix.cols()); this->m_OriginalMatrixR.update(matrix); this->m_MatrixR.update(matrix); }
  MatrixType GetMatrixR(  ) { return this->m_MatrixR; }

  RealType RunSCCAN2multiple( unsigned int n_vecs );
  RealType RunSCCAN2( );
  RealType RunSCCAN3();
 
  VectorType SoftThreshold( VectorType v_in, RealType fractional_goal , bool allow_negative_weights );
  VectorType InitializeV( MatrixType p );
  MatrixType NormalizeMatrix(MatrixType p);
  /** needed for partial scca */
  MatrixType WhitenMatrixOrGetInverseCovarianceMatrix(MatrixType p , bool white_else_invcov=true ); 
  MatrixType InverseCovarianceMatrix(MatrixType p) { return this->WhitenMatrixOrGetInverseCovarianceMatrix(p, false); }
  MatrixType WhitenMatrix(MatrixType p) { return this->WhitenMatrixOrGetInverseCovarianceMatrix(p); }
  VectorType TrueCCAPowerUpdate(RealType penaltyP, MatrixType p , VectorType w_q , MatrixType q, bool keep_pos, VectorType covar, bool factorOutR);
  MatrixType PartialOutZ( MatrixType X, MatrixType Y, MatrixType Z ) {
    /** compute the effect of Z and store it for later use */
  }

  VectorType GetPWeights() { return this->m_WeightsP; }
  VectorType GetQWeights() { return this->m_WeightsQ; }
  VectorType GetRWeights() { return this->m_WeightsR; }
  RealType GetCorrelationForSignificanceTest() { return this->CorrelationForSignificanceTest; }

  void SetCovariatesP(VectorType vec) { this->m_CovariatesP=vec; }
  void SetCovariatesQ(VectorType vec) { this->m_CovariatesQ=vec; }

protected:

// for pscca 
  void UpdatePandQbyR( );

  RealType CountNonZero( VectorType v ) 
  {
    unsigned long ct=0;
    for ( unsigned int i=0; i<v.size(); i++) 
      if ( v[i] != 0 ) ct++;
    return (RealType)ct/(RealType)v.size();
  }

  RealType PearsonCorr(VectorType v1, VectorType v2 )
  {
  double xysum=0;
  for ( unsigned int i=0; i<v1.size(); i++) xysum+=v1(i)*v2(i);
  double frac=1.0/(double)v1.size();
  double xsum=v1.sum(),ysum=v2.sum();
  double xsqr=v1.squared_magnitude();
  double ysqr=v2.squared_magnitude();
  double numer=xysum - frac*xsum*ysum;
  double denom=sqrt( ( xsqr - frac*xsum*xsum)*( ysqr - frac*ysum*ysum) );
  if ( denom <= 0 ) return 0;
  return numer/denom;
  }

  void SoftThresholdByRMask()
  {
    std::cout <<" enter r mask for special cognitive domain processing " << std::endl;
    unsigned long ct=0;
    unsigned int cols=this->m_MatrixR.columns();
    unsigned int vox=this->m_MaskImageR->GetLargestPossibleRegion().GetNumberOfPixels();
    std::cout << " got vox  " << vox << " cols " << cols << std::endl;
    typedef itk::ImageRegionIteratorWithIndex<ImageType> Iterator;
    if ( cols != vox ) 
      {
	std::cout << " exit r mask " <<std::endl;
        return;
      }
    std::cout <<" begin r mask " << std::endl;
    // zero out R weights that correspond to entries with value 0 in the mask 
    Iterator vfIter(this->m_MaskImageR,this->m_MaskImageR->GetLargestPossibleRegion() );
    for(  vfIter.GoToBegin(); !vfIter.IsAtEnd(); ++vfIter ){
      if ( vfIter.Get() == 0 ) {
        VectorType zero( this->m_MatrixR.rows() );  zero.fill(0);
	this->m_MatrixR.set_column(ct,zero);
      }
      ct++;
    }
  }

  void FactorOutCovariates()
  {
    typedef itk::ImageRegionIteratorWithIndex<ImageType> Iterator;
    unsigned long ct=0;
    VectorType WP(this->m_WeightsP.size());  WP.fill(0);
    {// zero out P weights that correspond to entries with value 2 in the mask 
      Iterator vfIter(this->m_MaskImageP,this->m_MaskImageP->GetLargestPossibleRegion() );
      for(  vfIter.GoToBegin(); !vfIter.IsAtEnd(); ++vfIter ){
        if ( vfIter.Get() == 2 ) { WP[ct]=this->m_WeightsP[ct]; this->m_WeightsP[ct]=0; }
        if ( vfIter.Get() != 0 ) ct++;
      }
    }
    ct=0;
    VectorType WQ(this->m_WeightsQ.size());  WQ.fill(0);
    {// zero out Q weights that correspond to entries with value 2 in the mask 
      Iterator vfIter(this->m_MaskImageQ,this->m_MaskImageQ->GetLargestPossibleRegion() );
      for(  vfIter.GoToBegin(); !vfIter.IsAtEnd(); ++vfIter ){
        if ( vfIter.Get() == 2 ) { WQ[ct]=this->m_WeightsQ[ct]; /*std::cout << " wt " << WQ[ct] << std::endl; */this->m_WeightsQ[ct]=0;  }
        if ( vfIter.Get() != 0 ) ct++;
      }
    }
    //    this->m_CovariatesQ=this->m_MatrixQ*WQ;
    // this->m_CovariatesP=this->m_MatrixP*WP;
    this->m_CovariatesQ=this->m_CovariatesQ+this->m_MatrixQ*WQ*0.2;
    this->m_CovariatesP=this->m_CovariatesP+this->m_MatrixP*WP*0.2;
  }

  RealType SpecializedCorrelation()
  {
    typedef itk::ImageRegionIteratorWithIndex<ImageType> Iterator;
    unsigned long ct=0;
    {// zero out P weights that correspond to entries with value 2 in the mask 
    Iterator vfIter(this->m_MaskImageP,this->m_MaskImageP->GetLargestPossibleRegion() );
    for(  vfIter.GoToBegin(); !vfIter.IsAtEnd(); ++vfIter ){
      if ( vfIter.Get() == 2 ) this->m_WeightsP[ct]=0;
      if ( vfIter.Get() != 0 ) ct++;
    }
    }
    ct=0;
    {// zero out Q weights that correspond to entries with value 2 in the mask 
    Iterator vfIter(this->m_MaskImageQ,this->m_MaskImageQ->GetLargestPossibleRegion() );
    for(  vfIter.GoToBegin(); !vfIter.IsAtEnd(); ++vfIter ){
      if ( vfIter.Get() == 2 ) { //std::cout << " zeroing " <<  this->m_WeightsQ[ct] << std::endl; 
	this->m_WeightsQ[ct]=0; }
      if ( vfIter.Get() != 0 ) ct++;
    }
    }
    ct=0;
    if ( this->m_MaskImageR && this->m_WeightsR.size() > 0 )
    {// zero out R weights that correspond to entries with value 2 in the mask 
    Iterator vfIter(this->m_MaskImageR,this->m_MaskImageR->GetLargestPossibleRegion() );
    for(  vfIter.GoToBegin(); !vfIter.IsAtEnd(); ++vfIter ){
      if ( vfIter.Get() == 2 ) {
	//	std::cout << " post-zeroing " << ct << " wt " << this->m_WeightsR[ct] <<std::endl;
        this->m_WeightsR[ct]=0;
      }
      if ( vfIter.Get() != 0 ) ct++;
    }
    }
    
    VectorType pvec=this->m_MatrixP*this->m_WeightsP;
    VectorType qvec=this->m_MatrixQ*this->m_WeightsQ;
    double corrpq=this->PearsonCorr( pvec , qvec );

    double corrpr=0;
    double corrqr=0;
    VectorType rvec;
    if ( this->m_MaskImageR && this->m_WeightsR.size() > 0 )
    {
      rvec=this->m_MatrixR*this->m_WeightsR;
      corrpr=this->PearsonCorr( pvec , rvec );
      corrqr=this->PearsonCorr( rvec , qvec );
    }
    else {
      std::cout <<" pre-specialization corr " << this->m_CorrelationForSignificanceTest << " post-specialization " <<  corrpq << std::endl;
      return corrpq;
    }
    /** FIXME!! this is biserial!!! */
    std::cout << "USING pr+qr correlation for significance: " <<corrpr+corrqr<<" not 3corr: "<<corrpr+corrqr+corrpq <<std::endl;
    return corrpr+corrqr;
    return corrpq+corrpr+corrqr;
  }

  antsSCCANObject(); 
  ~antsSCCANObject() {  }
 
  void PrintSelf( std::ostream& os, Indent indent ) const
  {
    if ( this->m_MaskImageP && this->m_MaskImageQ && this->m_MaskImageR ) std::cout << " 3 matrices " << std::endl;
    else if ( this->m_MaskImageP && this->m_MaskImageQ  ) std::cout << " 2 matrices " << std::endl;
    else std::cout << " fewer than 2 matrices " << std::endl;
  }

private:
  MatrixType m_OriginalMatrixP;
  MatrixType m_OriginalMatrixQ;
  MatrixType m_OriginalMatrixR;

  antsSCCANObject(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

  unsigned int                                   m_ElapsedIterations;
  unsigned int                                   m_MaximumNumberOfIterations;
  RealType                                       m_CurrentConvergenceMeasurement;
  RealType                                       m_ConvergenceThreshold;

  SCCANFormulationType            m_SCCANFormulation;
  RealType m_PinvTolerance;
  RealType m_PercentVarianceForPseudoInverse;

  VectorType m_WeightsP;
  MatrixType m_MatrixP;
  ImagePointer m_MaskImageP;
  RealType   m_FractionNonZeroP;
  bool       m_KeepPositiveP;
  VectorType m_CovariatesP;

  VectorType m_WeightsQ;
  MatrixType m_MatrixQ;
  ImagePointer m_MaskImageQ;
  RealType   m_FractionNonZeroQ;
  bool       m_KeepPositiveQ;
  VectorType m_CovariatesQ;

  VariateType m_VariatesP;
  VariateType m_VariatesQ;

  VectorType m_WeightsR;
  MatrixType m_MatrixR;
  ImagePointer m_MaskImageR;
  RealType   m_FractionNonZeroR;
  bool       m_KeepPositiveR;
  VectorType m_CovariatesR;
/** a special variable for pscca, holds R^T R */
  MatrixType m_MatrixRRt;


  bool m_AlreadyWhitened;
  bool m_SpecializationForHBM2011;
  RealType m_CorrelationForSignificanceTest;

};

} // namespace ants
} // namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "antsSCCANObject.txx"
#endif

#endif


