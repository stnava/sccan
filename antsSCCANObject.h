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

  MatrixType PseudoInverse( MatrixType );

  VectorType Orthogonalize(VectorType Mvec, VectorType V , MatrixType* projecter = NULL  )
  {
    if ( ! projecter ) {
      double ratio=inner_product(Mvec,V)/inner_product(V,V);
      VectorType  ortho=Mvec-V*ratio;
      return ortho;
    } else {
      double ratio=inner_product(*projecter*Mvec,*projecter*V)/inner_product(*projecter*V,*projecter*V);
      VectorType  ortho=Mvec-V*ratio;
      return ortho;
    }
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

  itkSetMacro( FractionNonZeroP, RealType );
  itkSetMacro( KeepPositiveP, bool );
  void SetMaskImageP( ImagePointer mask ) { this->m_MaskImageP=mask; }
  void SetMatrixP(  MatrixType matrix ) { this->m_OriginalMatrixP.set_size(matrix.rows(),matrix.cols());  this->m_MatrixP.set_size(matrix.rows(),matrix.cols()); this->m_OriginalMatrixP.update(matrix); this->m_MatrixP.update(matrix); }
  
  itkSetMacro( FractionNonZeroQ, RealType );
  itkSetMacro( KeepPositiveQ, bool );
  void SetMaskImageQ( ImagePointer mask ) { this->m_MaskImageQ=mask; }
  void SetMatrixQ(  MatrixType  matrix ) {  this->m_OriginalMatrixQ.set_size(matrix.rows(),matrix.cols());  this->m_MatrixQ.set_size(matrix.rows(),matrix.cols()); this->m_OriginalMatrixQ.update(matrix); this->m_MatrixQ.update(matrix);}

  itkSetMacro( FractionNonZeroR, RealType );
  itkSetMacro( KeepPositiveR, bool );
  void SetMaskImageR( ImagePointer mask ) { this->m_MaskImageR=mask; }
  void SetMatrixR(  MatrixType matrix ) {  this->m_OriginalMatrixR.set_size(matrix.rows(),matrix.cols());  this->m_MatrixR.set_size(matrix.rows(),matrix.cols()); this->m_OriginalMatrixR.update(matrix); this->m_MatrixR.update(matrix); }

  MatrixType GetMatrixP(  ) { return this->m_MatrixP; }
  MatrixType GetMatrixQ(  ) { return this->m_MatrixQ; }
  MatrixType GetMatrixR(  ) { return this->m_MatrixR; }
  MatrixType GetOriginalMatrixP(  ) { return this->m_OriginalMatrixP; }
  MatrixType GetOriginalMatrixQ(  ) { return this->m_OriginalMatrixQ; }
  MatrixType GetOriginalMatrixR(  ) { return this->m_OriginalMatrixR; }

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
  VectorType TrueCCAPowerUpdate(RealType penaltyP, MatrixType p , VectorType w_q , MatrixType q, bool keep_pos, bool factorOutR);
  MatrixType PartialOutZ( MatrixType X, MatrixType Y, MatrixType Z ) {
    /** compute the effect of Z and store it for later use */
  }

  VectorType GetPWeights() { return this->m_WeightsP; }
  VectorType GetQWeights() { return this->m_WeightsQ; }
  VectorType GetRWeights() { return this->m_WeightsR; }
  RealType GetCorrelationForSignificanceTest() { return this->CorrelationForSignificanceTest; }

  VectorType GetCanonicalCorrelations( ) 
  { 
    return this->m_CanonicalCorrelations;
  }

  VectorType GetVariateP( unsigned int i = 0 ) 
  { 
    return this->m_VariatesP.get_column(i); 
  }
  VectorType GetVariateQ( unsigned int i = 0 ) 
  { 
    return this->m_VariatesQ.get_column(i); 
  }
  MatrixType GetVariatesP() 
  { 
    return this->m_VariatesP; 
  }
  MatrixType GetVariatesQ() 
  { 
    return this->m_VariatesQ; 
  }

protected:

// for pscca 
  void UpdatePandQbyR( );

  MatrixType  DeleteCol( MatrixType p_in , unsigned int col)
  {
  unsigned int ncols=p_in.cols()-1;
  if ( col >= ncols ) ncols=p_in.cols();
  MatrixType p(p_in.rows(),ncols);      
  unsigned int colct=0;
  for ( long i=0; i<p.cols(); ++i) { // loop over cols
    if ( i != col ) {
      p.set_column(colct,p_in.get_column(i));
      colct++;
    }
  }
  return p;
  } 

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

  antsSCCANObject(); 
  ~antsSCCANObject() {  }
 
  void PrintSelf( std::ostream& os, Indent indent ) const
  {
    if ( this->m_MaskImageP && this->m_MaskImageQ && this->m_MaskImageR ) std::cout << " 3 matrices " << std::endl;
    else if ( this->m_MaskImageP && this->m_MaskImageQ  ) std::cout << " 2 matrices " << std::endl;
    else std::cout << " fewer than 2 matrices " << std::endl;
  }
 
  void RunDiagnostics(unsigned int);

private:
  bool m_Debug;
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

  VectorType m_WeightsQ;
  MatrixType m_MatrixQ;
  ImagePointer m_MaskImageQ;
  RealType   m_FractionNonZeroQ;
  bool       m_KeepPositiveQ;


  VectorType  m_CanonicalCorrelations;
  VariateType m_VariatesP;
  VariateType m_VariatesQ;

  VectorType m_WeightsR;
  MatrixType m_MatrixR;
  ImagePointer m_MaskImageR;
  RealType   m_FractionNonZeroR;
  bool       m_KeepPositiveR;
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


