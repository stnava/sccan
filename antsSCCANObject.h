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
  typedef vnl_diag_matrix<RealType>   DiagonalMatrixType;

  enum SCCANFormulationType{ biSCCA, triSCCA };

  /** ivars Set/Get functionality */
  itkSetMacro( MaximumNumberOfIterations, unsigned int );
  itkGetConstMacro( MaximumNumberOfIterations, unsigned int );
  itkSetMacro( ConvergenceThreshold, RealType );
  itkGetConstMacro( ConvergenceThreshold, RealType );
  itkGetConstMacro( CurrentConvergenceMeasurement, RealType );
  itkGetConstMacro( ElapsedIterations, unsigned int );
  itkSetMacro( SCCANFormulation, SCCANFormulationType );
  itkGetConstMacro( SCCANFormulation, SCCANFormulationType );
  itkGetConstMacro( NumberOfInputMatrices, unsigned int );

  void SetPseudoInversePercentVariance( RealType p ) { this->m_PercentVarianceForPseudoInverse=p; }

  itkSetMacro( FractionNonZeroP, RealType );
  itkSetMacro( KeepPositiveP, bool );
  void SetMaskImage1( ImagePointer mask ) { this->m_Mask1=mask; }
  void SetMatrixP(  MatrixType matrix ) { this->m_MatrixP=matrix; }

  itkSetMacro( FractionNonZeroQ, RealType );
  itkSetMacro( KeepPositiveQ, bool );
  void SetMaskImageQ( ImagePointer mask ) { this->m_MaskQ=mask; }
  void SetMatrixQ(  MatrixType  matrix ) { this->m_MatrixQ=matrix; }

  itkSetMacro( FractionNonZeroR, RealType );
  itkSetMacro( KeepPositiveR, bool );
  void SetMaskImageR( ImagePointer mask ) { this->m_MaskR=mask; }
  void SetMatrixR(  MatrixType matrix ) { this->m_MatrixR=matrix; }

  RealType RunSCCAN2();
  RealType RunSCCAN3();
 
  VectorType SoftThreshold( VectorType v_in, RealType fractional_goal , bool allow_negative_weights );
  VectorType InitializeV( MatrixType p );
  MatrixType NormalizeMatrix(MatrixType p);
  MatrixType WhitenMatrix(MatrixType p); 
  VectorType TrueCCAPowerUpdate(RealType penaltyP, MatrixType p , VectorType w_q , MatrixType q, bool keep_pos);

  VectorType GetPWeights() { return this->m_WeightsP; }
  VectorType GetQWeights() { return this->m_WeightsQ; }
  VectorType GetRWeights() { return this->m_WeightsR; }

protected:

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

  antsSCCANObject() 
  {
    this->m_PinvTolerance=1.e-6;
    this->m_PercentVarianceForPseudoInverse=0.99;
    this->m_MaximumNumberOfIterations=15;
    this->m_MaskImageP=NULL;
    this->m_MaskImageQ=NULL;
    this->m_MaskImageR=NULL;
    this->m_KeepPositiveP=true;
    this->m_KeepPositiveQ=true;
    this->m_KeepPositiveR=true;
    this->m_FractionNonZeroP=0.5;
    this->m_FractionNonZeroQ=0.5;
    this->m_FractionNonZeroR=0.5;
    this->m_NumberOfInputMatrices=0;
    this->m_ConvergenceThreshold=1.e-6;
  } 
  ~antsSCCANObject() {  }
 
  void PrintSelf( std::ostream& os, Indent indent ) const
  {
    if ( this->m_MaskImageP && this->m_MaskImageQ && this->m_MaskImageR ) std::cout << " 3 matrices " << std::endl;
    else if ( this->m_MaskImageP && this->m_MaskImageQ  ) std::cout << " 2 matrices " << std::endl;
    else std::cout << " fewer than 2 matrices " << std::endl;
  }

private:
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

  VectorType m_WeightsR;
  MatrixType m_MatrixR;
  ImagePointer m_MaskImageR;
  RealType   m_FractionNonZeroR;
  bool       m_KeepPositiveR;
  unsigned int m_NumberOfInputMatrices;
};

} // namespace ants
} // namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "antsSCCANObject.txx"
#endif

#endif


