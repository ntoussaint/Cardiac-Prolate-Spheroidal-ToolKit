/*=========================================================================

Program:   ImagingSciences
Module:    $Id: itkNewUoaOptimizer.h 1 2010-05-21 14:00:33Z nt08 $
Language:  C++
Author:    $Author: nt08 $
Date:      $Date: 2010-05-21 14:00:33 +0000 (Fri, 21 May 2010) $
Version:   $Revision: 1 $

Copyright (c) 2010 King's College London - Division of Imaging Sciences. All rights reserved.
See Copyright.txt for details.

This software is distributed WITHOUT ANY WARRANTY; without even
the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkNewUoaOptimizer_h
#define __itkNewUoaOptimizer_h


#include "itkSingleValuedNonLinearOptimizer.h"
#include "itkVariableLengthVector.h"


namespace itk
{
  class ITK_EXPORT NewUoaOptimizer :
  public SingleValuedNonLinearOptimizer
  {
  public:

    /** Standard typedefs. */
    typedef NewUoaOptimizer                 Self;
    typedef SingleValuedNonLinearOptimizer  SuperClass;
    typedef SmartPointer<Self>              Pointer;
    typedef SmartPointer<const Self>        ConstPointer;
    typedef SuperClass::ParametersType      ParametersType;
    typedef double                          ValueType;
    typedef VariableLengthVector<ValueType> ParameterVectorType;
    
    /** Method for creation through the object factory */
    itkNewMacro(Self) ;
    /**  Run-time type information (and related methods)  */
    itkTypeMacro( NewUoaOptimizer, SingleValuedNonLinearOptimizer);
    /** Set/Get space dimension */
    itkSetMacro( SpaceDimension, unsigned int);
    itkGetConstReferenceMacro( SpaceDimension, unsigned int) ;
    /** Set/Get the maximum number of function calls */
    itkSetMacro( MaxFunctionCalls, unsigned int);
    itkGetConstReferenceMacro( MaxFunctionCalls, unsigned int) ;
    /** Get the actual number of function calls */
    itkGetConstReferenceMacro( NbFunctionCalls, unsigned int) ;
    /** Get the best value of the objective function */
    itkGetConstReferenceMacro( BestValue, double) ;
    /** Set/Get rhobeg and rhoend */
    itkSetMacro( RhoBeg, double);
    itkGetConstReferenceMacro( RhoBeg, double);
    /** Set/Get rhobeg and rhoend */  
    itkSetMacro( RhoEnd, double);
    itkGetConstReferenceMacro( RhoEnd, double);
    /** Set/Get the number of interpolation conditions */
    itkSetMacro( NbInterp, long int);
    itkGetConstReferenceMacro( NbInterp, long int);
    /** Start the optimizer */
    void StartOptimization();
    /** Set verbosity of the process 0/1/2/3 */
    itkSetClampMacro (Verbosity, unsigned int, 0, 3);
    itkGetMacro (Verbosity, unsigned int);
    itkBooleanMacro (Verbosity);
    /** set the bounds of the parameters */
    void SetBounds (ParameterVectorType* bds)
    { m_Bounds = bds; }
    ParameterVectorType* GetBounds (void)
    { return m_Bounds; }
    
    
  protected:
    NewUoaOptimizer();
    virtual ~NewUoaOptimizer();
    void PrintSelf(std::ostream& os, Indent indent) const ;
  
    // NewUOA optimization section
    int newuoa(double *w, long int *n, long int *npt, double *x, 
	       double *rhobeg, double *rhoend, unsigned int *iprint, long int * maxfun);

    int newuob(long int *n, long int *npt, double *x, 
	       double *rhobeg, double *rhoend, unsigned int *iprint, long int * maxfun, 
	       double *xbase, double *xopt, double *xnew, 
	       double *xpt, double *fval, double *gq, double *hq, 
	       double *pq, double *bmat, double *zmat, long int *ndim, 
	       double *d__, double *vlag, double *w);
	
    int biglag(long int *, long int *, double *, 
	       double *, double *, double *, long int *, long int *, 
	       long int *, double *, double *, double *, double *,
	       double *, double *, double *, double *);
		
    int bigden(long int *, long int *, double *, double *, double *, 
	       double *, long int *, long int *, long int *, long int *, 
	       double *, double *, double *, double *, 
	       double *, double *, double *);
	
    int update(long int *, long int *, 
	       double *, double *, long int *, long int *, double *, 
	       double *, long int *, double *);
	
    int trsapp(long int *, long int *, double *, 
	       double *, double *, double *, double *, 
	       double *, double *, double *, double *, 
	       double *, double *, double *);

  private:
    NewUoaOptimizer(const Self&);  // Purposely not implemented
    void operator=(const Self&)  ; // Purposely not implemented

    /** Space dimension */
    unsigned int m_SpaceDimension;
    /** Maximum number of function calls */
    unsigned int m_MaxFunctionCalls;
    /** Number of function calls. */
    unsigned int m_NbFunctionCalls;
    /** Best function value. */
    double       m_BestValue;
    /** rhobeg and rhoend. */
    double       m_RhoBeg;
    double       m_RhoEnd;
    /** Number of interpolation conditions. */
    long int     m_NbInterp;
    /** Verbosity of the optimization */
    unsigned int m_Verbosity;

    ParameterVectorType* m_Bounds;
    

  } ; // end of class

} // end of namespace itk


#endif
