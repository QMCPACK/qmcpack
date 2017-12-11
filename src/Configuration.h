//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    


#ifndef QMCPLUSPLUS_TRAITS_H
#define QMCPLUSPLUS_TRAITS_H

#include <config.h>
#include <string>
#include <vector>
#include <map>
#include <complex>
#include <OhmmsPETE/TinyVector.h>
#include <OhmmsPETE/Tensor.h>
#include <OhmmsData/OhmmsElementBase.h>
#include <OhmmsData/RecordProperty.h>
#if OHMMS_DIM==3
#include <Lattice/Uniform3DGridLayout.h>
#elif OHMMS_DIM==2
#include <Lattice/Uniform2DGridLayout.h>
#else
#error "Only 2D and 3D are implemented.\n"
#endif
#include <ParticleBase/ParticleAttrib.h>
#include <Utilities/OutputManager.h>
#include <Message/Communicate.h>

//define empty DEBUG_MEMORY
#define DEBUG_MEMORY(msg)
//uncomment this out to trace the call tree of destructors
//#define DEBUG_MEMORY(msg) std::cerr << "<<<< " << msg << std::endl;

#if defined(DEBUG_PSIBUFFER_ON)
#define DEBUG_PSIBUFFER(who,msg) std::cerr << "PSIBUFFER " << who << " " << msg << std::endl; std::cerr.flush();
#else
#define DEBUG_PSIBUFFER(who,msg)
#endif

namespace qmcplusplus
{

/** traits for the common particle attributes
 *
 *This is an alternative to the global typedefs.
 */
struct PtclAttribTraits
{
  typedef int                                                     Index_t;
  typedef ParticleAttrib<Index_t>                                 ParticleIndex_t;
  typedef ParticleAttrib<OHMMS_PRECISION>                         ParticleScalar_t;
  typedef ParticleAttrib<TinyVector<OHMMS_PRECISION, OHMMS_DIM> > ParticlePos_t;
  typedef ParticleAttrib<Tensor<OHMMS_PRECISION, OHMMS_DIM> >     ParticleTensor_t;

};


/** traits for QMC variables
 *
 *typedefs for the QMC data types
 */
struct QMCTraits
{
  enum {DIM = OHMMS_DIM};
  typedef OHMMS_INDEXTYPE                IndexType;
  typedef OHMMS_PRECISION                RealType;
  typedef OHMMS_PRECISION_FULL           EstimatorRealType;
#if defined(QMC_COMPLEX)
  typedef std::complex<OHMMS_PRECISION>  ValueType;
#ifdef QMC_CUDA
  typedef std::complex<CUDA_PRECISION>   CudaValueType;
#endif
#else
  typedef OHMMS_PRECISION                ValueType;
#ifdef QMC_CUDA
  typedef CUDA_PRECISION                 CudaValueType;
#endif
#endif
  typedef std::complex<RealType>         ComplexType;
  typedef TinyVector<RealType,DIM>       PosType;
  typedef TinyVector<ValueType,DIM>      GradType;
  typedef Tensor<RealType,DIM>           TensorType;
  ///define PropertyList_t
  typedef RecordNamedProperty<EstimatorRealType>    PropertySetType;
#ifdef QMC_CUDA
  typedef CUDA_PRECISION                 CudaRealType;
  typedef TinyVector<CudaValueType,DIM>  CudaGradType;
  typedef TinyVector<CudaRealType,DIM>   CudaPosType;
  typedef std::complex<CUDA_PRECISION>   CudaComplexType;
#endif
};

/** Particle traits to use UniformGridLayout for the ParticleLayout.
 */
struct PtclOnLatticeTraits
{
#if OHMMS_DIM==3
  typedef Uniform3DGridLayout                          ParticleLayout_t;
#elif OHMMS_DIM==2
  typedef Uniform2DGridLayout                          ParticleLayout_t;
#else
  typedef UniformGridLayout<OHMMS_PRECISION,OHMMS_DIM> ParticleLayout_t;
#endif

  typedef int                                          Index_t;
  typedef OHMMS_PRECISION_FULL                         Scalar_t;
  typedef std::complex<Scalar_t>                       Complex_t;

  typedef ParticleLayout_t::SingleParticleIndex_t      SingleParticleIndex_t;
  typedef ParticleLayout_t::SingleParticlePos_t        SingleParticlePos_t;
  typedef ParticleLayout_t::Tensor_t                   Tensor_t;

  typedef ParticleAttrib<Index_t>                      ParticleIndex_t;
  typedef ParticleAttrib<Scalar_t>                     ParticleScalar_t;
  typedef ParticleAttrib<SingleParticlePos_t>          ParticlePos_t;
  typedef ParticleAttrib<Tensor_t>                     ParticleTensor_t;

#if defined(QMC_COMPLEX)
  typedef ParticleAttrib<TinyVector<Complex_t,OHMMS_DIM> > ParticleGradient_t;
  typedef ParticleAttrib<Complex_t>                      ParticleLaplacian_t;
  typedef Complex_t                                      ParticleValue_t;
#else
  typedef ParticleAttrib<TinyVector<Scalar_t,OHMMS_DIM> > ParticleGradient_t;
  typedef ParticleAttrib<Scalar_t>                       ParticleLaplacian_t;
  typedef Scalar_t                                       ParticleValue_t;
#endif
};


// For unit tests
//  Check if we are compiling with Catch defined.  Could use other symbols if needed.
#ifdef TEST_CASE
#ifdef QMC_COMPLEX
typedef ComplexApprox ValueApprox;
#else
typedef Approx ValueApprox;
#endif
#endif

}

#endif
