 //////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2018 QMCPACK developers
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Lab
//
// File created by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Lab
//////////////////////////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_CUDA_TYPE_ALIASES_H
#define QMCPLUSPLUS_CUDA_TYPE_ALIASES_H

#include "config.h"
#include <OhmmsPETE/TinyVector.h>
#include <OhmmsPETE/Tensor.h>

namespace qmcplusplus
{

/* Facilitates use of CUDA real/complex full/mixed precision 
 * without rebuilding the code
 *
 * a CUDA device implmentation templated on precision (P) 
 * and value type (V)
 * would write this template<typename P, V>
 * class FooDeviveCUDA
 * {
 *   using CTA = CUDATypes<P, V>;
 *   ...
 * };    
 * Then instead of CudaGradType CTA::GradType.
 *
 * on P and V different in complexity is supported
 */
template<typename P, typename V, int DIM>
class CUDATypes;

// Partial Specializations

// Real
template<typename P, int DIM>
class CUDATypes<P, P, DIM>  final
{
public:
  using RealType = P;
  using ComplexType = std::complex<P>;
  using GradType = TinyVector<P,DIM>;
  using ValueType = P;
  using PosType = TinyVector<P,DIM>;
};

// Complex
template<typename P, int DIM>
class CUDATypes<P, std::complex<P>, DIM>  final
{
public:
  using RealType = P;
  using ComplexType = std::complex<P>;
  using GradType = TinyVector<std::complex<P>,DIM>;
  using ValueType = std::complex<P>;
  using PosType = TinyVector<P,DIM>;
};

// Right now we build only one CUDA precision and complexity at a time.
// This was set once in QMCTraits then inherited everywhere
// The CUDAGlobalTypeAliases reproduces this old behavior for now
#ifdef QMC_COMPLEX
using CUDAGlobalTypeAliases = CUDATypes<CUDA_PRECISION, std::complex<CUDA_PRECISION>, OHMMS_DIM>;
#else
using CUDAGlobalTypeAliases = CUDATypes<CUDA_PRECISION, CUDA_PRECISION, OHMMS_DIM>;
#endif  
 
}
#endif
