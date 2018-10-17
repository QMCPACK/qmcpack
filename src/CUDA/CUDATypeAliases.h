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

/* Facilitates use of CUDA full and mixed precision without rebuilding the code
 *
 * a CUDA device implmentation templated on precision (P) 
 * would write this template<typename P>
 * class FooDeviveCUDA
 * {
 *   using CTA = CUDATypeAliases<P>;
 *   ...
 * };    
 * Then instead of CudaGradType CTA::GradType.
 */
template<typename P, int DIM>
class CUDATypeAliases final
{
public:
  using RealType = P;
  using ComplexType = std::complex<P>;
#ifdef QMC_COMPLEX
  using GradType = TinyVector<std::complex<P>,DIM>;
  using ValueType = std::complex<P>;
#else
  using GradType = TinyVector<P,DIM>;
  using ValueType = P;
#endif
  using PosType = TinyVector<P,DIM>;
};

// Right now we build only one CUDA precision at a time.
// This was set once in QMCTraits then inherited everywhere
using CUDAGlobalTypeAliases = CUDATypeAliases<CUDA_PRECISION, OHMMS_DIM>;
 
}
#endif
