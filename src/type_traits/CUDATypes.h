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

#ifndef QMCPLUSPLUS_CUDA_TYPES_H
#define QMCPLUSPLUS_CUDA_TYPES_H

#include "config.h"
#include "OhmmsPETE/TinyVector.h"
#include "OhmmsPETE/Tensor.h"
#include "type_traits/QMCTypes.h"

namespace qmcplusplus
{
/** \brief Provides consistent set of CUDA types 
 * 
 *  Multiple classes need to share the same set of types
 *  for CUDA implementation depending on full/mixed precision
 *  and complex or real values. 
 *  Currently we build only one combination of these based on
 *  QMC_COMPLEX and CUDA_PRECISION.
 *  These types are available via: \c CUDAGlobalTypes \c
 *
 *  usage:
\code{.cpp}
class FooCUDA
{
using CTS = CUDAGlobalTypes
};
\endcode
 *  then CudaGradType is replaced with CTS::GradType.
 *  Do not make local type alias or typedefs use CTS::*Type
 *
 *  A CUDA device implmentation class template 
 *  on precision (P) and value type (V)
 *  would write this. 
\code{.cpp}
template<typename P, V, DIM>
class FooDeviveCUDA
{
  using CTS = CUDATypes<P, V, DIM>;
  ...
};
\endcode
 *
 *  P and V are assumed to match in precision
 *  this is enforced by partial template 
 *  specialization
 */
template<typename P, typename V, int DIM>
class CUDATypes;

// Partial Specializations
// Real
template<typename P, int DIM>
class CUDATypes<P, P, DIM> final
{
public:
  using RealType    = P;
  using ComplexType = std::complex<P>;
  using GradType    = TinyVector<P, DIM>;
  using ValueType   = P;
  using PosType     = TinyVector<P, DIM>;
};

// Complex
template<typename P, int DIM>
class CUDATypes<P, std::complex<P>, DIM> final
{
public:
  using RealType    = P;
  using ComplexType = std::complex<P>;
  using GradType    = TinyVector<std::complex<P>, DIM>;
  using ValueType   = std::complex<P>;
  using PosType     = TinyVector<P, DIM>;
};

// Right now we build only one CUDA precision at a time.
// This was set once in QMCTraits then inherited everywhere
// The CUDAGlobalTypes reproduces this behavior
using CUDAGlobalTypes = QMCTypes<CUDA_PRECISION, OHMMS_DIM>;

} // namespace qmcplusplus
#endif
