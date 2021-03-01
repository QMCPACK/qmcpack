//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2018 QMCPACK developers
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Lab
//                    Ye Luo, yeluo@anl.gov, Argonne National Lab
//
// File created by: Ye Luo, yeluo@anl.gov, Argonne National Lab
//////////////////////////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_QMC_TYPES_H
#define QMCPLUSPLUS_QMC_TYPES_H

#include "OhmmsPETE/TinyVector.h"
#include "OhmmsPETE/Tensor.h"
#include <complex>

namespace qmcplusplus
{
/* Facilitates use of full/mixed precision without rebuilding the code
 *
 * a template class may define its local types
 * template<typename Prec>
 * class Foo
 * {
 *   using FooTypes = QMCTypes<Prec>;
 *   ...
 * };
 *
 */
template<typename Precision, int DIM>
class QMCTypes final
{
public:
  using RealType    = Precision;
  using ComplexType = std::complex<Precision>;
#ifdef QMC_COMPLEX
  using ValueType = ComplexType;
#else
  using ValueType = RealType;
#endif
  using GradType   = TinyVector<ValueType, DIM>;
  using PosType    = TinyVector<RealType, DIM>;
  using TensorType = Tensor<RealType, DIM>;
};

} // namespace qmcplusplus
#endif
