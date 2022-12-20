//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2021 QMCPACK developers
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Lab
//
// File created by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Lab
//////////////////////////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_QMC_WAVEFUNCTION_TYPES_HPP
#define QMCPLUSPLUS_QMC_WAVEFUNCTION_TYPES_HPP

#include "type_traits/complex_help.hpp"
#include "OhmmsPETE/OhmmsMatrix.h"

namespace qmcplusplus
{

/** Consistent way to get the set of types used in the QMCWaveFunction module, without resorting
 *  to ifdef'd type aliases accessed through inheritance.
 *
 *  You can at least test all flavors without tricky scoping or recompiling.
 *  Determines type set for class through composition instead of inheritance.
 *
 *  Its somewhat question whether it is even worth have all these shortened type aliases
 *  Defined differently in different modules of the code, but it needs further study.
 *
 *  I would have liked to use QMCTypes but they do something dumb with handling complex
 *  which basically bakes in multiple builds for Real and Complex.
 */
template<typename VALUE, typename FP_VALUE>
struct WaveFunctionTypes final
{
  // Should be factored up into a non broken QMCTypes
  using Value         = VALUE;
  using FullPrecValue = FP_VALUE;
  using Real          = RealAlias<Value>;
  using FullPrecReal  = RealAlias<FullPrecValue>;
  using Complex       = std::complex<Real>;
  // This is all that belongs here so far.
  using PsiValue = FP_VALUE;
  using Grad     = TinyVector<Value, OHMMS_DIM>;
  using Hess     = Tensor<Value, OHMMS_DIM>;
  using LogValue = std::complex<FullPrecReal>;
};

} // namespace qmcplusplus
#endif
