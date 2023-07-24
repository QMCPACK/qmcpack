//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by:
//
// File created by: Jeongnim Kim, jeongnim.kim@intel.com, Intel Corp.
//////////////////////////////////////////////////////////////////////////////////////

/** @file SoaAtomicBasisSet.h
 */
#ifndef QMCPLUSPLUS_SOA_SPHERICALORBITAL_BASISSET_H
#define QMCPLUSPLUS_SOA_SPHERICALORBITAL_BASISSET_H

#include "Configuration.h"
#include "SoaAtomicBasisSetT.h"

namespace qmcplusplus
{
/* A basis set for a center type 
   *
   * @tparam ROT : radial function type, e.g.,NGFunctor<T>
   * @tparam SH : spherical or carteisan Harmonics for (l,m) expansion
   *
   * \f$ \phi_{n,l,m}({\bf r})=R_{n,l}(r) Y_{l,m}(\theta) \f$
   */
template<typename ROT, typename SH>
using SoaAtomicBasisSet = SoaAtomicBasisSetT<ROT, SH, QMCTraits::ValueType>;

} // namespace qmcplusplus
#endif
