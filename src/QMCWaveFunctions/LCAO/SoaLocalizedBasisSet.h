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


/** @file SoaLocalizedBasisSet.h
 * @brief A derived class from BasisSetBase
 *
 * This is intended as a replacement for MolecularWaveFunctionComponent and
 * any other localized basis set.
 */
#ifndef QMCPLUSPLUS_SOA_LOCALIZEDBASISSET_H
#define QMCPLUSPLUS_SOA_LOCALIZEDBASISSET_H

#include "Configuration.h"
#include "QMCWaveFunctions/LCAO/SoaLocalizedBasisSetT.h"

namespace qmcplusplus
{
template<class COT, typename ORBT>
using SoaLocalizedBasisSet = SoaLocalizedBasisSetT<COT, ORBT>;

} // namespace qmcplusplus
#endif
