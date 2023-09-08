//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by:   Raymond Clay III, rclay@sandia.gov, Sandia National Laboratories
//
// File created by:    Raymond Clay III, rclay@sandia.gov, Sandia National Laboratories
//////////////////////////////////////////////////////////////////////////////////////


/** @file EinsplineSpinorSetBuilder.h
 * 
 * Derives EinsplineSetBuilder.  Overrides the createSPOSetFromXML method to read an up and down channel from hdf5
 *   and then construct an appropriate einspline spinor set object.
 *
 */
#ifndef QMCPLUSPLUS_EINSPLINE_SPINORSET_BUILDER_H
#define QMCPLUSPLUS_EINSPLINE_SPINORSET_BUILDER_H

#include "Configuration.h"
#include "QMCWaveFunctions/EinsplineSpinorSetBuilderT.h"

namespace qmcplusplus
{
using EinsplineSpinorSetBuilder = EinsplineSpinorSetBuilderT<QMCTraits::ValueType>;

} // namespace qmcplusplus

#endif
