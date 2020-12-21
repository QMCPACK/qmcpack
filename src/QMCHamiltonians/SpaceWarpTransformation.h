//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Raymond Clay
//
// File created by: Raymond Clay, Sandia National Laboratory, rclay@sandia.gov
//////////////////////////////////////////////////////////////////////////////////////


#ifndef QMCPLUSPLUS_SPACEWARP_H
#define QMCPLUSPLUS_SPACEWARP_H

#include "Configuration.h"

namespace qmcplusplus
{
/** @ingroup hamiltonian
 *\brief Calculates the AA Coulomb potential using PBCs
 *
 * Functionally identical to CoulombPBCAB but uses a templated version of
 * LRHandler.
 */
struct SpaceWarpTransformation : public QMCTraits
{
  SpaceWarpTransformation(){};
  ~SpaceWarpTransformation(){};

  RealType f(RealType r); 
};
}
#endif


