//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2021 QMCPACK developers.
//
// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////


/**
 * @file
 * Driver level walker (DriverWalker) related data structures.
 * Unlike MCWalkerConfiguration which only holds electron positions and weights.
 * Driver level walker includes all the per-walker data structures which depends on the type of driver
 */

#ifndef QMCPLUSPLUS_DRIVERWALKERTYPES_H
#define QMCPLUSPLUS_DRIVERWALKERTYPES_H

#include <ResourceCollection.h>

namespace qmcplusplus
{
/** DriverWalker multi walker resource collections
 * It currently supports VMC and DMC only
 */
struct DriverWalkerResourceCollection
{
  ResourceCollection pset_res;
  ResourceCollection twf_res;
  ResourceCollection ham_res;

  DriverWalkerResourceCollection() : pset_res("ParticleSet"), twf_res("TrialWaveFunction"), ham_res("Hamiltonian") {}
};
} // namespace qmcplusplus
#endif
