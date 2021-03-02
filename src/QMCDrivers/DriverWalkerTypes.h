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

/** DriverWalkerResourceCollection locks
 * Helper class for acquiring and releasing multi walker resources
 */
class DriverWalkerResourceCollectionLock
{
public:
  DriverWalkerResourceCollectionLock(DriverWalkerResourceCollection& driverwalker_res, ParticleSet& pset, TrialWaveFunction& twf, QMCHamiltonian& ham)
      : pset_res_lock_(driverwalker_res.pset_res, pset), twf_res_lock_(driverwalker_res.twf_res, twf), ham_res_lock_(driverwalker_res.ham_res, ham)
  {}

private:
  ResourceCollectionLock<ParticleSet> pset_res_lock_;
  ResourceCollectionLock<TrialWaveFunction> twf_res_lock_;
  ResourceCollectionLock<QMCHamiltonian> ham_res_lock_;
};

/** DriverWalkerResourceCollection locks
 * Helper class for acquiring and releasing multi walker resources.
 * Only ParticleSet and TrialWaveFunction resources are engaged.
 */
class DriverWalkerResourceCollection_PsetTWF_Lock
{
public:
  DriverWalkerResourceCollection_PsetTWF_Lock(DriverWalkerResourceCollection& driverwalker_res, ParticleSet& pset, TrialWaveFunction& twf)
      : pset_res_lock_(driverwalker_res.pset_res, pset), twf_res_lock_(driverwalker_res.twf_res, twf)
  {}

private:
  ResourceCollectionLock<ParticleSet> pset_res_lock_;
  ResourceCollectionLock<TrialWaveFunction> twf_res_lock_;
};
} // namespace qmcplusplus
#endif
