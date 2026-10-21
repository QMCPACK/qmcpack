//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Cynthia Gu, zg1@ornl.gov, Oak Ridge National Laboratory
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


#ifndef QMCPLUSPLUS_WALKER_OUTPUT_H
#define QMCPLUSPLUS_WALKER_OUTPUT_H

#include "Particle/WalkerConfigurations.h"
#include <utility>
#include "hdf/hdf_archive.h"

namespace qmcplusplus
{
/** Writes a set of walker configurations to an HDF5 file. */
class HDFWalkerOutput
{
  /** number of walkers when a state is dumped
   *
   * When the number of walkers per state has changed, NumOfWalkers is used
   * to reallocate the hdf5 group.
   */
  size_t number_of_walkers_;
  /** number of particles */
  const size_t number_of_particles_;
  ///communicator
  Communicate& myComm;

public:
  ///constructor
  HDFWalkerOutput(size_t num_ptcls, Communicate& c);
  ///destructor
  ~HDFWalkerOutput();

  /** dump configurations
   * @param w walkers
   * @param fname file root name
   * @param block current block
   */
  bool dump(const WalkerConfigurations& w, const std::string& fname, int block);

private:
  ///PooledData<T> is used to define the shape of multi-dimensional array
  using BufferType = PooledData<OHMMS_PRECISION>;
  std::array<BufferType, 2> RemoteData;
  std::array<std::vector<QMCTraits::FullPrecRealType>, 2> RemoteDataW;
  int block;
  void write_configuration(const WalkerConfigurations& W, hdf_archive& hout, int block);
};

} // namespace qmcplusplus
#endif
