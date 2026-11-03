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
#include <filesystem>
#include "hdf/hdf_archive.h"

namespace qmcplusplus
{
/** Writes a set of walker configurations and weights to an HDF5 restart file (<root>.config.h5).
 *
 *  The output file contains:
 *  - `/version`: current HDF version
 *  - `/state_0/block`: integer block index
 *  - `/state_0/number_of_walkers`: total number of walkers across all ranks
 *  - `/state_0/walker_partition`: rank offsets partition array
 *  - `/state_0/walkers`: 3D array of walker configurations `(total_walkers, num_particles, OHMMS_DIM)`
 *  - `/state_0/walker_weights`: 1D array of walker weights `(total_walkers)`
 *
 *  Supports two writing modes:
 *  - Parallel collective I/O (`write_hdf_parallel`): enabled when `ENABLE_PHDF5` is defined.
 *    Each rank writes its walker subset directly using parallel HDF5 hyperslabs.
 *  - Rank 0 gather I/O (`write_hdf_rank0`): used when `ENABLE_PHDF5` is not defined.
 *    Walkers and weights are gathered to rank 0 using MPI `gatherv`, and rank 0 writes the file.
 */
class HDFWalkerOutput
{
  /// Total number of walkers across all ranks when a state is dumped
  size_t number_of_walkers_;
  /// Number of particles per walker
  const size_t number_of_particles_;
  /// Communicator coordinating parallel execution
  Communicate& myComm;

public:
  /** Constructor
   * @param num_ptcls Number of particles per walker configuration
   * @param c Communicator
   */
  HDFWalkerOutput(size_t num_ptcls, Communicate& c);
  /// Destructor
  ~HDFWalkerOutput();

  /** Dump walker configurations to an HDF5 file (<fname>.config.h5).
   * @param w Walker configurations to dump
   * @param fname Base filename (without .config.h5 extension)
   * @param block Current block number
   * @return true if dump succeeded
   */
  bool dump(const WalkerConfigurations& w, const std::string& fname, int block);

private:
  /// PooledData buffer type for serialized walker particle positions
  using BufferType = PooledData<OHMMS_PRECISION>;
  /// Buffers for walker position data (index 0: local, index 1: gathered on rank 0)
  std::array<BufferType, 2> RemoteData;
  /// Buffers for walker weight data (index 0: local, index 1: gathered on rank 0)
  std::array<std::vector<QMCTraits::FullPrecRealType>, 2> RemoteDataW;
  /// Cached block number to prevent redundant buffer serialization within the same block
  int block;

  /** Write configurations using parallel HDF5 collective I/O and hyperslabs.
   * @param W Walker configurations
   * @param file_name Path to destination HDF5 file
   * @param nblock Current block number
   */
  void write_hdf_parallel(const WalkerConfigurations& W, const std::filesystem::path& file_name, int nblock);

  /** Gather walker configurations to rank 0 and write serially.
   * @param W Walker configurations
   * @param file_name Path to destination HDF5 file
   * @param nblock Current block number
   */
  void write_hdf_rank0(const WalkerConfigurations& W, const std::filesystem::path& file_name, int nblock);
};

} // namespace qmcplusplus
#endif
