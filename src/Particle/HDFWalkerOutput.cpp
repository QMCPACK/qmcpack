//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: D.C. Yang, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Cynthia Gu, zg1@ornl.gov, Oak Ridge National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


/** @file
 * @brief definition  of HDFWalkerOuput  and other support class and functions
 */
#include "HDFWalkerOutput.h"
#include "Utilities/IteratorUtility.h"
#include "OhmmsData/FileUtility.h"
#include "hdf/HDFVersion.h"
#include <numeric>
#include <iostream>
#include <sstream>
#include "Message/Communicate.h"
#include "Message/CommOperators.h"
#include "hdf/hdf_hyperslab.h"

namespace qmcplusplus
{
HDFWalkerOutput::HDFWalkerOutput(size_t num_ptcls, Communicate& c)
    : number_of_walkers_(0),
      number_of_particles_(num_ptcls),
      myComm(c)
{ block = -1; }

HDFWalkerOutput::~HDFWalkerOutput() = default;

bool HDFWalkerOutput::dump(const WalkerConfigurations& W, const std::string& fname, int nblock)
{
  const int wb = OHMMS_DIM * number_of_particles_;
  if (nblock > block)
  {
    RemoteData[0].resize(wb * W.getActiveWalkers());
    RemoteDataW[0].resize(W.getActiveWalkers());
    W.putConfigurations(RemoteData[0].data(), RemoteDataW[0].data());
    block = nblock;
  }

  auto& walker_offsets = W.getWalkerOffsets();
  number_of_walkers_   = walker_offsets[myComm.size()];

  std::filesystem::path file_name = fname;
  file_name.concat(hdf::config_ext);

#if defined(ENABLE_PHDF5)
  write_hdf_parallel(W, file_name, nblock);
#else
  write_hdf_rank0(W, file_name, nblock);
#endif

  return true;
}

void HDFWalkerOutput::write_hdf_parallel(const WalkerConfigurations& W,
                                         const std::filesystem::path& file_name,
                                         int nblock)
{
  hdf_archive hout(myComm, true);
  hout.create(file_name);
  HDFVersion cur_version;
  hout.write(cur_version.version, hdf::version);
  hout.push(hdf::main_state);
  hout.write(nblock, "block");

  auto& walker_offsets = W.getWalkerOffsets();

  hout.write(number_of_walkers_, hdf::num_walkers);
  { // write walker offset.
    // Though it is a small array, it needs to be written collectively in large scale runs.
    std::array<size_t, 1> gcounts{static_cast<size_t>(myComm.size()) + 1};
    std::array<size_t, 1> counts{0};
    std::array<size_t, 1> offsets{static_cast<size_t>(myComm.rank())};
    std::vector<int> myWalkerOffset;
    if (myComm.size() - 1 == myComm.rank())
    {
      counts[0] = 2;
      myWalkerOffset.push_back(walker_offsets[myComm.rank()]);
      myWalkerOffset.push_back(walker_offsets[myComm.size()]);
    }
    else
    {
      counts[0] = 1;
      myWalkerOffset.push_back(walker_offsets[myComm.rank()]);
    }
    hyperslab_proxy<std::vector<int>, 1> slab(myWalkerOffset, gcounts, counts, offsets);
    hout.write(slab, "walker_partition");
  }
  { // write walker configuration
    std::array<size_t, 3> gcounts{number_of_walkers_, number_of_particles_, OHMMS_DIM};
    std::array<size_t, 3> counts{W.getActiveWalkers(), number_of_particles_, OHMMS_DIM};
    std::array<size_t, 3> offsets{static_cast<size_t>(walker_offsets[myComm.rank()]), 0, 0};
    hyperslab_proxy<BufferType, 3> slab(RemoteData[0], gcounts, counts, offsets);
    hout.write(slab, hdf::walkers);
  }
  {
    std::array<size_t, 1> gcounts{number_of_walkers_};
    std::array<size_t, 1> counts{W.getActiveWalkers()};
    std::array<size_t, 1> offsets{static_cast<size_t>(walker_offsets[myComm.rank()])};
    hyperslab_proxy<std::vector<QMCTraits::FullPrecRealType>, 1> slab(RemoteDataW[0], gcounts, counts, offsets);
    hout.write(slab, hdf::walker_weights);
  }
  hout.close();
}

void HDFWalkerOutput::write_hdf_rank0(const WalkerConfigurations& W,
                                      const std::filesystem::path& file_name,
                                      int nblock)
{
  auto& walker_offsets = W.getWalkerOffsets();
  const int wb         = OHMMS_DIM * number_of_particles_;

  //gatherv to the master and master writes it, could use isend/irecv
  if (myComm.size() > 1)
  {
    std::vector<int> displ(myComm.size()), counts(myComm.size());
    for (int i = 0; i < myComm.size(); ++i)
    {
      counts[i] = wb * (walker_offsets[i + 1] - walker_offsets[i]);
      displ[i]  = wb * walker_offsets[i];
    }
    if (!myComm.rank())
      RemoteData[1].resize(wb * walker_offsets[myComm.size()]);
    myComm.gatherv(RemoteData[0], RemoteData[1], counts, displ);
    // update counts and displ for gathering walker weights
    for (int i = 0; i < myComm.size(); ++i)
    {
      counts[i] = (walker_offsets[i + 1] - walker_offsets[i]);
      displ[i]  = walker_offsets[i];
    }
    if (!myComm.rank())
      RemoteDataW[1].resize(walker_offsets[myComm.size()]);
    myComm.gatherv(RemoteDataW[0], RemoteDataW[1], counts, displ);
  }
  if (!myComm.rank())
  {
    hdf_archive hout;
    hout.create(file_name);
    HDFVersion cur_version;
    hout.write(cur_version.version, hdf::version);
    hout.push(hdf::main_state);
    hout.write(nblock, "block");

    hout.write(number_of_walkers_, hdf::num_walkers);
    hout.write(walker_offsets, "walker_partition");
    int buffer_id = (myComm.size() > 1) ? 1 : 0;
    {
      std::array<size_t, 3> gcounts{number_of_walkers_, number_of_particles_, OHMMS_DIM};
      hout.writeSlabReshaped(RemoteData[buffer_id], gcounts, hdf::walkers);
    }
    {
      std::array<size_t, 1> gcounts{number_of_walkers_};
      hout.writeSlabReshaped(RemoteDataW[buffer_id], gcounts, hdf::walker_weights);
    }
    hout.close();
  }
}

} // namespace qmcplusplus
