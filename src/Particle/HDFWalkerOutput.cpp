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
#include "mpi/collectives.h"
#include "hdf/hdf_hyperslab.h"

namespace qmcplusplus
{
/* constructor
 * @param W walkers to operate on
 * @param aroot the root file name
 * @param c communicator
 *
 * Create the HDF5 file "aroot.config.h5" for output.
 * Constructor is responsible to create a hdf5 file and create basic groups
 * so that subsequent write can utilize existing dataspaces.
 * HDF5 contains
 * - state_0
 *   -- block
 *   -- number of walkers
 *   -- walker_partition
 *   -- walker_weights
 *   -- walkers
 * - config_collection
 *   -- NumOfConfigurations current count of the configurations
 *   -- config0000 current configuration
 *   -- config#### configuration for each block
 * Other classes can add datasets as needed.
 * open/close functions are provided so that the hdf5 file is closed during
 * the life time of this object. This is necessary so that failures do not lead
 * to unclosed hdf5.
 */
HDFWalkerOutput::HDFWalkerOutput(size_t num_ptcls, const std::string& aroot, Communicate* c)
    : appended_blocks(0),
      number_of_walkers_(0),
      number_of_particles_(num_ptcls),
      myComm(c),
      currentConfigNumber(0),
      RootName(aroot)
{
  block = -1;
}

/** Destructor writes the state of random numbers and close the file */
HDFWalkerOutput::~HDFWalkerOutput() = default;

/** Write the set of walker configurations to the HDF5 file.
 * @param W set of walker configurations
 *
 * Dump is for restart and do not preserve the file
 * - version
 * - state_0
 *  - block (int)
 *  - number_of_walkers (int)
 *  - walker_partition (int array)
 *  - walkers (nw,np,3)
 */
bool HDFWalkerOutput::dump(const WalkerConfigurations& W, int nblock)
{
  std::filesystem::path FileName = myComm->getName();
  FileName.concat(hdf::config_ext);
  //rotate files
  //if(!myComm->rank() && currentConfigNumber)
  //{
  //  std::ostringstream o;
  //  o<<myComm->getName()<<".b"<<(currentConfigNumber%5) << hdf::config_ext;
  //  rename(prevFile.c_str(),o.str().c_str());
  //}

  //try to use collective
  hdf_archive dump_file(myComm, true);
  dump_file.create(FileName);
  HDFVersion cur_version;
  dump_file.write(cur_version.version, hdf::version);
  dump_file.push(hdf::main_state);
  dump_file.write(nblock, "block");

  write_configuration(W, dump_file, nblock);
  dump_file.close();

  currentConfigNumber++;
  prevFile = FileName;
  return true;
}

void HDFWalkerOutput::write_configuration(const WalkerConfigurations& W, hdf_archive& hout, int nblock)
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
  number_of_walkers_ = walker_offsets[myComm->size()];
  hout.write(number_of_walkers_, hdf::num_walkers);

  if (hout.is_parallel())
  {
    { // write walker offset.
      // Though it is a small array, it needs to be written collectively in large scale runs.
      std::array<size_t, 1> gcounts{static_cast<size_t>(myComm->size()) + 1};
      std::array<size_t, 1> counts{0};
      std::array<size_t, 1> offsets{static_cast<size_t>(myComm->rank())};
      std::vector<int> myWalkerOffset;
      if (myComm->size() - 1 == myComm->rank())
      {
        counts[0] = 2;
        myWalkerOffset.push_back(walker_offsets[myComm->rank()]);
        myWalkerOffset.push_back(walker_offsets[myComm->size()]);
      }
      else
      {
        counts[0] = 1;
        myWalkerOffset.push_back(walker_offsets[myComm->rank()]);
      }
      hyperslab_proxy<std::vector<int>, 1> slab(myWalkerOffset, gcounts, counts, offsets);
      hout.write(slab, "walker_partition");
    }
    { // write walker configuration
      std::array<size_t, 3> gcounts{number_of_walkers_, number_of_particles_, OHMMS_DIM};
      std::array<size_t, 3> counts{W.getActiveWalkers(), number_of_particles_, OHMMS_DIM};
      std::array<size_t, 3> offsets{static_cast<size_t>(walker_offsets[myComm->rank()]), 0, 0};
      hyperslab_proxy<BufferType, 3> slab(RemoteData[0], gcounts, counts, offsets);
      hout.write(slab, hdf::walkers);
    }
    {
      std::array<size_t, 1> gcounts{number_of_walkers_};
      std::array<size_t, 1> counts{W.getActiveWalkers()};
      std::array<size_t, 1> offsets{static_cast<size_t>(walker_offsets[myComm->rank()])};
      hyperslab_proxy<std::vector<QMCTraits::FullPrecRealType>, 1> slab(RemoteDataW[0], gcounts, counts, offsets);
      hout.write(slab, hdf::walker_weights);
    }
  }
  else
  { //gaterv to the master and master writes it, could use isend/irecv
    hout.write(walker_offsets, "walker_partition");
    if (myComm->size() > 1)
    {
      std::vector<int> displ(myComm->size()), counts(myComm->size());
      for (int i = 0; i < myComm->size(); ++i)
      {
        counts[i] = wb * (walker_offsets[i + 1] - walker_offsets[i]);
        displ[i]  = wb * walker_offsets[i];
      }
      if (!myComm->rank())
        RemoteData[1].resize(wb * walker_offsets[myComm->size()]);
      mpi::gatherv(*myComm, RemoteData[0], RemoteData[1], counts, displ);
      // update counts and displ for gathering walker weights
      for (int i = 0; i < myComm->size(); ++i)
      {
        counts[i] = (walker_offsets[i + 1] - walker_offsets[i]);
        displ[i]  = walker_offsets[i];
      }
      if (!myComm->rank())
        RemoteDataW[1].resize(walker_offsets[myComm->size()]);
      mpi::gatherv(*myComm, RemoteDataW[0], RemoteDataW[1], counts, displ);
    }
    int buffer_id = (myComm->size() > 1) ? 1 : 0;
    {
      std::array<size_t, 3> gcounts{number_of_walkers_, number_of_particles_, OHMMS_DIM};
      hout.writeSlabReshaped(RemoteData[buffer_id], gcounts, hdf::walkers);
    }
    {
      std::array<size_t, 1> gcounts{number_of_walkers_};
      hout.writeSlabReshaped(RemoteDataW[buffer_id], gcounts, hdf::walker_weights);
    }
  }
}
} // namespace qmcplusplus
