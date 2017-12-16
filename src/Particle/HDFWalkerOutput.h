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

#include <Particle/MCWalkerConfiguration.h>
// #include <QMCDrivers/ForwardWalking/ForwardWalkingStructure.h>
#include <utility>
#include <io/hdf_archive.h>
#ifdef HAVE_ADIOS
#include <adios.h>
#ifdef ADIOS_VERIFY
#include "adios_read.h"
#include "adios_error.h"
#endif
#endif

namespace qmcplusplus
{


/** Writes a set of walker configurations to an HDF5 file. */
class HDFWalkerOutput
{
  ///if true, keep it in memory
  bool DoNotAppend;
  ///number of blocks for append
  int appended_blocks;
  /** number of walkers when a state is dumped
   *
   * When the number of walkers per state has changed, NumOfWalkers is used
   * to reallocate the hdf5 group.
   */
  int number_of_walkers;
  /** number of particles */
  int number_of_particles;
  ///current number of backups
  int number_of_backups;
  ///current number of backups
  int max_number_of_backups;
  ///communicator
  Communicate* myComm;
  int currentConfigNumber;
  ///rootname
  std::string RootName;
  std::string prevFile;
//     ///handle for the storeConfig.h5
//     hdf_archive fw_out;
public:
  ///constructor
  HDFWalkerOutput(MCWalkerConfiguration& W, const std::string& fname, Communicate* c);
  ///destructor
  ~HDFWalkerOutput();

#ifdef HAVE_ADIOS
  uint64_t get_group_size(MCWalkerConfiguration& W);

  bool adios_checkpoint(MCWalkerConfiguration& W, int64_t adios_handle, int block);
#ifdef ADIOS_VERIFY
  void adios_checkpoint_verify_variables(ADIOS_FILE* fp, const char* name, OHMMS_PRECISION* origin);
  void adios_checkpoint_verify_variables(ADIOS_FILE* fp, const char* name, int origin);
  void adios_checkpoint_verify(MCWalkerConfiguration& W,ADIOS_FILE* fp);
#endif
#endif

  /** dump configurations
   * @param w walkers
   */
  bool dump(MCWalkerConfiguration& w, int block);
//     bool dump(ForwardWalkingHistoryObject& FWO);

private:

  ///PooledData<T> is used to define the shape of multi-dimensional array
  typedef PooledData<OHMMS_PRECISION> BufferType;
  std::vector<Communicate::request> myRequest;
  std::vector<BufferType*> RemoteData;
  int block;

//     //define some types for the FW collection
//     typedef std::vector<ForwardWalkingData> FWBufferType;
//     std::vector<FWBufferType*> FWData;
//     std::vector<std::vector<int> > FWCountData;

  void write_configuration(MCWalkerConfiguration& W, hdf_archive& hout, int block);
};

}
#endif
