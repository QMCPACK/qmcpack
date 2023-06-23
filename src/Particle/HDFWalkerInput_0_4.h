//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


#ifndef QMCPLUSPLUS_WALKER_INPUT_V0_4_H
#define QMCPLUSPLUS_WALKER_INPUT_V0_4_H

#include "hdf/HDFVersion.h"
#include "OhmmsData/AttributeSet.h"
#include "WalkerConfigurations.h"
#include <filesystem>
#include <stack>

class Communicate;

namespace qmcplusplus
{
struct HDFWalkerInput_0_4
{
  struct IOInfo
  {
    bool parallel;
    bool collected;
    int rank;
    int nprocs;
    hid_t h_plist;
    HDFVersion version;
    inline IOInfo() { reset(); }
    inline void reset()
    {
      parallel  = false;
      collected = false;
      rank      = 0;
      nprocs    = 1;
    }
  };

  /// reference to the list of walker configurations to be read from file
  WalkerConfigurations& wc_list_;
  /// number of particles
  const size_t num_ptcls_;
  //pointer to the communicator
  Communicate* myComm;
  //current version this class supports
  HDFVersion cur_version;
  //information of the input files
  IOInfo i_info;
  //propoery list to handle parallel io
  hid_t h_plist;
  //the last file which was processed, extension removed
  std::filesystem::path FileName_noext;
  //list of files to process, extensions removed
  std::stack<std::filesystem::path> FileStack;

  /** constructor
   * @param wc_list target walker configurations
   * @param num_ptcls the number of particles in each walker
   * @param c communicator
   * @param v version
   */
  HDFWalkerInput_0_4(WalkerConfigurations& wc_list, size_t num_ptcls, Communicate* c, const HDFVersion& v);
  ~HDFWalkerInput_0_4();

  /** read walkers
   * @param W walker set to which new walkers will be added
   * @param rollback number of blocks to read
   * @return true if successful
   */
  bool put(xmlNodePtr cur);

  /** check options from xml */
  void checkOptions(xmlNodePtr cur);

  /** read walkers for small number of MPI tasks */
  bool read_hdf5(const std::filesystem::path& h5name);
  /** read walkers. Master reads and scatter the walkers */
  bool read_hdf5_scatter(const std::filesystem::path& h5name);
  /** read walkers using PHDF5 */
  bool read_phdf5(const std::filesystem::path& h5name);
};

} // namespace qmcplusplus
#endif
