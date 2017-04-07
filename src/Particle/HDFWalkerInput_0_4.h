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

#include "HDFVersion.h"
#include "OhmmsData/AttributeSet.h"
#include <stack>

class Communicate;

namespace qmcplusplus
{

class MCWalkerConfiguration;

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
    inline IOInfo()
    {
      reset();
    }
    inline void reset()
    {
      parallel=false;
      collected=false;
      rank=0;
      nprocs=1;
    }
  };

  //reference to target walker configuration
  MCWalkerConfiguration& targetW;
  //pointer to the communicator
  Communicate* myComm;
  //current version this class supports
  HDFVersion cur_version;
  //information of the input files
  IOInfo i_info;
  //propoery list to handle parallel io
  hid_t h_plist;
  //the last file which was processed
  std::string FileName;
  //list of files to process
  std::stack<std::string> FileStack;

  /** constructor
   * @param W target MCWalkerConfiguration
   * @param c communicator
   * @param v version
   */
  HDFWalkerInput_0_4(MCWalkerConfiguration& W, Communicate* c, const HDFVersion& v);
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
  bool read_hdf5( std::string h5name);
  /** read walkers. Master reads and scatter the walkers */
  bool read_hdf5_scatter( std::string h5name);
  /** read walkers using PHDF5 */
  bool read_phdf5( std::string h5name);
  /**  for adios validation */
  bool read_adios(xmlNodePtr cur);

};

}
#endif
