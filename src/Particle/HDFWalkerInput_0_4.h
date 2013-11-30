//////////////////////////////////////////////////////////////////
// (c) Copyright 2005- by Jeongnim Kim
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//   Jeongnim Kim
//   National Center for Supercomputing Applications &
//   Materials Computation Center
//   University of Illinois, Urbana-Champaign
//   Urbana, IL 61801
//   e-mail: jnkim@ncsa.uiuc.edu
//
// Supported by
//   National Center for Supercomputing Applications, UIUC
//   Materials Computation Center, UIUC
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
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
  string FileName;
  //list of files to process
  stack<string> FileStack;

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
  bool read_hdf5(string h5name);
  /** read walkers. Master reads and scatter the walkers */
  bool read_hdf5_scatter(string h5name);
  /** read walkers using PHDF5 */
  bool read_phdf5(string h5name);
  /**  for adios validation */
  bool read_adios(xmlNodePtr cur);

};

}
#endif
/***************************************************************************
 * $RCSfile$   $Author: jnkim $
 * $Revision: 1575 $   $Date: 2007-01-04 09:36:37 -0600 (Thu, 04 Jan 2007) $
 * $Id: HDFWalkerInput_0_4.h 1575 2007-01-04 15:36:37Z jnkim $
 ***************************************************************************/
