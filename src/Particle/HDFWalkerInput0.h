//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    



#ifndef QMCPLUSPLUS_WALKER_INPUT_FIRST_H
#define QMCPLUSPLUS_WALKER_INPUT_FIRST_H

#include "OhmmsData/HDFAttribIO.h"

namespace qmcplusplus
{

/** Reads a set of walker configurations from an HDF5 file. */

class HDFWalkerInput0
{

  ///id of the first set
  int FirstSet;
  ///id of the last set
  int LastSet;
  ///number of times file has been accesed
  int Counter;
  ///number of sets of walker configurations
  int NumSets;
  ///id for HDF5 file
  hid_t h_file;
  ///id for HDF5 main group
  hid_t h_config;

public:

  HDFWalkerInput0(const std::string& fname, int ipart=0, int nparts=1);
  ~HDFWalkerInput0();
  //int put(MCWalkerConfiguration&);
  bool put(MCWalkerConfiguration&, int ic);

  //bool append(MCWalkerConfiguration& w);
  bool append(MCWalkerConfiguration& w, int blocks);

  template<class CT>
  void read(CT& anything)
  {
    anything.read(h_config);
  }


  /** read the state of the number generator
   * @param restart if true, read the state
   */
  void getRandomState(bool restart);

};

}
#endif
