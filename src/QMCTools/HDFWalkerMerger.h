//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    



#ifndef QMCPLUSPLUS_WALKER_MERGE_IO_H
#define QMCPLUSPLUS_WALKER_MERGE_IO_H

#include <string>
#include "OhmmsPETE/OhmmsMatrix.h"
#include "Numerics/HDFSTLAttrib.h"

class Communicate;

namespace qmcplusplus
{

/** Writes a set of walker configurations to an HDF5 file. */

class HDFWalkerMerger
{

  ///the physical dimension
  int Dimension;
  ///the number of separate files
  int NumCopy;
  ///the number of configuration
  int NumConfig;
  ///the number of particles
  int NumPtcl;
  ///communicator
  Communicate* myComm;

  ///maxmimum number of walkers for any config
  hsize_t MaxNumWalkers;
  ///file root
  std::string FileRoot;

  ///numWalkerIn[node]->operator[](iconfig) returns the number of walkers
  std::vector<std::vector<hsize_t>*> numWalkersIn;
  /** offset of the walkers
   *
   * OffSet(NumCopy,*) is the total number of walkers for the merged
   * configuration.
   */
  Matrix<hsize_t>                    OffSet;

  /** summary is accumulated */
  std::vector<double>                Summary;

  /** write header information */
  void writeHeader(hid_t git);

  /** initialize sizes */
  void init();


public:

  HDFWalkerMerger(const std::string& froot, int ncopy);
  ~HDFWalkerMerger();

  /** set communicator */
  void setCommunicator(Communicate* c);

  void merge();
};

}
#endif
