//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    


#ifndef QMCPLUSPLUS_MULTICHAIN_HDF_IO_H
#define QMCPLUSPLUS_MULTICHAIN_HDF_IO_H

#include "QMCDrivers/MultiChain.h"

namespace qmcplusplus
{

/** Writes a set of walker configurations to an HDF5 file. */
class HDFMultiChainOutput
{

  bool AppendMode;
  ///number of times file has been written to
  int Counter;
  ///version number
  int m_version;
  ///string file name
  std::string h5file;

public:

  HDFMultiChainOutput(const std::string& fname, int count);
  ~HDFMultiChainOutput();

  bool get(MultiChain&);

  //template<class CT>
  //bool write(CT& anything) {
  //  anything.write(h_config);
  //  return true;
  //}

};

/** Reads a set of walker configurations from an HDF5 file. */

class HDFMultiChainInput
{

  ///id of the first set
  int FirstSet;
  ///id of the last set
  int LastSet;
  ///number of times file has been accesed
  int Counter;
  ///number of sets of walker configurations
  hsize_t NumSets;
  ///id for HDF5 file
  hid_t h_file;
  ///id for HDF5 main group
  hid_t h_config;

public:

  HDFMultiChainInput(const std::string& fname);
  ~HDFMultiChainInput();

  ///version number
  int m_version;
  ///string file name
  std::string h5file;

  bool put(MultiChain&);

  template<class CT>
  bool read(CT& anything)
  {
    anything.read(h_config);
    return true;
  }

};

}
#endif
