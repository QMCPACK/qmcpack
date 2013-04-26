//////////////////////////////////////////////////////////////////
// (c) Copyright 2003- by Jeongnim Kim
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//   Jeongnim Kim
//   National Center for Supercomputing Applications &
//   Materials Computation Center
//   University of Illinois, Urbana-Champaign
//   Urbana, IL 61801
//   e-mail: jnkim@ncsa.uiuc.edu
//   Tel:    217-244-6319 (NCSA) 217-333-3324 (MCC)
//
// Supported by
//   National Center for Supercomputing Applications, UIUC
//   Materials Computation Center, UIUC
//   Department of Physics, Ohio State University
//   Ohio Supercomputer Center
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
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
  string h5file;

public:

  HDFMultiChainOutput(const string& fname, int count);
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

  HDFMultiChainInput(const string& fname);
  ~HDFMultiChainInput();

  ///version number
  int m_version;
  ///string file name
  string h5file;

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
/***************************************************************************
 * $RCSfile$   $Author: jnkim $
 * $Revision: 760 $   $Date: 2005-11-02 16:23:38 -0600 (Wed, 02 Nov 2005) $
 * $Id: MultiChainIO.h 760 2005-11-02 22:23:38Z jnkim $
 ***************************************************************************/
