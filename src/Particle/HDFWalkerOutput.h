//////////////////////////////////////////////////////////////////
// (c) Copyright 2004- by Jeongnim Kim and Jordan Vincent
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
#ifndef QMCPLUSPLUS_WALKER_OUTPUT_H
#define QMCPLUSPLUS_WALKER_OUTPUT_H

#include "OhmmsData/HDFAttribIO.h"

namespace qmcplusplus {

  /** Writes a set of walker configurations to an HDF5 file. */

  class HDFWalkerOutput {

    bool AppendMode;
    ///number of times file has been written to
    int Counter;
    ///id for HDF5 file 
    hid_t h_file;
    ///id for HDF5 main group 
    hid_t h_config;
    ///id for the random number generator
    hid_t h_random;
  public:

    HDFWalkerOutput(const string& fname, bool append=false, int count=0);
    ~HDFWalkerOutput();
    bool get(MCWalkerConfiguration&);

    template<class CT>
    void write(CT& anything) {
      anything.write(h_config,AppendMode);
    }

    /** return the file ID **/
    hid_t getFileID() { return h_file;}

    /** return the config_collection file ID **/
    hid_t getConfigID() { return h_config;}
  };

}
#endif
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
