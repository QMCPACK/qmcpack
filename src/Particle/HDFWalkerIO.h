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
#ifndef OHMMS_QMC_WALKER_IO_H
#define OHMMS_QMC_WALKER_IO_H

#include "OhmmsData/HDFAttribIO.h"

namespace ohmmsqmc {

  /** Writes a set of walker configurations to an HDF5 file. */

  class HDFWalkerOutput {

    bool AppendMode;
    ///number of times file has been written to
    int Counter;
    ///id for HDF5 file 
    hid_t h_file;
    ///id for HDF5 main group 
    hid_t h_config;

  public:

    HDFWalkerOutput(const string& fname, bool append=true);
    ~HDFWalkerOutput();
    bool get(MCWalkerConfiguration&);

  };

  /** Reads a set of walker configurations from an HDF5 file. */

  class HDFWalkerInput {

    ///number of times file has been accesed
    int Counter;
    ///number of sets of walker configurations
    hsize_t NumSets;
    ///id for HDF5 file 
    hid_t h_file;
    ///id for HDF5 main group 
    hid_t h_config;

  public:

    HDFWalkerInput(const string& fname);
    ~HDFWalkerInput();
    bool put(MCWalkerConfiguration&);
    bool put(MCWalkerConfiguration&, int ic);
  };

}
#endif
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
