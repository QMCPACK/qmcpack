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
//
// Supported by 
//   National Center for Supercomputing Applications, UIUC
//   Materials Computation Center, UIUC
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
#ifndef QMCPLUSPLUS_WALKER_OUTPUT_H
#define QMCPLUSPLUS_WALKER_OUTPUT_H

#include "OhmmsData/HDFAttribIO.h"
#include "Particle/MCWalkerConfiguration.h"

namespace qmcplusplus {

  /** Writes a set of walker configurations to an HDF5 file. */
  class HDFWalkerOutput {
    ///number of times file has been written to
    int NumOfWalkers;
    ///id for HDF5 file 
    hid_t h_file;
    ///id for HDF5 main group 
    hid_t h_config;
    ///save file name
    string h5FileName;
  public:

    HDFWalkerOutput(MCWalkerConfiguration& W, const string& fname);
    ~HDFWalkerOutput();
    //bool get(MCWalkerConfiguration&);

    /** dump configurations
     * @param w walkers
     */
    bool dump(MCWalkerConfiguration& w);

    /** append configurations
     * @param w walkers
     * @param c counter
     */
    bool append(MCWalkerConfiguration& w, int c);

    template<class CT>
    void write(CT& anything) {
      anything.write(h_config,false);
    }

    /** return the file ID **/
    hid_t getFileID() { return h_file;}

    /** return the config_collection file ID **/
    hid_t getConfigID() { return h_config;}

#if defined(HAVE_LIBHDF5)
    inline void open()
    {
      h_file =  H5Fopen(h5FileName.c_str(),H5F_ACC_RDWR,H5P_DEFAULT);
      h_config = H5Gopen(h_file,"config_collection");
    }
    inline void close()
    {
      H5Gclose(h_config); h_config=-1;
      H5Fclose(h_file); h_file=-1;
    }
#else
    inline void open(){}
    inline void close(){}
#endif
  };

}
#endif
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
