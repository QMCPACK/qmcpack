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
#ifndef QMCPLUSPLUS_WALKER_MERGE_IO_H
#define QMCPLUSPLUS_WALKER_MERGE_IO_H

#include <string>
#include "OhmmsPETE/OhmmsMatrix.h"
#include "Numerics/HDFSTLAttrib.h"

class Communicate;

namespace qmcplusplus {

  /** Writes a set of walker configurations to an HDF5 file. */

  class HDFWalkerMerger {

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
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
