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
#include <utility>

namespace qmcplusplus {

  /** Writes a set of walker configurations to an HDF5 file. */
  class HDFWalkerOutput {
    ///if true, keep it in memory
    bool DoNotAppend;
    ///number of blocks for append
    int appended_blocks;
    /** number of walkers when a state is dumped
     *
     * When the number of walkers per state has changed, NumOfWalkers is used
     * to reallocate the hdf5 group.
     */
    int number_of_walkers;
    ///current number of backups
    int number_of_backups;
    ///current number of backups
    int max_number_of_backups;
    ///id for HDF5 file 
    hid_t h_file;
    ///handler for parallel I/O
    hid_t h_plist;
    ///id for HDF5 main group 
    hid_t h_state;
    ///id for debug HDF5 file 
    hid_t h_debug_file;
    ///communicator
    Communicate* myComm;
  public:
    ///save file name
    string FileName;

    ///constructor
    HDFWalkerOutput(MCWalkerConfiguration& W, const string& fname, Communicate* c);
    ///destructor
    ~HDFWalkerOutput();

    /** dump configurations
     * @param w walkers
     */
    bool dump(MCWalkerConfiguration& w);

    /** append configurations
     * @param w walkers
     * @param c counter
     */
    bool append(MCWalkerConfiguration& w);

    template<class CT>
    void write(CT& anything, bool doappend) 
    {
      anything.write(h_state,doappend);
    }

    inline hid_t getFileNode() const { return h_file;}
    inline hid_t getStateNode() const { return h_state;}

#if defined(HAVE_LIBHDF5)
    void open();
    inline void close()
    {
      H5Gclose(h_state); h_state=-1;
      H5Fclose(h_file); h_file=-1;
    }
#else
    inline void open(){}
    inline void close(){}
#endif

  private:
    /** create hdf5 file for a node and return hanlder for a file and configuration
     * @param fname name of a hdf5 file
     * @param w wakler ensemble
     */
    std::pair<hid_t,hid_t> createH5FileSingle(const string& fname, MCWalkerConfiguration& w);
    ///dump content to a group
    bool dumpSingle(hid_t gid, MCWalkerConfiguration& w);
  };

}
#endif
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
