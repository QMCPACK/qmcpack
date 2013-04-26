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
#ifndef QMCPLUSPLUS_WALKER_INPUT_COLLECT_H
#define QMCPLUSPLUS_WALKER_INPUT_COLLECT_H

#include "OhmmsData/HDFAttribIO.h"
#include "Particle/MCWalkerConfiguration.h"

namespace qmcplusplus
{

/** Reads a set of walker configurations from an HDF5 file. */

class HDFWalkerInputCollect
{

  hid_t fileID;
  int prevNContexts;
  int curNContexts;
  bool CollectMode;
  bool RewindMode;

  vector<hsize_t> OffSet;

  /** distribute nw among the processors */
  void distribute(hsize_t nw);

public:

  HDFWalkerInputCollect(const std::string& aroot);
  ~HDFWalkerInputCollect();

  /** open a hdf5 file to read in a collective mode
   * @param arrot root name
   * @return true if the file can be open
   */
  bool open(const string& aroot);
  /** read walkers
   * @param W walker set to which new walkers will be added
   * @param rollback number of blocks to read
   * @return true if successful
   */
  bool put(MCWalkerConfiguration& W, int rollback=1);

  /** read all the walkers
   *
   * A special put function to read entire configurations without sharing
   */
  bool putSingle(MCWalkerConfiguration& W);

  /** read walkers for rollback blocks
   * @param W walker set to which new walkers will be added
   * @param rollback number of blocks to read
   */
  bool rewind(MCWalkerConfiguration& W, int rollback);

  /** read walkers within the blocks [firstConf,lastConf)
   * @param W walker set to which new walkers will be added
   * @param firstConf first block to read
   * @param lastConf last block to read
   */
  bool read(MCWalkerConfiguration& W, int firstConf, int lastConf);


  /** read random state when RewindMode=false
   */
  void readRandomState();

  bool close();
};

}
#endif
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$
 ***************************************************************************/
