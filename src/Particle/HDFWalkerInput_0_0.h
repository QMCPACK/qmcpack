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
/** @file HDFWalkerInput_0_0.h
 * @brief definition of an input engine for walkers
 *
 * Supports hdf5 files with unspecified version.
 */
#ifndef QMCPLUSPLUS_WALKER_INPUT_0_0_H
#define QMCPLUSPLUS_WALKER_INPUT_0_0_H

#include "HDFVersion.h"
#include "OhmmsData/AttributeSet.h"

namespace qmcplusplus
{

class MCWalkerConfiguration;

/** Reads a set of walker configurations from an HDF5 file. */

class HDFWalkerInput_0_0
{
  ///supported version
  HDFVersion version;
  ///target MCWalkerConfiguration
  MCWalkerConfiguration& targetW;
  ///name of hdf file
  string FileName;
public:
  /** constructor
   * @param W reference configuration to fill in
   * @param c communicator
   * @param aroot fileroot
   */
  HDFWalkerInput_0_0(MCWalkerConfiguration& W, const string& aroot);

  ~HDFWalkerInput_0_0();
  /** process xml node
   *
   * Append walkers based on the options
   */
  bool put(xmlNodePtr cur);

  //template<class CT>
  //void read(CT& anything) {
  //  anything.read(h_config);
  //}

  ///** read the state of the number generator
  // * @param restart if true, read the state
  // */
  //void getRandomState(bool restart);

};

}
#endif
/***************************************************************************
 * $RCSfile$   $Author: jnkim $
 * $Revision: 1575 $   $Date: 2007-01-04 09:36:37 -0600 (Thu, 04 Jan 2007) $
 * $Id: HDFWalkerInput_0_0.h 1575 2007-01-04 15:36:37Z jnkim $
 ***************************************************************************/
