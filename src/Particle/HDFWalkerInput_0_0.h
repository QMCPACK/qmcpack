//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    



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
  std::string FileName;
public:
  /** constructor
   * @param W reference configuration to fill in
   * @param c communicator
   * @param aroot fileroot
   */
  HDFWalkerInput_0_0(MCWalkerConfiguration& W, const std::string& aroot);

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
