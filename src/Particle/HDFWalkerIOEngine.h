//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@intel.com, Intel Inc.
//                    Jeremy McMinnis, jmcminis@gmail.com, Navar Inc.
//
// File created by: Jeongnim Kim, jeongnim.kim@intel.com, Intel Inc.
//////////////////////////////////////////////////////////////////////////////////////
    
    



/** @file HDFWalkerIOEngine.h
 * @brief definition  of reader/writer of MCWalkerConfiguration
 */
#ifndef QMCPLUSPLUS_WALKERPACKEDIO_H
#define QMCPLUSPLUS_WALKERPACKEDIO_H

#include <Configuration.h>
#include <io/hdf_archive.h>

class Communicate;
namespace qmcplusplus
{
class MCWalkerConfiguration;
/** IO engine for Walker :  read/write the positions of all the walkers.
 *
 * This class is to optimize writing walker configurations for restart and
 * for variational optimizations.
 */
struct HDFWalkerIOEngine
{

  ///reference to the walkers
  MCWalkerConfiguration& W;
  ///if true, the content is replaced
  bool replace;
  ///transfer mode
  hid_t xfer_plist;

  HDFWalkerIOEngine(MCWalkerConfiguration& a, bool reuse=false);

  void read(hid_t grp, const char* name);

  //read collectively
  void readAll(hid_t grp, const char* name, Communicate* myComm);
};
}
#endif
/***************************************************************************
 * $RCSfile$   $Author: jnkim $
 * $Revision: 2303 $   $Date: 2007-11-19 13:25:50 -0600 (Mon, 19 Nov 2007) $
 * $Id: HDFWalkerOutput.cpp 2303 2007-11-19 19:25:50Z jnkim $
 ***************************************************************************/
