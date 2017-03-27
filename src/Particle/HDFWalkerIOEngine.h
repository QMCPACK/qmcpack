//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
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
