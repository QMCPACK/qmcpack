//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    


#include "QMCDrivers/SimpleFixedNodeBranch.h"
//#include <boost/archive/text_oarchive.hpp>

#ifdef HAVE_ADIOS
#include "adios_read.h"
#include "adios_error.h"
#endif

#ifndef QMCPLUSPLUS_BRANCHIO_H
#define QMCPLUSPLUS_BRANCHIO_H
namespace qmcplusplus
{
struct BranchIO
{
  typedef SimpleFixedNodeBranch::RealType RealType;
  typedef SimpleFixedNodeBranch::BranchModeType BranchModeType;
  typedef SimpleFixedNodeBranch::IParamType IParamType;
  typedef SimpleFixedNodeBranch::VParamType VParamType;

  SimpleFixedNodeBranch& ref;
  Communicate* myComm;
  BranchIO(SimpleFixedNodeBranch& source, Communicate* c): ref(source),myComm(c) {}
#ifdef HAVE_ADIOS
  int64_t get_Checkpoint_size();
  void adios_checkpoint(int64_t adios_handle);
#ifdef ADIOS_VERIFY
  void adios_checkpoint_verify(ADIOS_FILE* fp);
#endif
#endif

  bool write(const std::string& fname);
  bool read(const std::string& fname);
  bool read_adios(const std::string& fname);
  void bcast_state();

  static std::vector<std::string> vParamName;
  static std::vector<std::string> iParamName;

  static void initAttributes();
};
}
#endif


