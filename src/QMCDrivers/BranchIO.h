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

#ifndef QMCPLUSPLUS_BRANCHIO_H
#define QMCPLUSPLUS_BRANCHIO_H
namespace qmcplusplus
{
struct BranchIO
{
  using SFNB = SimpleFixedNodeBranch;
  using RealType = SFNB::RealType;
  using BranchModeType = SFNB::BranchModeType;
  using IParamType = SFNB::IParamType;
  using VParamType = SFNB::VParamType;

  SimpleFixedNodeBranch& ref;
  Communicate* myComm;
  BranchIO(SimpleFixedNodeBranch& source, Communicate* c) : ref(source), myComm(c) {}

  bool write(const std::string& fname);
  bool read(const std::string& fname);
  void bcast_state();

  // TODO: names should be defined with the enum in SimpleFixedNodeBranch.h
  static std::vector<std::string> vParamName;
  static std::vector<std::string> iParamName;

  static void initAttributes();
};
} // namespace qmcplusplus
#endif
