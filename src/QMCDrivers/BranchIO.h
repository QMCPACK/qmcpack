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

#ifndef QMCPLUSPLUS_BRANCHIO_H
#define QMCPLUSPLUS_BRANCHIO_H
namespace qmcplusplus
{
template<class SFNB>
class BranchIO
{
public:
  using RealType = typename SFNB::RealType;
  using BranchModeType = typename SFNB::BranchModeType;
  using IParamType = typename SFNB::IParamType;
  using VParamType = typename SFNB::VParamType;

  SFNB& ref;
  Communicate* myComm;
  BranchIO(SFNB& source, Communicate* c) : ref(source), myComm(c) {}

  bool write(const std::string& fname);
  bool read(const std::string& fname);
  void bcast_state();

  // TODO: names should be defined with the enum in SimpleFixedNodeBranch.h
  static std::vector<std::string> vParamName;
  static std::vector<std::string> iParamName;

  static void initAttributes();
};

extern template class BranchIO<SimpleFixedNodeBranch>;

} // namespace qmcplusplus
#endif
