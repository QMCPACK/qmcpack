//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Miguel Morales, moralessilva2@llnl.gov, Lawrence Livermore National Laboratory
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Miguel Morales, moralessilva2@llnl.gov, Lawrence Livermore National Laboratory
//////////////////////////////////////////////////////////////////////////////////////


#ifndef QMCPLUSPLUS_BACKFLOW_BUILDER_H
#define QMCPLUSPLUS_BACKFLOW_BUILDER_H

#include <map>
#include <cmath>
#include "Configuration.h"
#include "Numerics/OneDimGridBase.h"
#include "QMCWaveFunctions/Fermion/BackflowFunctionBase.h"
#include "LongRange/LRHandlerBase.h"

namespace qmcplusplus
{
class BackflowTransformation;
class Backflow_ee_kSpace;
template<class T>
struct BsplineFunctor;
template<class FT>
class Backflow_ee;

class BackflowBuilder
{
  using RealType    = BackflowFunctionBase::RealType;
  using HandlerType = LRHandlerBase;
  using GridType    = LinearGrid<RealType>;
  using PSetMap     = std::map<std::string, const std::unique_ptr<ParticleSet>>;

public:
  BackflowBuilder(ParticleSet& p, const PSetMap& pool);

  std::unique_ptr<BackflowTransformation> buildBackflowTransformation(xmlNodePtr cur);

  RealType cutOff;

private:
  ParticleSet& targetPtcl;
  const PSetMap& ptclPool;
  bool IgnoreSpin;
  RealType Rs;
  RealType Kc;
  RealType Rcut;
  bool OneBody;
  bool TwoBody;

  HandlerType* myHandler;

  std::unique_ptr<BackflowFunctionBase> addOneBody(xmlNodePtr cur);

  std::unique_ptr<BackflowFunctionBase> addTwoBody(xmlNodePtr cur);

  std::unique_ptr<BackflowFunctionBase> addRPA(xmlNodePtr cur);

  void makeShortRange_oneBody();

  void makeLongRange_oneBody();

  void makeShortRange_twoBody(xmlNodePtr cur, Backflow_ee<BsplineFunctor<RealType>>* tbf, std::vector<int>& offsets);

  void makeLongRange_twoBody(xmlNodePtr cur, Backflow_ee_kSpace* tbf, std::vector<int>& offsets);
};

} // namespace qmcplusplus

#endif
