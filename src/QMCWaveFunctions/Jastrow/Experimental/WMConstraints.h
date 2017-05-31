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
    
    
#ifndef QMCPLUSPLUS_WMFUNCTOR_SM_CONSTRAINTS_H
#define QMCPLUSPLUS_WMFUNCTOR_SM_CONSTRAINTS_H
#include "QMCWaveFunctions/OrbitalConstraintsBase.h"
#include "QMCWaveFunctions/Jastrow/WMFunctor.h"
#include "QMCWaveFunctions/Jastrow/SplineFunctors.h"
#include "QMCWaveFunctions/ComboOrbital.h"

namespace qmcplusplus
{

struct WMConstraints: public OrbitalConstraintsBase
{
  ///basis type
  typedef OptimizableFunctorBase<RealType>  BasisType;
  ///basis set type
  typedef std::vector<BasisType*> BasisSetType;
  ///analytic functor
  typedef ComboFunctor<RealType> InFuncType;
  ///numerical functor
  typedef CubicBsplineSingle<RealType> FuncType;
  ///flag to tunr on/off spin-dependent term, always off
  bool IgnoreSpin;
  ///cutoff radius
  RealType Rcut;

  std::map<std::string,BasisSetType*> myBasisSet;
  std::vector<InFuncType*> InFuncList;
  std::vector<FuncType*> FuncList;

  ~WMConstraints();
  WMConstraints(bool nospin=true):IgnoreSpin(nospin) {}
  void apply();
  void addOptimizables(VarRegistry<RealType>& outVars);
  OrbitalBase* createTwoBody(ParticleSet& target);
  OrbitalBase* createOneBody(ParticleSet& target, ParticleSet& source);
  inline void addTwoBodyPart(ParticleSet& target, ComboOrbital* jcombo)
  {
    OrbitalBase* twobody = createTwoBody(target);
    if (twobody)
      jcombo->Psi.push_back(twobody);
  }
  bool put(xmlNodePtr cur);
  void addBasisGroup(xmlNodePtr cur);
  InFuncType* createCorrelation(xmlNodePtr cur, BasisSetType* basis);
  InFuncType* createDefaultTwoBody(xmlNodePtr cur, const std::string& tname);
};

}
#endif
