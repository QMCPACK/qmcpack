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
    
    
#ifndef QMCPLUSPLUS_PADE_COMBO_CONSTRAINTS_H
#define QMCPLUSPLUS_PADE_COMBO_CONSTRAINTS_H
#include "QMCWaveFunctions/Jastrow/PadeFunctors.h"
#include "QMCWaveFunctions/OrbitalConstraintsBase.h"
#include "QMCWaveFunctions/ComboOrbital.h"
#include "QMCWaveFunctions/Jastrow/NumericalJastrowFunctor.h"

namespace qmcplusplus
{

struct PadeConstraints: public OrbitalConstraintsBase
{
  typedef PadeFunctor<RealType> FuncType;
  bool IgnoreSpin;
  RealType B;
  std::string ID_B;
  xmlNodePtr node_B;
  std::vector<FuncType*> FuncList;

  ~PadeConstraints();

  PadeConstraints(ParticleSet& p, TrialWaveFunction& psi, bool nospin=true);

  void resetParameters(const opt_variables_type& active);

  OrbitalBase* createTwoBody();
  OrbitalBase* createOneBody(ParticleSet& source);
  void addExtra2ComboOrbital(ComboOrbital* jcombo) {}
  bool put(xmlNodePtr cur);
};

struct ScaledPadeConstraints: public OrbitalConstraintsBase
{
  typedef ScaledPadeFunctor<RealType> FuncType;
  bool IgnoreSpin;
  RealType B;
  RealType C;
  std::string ID_B;
  std::string ID_C;
  std::vector<FuncType*> FuncList;

  ScaledPadeConstraints(ParticleSet& p, TrialWaveFunction& psi, bool nospin=true);

  ~ScaledPadeConstraints();

  void resetParameters(const opt_variables_type& active);

  OrbitalBase* createTwoBody();
  OrbitalBase* createOneBody(ParticleSet& source);
  void addExtra2ComboOrbital(ComboOrbital* jcombo) {}
  bool put(xmlNodePtr cur);
};

//struct PadeOnGridConstraints: public OrbitalConstraintsBase {
//  ///analytic functor
//  typedef PadeFunctor<RealType> InFuncType;
//  ///numerical functor
//  typedef NumericalJastrow<RealType> FuncType;
//  bool IgnoreSpin;
//  RealType B;
//  std::string ID;
//  std::vector<InFuncType*> InFuncList;
//  std::vector<FuncType*> FuncList;

//  ~PadeOnGridConstraints();
//  PadeOnGridConstraints(bool nospin=true):IgnoreSpin(nospin) {}
//  void apply();
//  void addOptimizables(VarRegistry<RealType>& outVars);
//  OrbitalBase* createTwoBody(ParticleSet& target);
//  OrbitalBase* createOneBody(ParticleSet& target, ParticleSet& source);
//  inline void addTwoBodyPart(ParticleSet& target, ComboOrbital* jcombo) {
//    OrbitalBase* j2 = createTwoBody(target);
//    if (j2) jcombo->Psi.push_back(j2);
//  }
//  bool put(xmlNodePtr cur);
//};

}
#endif
