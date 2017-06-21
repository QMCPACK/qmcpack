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
    
    
#ifndef QMCPLUSPLUS_POLYNOMIAL_CONSTRAINTS_H
#define QMCPLUSPLUS_POLYNOMIAL_CONSTRAINTS_H
#include "QMCWaveFunctions/OrbitalConstraintsBase.h"
#include "QMCWaveFunctions/Jastrow/PolyFunctor.h"
#include "QMCWaveFunctions/Jastrow/SplineFunctors.h"
#include "QMCWaveFunctions/ComboOrbital.h"

namespace qmcplusplus
{

struct PolyConstraints: public OrbitalConstraintsBase
{
  ///analytic functor
  typedef PolyFunctor<RealType> BasisGroupType;
  ///typedef for BasisGroupContainer
  typedef std::map<std::string,BasisGroupType*> BGContainerType;
  ///flag to tunr on/off spin-dependent term, always off
  bool IgnoreSpin;
  ///has all the basis groups
  BGContainerType BasisGroups;

  PolyConstraints(ParticleSet& p, TrialWaveFunction& psi, bool spinInd=true);
  ~PolyConstraints();

  void addOptimizables(OptimizableSetType& vlist);
  void resetParameters(OptimizableSetType& optVariables);

  OrbitalBase* createTwoBody();
  OrbitalBase* createOneBody(ParticleSet& source);
  void addExtra2ComboOrbital(ComboOrbital* jcombo) {}
  bool put(xmlNodePtr cur);
  void addSingleBasisPerSpecies(xmlNodePtr cur);
  void createBasisGroup(xmlNodePtr cur, const std::string& elementType, RealType rcut);
};

}
#endif
