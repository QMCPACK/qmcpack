//////////////////////////////////////////////////////////////////
// (c) Copyright 2003-  by Jeongnim Kim
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//   National Center for Supercomputing Applications &
//   Materials Computation Center
//   University of Illinois, Urbana-Champaign
//   Urbana, IL 61801
//   e-mail: jnkim@ncsa.uiuc.edu
//
// Supported by
//   National Center for Supercomputing Applications, UIUC
//   Materials Computation Center, UIUC
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
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
  typedef map<string,BasisGroupType*> BGContainerType;
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
  void createBasisGroup(xmlNodePtr cur, const string& elementType, RealType rcut);
};

}
#endif
/***************************************************************************
 * $RCSfile$   $Author: jnkim $
 * $Revision: 1693 $   $Date: 2007-02-02 12:13:11 -0600 (Fri, 02 Feb 2007) $
 * $Id: AnyConstratints.h 1693 2007-02-02 18:13:11Z jnkim $
 ***************************************************************************/
