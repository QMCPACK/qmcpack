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
#ifndef QMCPLUSPLUS_WM_JASTROW_BUILDER_H
#define QMCPLUSPLUS_WM_JASTROW_BUILDER_H
#include "QMCWaveFunctions/OrbitalBuilderBase.h"
#include "QMCWaveFunctions/Jastrow/WMFunctor.h"
#include "QMCWaveFunctions/Jastrow/SplineFunctors.h"

namespace qmcplusplus
{

/** JastrowBuilder using Pade functor
 *
 * To replace PadeConstraints
 */
class WMJastrowBuilder: public OrbitalBuilderBase
{

public:
  WMJastrowBuilder(ParticleSet& target, TrialWaveFunction& psi, ParticleSet* source=0);
  bool put(xmlNodePtr cur);

private:
  ParticleSet* sourcePtcl;
  bool addOneBody(xmlNodePtr cur);
  bool addTwoBody(xmlNodePtr cur);

  typedef OptimizableFunctorBase InFuncType;
  ///spline engine
  typedef CubicBspline<RealType,LINEAR_1DGRID,FIRSTDERIV_CONSTRAINTS> SplineEngineType;
  ///numerical functor
  typedef CubicSplineSingle<RealType,SplineEngineType> OutFuncType;
  ///derivative functor type
  typedef WMFunctorSet<RealType> DerivFuncType;
  /** class to define a basisGroup to represent a radial function
  */
  struct BasisGroupType
  {
    InFuncType* In_;
    DerivFuncType* Deriv_;
    OutFuncType* Out_;
    int NumGridPoints;
    RealType Rcut;
    RealType Cusp;
    BasisGroupType():In_(0),Deriv_(0),Out_(0),NumGridPoints(81),Rcut(8.0),Cusp(0.0) {}
    void setGrid(RealType rc, int npts)
    {
      Rcut=rc;
      NumGridPoints=npts;
    }
  };
  ///typedef for BasisGroupContainer
  typedef map<string,BasisGroupType*> BGContainerType;
  ///flag to tunr on/off spin-dependent term, always off
  bool IgnoreSpin;
  ///has all the basis groups
  BGContainerType BasisGroups;

  void addSingleBasisPerSpecies(xmlNodePtr cur);
  BasisGroupType* createBasisGroup(xmlNodePtr cur, const string& elementType);
  void add2BasisGroup(BasisGroupType* curBG, xmlNodePtr cur);
};
}

#endif
/***************************************************************************
 * $RCSfile$   $Author: jnkim $
 * $Revision: 2930 $   $Date: 2008-07-31 10:30:42 -0500 (Thu, 31 Jul 2008) $
 * $Id: PadeConstraints.h 2930 2008-07-31 15:30:42Z jnkim $
 ***************************************************************************/
