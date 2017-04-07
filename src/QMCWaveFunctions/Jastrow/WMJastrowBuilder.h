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
  typedef std::map<std::string,BasisGroupType*> BGContainerType;
  ///flag to tunr on/off spin-dependent term, always off
  bool IgnoreSpin;
  ///has all the basis groups
  BGContainerType BasisGroups;

  void addSingleBasisPerSpecies(xmlNodePtr cur);
  BasisGroupType* createBasisGroup(xmlNodePtr cur, const std::string& elementType);
  void add2BasisGroup(BasisGroupType* curBG, xmlNodePtr cur);
};
}

#endif
