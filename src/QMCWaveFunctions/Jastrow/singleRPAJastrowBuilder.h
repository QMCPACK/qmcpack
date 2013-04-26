//////////////////////////////////////////////////////////////////
// (c) Copyright 2008-
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//   National Center for Supercomputing Applications &
//   Materials Computation Center
//   University of Illinois, Urbana-Champaign
//   Urbana, IL 61801
//   e-mail:
//   Tel:
//
// Supported by
//   National Center for Supercomputing Applications, UIUC
//   Materials Computation Center, UIUC
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
#ifndef QMCPLUSPLUS_SINGLE_RPA_JASTROW_H
#define QMCPLUSPLUS_SINGLE_RPA_JASTROW_H
#include "QMCWaveFunctions/OrbitalBuilderBase.h"
#include "LongRange/LRHandlerBase.h"
#include "QMCWaveFunctions/Jastrow/SplineFunctors.h"
#include "QMCWaveFunctions/Jastrow/LRBreakupUtilities.h"
#include "QMCWaveFunctions/Jastrow/OneBodyJastrowOrbital.h"
#include "Particle/DistanceTableData.h"
#include "Particle/DistanceTable.h"


namespace qmcplusplus
{
class ParticleSet;

struct singleRPAJastrowBuilder: public OrbitalBuilderBase
{
  typedef LRHandlerBase HandlerType;
  typedef CubicBspline<RealType,LINEAR_1DGRID,FIRSTDERIV_CONSTRAINTS> SplineEngineType;
  typedef CubicSplineSingle<RealType,SplineEngineType> FuncType;
  typedef LinearGrid<RealType> GridType;
  typedef OneBodyJastrowOrbital<FuncType> JneType;

  RealType Rs;
  RealType tlen;
  RealType Kc;
  double Kc_max;
  int indx,ng;
  RealType Rcut;
  string ID_Rs;
  string MyName;
  HandlerType* myHandler;
  FuncType* nfunc;
  ShortRangePartAdapter<RealType>* SRA;
  JneType* J1s;
  ParticleSet* sourcePtcl;

  singleRPAJastrowBuilder(ParticleSet& target, TrialWaveFunction& psi,
                          ParticleSet& source) :
    OrbitalBuilderBase(target,psi), sourcePtcl(&source), myHandler(0),J1s(0)
  {
    tlen = std::pow(3.0/4.0/M_PI*target.Lattice.Volume/ static_cast<RealType>(target.getTotalNum()) ,1.0/3.0);
//         indx = target.SK->KLists.ksq.size()-1;
//         Kc_max=std::pow(target.SK->KLists.ksq[indx],0.5);
//         Kc=target.Lattice.LR_kc;
//         Rcut=target.Lattice.LR_rc;
//         ng = 0;
    ng=source.getSpeciesSet().getTotalNum();
//         J1s = new JneType (source,target);
  }
  OrbitalBase* getOrbital();
  bool put(xmlNodePtr cur);
  bool put(xmlNodePtr cur, int addOrbital);


};

}
#endif
