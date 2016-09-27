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
  std::string ID_Rs;
  std::string MyName;
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
