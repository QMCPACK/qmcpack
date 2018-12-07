//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Luke Shulenburger, lshulen@sandia.gov, Sandia National Laboratories
//
// File created by: Luke Shulenburger, lshulen@sandia.gov, Sandia National Laboratories
//////////////////////////////////////////////////////////////////////////////////////

#include "QMCWaveFunctions/Jastrow/CountingJastrowBuilder.h"
#include "QMCWaveFunctions/Jastrow/CountingJastrowOrbital.h"
#include "QMCWaveFunctions/Jastrow/DiffCountingJastrowOrbital.h"
#include "Utilities/ProgressReportEngine.h"
#include <iostream>

namespace qmcplusplus
{

CountingJastrowBuilder::CountingJastrowBuilder(ParticleSet& target, TrialWaveFunction& psi):
  WaveFunctionComponentBuilder(target, psi)
  {
  ClassName="CountingJastrowBuilder";
  NameOpt="0";
  TypeOpt="unknown";
  RegionOpt="unknown";
  SourceOpt=targetPtcl.getName();
}

template<class precision, template<class> class CountingRegionType>
class CountingJastrowTypeHelper
{
  public:
    using rft = CountingRegionType<precision>;
    using CJOrbitalType = CountingJastrowOrbital<rft>;
    using DiffCJOrbitalType = DiffCountingJastrowOrbital<rft>;
};



template<template<class> class CountingRegionType>
bool CountingJastrowBuilder::createCJ(xmlNodePtr cur)
{
  ReportEngine PRE(ClassName,"createCJ(xmlNodePtr)");
  using RegionType = CountingRegionType<RT>;
  using FunctorType = typename RegionType::FunctorType;
  using CJOrbitalType = typename CountingJastrowTypeHelper<RT,CountingRegionType>::CJOrbitalType;
  using DiffCJOrbitalType = typename CountingJastrowTypeHelper<RT,CountingRegionType>::DiffCJOrbitalType;

  SpeciesSet& species(targetPtcl.getSpeciesSet());

  auto *CJ = new CJOrbitalType(targetPtcl);
  auto *dCJ = new DiffCJOrbitalType(targetPtcl);
  auto *CR = new RegionType(targetPtcl);

  Matrix<RealType>* F = new Matrix<RealType>();
  std::vector<RealType>* G = new std::vector<RealType>();
  bool opt_R = true, opt_F = true, opt_G = true;

  // standard child loop
  cur = cur->xmlChildrenNode;
  while(cur != NULL)
  {
    std::string cname((const char*) cur->name);
    // create counting region, populate with functions
    if(cname == "region")
    {
      // read in opt_R option
      OhmmsAttributeSet oAttrib;
      oAttrib.add(opt_R,"opt");
      oAttrib.put(cur);
      // add functions
      xmlNodePtr cur2 = cur->xmlChildrenNode;
      std::string cname2((const char*) cur2->name);
      if(cname2 == "function")
      {
        // read id
        std::string fid;
        OhmmsAttributeSet oAttrib2;
        oAttrib2.add(fid,"id");
        oAttrib2.put(cur2);

        // get functor, add to function
        FunctorType* func = new FunctorType();
        func->put(cur2);
        CR->addFunc(func, fid);
      }
    }
    if(cname == "var")
    {
      xmlNodePtr cur2 = cur->xmlChildrenNode;
      std::string cname2((const char*) cur2->name);
      if(cname2 == "F")
      {
        // input F matrix
      }
      if(cname2 == "G")
      {
        // input G vector
      }
      if(cname2 == "debug")
      {
        int period = 10000, seqlen = 10;
        OhmmsAttributeSet oAttrib;
        oAttrib.add(period,"period");
        oAttrib.add(seqlen,"seqlen");
        oAttrib.put(cur);
        CJ->addDebug(seqlen, period);
        dCJ->addDebug(seqlen, period);
      }
    }
    cur = cur->next;
  }

  CJ->addRegion(CR,F,G, opt_R,opt_G,opt_F);
  dCJ->addRegion(CR);

  CJ->setOptimizable(opt_R || opt_G || opt_F);
  CJ->dPsi = dCJ;

  std::string cjname = "CJ_"+RegionOpt;
  targetPsi.addOrbital(CJ,cjname.c_str());
  return true;
  
}

bool CountingJastrowBuilder::put(xmlNodePtr cur)
{
  // typedefs
//  typedef NormalizedGaussianRegion<RT> NormGaussRegionType;
//  typedef SigmoidRegion<RT> SigmoidRegionType;

  OhmmsAttributeSet oAttrib;
  oAttrib.add(RegionOpt,"region");
  oAttrib.add(TypeOpt,"type");
  oAttrib.add(NameOpt,"name");
  oAttrib.put(cur);
  if(RegionOpt.find("normalized_gaussian") < RegionOpt.size())
  {
    createCJ<NormalizedGaussianRegion>(cur);
  }
  if(RegionOpt.find("sigmoid") < RegionOpt.size())
  {
    createCJ<SigmoidRegion>(cur);
  }

  

  return true;
}

}
