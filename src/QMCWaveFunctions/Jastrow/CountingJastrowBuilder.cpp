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

#include "CountingJastrowBuilder.h"
#include "QMCWaveFunctions/Jastrow/CountingJastrowOrbital.h"
#include "QMCWaveFunctions/Jastrow/DiffCountingJastrowOrbital.h"
#include <iostream>

namespace qmcplusplus
{

CountingJastrowBuilder::CountingJastrowBuilder(ParticleSet& target, TrialWaveFunction& psi):
  WaveFunctionComponentBuilder(target, psi)
{
  ClassName="CountingJastrowBuilder";
  NameOpt="0";
  TypeOpt="unknown";
  Jastfunction="unknown";
  SourceOpt=targetPtcl.getName();
  SpinOpt="no";
}

template<class precision, template<class> class RegionType>
class JastrowTypeHelper
{
  public:
    using rft = RegionType<precision>;
    using CJOrbitalType = CountingJastrowOrbital<rft>;
    using DiffCJOrbitalType = DiffCountingJastrowOrbital<rft>;
}



template<template<class> class RegionType>
bool CountingJastrowBuilder::createCJ(xmlNodePtr cur)
{
  ReportEngine PRE(ClassName,"createCJ(xmlNodePtr)");
  using RegionType = RegionType<RT>;
  using FunctorType = typename RegionType::FunctorType;
  using CJOrbitalType = typename JastrowTypeHelper<RT, RegionType>::CJOrbitalType;
  using DiffCJOrbitalType = typename JastrowTypeHelper<RT, RegionType>::DiffCJOrbitalType;

  SpeciesSet& species(targetPtcl.getSpeciesSet());

  auto *CJ = new CJOrbitalType(targetPtcl);
  auto *dCJ = new DiffCJOrbitalType(targetPtcl);
  auto *CR = new RegionType();

  // standard child loop
  cur = cur->xmlChildrenNode;
  while(cur != NULL)
  {
    std::string cname((const char*) cur->name);
    // create counting region, populate with functions
    bool opt_R = true, opt_F = true, opt_G = true;
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
        oAttrib2.add(id,"id");
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
      }
    }
    cur = cur->next;
  }

  CJ->addRegion(CR,F,G, opt_R, opt_G, opt_F)
  dCJ->addRegion(CR);

  CJ->setOptimizable(opt_R || opt_G || opt_F);
  CJ->dPsi = dCJ;

  std::string cjname = "CJ_"+Jastfunction;
  targetPsi.addOrbital(CJ,cjname.c_str());
  return true;
  
}

bool CountingJastrowBuilder::put(xmlNodePtr cur)
{
  typedef RealType RT;
  // typedefs
  typedef NormalizedGaussianRegion<RT> NormGaussRegionType;
  typedef SigmoidRegion<RT> SigmoidRegionType;

  OhmmsAttributeSet oAttrib;
  oAttrib.add(regionOpt,"region");
  oAttrib.put(cur)
  if(regionOpt.find("normalized_gaussian") < regionOpt.size())
  {
    createCJ<NormGaussRegionType>(cur);
  }
  if(regionOpt.find("sigmoid") < regionOpt.size())
  {
    createCJ<SigmoidRegionType>(cur);
  }

  

  return true;
}

}
