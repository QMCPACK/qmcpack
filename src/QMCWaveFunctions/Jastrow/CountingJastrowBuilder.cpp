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
  using CountingFunctorType = RegionType<RT>;
  using CJOrbitalType = typename JastrowTypeHelper<RT, RegionType>::CJOrbitalType;
  using DiffCJOrbitalType = typename JastrowTypeHelper<RT, RegionType>::DiffCJOrbitalType;

  SpeciesSet& species(targetPtcl.getSpeciesSet());



  // not sure what this is for
  //std::string init_mode("0");
  //{
  //  OhmmsAttributeSet hAttrib;
  //  hAttrib.add(init_mode,"init");
  //  hAttrib.put(cur);
  //}

  // standard child loop
  cur = cur->xmlChildrenNode;
  while(cur != NULL)
  {
    std::string cname((const char*) cur->name);
    
    // put() logic to instantiate and connect the structure of 
    // orbital/region/functor classes together

    // read debug options <var="debug" seqlen="" period="" /> under jastrow
    // if(...)
    //   CJ->addDebug(seqlen, period);
    // put these into the wavefunction somehow?
  
    // auto *CR = new CountingRegionType(targetPtcl)
    // iterate through and initialize the functors
    // for cname = ...
    //   if cname=="functor" "function"
    //     auto *CF = new CountingFunctorType(params, opt_flags)
           // relegate put to functors because these are really short and numerical 
           // and depends on the functor and this is otherwise clutter
    //     CF->put(cur);
    //     CR->addFunc(CF,...) 

    cur = cur->next;
  }

  // auto *CJ = new CJOrbitalType(targetPtcl);
  // CJ->addRegion(CR,F,G, opt_R, opt_G, opt_F)
  // CJ->setOptimizable(opt_R || opt_G || opt_F);
  // auto *dCJ = new DiffCJOrbitalType(targetPtcl);
  // dCJ->addRegion(CR);

  CJ->dPsi = dCJ;
  std::string cjname = "CJ_"+Jastfunction;
  targetPsi.addOrbital(CJ,cjname.c_str());
  return true;
  
}

bool CountingJastrowBuilder::put(xmlNodePtr cur)
{
  // template class by region type
  OhmmsAttributeSet oAttrib;
  oAttrib.add(regionOpt,"region");
  oAttrib.put(cur)
  if(regionOpt.find("normalized_gaussian") < regionOpt.size())
  {
    // initialize a CountingRegion<...> region using hardcoded functor types
    // NormalizedGaussianRegion<type> CR = ...
    // instantiate template of CountingRegion based on region node pointer
    // instantiate template based on NormalizedCountingRegion<GaussianFunctor>
    //   createCJ< >
  }
  if(regionOpt.find("sigmoid") < regionOpt.size())
  {
    // initialize a CountingRegion<...> region using hardcoded functor types
    // instantiate template of CountingRegion after finding <region /> node pointer
    // instantiate template based on the CountingRegion<FermiDiracFunctor>
    //   createCJ< >
  }

  

  return true;
}

}
