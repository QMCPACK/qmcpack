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
#include "QMCWaveFucntions/Jastrow/CountingJastrowOrbital.h"
#include "QMCWaveFucntions/Jastrow/DiffCountingJastrowOrbital.h"
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
  int taskid = (targetPsi.is_manager())?targetPsi.getGroupID():-1;



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
    // put these into the wavefunction somehow?
  
    // auto *CR = new CountingRegionType(num_els)
    // iterate through and initialize the functors
    // for cname = ...
    //   if cname=="functor" "function"
    //     auto *CF = new CountingFunctorType(params, opt_flags)
    //     CF->put(cur); // because these are really short and numerical and this would otherwise be clutter
    //     CR->addFunc(CF,...) 

    // Not sure what this is for.
    //if(qmc_common.io_node)
    //{
    //	char fname[32];
    //	if(qmc_common.mpi_groups>1)
    //	  sprintf(fname,"J2.%s.g%03d.dat",pairType.c_str(),taskid);
    //	else
    //	  sprintf(fname,"J2.%s.dat",pairType.c_str());
    //	std::ofstream os(fname);
    //	print(*functor, os);	
    //}
    cur = cur->next;
  }

  // auto *CJ = new CJOrbitalType(targetPtcl, taskid);
  // CJ->add_linear_parameters(F,G)
  // CJ->add_debug_parameters(seqlen, period);
  // CJ->addRegion(CR);
  // auto *dCJ = new DiffCJOrbitalType(targetPtcl);
  // dCJ->add_opt_parameters(F, G, opt_G, opt_F) // need these for Counting function derivatives
  // dCJ->add_debug_parameters(seqlen, period);
  // dCJ->addRegion(CR);

  CJ->dPsi = dCJ;
  std::string cjname = "CJ_"+Jastfunction;
  targetPsi.addOrbital(CJ,cjname.c_str());
  CJ->setOptimizable(true);
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
