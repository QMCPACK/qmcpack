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
  RegionOpt="normalized_gaussian";
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
  bool opt_C = true, opt_F = true, opt_G = true;

  // standard child loop
  cur = cur->xmlChildrenNode;
  while(cur != NULL)
  {
    std::string cname((const char*) cur->name);
    // create counting region, populate with functions
    if(cname == "region")
    {
      // read in opt_C option
      OhmmsAttributeSet oAttrib;
      std::string opt = "true";
      oAttrib.add(opt, "opt");
      oAttrib.put(cur);
      opt_C = (opt == "true" || opt == "yes");
      // call RegionType->put
      CR->put(cur);
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
        FunctorType* func = new FunctorType(fid);
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
        // read in opt_F option and form
        std::string form = "upper_triang";
        std::string opt = "true";
        OhmmsAttributeSet rAttrib2;
        rAttrib2.add(opt,"opt");
        rAttrib2.add(form,"form");
        rAttrib2.put(cur2);
        opt_F = (opt == "yes" || opt == "true");
        // read in F matrix
        if(form == "upper_triang")
        {
          // read in upper triangle, set dimension
          std::vector<RealType> F_utri;
          putContent(F_utri,cur);
          int Fdim = (std::sqrt(8*F_utri.size() + 1) - 1)/2;

          if(!(F_utri.size() == Fdim*(Fdim+1)/2))
          {
            std::ostringstream err;
            err << "CountingJastrowOrbital::put: F cannot be the upper-triangular component of a square matrix: " << F_utri.size() << " != " << Fdim*(Fdim+1)/2 << std::endl;
            APP_ABORT(err.str());
          }
          // set F from upper triangular elements
          F->resize(Fdim, Fdim);
          auto it = F_utri.begin();
          for(int I = 0; I < Fdim; ++I)
            for(int J = I; J < Fdim; ++J, ++it)
              (*F)(I,J) = (*F)(J,I) = (*it);
        }
        else if (form == "full_matrix")
          putContent(*F, cur);
      }
      if(cname2 == "G")
      {
        // read in opt_G
        OhmmsAttributeSet rAttrib2;
        std::string opt = "true";
        rAttrib2.add(opt,"opt");
        rAttrib2.put(cur2);
        opt = (opt == "yes" || opt == "true");
        // read in G vector
        putContent(*G, cur2);
      }
      if(cname2 == "debug")
      {
        // read in debug options
        int period = 10000, seqlen = 10;
        OhmmsAttributeSet rAttrib2;
        rAttrib2.add(period,"period");
        rAttrib2.add(seqlen,"seqlen");
        rAttrib2.put(cur);
        // set debug options
        CJ->addDebug(seqlen, period);
        dCJ->addDebug(seqlen, period);
      }
    }
    cur = cur->next;
  }

  CJ->addRegion(CR,F,G, opt_C,opt_G,opt_F);
  dCJ->addRegion(CR);

  CJ->setOptimizable(opt_C || opt_G || opt_F);
  CJ->dPsi = dCJ;

  std::string cjname = "CJ_"+RegionOpt;
  targetPsi.addOrbital(CJ,cjname.c_str());
  return true;
  
}

bool CountingJastrowBuilder::put(xmlNodePtr cur)
{
  OhmmsAttributeSet oAttrib;
  oAttrib.add(RegionOpt,"region");
  oAttrib.add(TypeOpt,"type");
  oAttrib.add(NameOpt,"name");
  oAttrib.put(cur);
  if(RegionOpt.find("normalized_gaussian") < RegionOpt.size())
  {
    createCJ<NormalizedGaussianRegion>(cur);
  }
//  if(RegionOpt.find("sigmoid") < RegionOpt.size())
//  {
//    createCJ<SigmoidRegion>(cur);
//  }

  

  return true;
}

}
