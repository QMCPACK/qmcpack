//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Brett Van Der Goetz, bvdg@berkeley.edu, University of California at Berkeley
//
// File created by: Brett Van Der Goetz, bvdg@berkeley.edu, University of California at Berkeley
//////////////////////////////////////////////////////////////////////////////////////

#include "QMCWaveFunctions/Jastrow/CountingJastrowBuilder.h"
#include "QMCWaveFunctions/Jastrow/CountingJastrow.h"
#include "Utilities/ProgressReportEngine.h"
#include "OhmmsData/AttributeSet.h"
#include <iostream>

namespace qmcplusplus
{
CountingJastrowBuilder::CountingJastrowBuilder(ParticleSet& target, TrialWaveFunction& psi)
    : WaveFunctionComponentBuilder(target, psi)
{
  ClassName = "CountingJastrowBuilder";
  NameOpt   = "0";
  TypeOpt   = "unknown";
  RegionOpt = "normalized_gaussian";
  SourceOpt = targetPtcl.getName();
}

bool CountingJastrowBuilder::createCJ(xmlNodePtr cur)
{
  ReportEngine PRE(ClassName, "createCJ(xmlNodePtr)");
  using RegionType    = NormalizedGaussianRegion<RealType>;
  using FunctorType   = typename RegionType::FunctorType;
  using CJOrbitalType = CountingJastrow<RegionType>;

  SpeciesSet& species(targetPtcl.getSpeciesSet());

  auto* C = new RegionType(targetPtcl);

  Matrix<RealType> F;
  std::vector<RealType> G;
  bool opt_C = true, opt_F = true, opt_G = true;
  bool debug_flag = false;
  int period = 0, seqlen = 0;

  // standard child loop
  cur = cur->xmlChildrenNode;
  while (cur != NULL)
  {
    std::string cname((const char*)cur->name);
    // create counting region, populate with functions
    if (cname == "region")
    {
      // read in opt_C option
      OhmmsAttributeSet oAttrib;
      std::string opt = "true";
      oAttrib.add(opt, "opt");
      oAttrib.put(cur);
      // add functions
      xmlNodePtr cur2 = cur->xmlChildrenNode;
      while (cur2 != NULL)
      {
        std::string cname2((const char*)cur2->name);
        if (cname2 == "function")
        {
          // read id
          std::string fid;
          OhmmsAttributeSet oAttrib2;
          oAttrib2.add(fid, "id");
          oAttrib2.put(cur2);
          // get functor, add to function
          FunctorType* func = new FunctorType(fid);
          func->put(cur2);
          C->addFunc(func, fid);
        }
        cur2 = cur2->next;
      }
      // call RegionType->put
      opt_C = (opt == "true" || opt == "yes");
      C->put(cur);
    }
    if (cname == "var")
    {
      OhmmsAttributeSet oAttrib;
      std::string namestr = "none";
      oAttrib.add(namestr, "name");
      oAttrib.put(cur);
      if (namestr == "F")
      {
        // read in opt_F option and form
        std::string form = "upper_triang";
        std::string opt  = "true";
        OhmmsAttributeSet rAttrib2;
        rAttrib2.add(opt, "opt");
        rAttrib2.add(form, "form");
        rAttrib2.put(cur);
        opt_F = (opt == "yes" || opt == "true");
        // read in F matrix
        if (form == "upper_triang")
        {
          // read in upper triangle, set dimension
          std::vector<RealType> F_utri;
          putContent(F_utri, cur);
          int Fdim = (std::sqrt(8 * F_utri.size() + 1) - 1) / 2;

          if (!(F_utri.size() == Fdim * (Fdim + 1) / 2))
          {
            std::ostringstream err;
            err << "CountingJastrow::put: F cannot be the upper-triangular component of a square matrix: "
                << F_utri.size() << " != " << Fdim * (Fdim + 1) / 2 << std::endl;
            APP_ABORT(err.str());
          }
          // set F from upper triangular elements
          F.resize(Fdim, Fdim);
          auto it = F_utri.begin();
          for (int I = 0; I < Fdim; ++I)
            for (int J = I; J < Fdim; ++J, ++it)
              F(I, J) = F(J, I) = (*it);
        }
        else if (form == "full_matrix")
          putContent(F, cur);
      }
      if (namestr == "G")
      {
        // read in opt_G
        OhmmsAttributeSet rAttrib2;
        std::string opt = "true";
        rAttrib2.add(opt, "opt");
        rAttrib2.put(cur);
        opt = (opt == "yes" || opt == "true");
        // read in G vector
        putContent(G, cur);
      }
      if (namestr == "debug")
      {
        // read in debug options
        OhmmsAttributeSet rAttrib2;
        debug_flag = true;
        rAttrib2.add(period, "period");
        rAttrib2.add(seqlen, "seqlen");
        rAttrib2.put(cur);
      }
    }
    cur = cur->next;
  }

  auto* CJ = new CJOrbitalType(targetPtcl, C, F, G);

  CJ->addDebug(debug_flag, seqlen, period);
  CJ->addOpt(opt_C, opt_G, opt_F);
  CJ->setOptimizable(opt_C || opt_G || opt_F);
  CJ->initialize();

  std::string cjname = "CJ_" + RegionOpt;
  targetPsi.addComponent(CJ, cjname.c_str());
  return true;
}

bool CountingJastrowBuilder::put(xmlNodePtr cur)
{
  OhmmsAttributeSet oAttrib;
  oAttrib.add(RegionOpt, "region");
  oAttrib.add(TypeOpt, "type");
  oAttrib.add(NameOpt, "name");
  oAttrib.put(cur);
  if (RegionOpt.find("normalized_gaussian") < RegionOpt.size())
  {
    createCJ(cur);
  }
  app_log() << "end CountingRegionOrbital::put" << std::endl;

  return true;
}

} // namespace qmcplusplus
