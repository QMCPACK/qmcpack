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

#include <iostream>
#include "CountingJastrowBuilder.h"
#include "Utilities/ProgressReportEngine.h"
#include "OhmmsData/AttributeSet.h"
#include "QMCWaveFunctions/Jastrow/CountingJastrow.h"

namespace qmcplusplus
{
CountingJastrowBuilder::CountingJastrowBuilder(Communicate* comm, ParticleSet& target, ParticleSet& source)
    : WaveFunctionComponentBuilder(comm, target), SourcePtcl(&source)
{
  ClassName = "CountingJastrowBuilder";
  NameOpt   = "0";
  TypeOpt   = "Counting";
  RegionOpt = "voronoi";
  SourceOpt = SourcePtcl->getName();
}

CountingJastrowBuilder::CountingJastrowBuilder(Communicate* comm, ParticleSet& target)
    : WaveFunctionComponentBuilder(comm, target)
{
  ClassName  = "CountingJastrowBuilder";
  NameOpt    = "0";
  TypeOpt    = "Counting";
  RegionOpt  = "normalized_gaussian";
  SourceOpt  = "none";
  SourcePtcl = NULL;
}

std::unique_ptr<WaveFunctionComponent> CountingJastrowBuilder::createCJ(xmlNodePtr cur)
{
  ReportEngine PRE(ClassName, "createCJ(xmlNodePtr)");
  using RegionType    = CountingGaussianRegion;
  using FunctorType   = CountingGaussian;
  using CJOrbitalType = CountingJastrow<RegionType>;

  SpeciesSet& species(targetPtcl.getSpeciesSet());

  auto C = std::make_unique<RegionType>(targetPtcl);

  Matrix<RealType> F;
  bool opt_C = true, opt_F = true;
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
      std::string opt  = "true";
      std::string type = RegionOpt;
      oAttrib.add(opt, "opt");
      oAttrib.add(type, "type");
      oAttrib.put(cur);
      opt_C = (opt == "true" || opt == "yes");
      // add based on <function /> tags
      xmlNodePtr cur2 = cur->xmlChildrenNode;
      if (type == "normalized_gaussian")
      {
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
            auto func = std::make_unique<FunctorType>(fid);
            func->put(cur2);
            C->addFunc(std::move(func), fid);
          }
          cur2 = cur2->next;
        }
      }
      // add functions based on source
      else if (type == "voronoi")
      {
        RealType alpha = 1.0;
        while (cur2 != NULL)
        {
          std::string cname2((const char*)cur2->name);
          if (cname2 == "var")
          {
            // read id
            std::string vname;
            OhmmsAttributeSet oAttrib2;
            oAttrib2.add(vname, "name");
            oAttrib2.put(cur2);
            if (vname == "alpha")
            {
              putContent(alpha, cur2);
            }
          }
          cur2 = cur2->next;
        }
        // read in source
        if (SourcePtcl == NULL)
        {
          // quit with error - need a valid source
        }
        std::ostringstream os;
        // add a function for each source particle
        for (int i = 0; i < SourcePtcl->R.size(); ++i)
        {
          PosType gr = SourcePtcl->R[i];
          bool opt_g = opt_C && (i != 0);
          os.str("");
          os << "g" << i;
          std::string fid = os.str();
          auto func       = std::make_unique<FunctorType>(fid, alpha, gr, opt_g);
          C->addFunc(std::move(func), fid);
        }
        // set default F value to all zeroes with appropriate dimension
        if (F.size() == 0)
        {
          int Fdim = SourcePtcl->R.size();
          F.resize(Fdim, Fdim);
          for (int I = 0; I < Fdim; ++I)
            for (int J = 0; J < Fdim; ++J)
              F(I, J) = 0;
        }
      }
      // read in the counting region
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
            err << "CountingJastrow::put: F needs to be the upper-triangular component of a square matrix: "
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

        // transform the F matrix to put a zero in the lower-right corner

        RealType x = F(F.rows() - 1, F.cols() - 1);
        for (int I = 0; I < F.rows(); ++I)
          for (int J = 0; J < F.cols(); ++J)
            F(I, J) -= x;
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
  auto CJ = std::make_unique<CJOrbitalType>(targetPtcl, std::move(C), F, opt_C, opt_F);

  CJ->addDebug(debug_flag, seqlen, period);
  CJ->initialize();
  CJ->reportStatus(app_log());

  return CJ;
}

std::unique_ptr<WaveFunctionComponent> CountingJastrowBuilder::buildComponent(xmlNodePtr cur)
{
  OhmmsAttributeSet oAttrib;
  oAttrib.add(RegionOpt, "region");
  oAttrib.add(TypeOpt, "type");
  oAttrib.add(NameOpt, "name");
  oAttrib.put(cur);
  return createCJ(cur);
}

} // namespace qmcplusplus
