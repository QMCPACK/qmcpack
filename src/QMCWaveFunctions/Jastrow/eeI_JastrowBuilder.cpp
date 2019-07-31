//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


#include "Particle/DistanceTableData.h"
#include "QMCWaveFunctions/Jastrow/eeI_JastrowBuilder.h"
#include "QMCWaveFunctions/Jastrow/eeI_JastrowOrbital.h"
#include "QMCWaveFunctions/Jastrow/JeeIOrbitalSoA.h"
#include "QMCWaveFunctions/Jastrow/DiffOneBodyJastrowOrbital.h"
#include "QMCWaveFunctions/Jastrow/TwoBodyJastrowOrbital.h"
#include "QMCWaveFunctions/Jastrow/DiffTwoBodyJastrowOrbital.h"
#include "Utilities/ProgressReportEngine.h"
#include "QMCWaveFunctions/Jastrow/PolynomialFunctor3D.h"

namespace qmcplusplus
{
template<typename J3type>
bool eeI_JastrowBuilder::putkids(xmlNodePtr kids, J3type& J3)
{
  std::string jname = "JeeI";
  SpeciesSet& iSet  = sourcePtcl->getSpeciesSet();
  SpeciesSet& eSet  = targetPtcl.getSpeciesSet();
  int numiSpecies   = iSet.getTotalNum();
  //read in xml
  while (kids != NULL)
  {
    std::string kidsname = (char*)kids->name;
    if (kidsname == "correlation")
    {
      RealType ee_cusp = 0.0;
      RealType eI_cusp = 0.0;
      std::string iSpecies, eSpecies1("u"), eSpecies2("u");
      OhmmsAttributeSet rAttrib;
      rAttrib.add(iSpecies, "ispecies");
      rAttrib.add(eSpecies1, "especies1");
      rAttrib.add(eSpecies2, "especies2");
      rAttrib.add(ee_cusp, "ecusp");
      rAttrib.add(eI_cusp, "icusp");
      rAttrib.put(kids);
      typedef typename J3type::FuncType FT;
      FT* functor        = new FT(ee_cusp, eI_cusp);
      functor->iSpecies  = iSpecies;
      functor->eSpecies1 = eSpecies1;
      functor->eSpecies2 = eSpecies2;
      int iNum           = iSet.findSpecies(iSpecies);
      int eNum1          = eSet.findSpecies(eSpecies1);
      int eNum2          = eSet.findSpecies(eSpecies2);
      if (iNum == iSet.size())
      {
        APP_ABORT("ion species " + iSpecies + " requested for Jastrow " + jname + " does not exist in ParticleSet " +
                  sourcePtcl->getName());
      }
      std::string illegal_eSpecies;
      if (eNum1 == eSet.size())
        illegal_eSpecies = eSpecies1;
      if (eNum2 == eSet.size())
      {
        if (illegal_eSpecies.size())
          illegal_eSpecies += " and ";
        illegal_eSpecies += eSpecies2;
      }
      if (illegal_eSpecies.size())
        APP_ABORT("electron species " + illegal_eSpecies + " requested for Jastrow " + jname +
                  " does not exist in ParticleSet " + targetPtcl.getName());
      functor->put(kids);
      if (sourcePtcl->Lattice.SuperCellEnum != SUPERCELL_OPEN)
      {
        const RealType WSRadius = sourcePtcl->Lattice.WignerSeitzRadius;
        if (functor->cutoff_radius > WSRadius)
        {
          if (functor->cutoff_radius - WSRadius > 1e-4)
          {
            APP_ABORT("  The eeI Jastrow cutoff specified should not be larger than Wigner-Seitz radius.");
          }
          else
          {
            app_log() << "  The eeI Jastrow cutoff specified is slightly larger than the Wigner-Seitz radius.";
            app_log() << "  Setting to Wigner-Seitz radius = " << WSRadius << ".\n";
            functor->cutoff_radius = WSRadius;
            functor->reset();
          }
        }
        if (functor->cutoff_radius < 1.0e-6)
        {
          app_log() << "  eeI functor rcut is currently zero.\n"
                    << "  Setting to Wigner-Seitz radius = " << WSRadius << std::endl;
          functor->cutoff_radius = WSRadius;
          functor->reset();
        }
      }
      else if (functor->cutoff_radius < 1.0e-6)
      {
        APP_ABORT("  eeI Jastrow cutoff unspecified.  Cutoff must be given when using open boundary conditions");
      }
      J3.addFunc(iNum, eNum1, eNum2, functor);
    }
    kids = kids->next;
  }
  //check that each ion species has up and down components
  J3.check_complete();
  targetPsi.addComponent(&J3, jname.c_str());
  J3.setOptimizable(true);
  return true;
}

bool eeI_JastrowBuilder::put(xmlNodePtr cur)
{
  ReportEngine PRE(ClassName, "put(xmlNodePtr)");
  xmlNodePtr kids = cur->xmlChildrenNode;
  // Create a three-body Jastrow
  if (sourcePtcl)
  {
    std::string ftype("polynomial");
    OhmmsAttributeSet tAttrib;
    tAttrib.add(ftype, "function");
    tAttrib.put(cur);
    SpeciesSet& iSet = sourcePtcl->getSpeciesSet();
    SpeciesSet& eSet = targetPtcl.getSpeciesSet();
    int numiSpecies  = iSet.getTotalNum();
    if (ftype == "polynomial")
    {
#ifdef ENABLE_SOA
      typedef JeeIOrbitalSoA<PolynomialFunctor3D> J3Type;
#else
      typedef eeI_JastrowOrbital<PolynomialFunctor3D> J3Type;
#endif
      J3Type& J3 = *(new J3Type(*sourcePtcl, targetPtcl, true));
      putkids(kids, J3);
    }
    else
    {
      app_error() << "Unknown function \"" << ftype << "\" in"
                  << " eeI_JastrowBuilder.  Aborting.\n";
      abort();
    }
  }
  else
    app_error() << "You must specify the \"source\" particleset for a three-body Jastrow.\n";
  return true;
}
} // namespace qmcplusplus
