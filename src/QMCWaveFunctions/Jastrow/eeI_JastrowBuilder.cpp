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


#include "Particle/DistanceTable.h"
#include "eeI_JastrowBuilder.h"
#include "QMCWaveFunctions/Jastrow/JeeIOrbitalSoA.h"
#include "Utilities/ProgressReportEngine.h"
#include "QMCWaveFunctions/Jastrow/PolynomialFunctor3D.h"

namespace qmcplusplus
{
eeI_JastrowBuilder::eeI_JastrowBuilder(Communicate* comm, ParticleSet& target, ParticleSet& source)
    : WaveFunctionComponentBuilder(comm, target), sourcePtcl(&source)
{
  ClassName = "eeI_JastroBuilder";
}


template<typename J3type>
bool eeI_JastrowBuilder::putkids(xmlNodePtr kids, J3type& J3)
{
  auto& jname      = J3.getName();
  SpeciesSet& iSet = sourcePtcl->getSpeciesSet();
  SpeciesSet& eSet = targetPtcl.getSpeciesSet();
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
      using FT           = typename J3type::FuncType;
      const auto coef_id = extractCoefficientsID(kids);
      auto functor =
          std::make_unique<FT>(coef_id.empty() ? jname + "_"  + iSpecies + eSpecies1 + eSpecies2 : coef_id, ee_cusp, eI_cusp);
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
      if (sourcePtcl->getLattice().SuperCellEnum != SUPERCELL_OPEN)
      {
        const RealType WSRadius = sourcePtcl->getLattice().WignerSeitzRadius;
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
      J3.addFunc(iNum, eNum1, eNum2, std::move(functor));
    }
    kids = kids->next;
  }
  //check that each ion species has up and down components
  J3.check_complete();
  return true;
}

std::unique_ptr<WaveFunctionComponent> eeI_JastrowBuilder::buildComponent(xmlNodePtr cur)
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

    std::string input_name(getXMLAttributeValue(cur, "name"));
    std::string jname = input_name.empty() ? "JeeI_" + ftype : input_name;
    SpeciesSet& iSet  = sourcePtcl->getSpeciesSet();
    if (ftype == "polynomial")
    {
      using J3Type = JeeIOrbitalSoA<PolynomialFunctor3D>;
      auto J3      = std::make_unique<J3Type>(jname, *sourcePtcl, targetPtcl);
      putkids(kids, *J3);
      return J3;
    }
    else
    {
      std::ostringstream err_msg;
      err_msg << "Unknown function \"" << ftype << "\" in"
              << " eeI_JastrowBuilder.  Aborting.\n";
      APP_ABORT(err_msg.str());
    }
  }
  else
    APP_ABORT("You must specify the \"source\" particleset for a three-body Jastrow.\n");
  return nullptr;
}

} // namespace qmcplusplus
