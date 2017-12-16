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
#include "Particle/DistanceTable.h"
#include "QMCWaveFunctions/Jastrow/eeI_JastrowBuilder.h"
#include "QMCWaveFunctions/Jastrow/BsplineFunctor.h"
#include "QMCWaveFunctions/Jastrow/eeI_JastrowOrbital.h"
#include "QMCWaveFunctions/Jastrow/JeeIOrbitalSoA.h"
#include "QMCWaveFunctions/Jastrow/DiffOneBodyJastrowOrbital.h"
#include "QMCWaveFunctions/Jastrow/TwoBodyJastrowOrbital.h"
#include "QMCWaveFunctions/Jastrow/DiffTwoBodyJastrowOrbital.h"
#include "Utilities/ProgressReportEngine.h"
#include "QMCWaveFunctions/Jastrow/BsplineFunctor3D.h"
#include "QMCWaveFunctions/Jastrow/PolynomialFunctor3D.h"

namespace qmcplusplus
{

template<typename J3type>
bool eeI_JastrowBuilder::putkids (xmlNodePtr kids, J3type &J3)
{
  SpeciesSet &iSet = sourcePtcl->getSpeciesSet();
  SpeciesSet &eSet = targetPtcl.getSpeciesSet();
  int numiSpecies = iSet.getTotalNum();
  bool success=false;
  //read in xml
  while (kids != NULL)
  {
    std::string kidsname = (char*)kids->name;
    if (kidsname == "correlation")
    {
      RealType ee_cusp=0.0;
      RealType eI_cusp=0.0;
      std::string iSpecies, eSpecies1("u"), eSpecies2("u");
      OhmmsAttributeSet rAttrib;
      rAttrib.add(iSpecies,"ispecies");
      rAttrib.add(eSpecies1,"especies1");
      rAttrib.add(eSpecies2,"especies2");
      rAttrib.add(ee_cusp,"ecusp");
      rAttrib.add(eI_cusp,"icusp");
      rAttrib.put(kids);
      typedef typename J3type::FuncType FT;
      FT *functor = new FT(ee_cusp, eI_cusp);
      functor->iSpecies = iSpecies;
      functor->eSpecies1 = eSpecies1;
      functor->eSpecies2 = eSpecies2;
      int iNum = iSet.findSpecies (iSpecies);
      int eNum1 = eSet.findSpecies (eSpecies1);
      int eNum2 = eSet.findSpecies (eSpecies2);
      functor->put (kids);
      if (sourcePtcl->Lattice.SuperCellEnum != SUPERCELL_OPEN)
      {
        const RealType WSRadius = sourcePtcl->Lattice.WignerSeitzRadius;
        if (functor->cutoff_radius > WSRadius)
        {
          if ( functor->cutoff_radius - WSRadius > 1e-4 )
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
          app_log()  << "  eeI functor rcut is currently zero.\n"
                     << "  Setting to Wigner-Seitz radius = "
                     << WSRadius << std::endl;
          functor->cutoff_radius = WSRadius;
          functor->reset();
        }
      }
      else
        if (functor->cutoff_radius < 1.0e-6)
        {
          APP_ABORT("  eeI Jastrow cutoff unspecified.  Cutoff must be given when using open boundary conditions");
        }
      J3.addFunc(iNum, eNum1, eNum2, functor);
    }
    kids = kids->next;
  }
  //check that each ion species has up and down components
  J3.check_complete();
  targetPsi.addOrbital(&J3,"JeeI");
  J3.setOptimizable(true);
  return true;
}

bool eeI_JastrowBuilder::put(xmlNodePtr cur)
{
  ReportEngine PRE(ClassName,"put(xmlNodePtr)");
  bool PrintTables=true;
  xmlNodePtr kids = cur->xmlChildrenNode;
  typedef BsplineFunctor3D FuncType;
  // Create a three-body Jastrow
  if (sourcePtcl)
  {
//       std::cerr << "sourcePtcl = " << sourcePtcl << std::endl;
    std::string ftype("Bspline");
    OhmmsAttributeSet tAttrib;
    tAttrib.add(ftype,"function");
    tAttrib.put (cur);
    SpeciesSet &iSet = sourcePtcl->getSpeciesSet();
    SpeciesSet &eSet = targetPtcl.getSpeciesSet();
    int numiSpecies = iSet.getTotalNum();
    if (ftype == "Bspline")
    {
#ifdef ENABLE_SOA
      typedef JeeIOrbitalSoA<BsplineFunctor3D> J3Type;
#else
      typedef eeI_JastrowOrbital<BsplineFunctor3D> J3Type;
#endif
      J3Type &J3 = *(new J3Type(*sourcePtcl, targetPtcl, true));
      putkids (kids, J3);
    }
    else if (ftype == "polynomial")
    {
#ifdef ENABLE_SOA
      typedef JeeIOrbitalSoA<PolynomialFunctor3D> J3Type;
#else
      typedef eeI_JastrowOrbital<PolynomialFunctor3D> J3Type;
#endif
      J3Type &J3 = *(new J3Type(*sourcePtcl, targetPtcl, true));
      putkids (kids, J3);
    }
    else
    {
      app_error() << "Unknown function \"" << ftype << "\" in"
                  << " eeI_JastrowBuilder.  Aborting.\n";
      abort();
    }
    // 	// Find the number of the source species
    // 	bool success=false;
    // 	while (kids != NULL) {
    // 	  std::string kidsname = (char*)kids->name;
    // 	  if (kidsname == "correlation") {
    // 	    RealType ee_cusp=0.0;
    // 	    RealType eI_cusp=0.0;
    // 	    std::string iSpecies, eSpecies1("u"), eSpecies2("u");
    // 	    OhmmsAttributeSet rAttrib;
    // 	    rAttrib.add(iSpecies,"ispecies");
    // 	    rAttrib.add(eSpecies1,"especies1");
    // 	    rAttrib.add(eSpecies2,"especies2");
    // 	    rAttrib.add(ee_cusp,"ecusp");
    // 	    rAttrib.add(eI_cusp,"icusp");
    // 	    rAttrib.put(kids);
    // 	    BsplineFunctor3D *functor = new BsplineFunctor3D(ee_cusp, eI_cusp);
    // 	    functor->iSpecies = iSpecies;
    // 	    functor->eSpecies1 = eSpecies1;
    // 	    functor->eSpecies2 = eSpecies2;
    // 	    int iNum = iSet.findSpecies (iSpecies);
    // 	    int eNum1 = eSet.findSpecies (eSpecies1);
    // 	    int eNum2 = eSet.findSpecies (eSpecies2);
    // 	    functor->put (kids);
    // 	    std::strstream aname;
    // 	    aname << iSpecies << "_" << eSpecies1 << "_" << eSpecies2;
    // 	    J3->addFunc(aname.str(), iNum, eNum1, eNum2, functor);
    // 	  }
    // 	  kids = kids->next;
    // 	}
    // 	targetPsi.addOrbital(J3,"eeI_bspline");
    // 	J3->setOptimizable(true);
    // }
  }
  else
    app_error() << "You must specify the \"source\" particleset for a three-body Jastrow.\n";
  return true;
  //   // Find the number of the source species
  //   SpeciesSet &sSet = sourcePtcl->getSpeciesSet();
  //   int numSpecies = sSet.getTotalNum();
  //   bool success=false;
  //   while (kids != NULL)
  //   {
  // 	std::string kidsname = (char*)kids->name;
  // 	if (kidsname == "correlation")
  //     {
  //       RealType cusp=0.0;
  //       std::string elementType;
  //       OhmmsAttributeSet rAttrib;
  //       rAttrib.add(elementType,"elementType");
  //       rAttrib.add(cusp,"cusp");
  //       rAttrib.put(kids);
  // 	  BsplineFunctor<double> *functor = new BsplineFunctor<double>(cusp);
  // 	  functor->elementType = elementType;
  // 	  int ig = sSet.findSpecies (elementType);
  // 	  if (ig < numSpecies)
  //       {//ignore
  //         functor->put (kids);
  // 	    if (functor->cutoff_radius < 1.0e-6) {
  // 	      app_log()  << "  BsplineFunction rcut is currently zero.\n"
  // 			 << "  Setting to Wigner-Seitz radius = "
  // 			 << sourcePtcl->Lattice.WignerSeitzRadius << std::endl;
  // 	      functor->cutoff_radius = sourcePtcl->Lattice.WignerSeitzRadius;
  // 	      functor->reset();
  // 	    }
  //         J1->addFunc (ig,functor);
  // 	    success = true;
  //         dJ1->addFunc(ig,functor);
  //         if(ReportLevel)
  //         {
  //           std::string fname="J1."+elementType+".dat";
  //           std::ofstream fout(fname.c_str());
  //           fout.setf(std::ios::scientific, std::ios::floatfield);
  //           fout << "# One-body Jastrow generated by BsplineJastrowBuilder" << std::endl;
  //           functor->print(fout);
  //         }
  //       }
  // 	}
  // 	kids = kids->next;
  //   }
  //   if(success)
  //   {
  //     //assign derivatives to J1
  //     //dJ1->initialize();
  //     //J1->setDiffOrbital(dJ1);
  //     J1->dPsi=dJ1;
  //     targetPsi.addOrbital(J1,"J1_bspline");
  //     J1->setOptimizable(true);
  //   }
  //   else
  //   {
  //     PRE.warning("eeI_JastrowBuilder failed to add an One-Body Jastrow.");
  //     delete J1;
  //     delete dJ1;
  //   }
  // }
  // // Create a two-body Jastrow
  // else
  // {
  //   typedef TwoBodyJastrowOrbital<BsplineFunctor<RealType> > J2Type;
  //   typedef DiffTwoBodyJastrowOrbital<BsplineFunctor<RealType> > dJ2Type;
  //   J2Type *J2 = new J2Type(targetPtcl,targetPsi.is_manager());
  //   dJ2Type *dJ2 = new dJ2Type(targetPtcl);
  //   SpeciesSet& species(targetPtcl.getSpeciesSet());
  //   RealType q=species(0,species.addAttribute("charge"));
  //   //std::map<std::string,RadFuncType*> functorMap;
  //   while (kids != NULL)
  //   {
  // 	std::string kidsname((const char*)kids->name);
  // 	if (kidsname == "correlation")
  //     {
  //       OhmmsAttributeSet rAttrib;
  //       RealType cusp=-1e10;
  //       std::string pairType("0");
  //       std::string spA(species.speciesName[0]);
  //       std::string spB(species.speciesName[0]);
  //       rAttrib.add(spA,"speciesA");
  //       rAttrib.add(spB,"speciesB");
  //       rAttrib.add(pairType,"pairType");
  //       rAttrib.add(cusp,"cusp");
  //       rAttrib.put(kids);
  //       if(pairType[0]=='0')
  //       {
  //         pairType=spA+spB;
  //       }
  //       else
  //       {
  //         PRE.warning("pairType is deprecated. Use speciesA/speciesB");
  //         //overwrite the species
  //         spA=pairType[0]; spB=pairType[1];
  //       }
  //       int ia = species.findSpecies(spA);
  //       int ib = species.findSpecies(spB);
  //       if(ia==species.size() || ib == species.size())
  //       {
  //         PRE.error("Failed. Species are incorrect.",true);
  //       }
  //       if(cusp<-1e6)
  //       {
  //         if(ia==ib)
  //           cusp=-0.25*q*q;
  //         else
  //           cusp=-0.5*q*q;
  //       }
  //       app_log() << "  eeI_JastrowBuilder adds a functor with cusp = " << cusp << std::endl;
  // 	  RadFuncType *functor = new RadFuncType(cusp);
  // 	  functor->put (kids);
  // 	  functor->elementType=pairType;
  // 	  if (functor->cutoff_radius < 1.0e-6) {
  // 	    app_log()  << "  BsplineFunction rcut is currently zero.\n"
  // 		       << "  Setting to Wigner-Seitz radius = "
  // 	    	       << targetPtcl.Lattice.WignerSeitzRadius << std::endl;
  // 	     functor->cutoff_radius = targetPtcl.Lattice.WignerSeitzRadius;
  // 	     functor->reset();
  // 	  }
  //       J2->addFunc(pairType,ia,ib,functor);
  //       dJ2->addFunc(pairType,ia,ib,functor);
  //       if(ReportLevel)
  //       {
  //         std::string fname="J2."+pairType+".dat";
  //         std::ofstream fout(fname.c_str());
  //         fout.setf(std::ios::scientific, std::ios::floatfield);
  //         fout << "# Two-body Jastrow generated by eeI_JastrowBuilder" << std::endl;
  //         functor->print(fout);
  //       }
  // 	}
  // 	kids = kids->next;
  //   }
  //   //dJ2->initialize();
  //   //J2->setDiffOrbital(dJ2);
  //   J2->dPsi=dJ2;
  //   targetPsi.addOrbital(J2,"J2_bspline");
  //   J2->setOptimizable(true);
  // }
  // return true;
}
}
