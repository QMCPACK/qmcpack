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
#include "QMCWaveFunctions/Jastrow/BsplineJastrowBuilder.h"
#include "QMCWaveFunctions/Jastrow/BsplineFunctor.h"
#include "QMCWaveFunctions/Jastrow/OneBodyJastrowOrbital.h"
#include "QMCWaveFunctions/Jastrow/DiffOneBodyJastrowOrbital.h"
#include "QMCWaveFunctions/Jastrow/OneBodySpinJastrowOrbital.h"
#include "QMCWaveFunctions/Jastrow/DiffOneBodySpinJastrowOrbital.h"
#include "QMCWaveFunctions/Jastrow/TwoBodyJastrowOrbital.h"
#include "QMCWaveFunctions/Jastrow/J1OrbitalSoA.h"
#include "QMCWaveFunctions/Jastrow/J2OrbitalSoA.h"
#include "QMCWaveFunctions/Jastrow/DiffTwoBodyJastrowOrbital.h"
#ifdef QMC_CUDA
#include "QMCWaveFunctions/Jastrow/OneBodyJastrowOrbitalBspline.h"
#include "QMCWaveFunctions/Jastrow/TwoBodyJastrowOrbitalBspline.h"
#endif
#include "LongRange/LRRPAHandlerTemp.h"
#include "QMCWaveFunctions/Jastrow/LRBreakupUtilities.h"
#include "Utilities/ProgressReportEngine.h"
#include "LongRange/LRJastrowSingleton.h"

namespace qmcplusplus
{

template<typename OBJT, typename DOBJT>
bool BsplineJastrowBuilder::createOneBodyJastrow(xmlNodePtr cur)
{
  ReportEngine PRE(ClassName,"createOneBodyJastrow(xmlNodePtr)");
  std::string j1name("J1");
  {
    OhmmsAttributeSet a;
    a.add(j1name,"name");
    a.put(cur);
  }
  int taskid=targetPsi.getGroupID();//(targetPsi.is_manager())?targetPsi.getGroupID():-1;
  OBJT* J1 =new OBJT(*sourcePtcl,targetPtcl);
  DOBJT *dJ1 = new DOBJT(*sourcePtcl, targetPtcl);
  xmlNodePtr kids = cur->xmlChildrenNode;
  // Find the number of the source species
  SpeciesSet &sSet = sourcePtcl->getSpeciesSet();
  SpeciesSet &tSet = targetPtcl.getSpeciesSet();
  int numSpecies = sSet.getTotalNum();
  bool success=false;
  bool Opt(false);
  while (kids != NULL)
  {
    std::string kidsname = (char*)kids->name;
    if (kidsname == "correlation")
    {
      RealType cusp=0.0;
      std::string speciesA;
      std::string speciesB;
      OhmmsAttributeSet rAttrib;
      rAttrib.add(speciesA,"elementType");
      rAttrib.add(speciesA,"speciesA");
      rAttrib.add(speciesB,"speciesB");
      rAttrib.add(cusp,"cusp");
      rAttrib.put(kids);
      BsplineFunctor<RealType> *functor = new BsplineFunctor<RealType>(cusp);
      functor->elementType = speciesA;
      int ig = sSet.findSpecies (speciesA);
      functor->periodic = sourcePtcl->Lattice.SuperCellEnum != SUPERCELL_OPEN;
      functor->cutoff_radius = sourcePtcl->Lattice.WignerSeitzRadius;
      int jg=-1;
      if(speciesB.size())
        jg=tSet.findSpecies(speciesB);
      if (ig < numSpecies)
      {
        //ignore
        functor->put (kids);
        J1->addFunc (ig,functor,jg);
        success = true;
        dJ1->addFunc(ig,functor,jg);
        Opt=(!functor->notOpt or Opt);
        if(qmc_common.io_node)
        {
          char fname[128];
          if(qmc_common.mpi_groups>1)
          {
            if(speciesB.size())
              sprintf(fname,"%s.%s%s.g%03d.dat",j1name.c_str(),speciesA.c_str(),speciesB.c_str(),taskid);
            else
              sprintf(fname,"%s.%s.g%03d.dat",j1name.c_str(),speciesA.c_str(),taskid);
          }
          else
          {
            if(speciesB.size())
              sprintf(fname,"%s.%s%s.dat",j1name.c_str(),speciesA.c_str(),speciesB.c_str());
            else
              sprintf(fname,"%s.%s.dat",j1name.c_str(),speciesA.c_str());
          }
          functor->setReportLevel(ReportLevel,fname);
          functor->print();
        }
      }
    }
    kids = kids->next;
  }
  if(success)
  {
    J1->dPsi=dJ1;
    targetPsi.addOrbital(J1,"J1_bspline");
    J1->setOptimizable(Opt);
    return true;
  }
  else
  {
    PRE.warning("BsplineJastrowBuilder failed to add an One-Body Jastrow.");
    delete J1;
    delete dJ1;
    return false;
  }
}

/** class to initialize bsplin functions
 */
template<typename T>
struct BsplineInitializer
{

  std::vector<T> rpaValues;

  /** initialize with RPA
   * @param ref particleset
   * @param bfunc bspline function to be initialized
   * @param nopbc true, if not periodic
   */
  inline void initWithRPA(ParticleSet& P,  BsplineFunctor<T>& bfunc, T fac)
  {
    if(P.Lattice.SuperCellEnum==SUPERCELL_OPEN) // for open systems, do nothing
    {
      return;
      //T vol=std::pow(bfunc.cutoff_radius,3);
      //rpa.reset(ref.getTotalNum(),vol*0.5);
    }
    int npts=bfunc.NumParams;
    if(rpaValues.empty())
    {
      rpaValues.resize(npts);
      LRRPAHandlerTemp<RPABreakup<T>,LPQHIBasis> rpa(P,-1.0);
      rpa.Breakup(P,-1.0);
      T dr=bfunc.cutoff_radius/static_cast<T>(npts);
      T r=0;
      for (int i=0; i<npts; i++)
      {
        rpaValues[i]=rpa.evaluate(r,1.0/r); //y[i]=fac*rpa.evaluate(r,1.0/r);
        r += dr;
      }
    }
    T last=rpaValues[npts-1];
    //vector<T>  y(npts);
    //for (int i=0; i<npts; i++) y[i]=fac*rpaValues[i];
    //T last=y[npts-1];
    for(int i=0; i<npts; i++)
      bfunc.Parameters[i]=fac*(rpaValues[i]-last);
    bfunc.reset();
  }
};

bool BsplineJastrowBuilder::put(xmlNodePtr cur)
{
  ReportEngine PRE(ClassName,"put(xmlNodePtr)");
  bool PrintTables=false;
  typedef BsplineFunctor<RealType> RadFuncType;
  // Create a one-body Jastrow
  if (sourcePtcl)
  {
    std::string j1spin("no");
    OhmmsAttributeSet jAttrib;
    jAttrib.add(j1spin,"spin");
    jAttrib.put(cur);
#ifdef QMC_CUDA
    return createOneBodyJastrow<OneBodyJastrowOrbitalBspline,DiffOneBodySpinJastrowOrbital<RadFuncType> >(cur);
#else
    //if(sourcePtcl->IsGrouped)
    //{
    //  app_log() << "Creating OneBodySpinJastrowOrbital<T> " << std::endl;
    //  return createOneBodyJastrow<OneBodySpinJastrowOrbital<RadFuncType>,DiffOneBodySpinJastrowOrbital<RadFuncType> >(cur);
    //}
    //else
    //{
    //  app_log() << "Creating OneBodyJastrowOrbital<T> " << std::endl;
    //  return createOneBodyJastrow<OneBodyJastrowOrbital<RadFuncType>,DiffOneBodyJastrowOrbital<RadFuncType> >(cur);
    //}
    if(j1spin=="yes")
      return createOneBodyJastrow<OneBodySpinJastrowOrbital<RadFuncType>,DiffOneBodySpinJastrowOrbital<RadFuncType> >(cur);
    else
#if defined(ENABLE_SOA)
      return createOneBodyJastrow<J1OrbitalSoA<RadFuncType>,DiffOneBodyJastrowOrbital<RadFuncType> >(cur);
#else
      return createOneBodyJastrow<OneBodyJastrowOrbital<RadFuncType>,DiffOneBodyJastrowOrbital<RadFuncType> >(cur);
#endif

#endif
  }
  else // Create a two-body Jastrow
  {
    std::string init_mode("0");
    {
      OhmmsAttributeSet hAttrib;
      hAttrib.add(init_mode,"init");
      hAttrib.put(cur);
    }
    BsplineInitializer<RealType> j2Initializer;
    xmlNodePtr kids = cur->xmlChildrenNode;
#ifdef QMC_CUDA
    typedef TwoBodyJastrowOrbitalBspline J2Type;
#else

#if defined(ENABLE_SOA)
    typedef J2OrbitalSoA<BsplineFunctor<RealType> > J2Type;
#else
    typedef TwoBodyJastrowOrbital<BsplineFunctor<RealType> > J2Type;
#endif

#endif
    typedef DiffTwoBodyJastrowOrbital<BsplineFunctor<RealType> > dJ2Type;
    int taskid=(targetPsi.is_manager())?targetPsi.getGroupID():-1;
    J2Type *J2 = new J2Type(targetPtcl,taskid);
    dJ2Type *dJ2 = new dJ2Type(targetPtcl);
    SpeciesSet& species(targetPtcl.getSpeciesSet());
    int chargeInd=species.addAttribute("charge");
    //std::map<std::string,RadFuncType*> functorMap;
    bool Opt(false);
    while (kids != NULL)
    {
      std::string kidsname((const char*)kids->name);
      if (kidsname == "correlation")
      {
        OhmmsAttributeSet rAttrib;
        RealType cusp=-1e10;
        std::string pairType("0");
        std::string spA(species.speciesName[0]);
        std::string spB(species.speciesName[0]);
        rAttrib.add(spA,"speciesA");
        rAttrib.add(spB,"speciesB");
        rAttrib.add(pairType,"pairType");
        rAttrib.add(cusp,"cusp");
        rAttrib.put(kids);
        if(pairType[0]=='0')
        {
          pairType=spA+spB;
        }
        else
        {
          PRE.warning("pairType is deprecated. Use speciesA/speciesB");
          //overwrite the species
          spA=pairType[0];
          spB=pairType[1];
        }
        int ia = species.findSpecies(spA);
        int ib = species.findSpecies(spB);
        if(ia==species.size() || ib == species.size())
        {
          PRE.error("Failed. Species are incorrect.",true);
        }
        // prevent adding uu/dd correlation if there is only 1 u/d electron.
        if(ia==ib && (targetPtcl.last(ia)-targetPtcl.first(ia)==1))
          PRE.error("Failed to add "+spA+spB+" correlation for only 1 "+spA+" particle. Please remove it from two-body Jastrow.",true);
        if(cusp<-1e6)
        {
          RealType qq=species(chargeInd,ia)*species(chargeInd,ib);
          cusp = (ia==ib)? -0.25*qq:-0.5*qq;
        }
        app_log() << "  BsplineJastrowBuilder adds a functor with cusp = " << cusp << std::endl;
        RadFuncType *functor = new RadFuncType(cusp);
        functor->periodic      = targetPtcl.Lattice.SuperCellEnum != SUPERCELL_OPEN;
        functor->cutoff_radius = targetPtcl.Lattice.WignerSeitzRadius;
        bool initialized_p=functor->put(kids);
        functor->elementType=pairType;
        //RPA INIT
        if(!initialized_p && init_mode =="rpa")
        {
          app_log() << "  Initializing Two-Body with RPA Jastrow " << std::endl;
          j2Initializer.initWithRPA(targetPtcl,*functor,-cusp/0.5);
        }
        J2->addFunc(ia,ib,functor);
        dJ2->addFunc(ia,ib,functor);
        Opt=(!functor->notOpt or Opt);
        if(qmc_common.io_node)
        {
          char fname[32];
          if(qmc_common.mpi_groups>1)
            sprintf(fname,"J2.%s.g%03d.dat",pairType.c_str(),taskid);
          else
            sprintf(fname,"J2.%s.dat",pairType.c_str());
          functor->setReportLevel(ReportLevel,fname);
          functor->print();
        }
      }
      kids = kids->next;
    }
    //dJ2->initialize();
    //J2->setDiffOrbital(dJ2);
    J2->dPsi=dJ2;
    targetPsi.addOrbital(J2,"J2_bspline");
    J2->setOptimizable(Opt);

  }
  return true;
}

}
