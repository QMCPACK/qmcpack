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

#include "RadialJastrowBuilder.h"

#include "Utilities/IteratorUtility.h"
#include "Utilities/ProgressReportEngine.h"
#include "OhmmsData/AttributeSet.h"
#include "OhmmsData/ParameterSet.h"

#include "Particle/DistanceTableData.h"
#include "Particle/DistanceTable.h"



#include "QMCWaveFunctions/Jastrow/DiffOneBodyJastrowOrbital.h"
#include "QMCWaveFunctions/Jastrow/DiffOneBodySpinJastrowOrbital.h"
#include "QMCWaveFunctions/Jastrow/DiffTwoBodyJastrowOrbital.h"
#include "QMCWaveFunctions/Jastrow/OneBodyJastrowSpinOrbital.h"

#if defined(QMC_CUDA)
#if defined(ENABLE_SOA)
#include "QMCWaveFunctions/Jastrow/OneBodyJastrowOrbitalBspline.h"
#include "QMCWaveFunctions/Jastrow/TwoBodyJastrowOrbitalBspline.h"
#else
#include "QMCWaveFunctions/Jastrow/OneBodyJastrowOrbitalBsplineAoS.h"
#include "QMCWaveFunctions/Jastrow/TwoBodyJastrowOrbitalBsplineAoS.h"
#endif
#else // defined(QMC_CUDA)
#if defined(ENABLE_SOA)
#include "QMCWaveFunctions/Jastrow/J1OrbitalSoA.h"
#include "QMCWaveFunctions/Jastrow/J2OrbitalSoA.h"
#else
#include "QMCWaveFunctions/Jastrow/OneBodyJastrowOrbital.h"
#include "QMCWaveFunctions/Jastrow/TwoBodyJastrowOrbital.h"
#endif
#endif


namespace qmcplusplus
{

RadialJastrowBuilder::RadialJastrowBuilder(PartcleSet& target, TrialWaveFunction& psi,
					   PtclPoolType& psets):
  WaveFunctionComponentBuilder(target, psi),ptclPool(psets)
{
  ClassName="RadialJastrowBuilder";
}

// helper class to simplify and localize ugly ifdef stuff for types
template<typename precision, template<class> class RadFuncType>
class JastrowTypeHelper 
{
 private:
 public:
  typedef RadFuncType<precision> rft;
#if defined(ENABLE_SOA)
  typedef J1OrbitalSoA<rft> J1OrbitalType;
  typedef J1OrbitalSoA<rft> J1OrbitalTypeSpin;
  typedef DiffOneBodyJastrowOrbital<rft> DiffJ1OrbitalType;
  typedef DiffOneBodyJastrowOrbital<rft> DiffJ1OrbitalTypeSpin;
  typedef J2OrbitalSoA<rft> J2OrbitalType;
  typedef DiffTwoBodyJastrowOrbital<rft> DiffJ2OrbitalType;
#else
  typedef OneBodyJastrowOrbital<rft> J1OrbitalType;
  typedef OneBodyJastrowOrbital<rft> J1OrbitalTypeSpin;
  typedef DiffOneBodyJastrowOrbital<rft> DiffJ1OrbitalType;
  typedef DiffOneBodyJastrowOrbital<rft> DiffJ1OrbitalTypeSpin;
  typedef TwoBodyJastrowOrbital<rft> J2OrbitalType;
  typedef DiffTwoBodyJastrowOrbital<rft> DiffJ2OrbitalType;
#endif
};


// specialization for bspline (does this need to be so complicated?)
// note that this supports CUDA whereas the general case does not
template<typename precision>
class JastrowTypeHelper<precision, BsplineFunctor>
{
 private:
 public:
  typedef BsplineFunctor<precision> rft;
#if defined(QMC_CUDA) and defined(ENABLE_SOA)
  typedef OneBodyJastrowOrbitalBspline<rft> J1OrbitalType;
  typedef OneBodyJastrowOrbitalBspline<rft> J1OrbitalTypeSpin;
  typedef DiffOneBodySpinJastrowOrbital<rft> DiffJ1OrbitalType;
  typedef DiffOneBodySpinJastrowOrbital<rft> DiffJ1OrbitalTypeSpin;
  typedef TwoBodyJastrowOrbitalBspline<rft> J2OrbitalType;
  typedef DiffTwoBodyJastrowOrbital<rft> DiffJ2OrbitalType;
#if defined(QMC_CUDA) and !defined(ENABLE_SOA)
  typedef OneBodyJastrowOrbitalBsplineAoS J1OrbitalType;
  typedef OneBodyJastrowOrbitalBsplineAoS J1OrbitalTypeSpin;
  typedef DiffOneBodySpinJastrowOrbital<rft> DiffJ1OrbitalType;
  typedef DiffOneBodySpinJastrowOrbital<rft> DiffJ1OrbitalTypeSpin;
  typedef TwoBodyJastrowOrbitalBsplineAoS J2OrbitalType;
  typedef DiffTwoBodyJastrowOrbital<rft> DiffJ2OrbitalType;
#if !defined(QMC_CUDA) and defined(ENABLE_SOA)
  typedef J1OrbitalSoA<rft> J1OrbitalType;
  typedef OneBodySpinJastrowOrbital<rft> J1OrbitalTypeSpin;
  typedef DiffOneBodyJastrowOrbital<rft> DiffJ1OrbitalType;
  typedef DiffOneBodySpinJastrowOrbital<rft> DiffJ1OrbitalTypeSpin;
  typedef J2OrbitalSoA<rft> J2OrbitalType;
  typedef DiffTwoBodyJastrowOrbital<rft> DiffJ2OrbitalType;
#if !defined(QMC_CUDA) and !defined(ENABLE_SOA)
  typedef OneBodyJastrowOrbital<rft> J1OrbitalType;
  typedef OneBodySpinJastrowOrbital<rft> J1OrbitalTypeSpin;
  typedef DiffOneBodyJastrowOrbital<rft> DiffJ1OrbitalType;
  typedef DiffOneBodySpinJastrowOrbital<rft> DiffJ1OrbitalTypeSpin;
  typedef TwoBodyJastrowOrbital<rft> J2OrbitalType;
  typedef DiffTwoBodyJastrowOrbital<rft> DiffJ2OrbitalType;
};


}





















bool RadialJastrowBuilder::put(xmlNodePtr cur)
{
  // not sure if this can hande spin dependence yet
  // also not sure whether I should be using ptclPool
  ReportEngine PRE(ClassName,"put(xmlNodePtr)");
  std::string jastfunction("unknown");
  std::string sourceOpt=targetPtcl.getName();
  std::string spin="yes";
  std::string id_b="jee_b";
  OhmmsAttributeSet aAttrib;
  aAttrib.add(jastfunction,"function");
  aAttrib.add(spin, "spin");
  aAttrib.add(sourceOpt, "source");
  aAttrib.put(cur);
  bool success=false;

  SpeciesSet& species(targetPtcl.getSpeciesSet());
  int chargeInd=species.addAttribute("charge");
  
  if (sourceOpt == targetPtcl.getName())
  {
    // it's a two body jastrow factor
    int taskid=(targetPsi.is_manager())?targetPsi.getGroupID():-1;

    cur= cur->xmlChildrenNode;
    while(cur!=NULL)
    {
      std::string cname((const char*)cur->name);
      if (cname == "correlation")
      {
        OhmmsAttributeSet rAttrib;
        RealType cusp=-1e10;
        std::string spA(species.speciesName[0]);
        std::string spB(species.speciesName[0]);
        rAttrib.add(spA,"speciesA");
        rAttrib.add(spB,"speciesB");
        rAttrib.add(cusp,"cusp");
        rAttrib.put(cur);
        int ia = species.findSpecies(spA);
        int ib = species.findSpecies(spB);
        if(ia==species.size() || ib == species.size())
        {
          PRE.error("Failed. Species are incorrect.",true);
        }
        if(cusp<-1e6)
        {
          RealType qq=species(chargeInd,ia)*species(chargeInd,ib);
          cusp = (ia==ib)? -0.25*qq:-0.5*qq;
        }
        std::ostringstream o;
        o<<"j2"<<ia<<ib;
	// now have a series of if's to pick up various available functors (pade1, pade2, casino, bspline?)
	// and then call createJ2 with that functor type
	// we should have some logic to deal with species dependence (also need to have guards a-la bspline
	// for the case where there is only a single (or zero) electron of a given species
      }

  } 
  else    
  {
    // it's a one body jastrow factor
  }
}

template<typename RadFuncType>
bool RadialJastrowBuilder::createJ2(xmlNodePtr cur, const std::string& jname) 
{
  SpeciesSet& species(targetPtcl.getSpeciesSet());
  int taskid=(targetPsi.is_manager())?targetPsi.getGroupID():-1;
#if defined(ENABLE_SOA)
  typedef J2OrbitalSoA<RadFuncType> J2Type;
#else
  typedef TwoBodyJastrowOrbital<RadFuncType> J2Type;
#endif
  typedef DiffTwoBodyJastrowOrbital<RadFuncType> dJ2Type;
  int taskid=(targetPsi.is_manager())?targetPsi.getGroupID():-1;
  J2Type *J2 = new J2Type(targetPtcl,taskid);
  dJ2Type *dJ2 = new dJ2Type(targetPtcl);
  cur= cur->xmlChildrenNode;
  while(cur!=NULL)
  {
    std::string cname((const char*)cur->name);
    if (cname == "correlation")
    {
      OhmmsAttributeSet rAttrib;
      RealType cusp=-1e10;
      std::string spA(species.speciesName[0]);
      std::string spB(species.speciesName[0]);
      rAttrib.add(spA,"speciesA");
      rAttrib.add(spB,"speciesB");
      rAttrib.add(cusp,"cusp");
      rAttrib.put(cur);
      int ia = species.findSpecies(spA);
      int ib = species.findSpecies(spB);
      if(ia==species.size() || ib == species.size())
      {
	PRE.error("Failed. Species are incorrect.",true);
      }
      if(cusp<-1e6)
      {
	RealType qq=species(chargeInd,ia)*species(chargeInd,ib);
	cusp = (ia==ib)? -0.25*qq:-0.5*qq;
      }
      std::ostringstream o;
      o<<"j2"<<ia<<ib;
      RadFuncType *functor = new RadFuncType(cusp,o.str());
      functor->put(cur);
      J2->addFunc(ia,ib,functor);
      dJ2->addFunc(ia,ib,functor);
      }
    cur=cur->next;
  }
  J2->dPsi=dJ2;
  std::string j2name="J2_"+jname;
  targetPsi.addOrbital(J2,j2name);
  J2->setOptimizable(true);
  return true;
}
