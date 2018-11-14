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

#if defined(ENABLE_SOA)
#include "QMCWaveFunctions/Jastrow/J1OrbitalSoA.h"
#include "QMCWaveFunctions/Jastrow/J2OrbitalSoA.h"
#else
#include "QMCWaveFunctions/Jastrow/OneBodyJastrowOrbital.h"
#include "QMCWaveFunctions/Jastrow/TwoBodyJastrowOrbital.h"
#endif

#if defined(QMC_CUDA) and defined(ENABLE_SOA)
#include "QMCWaveFunctions/Jastrow/OneBodyJastrowOrbitalBspline.h"
#include "QMCWaveFunctions/Jastrow/TwoBodyJastrowOrbitalBspline.h"
#endif
#if defined(QMC_CUDA) and !defined(ENABLE_SOA)
#include "QMCWaveFunctions/Jastrow/OneBodyJastrowOrbitalBsplineAoS.h"
#include "QMCWaveFunctions/Jastrow/TwoBodyJastrowOrbitalBsplineAoS.h"
#endif
#if !defined(QMC_CUDA) and defined(ENABLE_SOA)
#include "QMCWaveFunctions/Jastrow/J1OrbitalSoA.h"
#include "QMCWaveFunctions/Jastrow/J2OrbitalSoA.h"
#endif
#if !defined(QMC_CUDA) and !defined(ENABLE_SOA)
#include "QMCWaveFunctions/Jastrow/OneBodyJastrowOrbital.h"
#include "QMCWaveFunctions/Jastrow/TwoBodyJastrowOrbital.h"
#endif

#include "QMCWaveFunctions/Jastrow/DiffOneBodyJastrowOrbital.h"
#include "QMCWaveFunctions/Jastrow/DiffOneBodySpinJastrowOrbital.h"
#include "QMCWaveFunctions/Jastrow/DiffTwoBodyJastrowOrbital.h"
#include "QMCWaveFunctions/Jastrow/OneBodyJastrowSpinOrbital.h"

#include "LongRange/LRHandlerBase.h"
#include "QMCWaveFunctions/Jastrow/SplineFunctors.h"
#include "QMCWaveFunctions/Jastrow/LRBreakupUtilities.h"
#include "LongRange/LRRPAHandlerTemp.h"
#include "QMCWaveFunctions/Jastrow/LRBreakupUtilities.h"



namespace qmcplusplus
{

RadialJastrowBuilder::RadialJastrowBuilder(ParticleSet& target, TrialWaveFunction& psi,
					   ParticleSet& source):
  WaveFunctionComponentBuilder(target, psi),SourcePtcl(source)
{
  ClassName="RadialJastrowBuilder";
  NameOpt="0";
  typeOpt="unknown";
  Jastfunction="unknown";
  SourceOpt=targetPtcl.getName();
  SpinOpt="no";
}

RadialJastrowBuilder::RadialJastrowBuilder(ParticleSet& target, TrialWaveFunction& psi):
  WaveFunctionComponentBuilder(target, psi),SourcePtcl(NULL)
{
  ClassName="RadialJastrowBuilder";
  NameOpt="0";
  typeOpt="unknown";
  Jastfunction="unknown";
  SourceOpt=targetPtcl.getName();
  SpinOpt="no";
}

// helper method for dealing with functor incompatible with Open Boundaries
void RadialJastrowBuilder::guardAgainstOBC()
{
  if (targetPtcl.Lattice.SuperCellEnum == SUPERCELL_OPEN) 
  {
    app_error() << Jastfunction << " relies on the total density for its form\n";
    app_error() << "but open boundary conditions are requested.  Please choose other forms of Jastrow\n";
  }
}
 
// helper method for dealing with functor incompatible with PBC
void RadialJastrowBuilder::guardAgainstPBC()
{
  if (targetPtcl.Lattice.SuperCellEnum != SUPERCELL_OPEN)
  {
    app_error() << Jastfunction << " does not support a cutoff, but is requested with\n";
    app_error() << "periodic boundary conditions, please choose other forms of Jastrow\n";
  }
}

// quick template helper to allow use of RPA
template <typename> 
class RPAFunctor { };

// helper class to simplify and localize ugly ifdef stuff for types
template<template<class> class RadFuncType>
class JastrowTypeHelper 
{
private:
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
template<>
class JastrowTypeHelper<BsplineFunctor>
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
#endif
#if defined(QMC_CUDA) and !defined(ENABLE_SOA)
  typedef OneBodyJastrowOrbitalBsplineAoS J1OrbitalType;
  typedef OneBodyJastrowOrbitalBsplineAoS J1OrbitalTypeSpin;
  typedef DiffOneBodySpinJastrowOrbital<rft> DiffJ1OrbitalType;
  typedef DiffOneBodySpinJastrowOrbital<rft> DiffJ1OrbitalTypeSpin;
  typedef TwoBodyJastrowOrbitalBsplineAoS J2OrbitalType;
  typedef DiffTwoBodyJastrowOrbital<rft> DiffJ2OrbitalType;
#endif
#if !defined(QMC_CUDA) and defined(ENABLE_SOA)
  typedef J1OrbitalSoA<rft> J1OrbitalType;
  typedef OneBodySpinJastrowOrbital<rft> J1OrbitalTypeSpin;
  typedef DiffOneBodyJastrowOrbital<rft> DiffJ1OrbitalType;
  typedef DiffOneBodySpinJastrowOrbital<rft> DiffJ1OrbitalTypeSpin;
  typedef J2OrbitalSoA<rft> J2OrbitalType;
  typedef DiffTwoBodyJastrowOrbital<rft> DiffJ2OrbitalType;
#endif
#if !defined(QMC_CUDA) and !defined(ENABLE_SOA)
  typedef OneBodyJastrowOrbital<rft> J1OrbitalType;
  typedef OneBodySpinJastrowOrbital<rft> J1OrbitalTypeSpin;
  typedef DiffOneBodyJastrowOrbital<rft> DiffJ1OrbitalType;
  typedef DiffOneBodySpinJastrowOrbital<rft> DiffJ1OrbitalTypeSpin;
  typedef TwoBodyJastrowOrbital<rft> J2OrbitalType;
  typedef DiffTwoBodyJastrowOrbital<rft> DiffJ2OrbitalType;
#endif
};


template<template<class> class RadFuncType>
void RadialJastrowBuilder::initTwoBodyFunctor(RadFuncType<RT>* functor, double fac) { }

template<>
void RadialJastrowBuilder::initTwoBodyFunctor(BsplineFunctor<RT>* bfunc, double fac) 
{
  if(P.Lattice.SuperCellEnum==SUPERCELL_OPEN) // for open systems, do nothing
  {
    return;
  }
  std::vector<RT> rpaValues;
  int npts=bfunc.NumParams;
  if(rpaValues.empty())
  {
    rpaValues.resize(npts);
    LRRPAHandlerTemp<RPABreakup<T>,LPQHIBasis> rpa(P,-1.0);
    rpa.Breakup(P,-1.0);
    RT dr=bfunc.cutoff_radius/static_cast<T>(npts);
    RT r=0;
    for (int i=0; i<npts; i++)
    {
      rpaValues[i]=rpa.evaluate(r,1.0/r); //y[i]=fac*rpa.evaluate(r,1.0/r);
      r += dr;
    }
  }
  RT last=rpaValues[npts-1];

  for(int i=0; i<npts; i++)
    bfunc.Parameters[i]=fac*(rpaValues[i]-last);
  bfunc.reset();
}



template<template<class> class RadFuncType>
bool RadialJastrowBuilder::createJ2(xmlNodePtr cur) 
{
  typedef RadFuncType<RT> RadFunctorType;
  typedef JastrowTypeHelper<RadFuncType>::J2OrbitalType J2OrbitalType;
  typedef JastrowTypeHelper<RadFuncType>::DiffJ2OrbitalType DiffJ2OrbitalType;
  
  SpeciesSet& species(targetPtcl.getSpeciesSet());
  int taskid=(targetPsi.is_manager())?targetPsi.getGroupID():-1;
  auto *J2 = new J2OrbitalType(targetPtcl,taskid);
  auto *dJ2 = new DiffJ2OrbitalType(targetPtcl);

  std::string init_mode("0");
  {
    OhmmsAttributeSet hAttrib;
    hAttrib.add(init_mode,"init");
    hAttrib.put(cur);
  }

  cur = cur->xmlChildrenNode;
  while(cur!=NULL)
  {
    std::string cname((const char*)cur->name);
    if (cname == "correlation")
    {
      OhmmsAttributeSet rAttrib;
      RealType cusp=-1e10;
      std::string spA(species.speciesName[0]);
      std::string spB(species.speciesName[0]);
      std::string pairType("0");
      rAttrib.add(spA,"speciesA");
      rAttrib.add(spB,"speciesB");
      rAttrib.add(pairType,"pairType");
      rAttrib.add(cusp,"cusp");
      rAttrib.put(cur);
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
      if(ia==ib && (targetPtcl.last(ia)-targetPtcl.first(ia)==1))
	PRE.error("Failed to add "+spA+spB+" correlation for only 1 "+spA+" particle. Please remove it from two-body Jastrow.",true);
      if(cusp<-1e6)
      {
	RealType qq=species(chargeInd,ia)*species(chargeInd,ib);
	cusp = (ia==ib)? -0.25*qq:-0.5*qq;
      }
      app_log() << "  RadialJastrowBuilder adds a functor with cusp = " << cusp << std::endl;

      auto *functor = new RadFunctorType();
      functor->setCusp(cusp);
      functor->setPeriodic(targetPtcl.Lattice.SuperCellEnum != SUPERCELL_OPEN);
      if (targetPtcl.Lattice.WignerSeitzRadius > 0) 
      {
	functor->cutoff_radius = targetPtcl.Lattice.WignerSeitzRadius;
      }
      else if (functor->cutoff_radius < 10.0)
      {
	functor->cutoff_radius = 10.0;
      }
      functor->put(cur);
      initTwoBodyFunctor(functor,-cusp/0.5);

      J2->addFunc(ia,ib,functor);
      dJ2->addFunc(ia,ib,functor);

      if(qmc_common.io_node)
      {
	char fname[32];
	if(qmc_common.mpi_groups>1)
	  sprintf(fname,"J2.%s.g%03d.dat",pairType.c_str(),taskid);
	else
	  sprintf(fname,"J2.%s.dat",pairType.c_str());
	ostream os(fname);
	print(*func, os);	
      }
    }
    cur=cur->next;
  }
  J2->dPsi=dJ2;
  std::string j2name="J2_"+Jastname;
  targetPsi.addOrbital(J2,j2name.c_str());
  J2->setOptimizable(true);
}

// specialiation for J2 RPA jastrow.  Note that the long range part is not implemented
template<>
bool RadialJastrowBuilder<RPAFunctor>::createJ2(xmlNodePtr cur) 
{
  RPAJastrow* rpajastrow = new RPAJastrow(targetPtcl,targetPsi.is_manager());
  rpajastrow->put(cur);
  targetPsi.addOrbital(rpajastrow,NameOpt);
  return true;
}

template<template<class> class RadFuncType>
bool RadialJastrowBuilder::createJ1(xmlNodePtr cur) 
{
  typedef RadFuncType<RT> RadFunctorType;
  typedef JastrowTypeHelper<RadFuncType>::J1OrbitalType J1OrbitalType;
  typedef JastrowTypeHelper<RadFuncType>::DiffJ1OrbitalType DiffJ1OrbitalType;
   
  int taskid=targetPsi.getGroupID();
  J1OrbitalType* J1= new J1OrbitalType(*sourcePtcl, targetPtcl);
  DiffJ1OrbitalType* dJ1= new DiffJ1OrbitalType(*sourcePtcl, targetPtcl);

  xmlNodePtr kids = cur->xmlChildrenNode;

  // Find the number of the source species
  SpeciesSet &sSet = sourcePtcl->getSpeciesSet();
  SpeciesSet &tSet = targetPtcl.getSpeciesSet();
  int numSpecies = sSet.getTotalNum();
  bool success=false;
  bool Opt(true);
  while (kids != NULL)
  {
    std::string kidsname = (char*)kids->name;
    tolower(kidsname);
    if (kidsname == "correlation")
    {
      RealType cusp=0.0;
      std::string speciesA;
      std::string speciesB;
      OhmmsAttributeSet rAttrib;
      rAttrib.add(speciesA,"elementType");
      rAttrib.add(speciesA,"speciesA");
      rAttrib.add(speciesB,"speciesB");
      rAttrib.put(kids);
      auto *functor = new RadFunctorType();
      int ig = sSet.findSpecies (speciesA);
      
      functor->setPeriodic(sourcePtcl->Lattice.SuperCellEnum != SUPERCELL_OPEN);
      if (targetPtcl.Lattice.WignerSeitzRadius > 0) 
      {
	functor->cutoff_radius = targetPtcl.Lattice.WignerSeitzRadius;
      }
      else if (functor->cutoff_radius < 10.0)
      {
	functor->cutoff_radius = 10.0;
      }
      int jg=-1;
      if(speciesB.size())
	jg=tSet.findSpecies(speciesB);
      if (ig < numSpecies)
      {
	functor->put(kids);
	J1->addFunc(ig,functor,jg);
	dJ1->addFunc(ig,functor,jg);
	success = true;
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
	  ostream os(fname);
	  print(*func, os);
	}
      }
    }
    kids = kids->next;
  }
  if(success)
  {
    J1->dPsi=dJ1;
    std::string jname = "J1_"+Jastfunction;
    targetPsi.addOrbital(J1,jname.c_str());
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

// specialiation for J1 RPA jastrow.  Note that the long range part is not implemented
template<>
bool RadialJastrowBuilder<RPAFunctor>::createJ1(xmlNodePtr cur) 
{
  typedef CubicBspline<RealType,LINEAR_1DGRID,FIRSTDERIV_CONSTRAINTS> SplineEngineType;
  typedef CubicSplineSingle<RT,SplineEngineType> RadFunctorType;
  typedef LinearGrid<RealType> GridType;
  typedef LRHandlerBase HandlerType;
  typedef JastrowTypeHelper<RadFuncType>::J1OrbitalType J1OrbitalType;
  typedef JastrowTypeHelper<RadFuncType>::DiffJ1OrbitalType DiffJ1OrbitalType;

  std::string MyName="Jep";
  std::string rpafunc="RPA";
  OhmmsAttributeSet a;
  a.add(MyName,"name");
  a.add(rpafunc,"function");
  a.put(cur);
  ParameterSet params;
  RealType Rs(-1.0);
  RealType Kc(-1.0);
  params.add(Rs,"rs","double");
  params.add(Kc,"kc","double");
  params.put(cur);

  HandlerType* myHandler;
  if(Rs<0)
  {
    Rs=std::pow(3.0/4.0/M_PI*targetPtcl.Lattice.Volume/ static_cast<RealType>(targetPtcl.getTotalNum()) ,1.0/3.0);
  }
  if(Kc<0)
  {
    Kc = 1e-6;
  }
  if (rpafunc=="RPA")
  {
    myHandler= new LRRPAHandlerTemp<EPRPABreakup<RealType>,LPQHIBasis>(targetPtcl,Kc);
    app_log()<<"  using e-p RPA"<< std::endl;
  }
  else if (rpafunc=="dRPA")
  {
    myHandler= new LRRPAHandlerTemp<derivEPRPABreakup<RealType>,LPQHIBasis>(targetPtcl,Kc);
    app_log()<<"  using e-p derivRPA"<< std::endl;
  }
  myHandler->Breakup(targetPtcl,Rs);
  
  Rcut = myHandler->get_rc()-0.1;
  GridType* myGrid = new GridType;
  int npts=static_cast<int>(Rcut/0.01)+1;
  myGrid->set(0,Rcut,npts);

  //create the numerical functor
  RadFunctorType* nfunc = new RadFunctorType;
  ShortRangePartAdapter<RT>* SRA = new ShortRangePartAdapter<RT>(myHandler);
  SRA->setRmax(Rcut);
  nfunc->initialize(SRA, myGrid);
  
  J1OrbitalType* J1= new J1OrbitalType(*sourcePtcl, targetPtcl);
  DiffJ1OrbitalType* dJ1= new DiffJ1OrbitalType(*sourcePtcl, targetPtcl);
  
  SpeciesSet &sSet = sourcePtcl->getSpeciesSet();
  for (int ig=0; ig< sSet.getTotalNum(); ig++) {
    J1->addFunc(ig,nfunc);
    dJ1->addFunc(ig,nfunc);
  }
  
  J1->dPsi=dJ1;
  std::string jname = "J1_"+Jastfunction;
  targetPsi.addOrbital(J1,jname.c_str());
  J1->setOptimizable(Opt);
  return true;
}


bool RadialJastrowBuilder::put(xmlNodePtr cur)
{
  ReportEngine PRE(ClassName,"put(xmlNodePtr)");
  OhmmsAttributeSet aAttrib;
  aAttrib.add(NameOpt,"name");
  aAttrib.add(TypeOpt,"type");
  aAttrib.add(Jastfunction,"function");
  aAttrib.add(SourceOpt, "source");
  aAttrib.add(SpinOpt, "spin");
  aAttrib.put(cur);
  tolower(NameOpt);
  tolower(TypeOpt);
  tolower(Jastfunction);
  tolower(SourceOpt);
  tolower(SpinOpt);

  bool success=false;

  SpeciesSet& species(targetPtcl.getSpeciesSet());
  int chargeInd=species.addAttribute("charge");  

  if (typeOpt.find("one") < typeOpt.size())
  {
    // it's a one body jastrow factor
    if (Jastfunction == "bspline") 
    {
      success = createJ1<BsplineFunctor>(cur);
    }
    else if (Jastfunction == "pade") 
    {
      guardAgainstPBC();
      success = createJ1<PadeFunctor>(cur);
    }
    else if (Jastfunction == "rpa") 
    {
#if !(OHMMS_DIM == 3)
      app_error() << "RPA for one-body jastrow is only available for 3D\n";
#endif
      guardAgainstOBC();
      success = createJ1<RPAFunctor>(cur);
    }
    else
    {
      app_error() << "Unknown one jastrow body function: " << Jastfunction << ".\n";
    }
  }
  else if (typeOpt.find("two") < typeOpt.size())
  {
    // it's a two body jastrow factor
    if (Jastfunction == "bspline") 
    {
      success = createJ2<BsplineFunctor>(cur);
    }
    else if (Jastfunction == "pade") 
    {
      guardAgainstPBC();
      success = createJ2<PadeFunctor>(cur);
    }
    else if (Jastfunction == "rpa" || Jastfunction = "yukawa") 
    {
#if !(OHMMS_DIM == 3)
      app_error() << "RPA for one-body jastrow is only available for 3D\n";
#endif
      guardAgainstOBC();
      success = createJ2<RPAFunctor>(cur);
    }
    else
    {
      app_error() << "Unknown two body jastrow function: " << Jastfunction << ".\n";
    }
    }
}      
