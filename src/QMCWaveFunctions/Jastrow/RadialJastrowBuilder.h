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

#ifndef QMCPLUSPLUS_RADIAL_JASTROW_BUILDER_H
#define QMCPLUSPLUS_RADIAL_JASTROW_BUILDER_H

#include "QMCWaveFunctions/WaveFunctionComponentBuilder.h"

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

#include "Utilities/IteratorUtility.h"
#include "Utilities/ProgressReportEngine.h"
#include "OhmmsData/AttributeSet.h"
#include "OhmmsData/ParameterSet.h"


#include <string>
#include <vector>


namespace qmcplusplus
{

/** JastrowBuilder using an analytic 1d functor
 * Should be able to eventually handle all one and two body jastrows
 * although spline based ones will come later
 */


struct RadialJastrowBuilder: public WaveFunctionComponentBuilder
{
public:
  // one body constructor
  RadialJastrowBuilder(ParticleSet& target, TrialWaveFunction& psi, ParticleSet& source);
  // two body constructor
  RadialJastrowBuilder(ParticleSet& target, TrialWaveFunction& psi);
  bool put(xmlNodePtr cur);
  
 private:
  typedef WaveFunctionComponent::RealType RT;
  
  ///jastrow/@name 
  std::string NameOpt;
  ///jastrow/@type
  std::string TypeOpt;
  ///jastrow/@function
  std::string Jastfunction;
  ///jastrow/@source
  std::string SourceOpt;
  ///jastrow/@spin
  std::string SpinOpt;
  ///particle set for source particle
  ParticleSet *SourcePtcl;

  // has a specialization for RPAFunctor in cpp file
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
	RadFunctorType *functor = new RadFunctorType();
	int ig = sSet.findSpecies (speciesA);

	/* Need to figure out how to either put this in the put 
	   function or some specialized initializer method...
	functor->periodic = sourcePtcl->Lattice.SuperCellEnum != SUPERCELL_OPEN;
	*/
	functor->cutoff_radius = sourcePtcl->Lattice.WignerSeitzRadius;
	int jg=-1;
	if(speciesB.size())
	  jg=tSet.findSpecies(speciesB);
	if (ig < numSpecies)
	{
	  functor->put(kids);
	  J1->addFunc(ig,functor,jg);
	  dJ1->addFunc(ig,functor,jg);
	  success = true;
	  // for debugging
	  /* to keep this, will need to make setReportLevel and print pure virtual
	     functions in the OptimizableFunctorBase class
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
	    functor->setReportLevel(rank()==0,fname);
	    functor->print();
	  }
	  */
        }
      }
      kids = kids->next
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

  // may have a specialization in the cpp file
  template<template<class> class RadFuncType>
  bool RadialJastrowBuilder::createJ2(xmlNodePtr cur) 
  {
    typedef RadFuncType<RT> RadFunctorType;
    typedef JastrowTypeHelper<RadFuncType>::J2OrbitalType J2OrbitalType;
    typedef JastrowTypeHelper<RadFuncType>::DiffJ2OrbitalType DiffJ2OrbitalType;

    SpeciesSet& species(targetPtcl.getSpeciesSet());
    int taskid=(targetPsi.is_manager())?targetPsi.getGroupID():-1;
    J2OrbitalType *J2 = new J2OrbitalType(targetPtcl,taskid);
    DiffJ2OrbitalType *dJ2 = new DiffJ2OrbitalType(targetPtcl);
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
	// need to figure out how to set the cusp
	// fairly sure now that setCusp needs to be in the base class
	RadFunctorType *functor = new RadFunctorType();
	functor->put(cur);
	J2->addFunc(ia,ib,functor);
	dJ2->addFunc(ia,ib,functor);
      }
      cur=cur->next;
    }
    J2->dPsi=dJ2;
    std::string j2name="J2_"+Jastname;
    targetPsi.addOrbital(J2,j2name.c_str());
    J2->setOptimizable(true);
  }


  void guardAgainstOBC();
  void guardAgainstPBC();

private:
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

};



#endif
