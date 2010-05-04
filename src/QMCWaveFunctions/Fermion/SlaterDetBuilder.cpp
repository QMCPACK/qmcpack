//////////////////////////////////////////////////////////////////
// (c) Copyright 2003-  by Jeongnim Kim
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//   National Center for Supercomputing Applications &
//   Materials Computation Center
//   University of Illinois, Urbana-Champaign
//   Urbana, IL 61801
//   e-mail: jnkim@ncsa.uiuc.edu
//
// Supported by
//   National Center for Supercomputing Applications, UIUC
//   Materials Computation Center, UIUC
//////////////////////////////////////////////////////////////////
#include "QMCWaveFunctions/BasisSetFactory.h"
#include "QMCWaveFunctions/Fermion/SlaterDetBuilder.h"
#include "Utilities/ProgressReportEngine.h"
#include "OhmmsData/AttributeSet.h"

#include "QMCWaveFunctions/MultiSlaterDeterminant.h"
//this is only for Bryan
#if defined(BRYAN_MULTIDET_TRIAL)
#include "QMCWaveFunctions/Fermion/DiracDeterminantIterative.h"
#include "QMCWaveFunctions/Fermion/DiracDeterminantTruncation.h"
#include "QMCWaveFunctions/Fermion/MultiDiracDeterminantBase.h"
#endif
//Cannot use complex and released node
#if !defined(QMC_COMPLEX)
#include "QMCWaveFunctions/Fermion/RNDiracDeterminantBase.h"
#include "QMCWaveFunctions/Fermion/RNDiracDeterminantBaseAlternate.h"
#endif
#ifdef QMC_CUDA
  #include "QMCWaveFunctions/Fermion/DiracDeterminantCUDA.h"
#endif
#if QMC_BUILD_LEVEL>2 && OHMMS_DIM==3
#include "QMCWaveFunctions/Fermion/SlaterDetWithBackflow.h"
#include "QMCWaveFunctions/Fermion/DiracDeterminantWithBackflow.h"
#endif
#include "QMCWaveFunctions/Fermion/DiracDeterminantOpt.h"


namespace qmcplusplus
{

  SlaterDetBuilder::SlaterDetBuilder(ParticleSet& els, TrialWaveFunction& psi,
      PtclPoolType& psets)
    : OrbitalBuilderBase(els,psi), ptclPool(psets)
      , myBasisSetFactory(0), slaterdet_0(0)
  {
    ClassName="SlaterDetBuilder";
#if QMC_BUILD_LEVEL>2 && OHMMS_DIM==3
    BFTrans=0;
    UseBackflow=false;
#endif
  }
  SlaterDetBuilder::~SlaterDetBuilder()
  {
    DEBUG_MEMORY("SlaterDetBuilder::~SlaterDetBuilder");
    if (myBasisSetFactory)
      {
        delete myBasisSetFactory;
      }
  }


  /** process <determinantset>
   *
   * determinantset 
   * - basiset 0..1
   * - sposet 0..* 
   * - slaterdeterminant
   *   - determinant 0..*
   * - ci
   */
  bool SlaterDetBuilder::put(xmlNodePtr cur)
  {

    ReportEngine PRE(ClassName,"put(xmlNodePtr)");

    ///save the current node
    xmlNodePtr curRoot=cur;
    xmlNodePtr BFnode;
    bool success=true;
    string cname, tname;
    std::map<string,SPOSetBasePtr> spomap;

    //check the basis set
    cur = curRoot->children;
    while (cur != NULL)//check the basis set
    {
      getNodeName(cname,cur);
      if (cname == basisset_tag)
      {
        if (myBasisSetFactory == 0)
        {
          myBasisSetFactory = new BasisSetFactory(targetPtcl,targetPsi, ptclPool);
          myBasisSetFactory->setReportLevel(ReportLevel);
        }
        myBasisSetFactory->createBasisSet(cur,curRoot);
      }
      else if (cname == sposet_tag) {
	app_log() << "Creating SPOSet in SlaterDetBuilder::put(xmlNodePtr cur).\n";
	string spo_name;
	OhmmsAttributeSet spoAttrib;
	spoAttrib.add (spo_name, "name");
	spoAttrib.put(cur);
	app_log() << "spo_name = " << spo_name << endl;

	if (myBasisSetFactory == 0)
        {
          myBasisSetFactory = new BasisSetFactory(targetPtcl,targetPsi, ptclPool);
          myBasisSetFactory->setReportLevel(ReportLevel);
        }
        myBasisSetFactory->createBasisSet(cur,cur);
	SPOSetBasePtr spo = myBasisSetFactory->createSPOSet(cur);
	spo->put(cur, spomap);
	if (spomap.find(spo_name) != spomap.end()) {
	  app_error() << "SPOSet name \"" << spo_name << "\" is already in use.\n";
	  abort();
	}
	spomap[spo_name] = spo;
	assert(spomap.find(spo_name) != spomap.end());
	//	slaterdet_0->add(spo,spo_name);
      }
#if QMC_BUILD_LEVEL>2 && OHMMS_DIM==3 
      else if(cname == backflow_tag) {
        app_log() <<"Creating Backflow transformation in SlaterDetBuilder::put(xmlNodePtr cur).\n";

        // to simplify the logic inside DiracDeterminantWithBackflow,
        // I'm requiring that only a single <backflow> block appears
        // in the xml file
        if(BFTrans != 0) {
           APP_ABORT("Only a single backflow block is allowed in the xml. Please collect all transformations into a single block. \n");
        }
        UseBackflow=true;
        // creating later due to problems with ParticleSets
        //BFTrans = new BackflowTransformation(targetPtcl,ptclPool);
        BFTrans = NULL;
        BFnode=cur;
// read xml later, in case some ParticleSets are read from hdf5 file.
        //BFTrans->put(cur);  
      }
#endif
      cur = cur->next;
    }

    //missing basiset, e.g. einspline
    if (myBasisSetFactory == 0)
    {
      myBasisSetFactory = new BasisSetFactory(targetPtcl,targetPsi, ptclPool);
      myBasisSetFactory->setReportLevel(ReportLevel);
      myBasisSetFactory->createBasisSet(curRoot,curRoot);
    }

    //add sposet
    
    cur = curRoot->children;
    while (cur != NULL)//check the basis set
    {
      getNodeName(cname,cur);
      if (cname == sd_tag)
      {
        if(slaterdet_0)
        {
          APP_ABORT("slaterdet is already instantiated.");
        }
#if QMC_BUILD_LEVEL>2 && OHMMS_DIM==3
        if(UseBackflow) 
          slaterdet_0 = new SlaterDetWithBackflow(targetPtcl,BFTrans);
        else 
#endif
          slaterdet_0 = new SlaterDeterminant_t(targetPtcl);
	// Copy any entries in sposetmap into slaterdet_0
	std::map<string,SPOSetBasePtr>::iterator iter;
	for (iter=spomap.begin(); iter!=spomap.end(); iter++) {
	  cerr << "Adding SPO \"" << iter->first << "\" to slaterdet_0.\n";
	  slaterdet_0->add(iter->second,iter->first);
	}


        int spin_group = 0;
        xmlNodePtr tcur = cur->children;
        while (tcur != NULL)
        {
          getNodeName(tname,tcur);
          if (tname == det_tag || tname == rn_tag)
          {
            if(putDeterminant(tcur, spin_group)) spin_group++;
          }
          tcur = tcur->next;
        }
      }
      cur = cur->next;
    }
    

    if (!slaterdet_0)
    {
      //fatal
      PRE.error("Failed to create a SlaterDeterminant.",true);
      return false;
    }

    // change DistanceTables if using backflow
#if QMC_BUILD_LEVEL>2 && OHMMS_DIM==3
    if(UseBackflow)   { 
       BFTrans = new BackflowTransformation(targetPtcl,ptclPool);
  // HACK HACK HACK, until I figure out a solution      
       SlaterDetWithBackflow* tmp = (SlaterDetWithBackflow*) slaterdet_0;
       tmp->BFTrans = BFTrans;
       for(int i=0; i<tmp->Dets.size(); i++) {
         DiracDeterminantWithBackflow* tmpD = (DiracDeterminantWithBackflow*) tmp->Dets[i]; 
         tmpD->BFTrans = BFTrans;
       }
       BFTrans->put(BFnode);
       tmp->resetTargetParticleSet(BFTrans->QP);
    }
#endif
    //only single slater determinant
    targetPsi.addOrbital(slaterdet_0,"SlaterDet");

    //buildMultiSlaterDetermiant();

    delete myBasisSetFactory;
    myBasisSetFactory=0;

    return success;
  }


  bool SlaterDetBuilder::putDeterminant(xmlNodePtr cur, int spin_group)
  {

    ReportEngine PRE(ClassName,"putDeterminant(xmlNodePtr,int)");

    string basisName("invalid");
    string detname("0"), refname("0");
    string s_detSize("0");
    string detMethod("");
    OhmmsAttributeSet aAttrib;
    aAttrib.add(basisName,basisset_tag);
    aAttrib.add(detname,"id");
    aAttrib.add(refname,"ref");
    aAttrib.add(detMethod,"DetMethod");
    aAttrib.add(s_detSize,"DetSize");

    string s_cutoff("0.0");
    string s_radius("0.0");
    int s_smallnumber(-999999);
    int rntype(0);
    aAttrib.add(s_cutoff,"Cutoff");
    aAttrib.add(s_radius,"Radius");
    aAttrib.add(s_smallnumber,"smallnumber");
    aAttrib.add(s_smallnumber,"eps");
    aAttrib.add(rntype,"primary");
    aAttrib.put(cur);


    map<string,SPOSetBasePtr>& spo_ref(slaterdet_0->mySPOSet);
    map<string,SPOSetBasePtr>::iterator lit(spo_ref.find(detname));
    SPOSetBasePtr psi;
    if (lit == spo_ref.end())
    {
      // cerr << "Didn't find sposet named \"" << detname << "\"\n";
#if defined(ENABLE_SMARTPOINTER)
      psi.reset(myBasisSetFactory->createSPOSet(cur));
#else
      psi = myBasisSetFactory->createSPOSet(cur);
#endif
      psi->put(cur);
      psi->checkObject();
      slaterdet_0->add(psi,detname);
      //SPOSet[detname]=psi;
      app_log() << "  Creating a new SPO set " << detname << endl;
    }
    else
    {
      app_log() << "  Reusing a SPO set " << detname << endl;
      psi = (*lit).second;
    }

    int firstIndex=targetPtcl.first(spin_group);
    int lastIndex=targetPtcl.last(spin_group);
    if(firstIndex==lastIndex) return true;

//    app_log() << "  Creating DiracDeterminant " << detname << " group=" << spin_group << " First Index = " << firstIndex << endl;
//    app_log() <<"   My det method is "<<detMethod<<endl;
//#if defined(BRYAN_MULTIDET_TRIAL)
//    if (detMethod=="Iterative")
//    {
//      //   string s_cutoff("0.0");
//      //   aAttrib.add(s_cutoff,"Cutoff");
//      app_log()<<"My cutoff is "<<s_cutoff<<endl;
//
//      double cutoff=std::atof(s_cutoff.c_str());
//      DiracDeterminantIterative *adet= new DiracDeterminantIterative(psi,firstIndex);
//      adet->set_iterative(firstIndex,psi->getOrbitalSetSize(),cutoff);
//      slaterdet_0->add(adet,spin_group);
//    }
//    else if (detMethod=="Truncation")
//    {
//      //   string s_cutoff("0.0");
//      //   aAttrib.add(s_cutoff,"Cutoff");
//      DiracDeterminantTruncation *adet= new DiracDeterminantTruncation(psi,firstIndex);
//      double cutoff=std::atof(s_cutoff.c_str());
//      double radius=std::atof(s_radius.c_str());
//      //   adet->set(firstIndex,psi->getOrbitalSetSize());
//      adet->set_truncation(firstIndex,psi->getOrbitalSetSize(),cutoff,radius);
//      slaterdet_0->add(adet,spin_group);
//    }
//    else if (detMethod=="Multi")
//    {
//      app_log()<<"BUILDING DIRAC DETERM "<<firstIndex<<endl;
//      MultiDiracDeterminantBase *adet = new MultiDiracDeterminantBase(psi,firstIndex);
//      int detSize=std::atof(s_detSize.c_str());
//      adet-> set_Multi(firstIndex,detSize,psi->getOrbitalSetSize());
//      slaterdet_0->add(adet,spin_group);
//    }
//    else
//      slaterdet_0->add(new Det_t(psi,firstIndex),spin_group);
//    }
//#else
    string dname;
    getNodeName(dname,cur);
    DiracDeterminantBase* adet=0;
#if !defined(QMC_COMPLEX)
    if (rn_tag == dname)
    {
      double bosonicEpsilon=s_smallnumber;
      app_log()<<"  BUILDING Released Node Determinant logepsilon="<<bosonicEpsilon<<endl;
      if (rntype==0)
        adet = new RNDiracDeterminantBase(psi,firstIndex);
      else
        adet = new RNDiracDeterminantBaseAlternate(psi,firstIndex);
      adet->setLogEpsilon(bosonicEpsilon);  
    }
    else
#endif
    {
#ifdef QMC_CUDA
      adet = new DiracDeterminantCUDA(psi,firstIndex);
#else
#if QMC_BUILD_LEVEL>2 && OHMMS_DIM==3
      if(UseBackflow) 
        adet = new DiracDeterminantWithBackflow(psi,BFTrans,firstIndex);
      else 
#endif
	if (psi->Optimizable)
	  adet = new DiracDeterminantOpt(targetPtcl, psi, firstIndex);
	else
	  adet = new DiracDeterminantBase(psi,firstIndex);
#endif
    }
    adet->set(firstIndex,lastIndex-firstIndex);
    slaterdet_0->add(adet,spin_group);
    if (psi->Optimizable)
      slaterdet_0->Optimizable = true;
    return true;
  }

 // void SlaterDetBuilder::buildMultiSlaterDetermiant()
 // {
 //   MultiSlaterDeterminant *multidet= new MultiSlaterDeterminant;
 //   for (int i=0; i<SlaterDetSet.size(); i++)
 //     {
 //       multidet->add(SlaterDetSet[i],sdet_coeff[i]);
 //     }
 //   // multidet->setOptimizable(true);
 //   //add a MultiDeterminant to the trial wavefuntion
 //   targetPsi.addOrbital(multidet,"MultiSlateDet");
 // }
}
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$
 ***************************************************************************/
