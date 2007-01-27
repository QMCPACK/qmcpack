//////////////////////////////////////////////////////////////////
// (c) Copyright 2003-  by Jeongnim Kim
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//   National Center for Supercomputing Applications &
//   Materials Computation Center
//   University of Illinois, Urbana-Champaign
//   Urbana, IL 61801
//   e-mail: jnkim@ncsa.uiuc.edu
//   Tel:    217-244-6319 (NCSA) 217-333-3324 (MCC)
//
// Supported by 
//   National Center for Supercomputing Applications, UIUC
//   Materials Computation Center, UIUC
//////////////////////////////////////////////////////////////////
#include "QMCWaveFunctions/Fermion/SlaterDetBuilder.h"
#include "QMCWaveFunctions/BasisSetFactory.h"
#include "OhmmsData/AttributeSet.h"

namespace qmcplusplus {

  SlaterDetBuilder::SlaterDetBuilder(ParticleSet& els, TrialWaveFunction& psi,
     PtclPoolType& psets):
    OrbitalBuilderBase(els,psi), 
    ptclPool(psets), 
    myBasisSetFactory(0)
  { 
  }   

  bool SlaterDetBuilder::put(xmlNodePtr cur){
    ///save the current node
    xmlNodePtr curRoot=cur;
    bool success=true;
    cur = cur->children;
    string cname, tname;
    while(cur != NULL) {
      getNodeName(cname,cur);
      if(cname == basisset_tag) {
        if(myBasisSetFactory == 0) {
          myBasisSetFactory = new BasisSetFactory(targetPtcl,targetPsi, ptclPool);
        }
        myBasisSetFactory->createBasisSet(cur,curRoot);
        //Children.push_back(bsFactory);
      } else if(cname == sd_tag) {
        int is=SlaterDetSet.size();
        //add a new SlaterDeterminant
        SlaterDetSet.push_back(new SlaterDeterminant_t);
        sdet_coeff.push_back(1.0);
        int firstIndex = 0;
        xmlNodePtr tcur = cur->children;
        while(tcur != NULL) {
          getNodeName(tname,tcur);
          if(tname == param_tag) {
            putContent(sdet_coeff[is],tcur);
          } else if(tname == det_tag) {
            firstIndex = putDeterminant(tcur, firstIndex);
          }
          tcur = tcur->next;
        }
      }
      cur = cur->next;
    }
    
    if(SlaterDetSet.empty()) {
      app_error() << "  Failed to create a SlaterDeterminant. Abort at SlaterDetBuilder::put " << endl;
      OHMMS::Controller->abort();
      return false;
    }

    if(SlaterDetSet.size()>1)
      buildMultiSlaterDetermiant();
    else
      buildSlaterDetermiant();
    return success;
  }


  int SlaterDetBuilder::putDeterminant(xmlNodePtr cur, int firstIndex) {

    string basisName("invalid");
    string detname("NONE"), refname("NONE");
    OhmmsAttributeSet aAttrib;
    aAttrib.add(basisName,basisset_tag);
    aAttrib.add(detname,"id");
    aAttrib.add(refname,"ref");
    aAttrib.put(cur);

    xmlNodePtr c_ptr = NULL, o_ptr=NULL;

    Det_t* adet=0;

    //index of the last SlaterDeterminant
    int dIndex=DetSet.size();

    if(refname == "NONE") { //create one and use detname
      if(detname =="NONE") { //no id is given, assign one
        char newname[8];
        sprintf(newname,"det%d",dIndex);
        detname=newname;
      }
    }

    map<string,SPOSetBasePtr>::iterator lit(SPOSet.find(detname));
    SPOSetBasePtr psi;
    if(lit == SPOSet.end()) {
#if defined(ENABLE_SMARTPOINTER)
      psi.reset(myBasisSetFactory->createSPOSet(cur)); 
#else
      psi = myBasisSetFactory->createSPOSet(cur); 
#endif
      psi->put(cur);
      psi->checkObject();
      SPOSet[detname]=psi;
    } else {
      psi = (*lit).second;
    }

    if(psi->getOrbitalSetSize()) {
      map<string,Det_t*>::iterator dit(DetSet.find(detname));
      if(dit == DetSet.end()) {
        adet = new Det_t(psi,firstIndex);
        adet->set(firstIndex,psi->getOrbitalSetSize());
        DetSet[detname]=adet;
      } else {
        adet = (*dit).second;
      }
      firstIndex += psi->getOrbitalSetSize();
    }

    //only if a determinant is not 0
    if(adet) SlaterDetSet.back()->add(adet);
    return firstIndex;
  }

  void SlaterDetBuilder::buildSlaterDetermiant() {
    if(SlaterDetSet.empty()) return;
    //add a SlaterDeterminant to the trial wavefuntion
    targetPsi.addOrbital(SlaterDetSet[0]);
  }

  void SlaterDetBuilder::buildMultiSlaterDetermiant() {
//    MultiSlaterDeterminant<LCOrbitalSet> *multidet= new MultiSlaterDeterminant<LCOrbitalSet>;
//    for(int i=0; i<SlaterDetSet.size(); i++) {
//      multidet->add(SlaterDetSet[i],sdet_coeff[i]);
//    }
//    multidet->setOptimizable(true);
//    //add a MultiDeterminant to the trial wavefuntion
//    targetPsi.addOrbital(multidet);
  }
}
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
