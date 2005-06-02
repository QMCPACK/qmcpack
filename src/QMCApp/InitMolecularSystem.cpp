//////////////////////////////////////////////////////////////////
// (c) Copyright 2003- by Jeongnim Kim
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//   Jeongnim Kim
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
//   Department of Physics, Ohio State University
//   Ohio Supercomputer Center
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
/**@file InitMolecuarSystem.cpp
 * @brief Implements InitMolecuarSystem operators.
 */
#include "QMCApp/InitMolecularSystem.h"
#include "QMCApp/ParticleSetPool.h"
using namespace std;
#include "OhmmsData/AttributeSet.h"
#include "Utilities/OhmmsInfo.h"
#include "Particle/DistanceTable.h"
#include "Particle/DistanceTableData.h"

namespace ohmmsqmc {

  InitMolecularSystem::InitMolecularSystem(ParticleSetPool* pset, 
      const char* aname):
    OhmmsElementBase(aname), ptclPool(pset){ }

  bool InitMolecularSystem::put(xmlNodePtr cur) {

    string target("e"), source("i");

    OhmmsAttributeSet hAttrib;
    hAttrib.add(target,"target");
    hAttrib.add(source,"source");
    hAttrib.put(cur);

    ParticleSet* els=ptclPool->getParticleSet(target);
    if(els == 0) {
      ERRORMSG("No target particle " << target << " exists.")
      return false;
    }

    ParticleSet* ions=ptclPool->getParticleSet(source);
    if(ions == 0) {
      ERRORMSG("No target particle " << target << " exists.")
      return false;
    }

    if(ions->getTotalNum()>1) {
      initMolecule(ions,els);
    } else {
      initAtom(els);
    }
     
    LOGMSG("Checking the new positions for " << target)
    els->get(cout);
    return true;
  }

  void InitMolecularSystem::initAtom(ParticleSet* els) {

    //3N-dimensional Gaussian
    ParticleSet::ParticlePos_t chi(els->getTotalNum());
    makeGaussRandom(chi);

    double q=sqrt(static_cast<double>(els->getTotalNum()))*0.5;
    int nel(els->getTotalNum()), items(0);
    while(nel) {
      els->R[items]=q*chi[items]; 
      --nel; ++items;
    }
  }

  void InitMolecularSystem::initMolecule(ParticleSet* ions, ParticleSet* els) {

    DistanceTableData* d_ii = DistanceTable::getTable(DistanceTable::add(*ions,*ions));
    d_ii->create(1);
    d_ii->evaluate(*ions);

    const ParticleSet::ParticleIndex_t& grID(ions->GroupID);
    SpeciesSet& Species(ions->getSpeciesSet());

    int Centers = ions->getTotalNum();
    vector<int> Qtot(Centers), Qcore(Centers), Qval(Centers,0);
    vector<double> Cut(Centers);

    //use charge as the core electrons first
    int icharge = Species.addAttribute("charge");
    for(int iat=0; iat<Centers; iat++) {
      Qtot[iat]=static_cast<int>(Species(icharge,grID[iat]));
    }

    //the number of core electrons
    int CoreTable[] ={0, /* index zero*/
      0,2,                    /* H He */
      2,2,2,2,2,2,2,10,       /*Li-Ne*/
      10,10,10,10,10,10,10,18 /*Na-Ar*/   
    };

    //Assign default core charge
    for(int iat=0; iat<Centers; iat++) Qcore[iat]=CoreTable[Qtot[grID[iat]]];

    //Overwrite the core charge
    int nattrib=Species.numAttributes();
    int icore=Species.addAttribute("CoreCharge");
    //store the max distance from atom
    if(icore < nattrib) {
      for(int iat=0; iat<Centers; iat++) 
        Qcore[iat]=static_cast<int>(Species(icore,grID[iat]));
    }

    //Add more: used Atomic Radius in AA
    //http://ccinfo.ims.ac.jp/periodic/periodic-main.html
    double CutTable[] = {1.0, /* index zero */
      0.37,1.5,                        /* H He */
      1.52,1.5,1.5,1.5,1.5,1.5,1.5,1.6 /* Li-Ne */
    };

    //search if AtomicRadius is given
    nattrib=Species.numAttributes();
    int icut = Species.addAttribute("AtomicRadius");
    if(icut<nattrib) {//overwrite this
      for(int iat=0; iat<Centers; iat++) Cut[iat] = Species(icut,grID[iat]);
    } else {
      //use default
      for(int iat=0; iat<Centers; iat++) Cut[iat] = CutTable[grID[iat]];
    }
 

    //3N-dimensional Gaussian
    ParticleSet::ParticlePos_t chi(els->getTotalNum());
    makeGaussRandom(chi);

    int numDown=els->getTotalNum()/2;
    int numUp= els->getTotalNum()-numDown;
 
    //assign the core
    int ncoreUp(0), items(0);
    for(int iat=0; iat<Centers; iat++) {
      double sep=0.8*Cut[iat];
      for(int iel=0; iel<Qcore[iat]/2; iel++) {
        els->R[ncoreUp]=ions->R[iat]+sep*chi[items++];
        els->R[ncoreUp+numUp]=ions->R[iat]+sep*chi[items++];
        ++ncoreUp;
      }
    }

    int vup=numUp-ncoreUp;
    int vdown=numDown-ncoreUp;
    int vtot=vup+vdown;
    int ic=0;
    while(vtot && ic<Centers) {
      for(int nn=d_ii->M[ic]; nn<d_ii->M[ic+1]; nn++){
        double bl = d_ii->r(nn);
        int jc = d_ii->J[nn];
        //only assign if the half bond-length < cutoff
        if(vtot && bl < Cut[ic]+Cut[jc]){
          ParticleSet::SingleParticlePos_t displ= ions->R[ic]+0.5*d_ii->dr(nn);
          bl*=0.1;
          if(vup) {
            els->R[ncoreUp] = displ+bl*chi[items];
            --vup;--vtot;++ncoreUp;++items;
          }
          if(vdown) {
            els->R[ncoreUp+numUp]=displ+bl*chi[items];
            --vdown;--vtot;++ncoreUp;++items;
          }
        }
      }
      ++ic;
    }
  }

  bool InitMolecularSystem::put(std::istream& is) {
    return true;
  }

  bool InitMolecularSystem::get(std::ostream& os) const {
    return true;
  }

  void InitMolecularSystem::reset() {
 
  }
}
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
