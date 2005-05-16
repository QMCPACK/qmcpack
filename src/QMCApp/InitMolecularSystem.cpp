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

    int nel(els->getTotalNum()), items(0);
    while(nel) {
      els->R[items]=0.5*chi[items]; 
      --nel; ++items;
    }
  }

  void InitMolecularSystem::initMolecule(ParticleSet* ions, ParticleSet* els) {

    DistanceTableData* d_ii = DistanceTable::getTable(DistanceTable::add(*ions,*ions));
    d_ii->create(1);
    d_ii->evaluate(*ions);

    vector<double> Cut, Core;
    int Centers = ions->getTotalNum();

    int nattrib=ions->Species.numAttributes();
    Cut.resize(Centers);

    //attribute id for cut
    int icut = ions->Species.addAttribute("cut");
    if(icut>nattrib) {
      for(int iat=0; iat<Centers; iat++) { Cut[iat] = 1.0; }
    } else {
      //store the max distance from atom
      for(int iat=0; iat<Centers; iat++) {
        Cut[iat] = ions->Species(icut,ions->GroupID[iat]);
      }
    }
 
    Core.resize(Centers);
    //use charge as the core electrons first
    int icore = ions->Species.addAttribute("charge");
    for(int iat=0; iat<Centers; iat++) {
      Core[iat]=ions->Species(icore,ions->GroupID[iat]);
    }

    nattrib=ions->Species.numAttributes();
    icore = ions->Species.addAttribute("core");
    //store the max distance from atom
    if(icore < nattrib) {
      for(int iat=0; iat<Centers; iat++) {
        Core[iat]=ions->Species(icore,ions->GroupID[iat]);
      }
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
      for(int iel=0; iel<Core[iat]/2; iel++) {
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
