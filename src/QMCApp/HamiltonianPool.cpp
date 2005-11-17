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
/**@file HamiltonianPool.cpp
 * @brief Implements HamiltonianPool operators.
 */
#include "QMCApp/HamiltonianPool.h"
#include "QMCApp/WaveFunctionPool.h"
#include "QMCApp/ParticleSetPool.h"
using namespace std;
#include "OhmmsData/AttributeSet.h"
#include "Utilities/OhmmsInfo.h"
#include "QMCHamiltonians/QMCHamiltonian.h"
#include "QMCHamiltonians/BareKineticEnergy.h"
#include "QMCHamiltonians/CoulombPotential.h"
#include "QMCHamiltonians/IonIonPotential.h"
#include "QMCHamiltonians/LocalPPotential.h"
#include "QMCHamiltonians/NonLocalPPotential.h"
#include "QMCHamiltonians/LocalCorePolPotential.h"
#if !defined(QMCPLUSPLUS_RELEASE)
#include "QMCHamiltonians/CoulombPBC.h"
#endif

namespace qmcplusplus {

  HamiltonianPool::HamiltonianPool(const char* aname):
    OhmmsElementBase(aname), curH(0), ptclPool(0), psiPool(0){ }

  bool HamiltonianPool::put(xmlNodePtr cur) {

    string id("h0"), htype("generic"), target("e"), source("i"),role("extra");

    OhmmsAttributeSet hAttrib;
    hAttrib.add(id,"id"); hAttrib.add(id,"name"); 
    hAttrib.add(htype,"type"); 
    hAttrib.add(target,"target");
    hAttrib.add(source,"source");
    hAttrib.add(role,"role");
    hAttrib.put(cur);

    curH = getHamiltonian(id);
    if(curH) {
      WARNMSG("hamiltonian with " << id << " is already created. Will add new terms.")
    } else {
      curH = new QMCHamiltonian;
      if(myPool.empty() || role == "primary" ) {
        primaryH = curH;
      }
      myPool[id]=curH;
    }

    ParticleSet* qp=ptclPool->getParticleSet(target);
    if(qp == 0) {
      ERRORMSG("No target particle " << target << " exists.")
      return false;
    }

    //important to add KineticEnergy first
    curH->add(new BareKineticEnergy,"Kinetic");

    //using hamiltonian/@type to populate QMCHamiltonian
    cur=cur->children;
    //check the child nodes of hamiltonian: pairpot 
    while(cur != NULL) {
      string cname((const char*)cur->name);
      const xmlChar* t = xmlGetProp(cur,(const xmlChar*)"type");
      if(t) { // accept only if it has type
        string pot_type((const char*)t);
        if(cname == "pairpot") {
          if(pot_type == "coulomb") {
	    addCoulombPotential(cur,qp);
          } else if(pot_type == "pseudo") {
            addPseudoPotential(cur,qp);
          } else if(pot_type == "cpp") {
            addCorePolPotential(cur,qp);
          }
        } else if(cname == "constant") { 
          if(pot_type == "coulomb") { //ugly!!!
            t = xmlGetProp(cur,(const xmlChar*)"source");
            string nuclei("i");
            if(t) nuclei=(const char*)t;
            ParticleSet* ion=ptclPool->getParticleSet(nuclei);

	    //CHECK PBC and create CoulombPBC for ion-ion
            if(ion) {
              if(ion->getTotalNum()>1) 
		if(ion->Lattice.BoxBConds[0]){
#if defined(QMCPLUSPLUS_RELEASE)
                  ERRORMSG("This version cannot handle PBC. The IonIon potential will be wrong.")
		  curH->add(new IonIonPotential(*ion),"IonIon");
#else
		  curH->add(new CoulombPBCAA(*ion),"IonIon");
#endif
		}
		else {
		  curH->add(new IonIonPotential(*ion),"IonIon");
		}
             }
          }
        }
      }
      cur = cur->next;
    }

    if(curH->size() == 1) {//no external potential is provided, use type

      WARNMSG("Using pre-determined hamiltonian for molecular systems.")
      ParticleSet* ion=ptclPool->getParticleSet(source);
      if(ion == 0) {
        ERRORMSG("No ionic system " << source << " exists.")
        return false;
      }

      if(qp->Lattice.BoxBConds[0]){
#if defined(QMCPLUSPLUS_RELEASE)
        ERRORMSG("This version cannot handle PBC. The ElecElec potential will be wrong.")
	curH->add(new CoulombPotentialAA(*qp),"ElecElec");
#else
	curH->add(new CoulombPBCAA(*qp),"ElecElec");
#endif
      }
      else{
	curH->add(new CoulombPotentialAA(*qp),"ElecElec");
      }

      if(htype == "molecule" || htype=="coulomb"){
	if(qp->Lattice.BoxBConds[0]){
#if defined(QMCPLUSPLUS_RELEASE)
          ERRORMSG("This version cannot handle PBC. The Coulomb potential will be wrong.")
	  curH->add(new CoulombPotentialAB(*ion,*qp),"Coulomb");
#else
	  curH->add(new CoulombPBCAB(*ion,*qp),"Coulomb");
#endif
	} else {
	  curH->add(new CoulombPotentialAB(*ion,*qp),"Coulomb");
	}
      } else if(htype == "siesta" || htype=="pseudo") {
        TrialWaveFunction* psi = psiPool->getPrimary();
        curH->add(new NonLocalPPotential(*ion,*qp,*psi),"NonLocal");
      //} else if(htype == "cpp") {
      //  xmlChar* att2=xmlGetProp(cur,(const xmlChar*)"species");
      //  string stype("Ge");
      //  if(att2) stype = (const char*)att2;
      //  curH->add(new LocalPPotential(*ion,*qp), "PseudoPot");
      //  curH->add(new LocalCorePolPotential(*ion,*qp), "GeCPP");
      } else {
        ERRORMSG(htype << " is diabled")
      }

      if(ion->getTotalNum()>1) 
	if(ion->Lattice.BoxBConds[0]){
#if defined(QMCPLUSPLUS_RELEASE)
          ERRORMSG("This version cannot handle PBC. The IonIon potential will be wrong.")
          curH->add(new IonIonPotential(*ion),"IonIon");
#else
          curH->add(new CoulombPBCAA(*ion),"IonIon");
#endif
	}
	else{
          curH->add(new IonIonPotential(*ion),"IonIon");
	}
    }


    return true;
  }

  void 
  HamiltonianPool::addCoulombPotential(xmlNodePtr cur, ParticleSet* target) {

    string a("e"),title("ElecElec");
    OhmmsAttributeSet hAttrib;
    hAttrib.add(title,"id"); hAttrib.add(title,"name"); 
    hAttrib.add(a,"source"); 
    hAttrib.put(cur);

    ParticleSet* source = ptclPool->getParticleSet(a);
    if(source ==0) {
      ERRORMSG("Missing source ParticleSet" << a)
      return;
    }

    if(source == target) {
      //CHECK PBC and create CoulombPBC for el-el
      if(source->getTotalNum()>1)  {
	  if(target->Lattice.BoxBConds[0]) {
#if defined(QMCPLUSPLUS_RELEASE)
            ERRORMSG("This version cannot handle PBC. The " << title << " potential will be wrong.")
	    curH->add(new CoulombPotentialAA(*target),title);
#else
	    curH->add(new CoulombPBCAA(*target),title);
#endif
	  }
	  else {
	    curH->add(new CoulombPotentialAA(*target),title);
	  }
      }
    } else {
      if(target->Lattice.BoxBConds[0]) {
#if defined(QMCPLUSPLUS_RELEASE)
        ERRORMSG("This version cannot handle PBC. The " << title << " potential will be wrong.")
	curH->add(new CoulombPotentialAB(*source,*target),title);
#else
	curH->add(new CoulombPBCAB(*source,*target),title);
#endif
      } else {
	curH->add(new CoulombPotentialAB(*source,*target),title);
      }
    }
  }
  
  void 
  HamiltonianPool::addPseudoPotential(xmlNodePtr cur, ParticleSet* target) {

    string src("i"),title("PseudoPot"),wfname("invalid");

    OhmmsAttributeSet pAttrib;
    pAttrib.add(title,"name");
    pAttrib.add(src,"source");
    pAttrib.add(wfname,"wavefunction");
    pAttrib.put(cur);

    ParticleSet* ion=ptclPool->getParticleSet(src);
    if(ion == 0) return;

    TrialWaveFunction* psi=0;

    if(wfname == "invalid") {
      psi = psiPool->getPrimary();
    } else {
      psi = psiPool->getWaveFunction(wfname);
    }

    if(psi == 0) return;

    curH->add(new NonLocalPPotential(*ion,*target,*psi), title);
  }

  void 
  HamiltonianPool::addCorePolPotential(xmlNodePtr cur, ParticleSet* target) {

    string src("i"),title("CorePol");

    OhmmsAttributeSet pAttrib;
    pAttrib.add(title,"name");
    pAttrib.add(src,"source");
    pAttrib.put(cur);

    ParticleSet* ion=ptclPool->getParticleSet(src);
    if(ion == 0) return;

    QMCHamiltonianBase* cpp=(new LocalCorePolPotential(*ion,*target));
    cpp->put(cur); 
    curH->add(cpp, title);
  }

  bool HamiltonianPool::put(std::istream& is) {
    return true;
  }

  bool HamiltonianPool::get(std::ostream& os) const {

    PoolType::const_iterator it(myPool.begin()), it_end(myPool.end());
    while(it != it_end) {
      os << "\n  Hamiltonian " << (*it).first << endl;;
      (*it).second->get(os);
      ++it;
    }
    return true;
  }

  void HamiltonianPool::reset() {
 
  }
}
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
