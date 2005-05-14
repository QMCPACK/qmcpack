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
#include "QMCApp/HamiltonianPool.h"
#include "QMCApp/WaveFunctionPool.h"
#include "QMCApp/ParticleSetPool.h"
using namespace std;
#include "OhmmsData/AttributeSet.h"
#include "Utilities/OhmmsInfo.h"
#include "QMCHamiltonians/QMCHamiltonian.h"
#include "QMCHamiltonians/CoulombPotential.h"
#include "QMCHamiltonians/IonIonPotential.h"
#include "QMCHamiltonians/LocalPPotential.h"
#include "QMCHamiltonians/NonLocalPPotential.h"
#include "QMCHamiltonians/BareKineticEnergy.h"

namespace ohmmsqmc {

  HamiltonianPool::HamiltonianPool(const char* aname):
    OhmmsElementBase(aname), curH(0), ptclPool(0), psiPool(0){ }

  bool HamiltonianPool::put(xmlNodePtr cur) {

    string id("primary"), htype("generic"), target("e");

    OhmmsAttributeSet hAttrib;
    hAttrib.add(id,"id"); hAttrib.add(id,"name"); 
    hAttrib.add(htype,"type"); 
    hAttrib.add(target,"target");

    hAttrib.put(cur);

    ParticleSet* qp=ptclPool->getParticleSet(target);

    if(qp == 0) {
      ERRORMSG("No target particle " << target << " exists.")
      return false;
    }

    curH = getHamiltonian(id);
    if(curH) {
      WARNMSG("hamiltonian with " << id << " is already created. Ignore the input")
      return true;
    }

    curH = new QMCHamiltonian;
    curH->add(new BareKineticEnergy,"Kinetic");
    myPool[id]=curH;
    cur=cur->children;
    while(cur != NULL) {
      string cname((const char*)cur->name);
      if(cname == "pairpot") {
         const xmlChar* t = xmlGetProp(cur,(const xmlChar*)"type");
         if(t == NULL) continue;
         string pot_type((const char*)t);
         if(pot_type == "coulomb") {
           addCoulombPotential(cur,qp);
        } else if(cname == "pseudo") {
           addPseudoPotential(cur,qp);
        }
      }
      cur = cur->next;
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

    if(source->getTotalNum()>1)  {
      if(source == target) {
        curH->add(new CoulombPotentialAA(*target),title);
      } else {
        curH->add(new CoulombPotentialAB(*source,*target),title);
      }
    }
  }
  
  void HamiltonianPool::addPseudoPotential(xmlNodePtr cur, ParticleSet* target) {

    string src("i"),title("PseudoPot"),wfname("primary");

    OhmmsAttributeSet pAttrib;
    pAttrib.add(title,"name");
    pAttrib.add(src,"source");
    pAttrib.add(wfname,"wavefunction");
    pAttrib.put(cur);

    ParticleSet* ion=ptclPool->getParticleSet(src);
    if(ion == 0) return;

    TrialWaveFunction* psi = psiPool->getWaveFunction(wfname);
    if(psi == 0) return;

    LOGMSG("Creating non-local pseudopotential nuclei = " << src)
    curH->add(new NonLocalPPotential(*ion,*target,*psi), title);
  }

  bool HamiltonianPool::put(std::istream& is) {
    return true;
  }

  bool HamiltonianPool::get(std::ostream& os) const {
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
