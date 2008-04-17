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
#include "Message/OpenMP.h"

namespace qmcplusplus {

  HamiltonianPool::HamiltonianPool(const char* aname):
    OhmmsElementBase(aname), curH(0), ptclPool(0), psiPool(0),
  curDoc(0){ }

  bool HamiltonianPool::put(xmlNodePtr cur) {

    string id("h0"), target("e"),role("extra");
    OhmmsAttributeSet hAttrib;
    hAttrib.add(id,"id"); hAttrib.add(id,"name"); 
    hAttrib.add(role,"role");
    hAttrib.add(target,"target");
    hAttrib.put(cur);

    ParticleSet* qp=ptclPool->getParticleSet(target);
    if(qp == 0) {
      ERRORMSG("No target particle " << target << " exists.")
      return false;
    }

    bool set2Primary=false;
    if(myPool.empty() || role == "primary" ) { set2Primary=true;}

    HamiltonianFactory *curH=0;
    PoolType::iterator hit(myPool.find(id));
    if(hit == myPool.end()) {
      curH= new HamiltonianFactory(qp, ptclPool->getPool(), psiPool->getPool());
      curH->setName(id);
      myPool[id]=curH;
    }
    else 
      curH=(*hit).second;

    bool success= curH->put(cur);
    
    if(set2Primary) {
      primaryH=curH->targetH;
    }

    return success;
  }

  void HamiltonianPool::clone(const MCWalkerConfiguration& qp, const TrialWaveFunction& psi, const QMCHamiltonian& h,
        vector<MCWalkerConfiguration*>& plist, vector<TrialWaveFunction*>& olist, vector<QMCHamiltonian*>& hlist) {

    //clone ParticleSet and TrialWaveFunction
    WaveFunctionFactory* psiFac=psiPool->getWaveFunctionFactory(psi.getName());
    int np=omp_get_max_threads();
    psiFac->setCloneSize(np);

    //capture cloned WaveFunctionFactory*
    vector<WaveFunctionFactory*> otemp;
    otemp.resize(np,0);

    //allocate the data on each thread
//#pragma omp parallel for
//    for(int ip=0; ip<np; ip++) {
#pragma omp parallel 
    {
      int ip=omp_get_thread_num();
#pragma omp critical 
      {
        if(ip) {
          char pname[16],oname[16];
          sprintf(pname,"%s.c%i",qp.getName().c_str(),ip);
          plist[ip]=new MCWalkerConfiguration(qp);
          plist[ip]->setName(pname);

          sprintf(oname,"%s.c%i",psi.getName().c_str(),ip);
          otemp[ip]= psiFac->clone(plist[ip],ip,oname);
        }
      }
    }
    
    //add the Clones to the pools
    for(int ip=1; ip<np; ip++) {
      ptclPool->addParticleSet(plist[ip]);
      psiPool->addFactory(otemp[ip]);
      olist[ip]=otemp[ip]->targetPsi;
      //if(ip%2==1) olist[ip]->reverse();
    }

    //find the HamiltonianFactory* to be cloned
    HamiltonianFactory* hFac=0;
    PoolType::iterator hit(myPool.find(h.getName()));
    if(hit == myPool.end()) {
      hFac=(*(myPool.begin())).second;
    } else {
      hFac=(*hit).second;
    }
    hFac->setCloneSize(np);

    vector<HamiltonianFactory*> htemp;
    htemp.resize(np,0);

//#pragma omp parallel for
//    for(int ip=0; ip<np; ip++) {
#pragma omp parallel 
    {
      int ip=omp_get_thread_num();
#pragma omp critical 
      {
        if(ip) {
          char hname[16];
          sprintf(hname,"%s.c%i",h.getName().c_str(),ip);
          htemp[ip]= hFac->clone(plist[ip],olist[ip],ip,hname);
        }
      }
    }

    for(int ip=1; ip<np; ip++) {
      myPool[htemp[ip]->getName()]=htemp[ip];
      hlist[ip]=htemp[ip]->targetH;
    }
  }

  bool HamiltonianPool::put(std::istream& is) {
    return true;
  }

  bool HamiltonianPool::get(std::ostream& os) const {

    PoolType::const_iterator it(myPool.begin()), it_end(myPool.end());
    while(it != it_end) {
      os << "\n  Hamiltonian " << (*it).first << endl;;
      (*it).second->targetH->get(os);
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
