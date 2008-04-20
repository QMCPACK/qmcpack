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
//
// Supported by 
//   National Center for Supercomputing Applications, UIUC
//   Materials Computation Center, UIUC
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
/**@file ParticleSetPool.cpp
 * @brief Implements ParticleSetPool operators.
 */
#include "QMCApp/ParticleSetPool.h"
#include "Particle/DistanceTable.h"
#include "OhmmsData/AttributeSet.h"
#include "ParticleIO/XMLParticleIO.h"
#include "Utilities/OhmmsInfo.h"

namespace qmcplusplus {
  
  ParticleSetPool::ParticleSetPool(Communicate* c, const char* aname): 
    MPIObjectBase(c), OhmmsElementBase(aname)
    { }

  void  ParticleSetPool::addParticleSet(ParticleSet* p) {
    PoolType::iterator pit(myPool.find(p->getName()));
    if(pit == myPool.end()) 
    {
      LOGMSG("  Adding " << p->getName() << " ParticleSet to the pool")
      myPool[p->getName()]=p;
    } else {
      WARNMSG("  " << p->getName() << " exists. Ignore addition")
    }
  }

  /** process an xml element
   * @param cur current xmlNodePtr
   * @return true, if successful.
   *
   * Creating MCWalkerConfiguration for all the ParticleSet
   * objects. 
   */
  bool ParticleSetPool::put(xmlNodePtr cur) {

    const ParticleSet::ParticleLayout_t* sc=DistanceTable::getSimulationCell();

    string id("e"), role("none");
    OhmmsAttributeSet pAttrib;
    pAttrib.add(id,"id"); pAttrib.add(id,"name"); 
    pAttrib.add(role,"role");
    pAttrib.put(cur);

    //backward compatibility
    if(id == "e" && role=="none") role="MC";

    ParticleSet* pTemp = getParticleSet(id);
    if(pTemp == 0) {
      app_log() << "  Creating " << id << " particleset" << endl;
      pTemp = new MCWalkerConfiguration;
      //if(role == "MC") 
      //  pTemp = new MCWalkerConfiguration;
      //else 
      //  pTemp = new ParticleSet;
      if(sc) {
        app_log() << "  Initializing the lattice of " << id << " by the global supercell" << endl;
        pTemp->Lattice.copy(*sc);
      }
      myPool[id] = pTemp;
      XMLParticleParser pread(*pTemp);
      bool success = pread.put(cur);
      pTemp->setName(id);
      return success;
    } else {
      app_warning() << "particleset " << id << " is already created. Ignore this" << endl;
    }

    return true;
  }

  bool ParticleSetPool::put(std::istream& is) {
    return true;
  }

  bool ParticleSetPool::get(std::ostream& os) const {
    PoolType::const_iterator it(myPool.begin()), it_end(myPool.end());
    while(it != it_end) {
      os << endl;
      (*it).second->get(os);
      ++it;
    }
    return true;
  }

  /** reset is used to initialize and evaluate the distance tables 
   */
  void ParticleSetPool::reset() {
    PoolType::iterator it(myPool.begin()), it_end(myPool.end());
    while(it != it_end) {
      ParticleSet* pt((*it).second);
      pt->update();
      ++it;
    }
    //DistanceTable::create(1);
    //PoolType::iterator it(myPool.begin()), it_end(myPool.end());
    //while(it != it_end) {
    //  DistanceTable::update(*((*it).second));
    //  ++it;
    //}
  }

// Experimental implementation of cloning ParticleSet*
// All taken care by HamiltonianPool
//  std::vector<ParticleSet*>
//  ParticleSetPool::clone(const string& pname, int np) {
//    ParticleSet* pRef=getParticleSet(pname);
//    vector<ParticleSet*> newPtclSets;
//    if(pRef == 0) return newPtclSets;
//
//    newPtclSets.resize(np,0);
//    newPtclSets[0]=pRef;
//    char pnameIP[128];
//
//    for(int ip=1; ip<np; ip++) {
//      sprintf(pnameIP,"%s.c%i",pname.c_str(),ip);
//      map<string,ParticleSet*>::iterator pit(myPool.find(pname));
//      if(pit == myPool.end())  {
//        myPool[pnameIP]=0;//add to the pool
//      } else {
//        newPtclSets[ip]=(*pit).second;
//      }
//    }
//
//#pragma omp parallel for
//    for(int ip=0; ip<np; ip++) {
//      if(newPtclSets[ip] == 0) {
//        newPtclSets[ip]=new ParticleSet(*pRef);
//      }
//    }
//
//    //add the cloned particle sets to myPool
//    for(int ip=1; ip<np; ip++) {
//      sprintf(pnameIP,"%s%i",pname.c_str(),ip);
//      myPool[pnameIP]=newPtclSets[ip];
//    }
//
//    return newPtclSets;
//  }
}
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
