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
/**@file ParticleSetPool.cpp
 * @brief Implements ParticleSetPool operators.
 */
#include "QMCApp/ParticleSetPool.h"
#include "Particle/DistanceTable.h"
#include "OhmmsData/AttributeSet.h"
#include "ParticleIO/XMLParticleIO.h"
#include "Utilities/OhmmsInfo.h"

namespace qmcplusplus {
  
  ParticleSetPool::ParticleSetPool(const char* aname):
    OhmmsElementBase(aname){ }

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
      LOGMSG("Creating " << id << " particleset")
      pTemp = new MCWalkerConfiguration;
      //if(role == "MC") 
      //  pTemp = new MCWalkerConfiguration;
      //else 
      //  pTemp = new ParticleSet;
      if(sc) {
        LOGMSG("Initializing the lattice of " << id << " by the global supercell")
        pTemp->Lattice.copy(*sc);
      }
      myPool[id] = pTemp;
      XMLParticleParser pread(*pTemp);
      bool success = pread.put(cur);
      pTemp->setName(id);
      return success;
    } else {
      WARNMSG("particleset " << id << " is already created. Ignore this")
    }

    return true;
  }

  bool ParticleSetPool::put(std::istream& is) {
    return true;
  }

  bool ParticleSetPool::get(std::ostream& os) const {
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
}
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
