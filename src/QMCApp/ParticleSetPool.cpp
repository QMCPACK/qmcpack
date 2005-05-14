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
#include "QMCApp/ParticleSetPool.h"
#include "OhmmsData/AttributeSet.h"
#include "ParticleIO/XMLParticleIO.h"

namespace ohmmsqmc {
  
  ParticleSetPool::ParticleSetPool(const char* aname):
    OhmmsElementBase(aname){ }

  bool ParticleSetPool::put(xmlNodePtr cur) {

    string id("e"), role("none");
    OhmmsAttributeSet pAttrib;
    pAttrib.add(id,"id"); pAttrib.add(id,"name"); 
    pAttrib.add(role,"role");

    //backward compatibility
    if(id == "e" && role=="none") role="MC";

    ParticleSet* pTemp = getParticleSet(id);
    if(pTemp == 0) {
      pTemp = new MCWalkerConfiguration;
      //if(role == "MC") 
      //  pTemp = new MCWalkerConfiguration;
      //else 
      //  pTemp = new ParticleSet;
      myPool[id] = pTemp;
      XMLParticleParser pread(*pTemp);
      return pread.put(cur);
    }//possible overwrite

    return true;
  }

  bool ParticleSetPool::put(std::istream& is) {
    return true;
  }

  bool ParticleSetPool::get(std::ostream& os) const {
    return true;
  }

  void ParticleSetPool::reset() {
 
  }
}
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
