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
#ifndef OHMMS_QMC_PARTICLESETPOOL_H
#define OHMMS_QMC_PARTICLESETPOOL_H

#include "OhmmsData/OhmmsElementBase.h"
#include "Particle/MCWalkerConfiguration.h"

namespace ohmmsqmc {

  /* A collection of ParticleSet
   */
  class ParticleSetPool : public OhmmsElementBase {

  public:

    typedef map<string,ParticleSet*> PoolType;

    ParticleSetPool(const char* aname = "particleset");

    bool get(std::ostream& os) const;
    bool put(std::istream& is);
    bool put(xmlNodePtr cur);
    void reset();
    inline bool empty() const { return myPool.empty();}

    ParticleSet* getParticleSet(const string& pname) {
      map<string,ParticleSet*>::iterator pit(myPool.find(pname));
      if(pit == myPool.end()) 
        return 0;
      else 
        (*pit).second;
    }

    MCWalkerConfiguration* getWalkerSet(const string& pname) {
      map<string,ParticleSet*>::iterator pit(myPool.find(pname));
      if(pit == myPool.end()) 
        return 0;
      else 
        return dynamic_cast<MCWalkerConfiguration*>((*pit).second);
    }

    inline PoolType& getPool() { return myPool;}

  private:

    map<string,ParticleSet*> myPool;
  };
}
#endif
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
