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
/**@file ParticleSetPool.h
 * @brief Declaration of ParticleSetPool
 */
#ifndef QMCPLUSPLUS_PARTICLESETPOOL_H
#define QMCPLUSPLUS_PARTICLESETPOOL_H

#include "OhmmsData/OhmmsElementBase.h"
#include "Particle/MCWalkerConfiguration.h"

namespace qmcplusplus {

  /** @ingroup qmcapp
   * @brief Manage a collection of ParticleSet objects
   *
   * This object handles \<particleset\> elements and
   * functions as a builder class for ParticleSet objects.
   */
  class ParticleSetPool : public OhmmsElementBase {

  public:

    typedef map<string,ParticleSet*> PoolType;

    /** constructor
     * @param aname xml tag
     */
    ParticleSetPool(const char* aname = "particleset");

    //implements virtual functions of OhmmsElementBase
    bool get(std::ostream& os) const;
    bool put(std::istream& is);
    bool put(xmlNodePtr cur);
    void reset();

    ///return true, if the pool is empty
    inline bool empty() const { return myPool.empty();}

    /** get a named ParticleSet
     * @param pname name of the ParticleSet
     * @return a ParticleSet object with pname
     *
     * When the named ParticleSet is not in this object, return 0.
     */
    ParticleSet* getParticleSet(const string& pname) {
      map<string,ParticleSet*>::iterator pit(myPool.find(pname));
      if(pit == myPool.end())  {
        return 0;
      }
      else  {
        return (*pit).second;
      }
    }

    /** get a named MCWalkerConfiguration
     * @param pname name of the MCWalkerConfiguration
     * @return a MCWalkerConfiguration object with pname
     *
     * When the named MCWalkerConfiguration is not in this object, return 0.
     */
    MCWalkerConfiguration* getWalkerSet(const string& pname) {
      map<string,ParticleSet*>::iterator pit(myPool.find(pname));
      if(pit == myPool.end()) 
        return 0;
      else 
        return dynamic_cast<MCWalkerConfiguration*>((*pit).second);
    }

    /** get the Pool object
     */
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
