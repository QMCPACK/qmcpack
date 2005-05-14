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
#ifndef OHMMS_QMC_QMCHAMILTONIANS_H
#define OHMMS_QMC_QMCHAMILTONIANS_H

#include "OhmmsData/OhmmsElementBase.h"
#include <map>

namespace ohmmsqmc {

  class ParticleSet;
  class QMCHamiltonian;
  class ParticleSetPool;
  class WaveFunctionPool;

  /* A collection of ParticleSet
   */
  class HamiltonianPool : public OhmmsElementBase {

  public:

    typedef std::map<std::string,QMCHamiltonian*> PoolType;

    HamiltonianPool(const char* aname = "hamiltonian");

    bool get(std::ostream& os) const;
    bool put(std::istream& is);
    bool put(xmlNodePtr cur);
    void reset();

    inline bool empty() const { return myPool.empty();}

    QMCHamiltonian* getHamiltonian(const std::string& pname) {
      PoolType::iterator hit(myPool.find(pname));
      if(hit == myPool.end()) 
        return 0;
      else 
        (*hit).second;
    }

    inline void setParticleSetPool(ParticleSetPool* pset) { ptclPool=pset;}
    inline void setWaveFunctionPool(WaveFunctionPool* pset) { psiPool=pset;}

    void addCoulombPotential(xmlNodePtr cur, ParticleSet* target);
    void addPseudoPotential(xmlNodePtr cur, ParticleSet* target);

  private:

    QMCHamiltonian* curH;
    ParticleSetPool* ptclPool;
    WaveFunctionPool* psiPool;
    PoolType myPool;
  };
}
#endif
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
