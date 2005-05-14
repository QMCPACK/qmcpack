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
#ifndef OHMMS_QMC_WAVEFUNCTIONPOOL_H
#define OHMMS_QMC_WAVEFUNCTIONPOOL_H

#include "OhmmsData/OhmmsElementBase.h"
#include <map>
#include <string>

namespace ohmmsqmc {

  class TrialWaveFunction;
  class ParticleSetPool;

  /* A collection of TrialWaveFunction s
   */
  class WaveFunctionPool : public OhmmsElementBase {

  public:

    WaveFunctionPool(const char* aname = "wavefunction");

    bool get(std::ostream& os) const;
    bool put(std::istream& is);
    bool put(xmlNodePtr cur);
    void reset();

    inline bool empty() const { return myPool.empty();}

    TrialWaveFunction* getWaveFunction(const std::string& pname) {
      std::map<std::string,TrialWaveFunction*>::iterator pit(myPool.find(pname));
      if(pit == myPool.end()) 
        return 0;
      else 
        (*pit).second;
    }

    inline void setParticleSetPool(ParticleSetPool* pset) { ptclPool=pset;}

  private:

    ParticleSetPool* ptclPool;
    std::map<std::string,TrialWaveFunction*> myPool;
  };
}
#endif
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
