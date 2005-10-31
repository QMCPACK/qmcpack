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
/**@file WaveFunctionPool.h
 * @brief Declaration of WaveFunctionPool
 */
#ifndef QMCPLUSPLUS_WAVEFUNCTIONPOOL_H
#define QMCPLUSPLUS_WAVEFUNCTIONPOOL_H

#include "OhmmsData/OhmmsElementBase.h"
#include <map>
#include <string>

namespace qmcplusplus {

  class TrialWaveFunction;
  class ParticleSetPool;
  class OrbitalBuilderBase;

  /** @ingroup qmcapp
   * @brief Manage a collection of TrialWaveFunction objects
   *
   * This object handles \<wavefunction\> elements and
   * functions as a builder class for TrialWaveFunction objects.
   */
  class WaveFunctionPool : public OhmmsElementBase {

  public:

    WaveFunctionPool(const char* aname = "wavefunction");

    bool get(std::ostream& os) const;
    bool put(std::istream& is);
    bool put(xmlNodePtr cur);
    void reset();

    inline bool empty() const { return myPool.empty();}

    TrialWaveFunction* getPrimary() {
      return primaryPsi;
    }

    TrialWaveFunction* getWaveFunction(const std::string& pname) {
      std::map<std::string,TrialWaveFunction*>::iterator pit(myPool.find(pname));
      if(pit == myPool.end()) 
        return 0;
      else 
        return (*pit).second;
    }

    /** assign a pointer of ParticleSetPool
     */
    inline void 
    setParticleSetPool(ParticleSetPool* pset) { ptclPool=pset;}

    /** return a xmlNode containing Jastrow
     * @param id name of the wave function
     *
     * If the wavefunction with id does not exist, return 0
     */
    xmlNodePtr getWaveFunctionNode(const std::string& id);

  private:

    TrialWaveFunction* primaryPsi;
    std::map<std::string,TrialWaveFunction*> myPool;

    /** pointer to ParticleSetPool
     *
     * TrialWaveFunction needs to know which ParticleSet object
     * is used as an input object for the evaluations.
     */
    ParticleSetPool* ptclPool;

    /** save the builder class to clean up at the end
     */
    std::vector<OrbitalBuilderBase*> orbitalBuilders;

    std::map<std::string,xmlNodePtr> m_wfsPtr;
  };
}
#endif
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
