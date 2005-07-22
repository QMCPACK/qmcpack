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
/**@file HamiltonianPool.h
 * @brief Declaration of HamiltonianPool
 */
#ifndef OHMMS_QMC_QMCHAMILTONIANS_H
#define OHMMS_QMC_QMCHAMILTONIANS_H

#include "OhmmsData/OhmmsElementBase.h"
#include <map>

namespace ohmmsqmc {

  class ParticleSet;
  class QMCHamiltonian;
  class ParticleSetPool;
  class WaveFunctionPool;

  /* A collection of QMCHamiltonian objects
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

    /** return the pointer to the primary QMCHamiltonian
     *
     * The first QMCHamiltonian is assigned to the primaryH.
     * The last QMCHamiltonian with role="primary" will be the primaryH.
     */
    inline QMCHamiltonian* getPrimary() {
      return primaryH;
    }

    /** return the pointer to a QMCHamiltonian with the name 
     * @param pname name of the QMCHamiltonian
     */
    inline QMCHamiltonian* getHamiltonian(const std::string& pname) {
      PoolType::iterator hit(myPool.find(pname));
      if(hit == myPool.end()) 
        return 0;
      else 
        return (*hit).second;
    }

    /** assign a pointer to a ParticleSetPool
     */
    inline void setParticleSetPool(ParticleSetPool* pset) { ptclPool=pset;}

    /** assign a pointer to a WaveFunctionPool
     */
    inline void setWaveFunctionPool(WaveFunctionPool* pset) { psiPool=pset;}

    void addCoulombPotential(xmlNodePtr cur, ParticleSet* target);
    void addPseudoPotential(xmlNodePtr cur, ParticleSet* target);
    void addCorePolPotential(xmlNodePtr cur, ParticleSet* target);

  private:

    /** pointer to the primary QMCHamiltonian
     */
    QMCHamiltonian* primaryH;

    /** pointer to a current QMCHamiltonian to be built.
     */
    QMCHamiltonian* curH;

    /** pointer to ParticleSetPool
     *
     * QMCHamiltonian needs to know which ParticleSet object
     * is used as an input object for the evaluations. 
     * Any number of ParticleSet can be used to describe
     * a QMCHamiltonian.
     */
    ParticleSetPool* ptclPool;

    /** pointer to WaveFunctionPool
     *
     * For those QMCHamiltonianBase that depends on TrialWaveFunction,
     * e.g., NonLocalPPotential.
     */
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
