//////////////////////////////////////////////////////////////////
// (c) Copyright 2005- by Jeongnim Kim and Simone Chiesa
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
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
#ifndef QMCPLUSPLUS_NONLOCAL_ECPOTENTIAL_H
#define QMCPLUSPLUS_NONLOCAL_ECPOTENTIAL_H
#include "QMCHamiltonians/NonLocalECPComponent.h"

namespace qmcplusplus {

  /** @ingroup hamiltonian
   * \brief Evaluate the semi local potentials
   */
  struct NonLocalECPotential: public QMCHamiltonianBase {

    int NumIons;
    ///the distance table containing electron-nuclei distances  
    DistanceTableData* d_table;
    ///the set of local-potentials (one for each ion)
    vector<NonLocalECPComponent*> PP;
    ///unique NonLocalECPComponent to remove
    map<int,NonLocalECPComponent*> PPset;
    ///reference to the center ion
    ParticleSet& IonConfig;
    ///target TrialWaveFunction
    TrialWaveFunction& Psi;

    NonLocalECPotential(ParticleSet& ions, ParticleSet& els, TrialWaveFunction& psi);

    ~NonLocalECPotential();

    void resetTargetParticleSet(ParticleSet& P);

    Return_t evaluate(ParticleSet& P);

    Return_t evaluate(ParticleSet& P, vector<NonLocalData>& Txy);

    /** Do nothing */
    bool put(xmlNodePtr cur) { return true; }

    bool get(std::ostream& os) const {
      os << "NonLocalECPotential: " << IonConfig.getName();
      return true;
    }

    void add(int groupID, NonLocalECPComponent* pp);

  };
}
#endif

/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/

