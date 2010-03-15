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
#include "QMCHamiltonians/ForceBase.h"

namespace qmcplusplus {

  /** @ingroup hamiltonian
   * \brief Evaluate the semi local potentials
   */
  struct NonLocalECPotential: public QMCHamiltonianBase,
  public ForceBase 
  {
    int NumIons;
    ///the distance table containing electron-nuclei distances  
    DistanceTableData* d_table;
    ///the set of local-potentials (one for each ion)
    vector<NonLocalECPComponent*> PP;
    ///unique NonLocalECPComponent to remove
    vector<NonLocalECPComponent*> PPset;
    ///reference to the center ion
    ParticleSet& IonConfig;
    ///target TrialWaveFunction
    TrialWaveFunction& Psi;
    ///true if we should compute forces
    bool ComputeForces;
    ParticleSet::ParticlePos_t PulayTerm;

    NonLocalECPotential(ParticleSet& ions, ParticleSet& els, 
			TrialWaveFunction& psi, bool computeForces=false); 

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

    QMCHamiltonianBase* makeClone(ParticleSet& qp, TrialWaveFunction& psi);

    void add(int groupID, NonLocalECPComponent* pp);

    void setRandomGenerator(RandomGenerator_t* rng);

    void addObservables(PropertySetType& plist, BufferType& collectables);

    void setObservables(PropertySetType& plist);

    void setParticlePropertyList(PropertySetType& plist, int offset);

    void registerObservables(vector<observable_helper*>& h5list,
			     hid_t gid) const;
  };
}
#endif

/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/

