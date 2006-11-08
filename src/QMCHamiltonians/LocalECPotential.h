//////////////////////////////////////////////////////////////////
// (c) Copyright 2006- by Jeongnim Kim
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
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
#ifndef QMCPLUSPLUS_LOCALECPPOTENTIAL_H
#define QMCPLUSPLUS_LOCALECPPOTENTIAL_H
#include "Particle/ParticleSet.h"
#include "QMCHamiltonians/QMCHamiltonianBase.h"
#include "Numerics/OneDimGridBase.h"
#include "Numerics/OneDimGridFunctor.h"
#include "Numerics/OneDimLinearSpline.h"

namespace qmcplusplus {

  /** @ingroup hamiltonian
   * \brief Evaluate the local potentials (either pseudo or full core) around each ion.
   */

  struct LocalECPotential: public QMCHamiltonianBase {

    typedef OneDimGridBase<RealType> GridType;
    typedef OneDimLinearSpline<RealType> RadialPotentialType;

    ///the number of ioncs
    int NumIons;
    ///the distance table containing electron-nuclei distances  
    DistanceTableData* d_table;
    ///unique set of local ECP to cleanup
    map<int,RadialPotentialType*> PPset;
    ///PP[iat] is the local potential for the iat-th particle
    vector<RadialPotentialType*> PP;
    ///effective charge per ion
    vector<RealType> Zeff;
    ///reference to the ionic configuration
    const ParticleSet& IonConfig;

    LocalECPotential(ParticleSet& ions, ParticleSet& els);
    
    ~LocalECPotential();

    void resetTargetParticleSet(ParticleSet& P);

    Return_t evaluate(ParticleSet& P);

    inline Return_t evaluate(ParticleSet& P, vector<NonLocalData>& Txy) {
      return evaluate(P);
    }

    bool put(xmlNodePtr cur) { return true;}

    bool get(std::ostream& os) const {
      os << "LocalECPotential: " << IonConfig.getName();
      return true;
    }

    /** Add a RadialPotentialType of a species
     * @param groupID index of the ion species
     * @param ppot local pseudopotential
     * @param z effective charge of groupID particle
     */
    void add(int groupID, RadialPotentialType* ppot, RealType z);
  };
}
#endif

/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/

