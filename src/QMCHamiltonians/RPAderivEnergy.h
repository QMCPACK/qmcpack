//////////////////////////////////////////////////////////////////
// (c) Copyright 2003  by Jeongnim Kim
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
#ifndef QMCPLUSPLUS_RPADERIVENERGY_H
#define QMCPLUSPLUS_RPADERIVENERGY_H
#include "Particle/ParticleSet.h"
#include "Particle/WalkerSetRef.h"
#include "QMCHamiltonians/QMCHamiltonianBase.h"
#include "ParticleBase/ParticleAttribOps.h"

namespace qmcplusplus {

  /** @ingroup hamiltonian
   @brief Evaluate the RPA deriv.
  **/

  struct RPADerivEnergy: public QMCHamiltonianBase {
    RPAPressureCorrection* bpcorr;
    QMCHamiltonian* Ham;

    /** constructor
     *
     * Pressure operators need to be re-evaluated during optimization.
     */
    RPADerivEnergy(ParticleSet& P) { 
      UpdateMode.set(OPTIMIZABLE,1);
    }
    ///destructor
    ~RPADerivEnergy() { }

    void resetTargetParticleSet(ParticleSet& P) { }

    inline Return_t 
    evaluate(ParticleSet& P) {
//       Value = (bpcorr->tValue)* (*Ham)[0];
      Value = P.PropertyList[LOCALENERGY];
      Value *= (bpcorr->tValue);
      return 0.0;
    }

    inline Return_t 
    evaluate(ParticleSet& P, vector<NonLocalData>& Txy) {
      return evaluate(P);
    }
    
    /** implements the virtual function.
     * 
     * Nothing is done but should check the mass
     */
    bool put(xmlNodePtr) {
      return true;
    }
    
    bool put(xmlNodePtr, RPAPressureCorrection* BP, QMCHamiltonian* H ) {
      bpcorr = BP;
      Ham = H;
      return true;
      
    }

    bool get(std::ostream& os) const {
      os << "E*dpsi";
      return true;
    }

    QMCHamiltonianBase* clone(ParticleSet& qp, TrialWaveFunction& psi)
    {
      return new RPADerivEnergy( qp);
    }

  };
}
#endif

/***************************************************************************
 * $RCSfile$   $Author: jnkim $
 * $Revision: 1581 $   $Date: 2007-01-04 10:02:14 -0600 (Thu, 04 Jan 2007) $
 * $Id: BareKineticEnergy.h 1581 2007-01-04 16:02:14Z jnkim $ 
 ***************************************************************************/

