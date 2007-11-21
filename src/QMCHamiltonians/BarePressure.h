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
#ifndef QMCPLUSPLUS_BAREPRESSURE_H
#define QMCPLUSPLUS_BAREPRESSURE_H
#include "Particle/ParticleSet.h"
#include "Particle/WalkerSetRef.h"
#include "QMCHamiltonians/QMCHamiltonianBase.h"
#include "ParticleBase/ParticleAttribOps.h"

namespace qmcplusplus {

  /** @ingroup hamiltonian
   @brief Evaluate the Bare Pressure.
   P=/frac{2T+V}{d* /Omega}
   where d is the dimension of space and /Omega is the volume.
  **/

  struct BarePressure: public QMCHamiltonianBase {
    double pNorm;

    /** constructor
     *
     * Pressure operators need to be re-evaluated during optimization.
     */
    BarePressure(ParticleSet& P) { 
      UpdateMode.set(OPTIMIZABLE,1);
      //for 3-D. For 2-D this should be 1.0/(2.0*P.getTotalNum())
      pNorm = 1.0/(3.0*P.Lattice.Volume);

    }
    ///destructor
    ~BarePressure() { }

    void resetTargetParticleSet(ParticleSet& P) { }

    inline Return_t 
    evaluate(ParticleSet& P) {
      Value=2.0*P.PropertyList[LOCALENERGY]-P.PropertyList[LOCALPOTENTIAL];
      //for 3-D
      Value*=pNorm;
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

    bool get(std::ostream& os) const {
      os << "BarePressure";
      return true;
    }

    QMCHamiltonianBase* clone(ParticleSet& qp, TrialWaveFunction& psi)
    {
      return new BarePressure( qp);
    }

  };
}
#endif

/***************************************************************************
 * $RCSfile$   $Author: jnkim $
 * $Revision: 1581 $   $Date: 2007-01-04 10:02:14 -0600 (Thu, 04 Jan 2007) $
 * $Id: BareKineticEnergy.h 1581 2007-01-04 16:02:14Z jnkim $ 
 ***************************************************************************/

