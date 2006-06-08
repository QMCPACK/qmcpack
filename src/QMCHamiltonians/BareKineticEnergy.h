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
#ifndef QMCPLUSPLUS_BAREKINETICENERGY_H
#define QMCPLUSPLUS_BAREKINETICENERGY_H
#include "Particle/ParticleSet.h"
#include "Particle/WalkerSetRef.h"
#include "QMCHamiltonians/QMCHamiltonianBase.h"
#include "ParticleBase/ParticleAttribOps.h"

namespace qmcplusplus {

  /** @ingroup hamiltonian
   @brief Evaluate the kinetic energy with a single mass
   
   The unit of the mass is AU, i.e., the electron mass \f$ m_e = 1 \f$.
   *
   * To evaluate the Bare Kinetic part of the local energy 
   \f$E_L({\bf R}) = \Psi^{-1}({\bf R})\hat{H}\Psi({\bf R}),\f$
   it is useful to use the following trick
   \f{eqnarray*}
   \nabla^2\Psi({\bf R}) &=& \nabla^2(\exp(\ln \Psi({\bf R})))\\
   &=&\nabla\cdot(\nabla\exp(\ln \Psi({\bf R}))) \\
   &=&\nabla\cdot(\nabla\ln \Psi({\bf R}))\exp(\ln \Psi({\bf R}))\\
   -\frac{1}{2}\frac{\nabla^2\Psi({\bf R})}{\Psi({\bf R})} &=& 
   -\frac{1}{2}\nabla^2\ln \Psi({\bf R})
   -\frac{1}{2}(\nabla\ln \Psi({\bf R}))^2
   \f}
  */
  struct BareKineticEnergy: public QMCHamiltonianBase {

    ///mass of the particle
    RealType M;
    ///\f$ 1/(2 m^*) \f$
    RealType OneOver2M;

    /** constructor
     *
     * Kinetic operators need to be re-evaluated during optimization.
     */
    BareKineticEnergy(RealType m=1.0): M(m),OneOver2M(0.5/m) { 
      UpdateMode.set(OPTIMIZABLE,1);
    }
    ///destructor
    ~BareKineticEnergy() { }

    void resetTargetParticleSet(ParticleSet& P) { }

    inline Return_t 
    evaluate(ParticleSet& P) {
      Value = Dot(P.G,P.G) + Sum(P.L); 
      //Value = 0.0;
      //Value=Dot(P.G,P.G)+Sum(P.L);
      return Value*=-OneOver2M;
    }
    
    /** implements the virtual function.
     * 
     * Nothing is done but should check the mass
     */
    bool put(xmlNodePtr) {
      return true;
    }

    bool get(std::ostream& os) const {
      os << "Kinetic energy";
      return true;
    }

    QMCHamiltonianBase* clone() {
      return new BareKineticEnergy(M);
    }
    //Not used anymore
    //void evaluate(WalkerSetRef& W, ValueVectorType& LE) {
    //  for(int iw=0; iw< W.walkers(); iw++) {
    //    RealType ke = 0.0;
    //    for(int iat=0; iat< W.particles(); iat++) {
    //      ke += dot(W.G(iw,iat),W.G(iw,iat)) + W.L(iw,iat);
    //    }
    //    LE[iw] -= M*ke;
    //  }
    //}
  };
}
#endif

/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/

