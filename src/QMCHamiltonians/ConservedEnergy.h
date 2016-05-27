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
#ifndef QMCPLUSPLUS_CONSERVEDENERGY_H
#define QMCPLUSPLUS_CONSERVEDENERGY_H

#include "Particle/ParticleSet.h"
#include "Particle/WalkerSetRef.h"
#include "QMCHamiltonians/QMCHamiltonianBase.h"
#include "ParticleBase/ParticleAttribOps.h"
#ifdef QMC_CUDA
#include "Particle/MCWalkerConfiguration.h"
#endif

namespace qmcplusplus
{

/** A fake Hamiltonian to check the sampling of the trial function.
 *
 * Integrating the expression
 \f[
 {\bf \nabla} \cdot (\Psi_T({\bf R}) {\bf \nabla}
 \Psi_T({\bf R})) = \Psi_T({\bf R}) \nabla^2
 \Psi_T({\bf R}) + ({\bf \nabla} \Psi_T({\bf R}))^2
 \f]
 leads to (using the divergence theorem)
 \f[
 0 = \int d^3 {\bf R} \: \left[\Psi_T({\bf R}) \nabla^2
 \Psi_T({\bf R}) + ({\bf \nabla} \Psi_T({\bf R}))^2\right]
 \f]
 or, written in terms of the probability distribution
 \f$ |\Psi_T({\bf R})|^2 \f$
 \f[ 0 = \int d^3 {\bf R} \: |\Psi_T({\bf R})|^2
 \left[\frac{\nabla^2 \Psi_T({\bf R})}{\Psi_T({\bf R})} +
 \left(\frac{{\bf \nabla} \Psi_T({\bf R})}{\Psi_T({\bf R})}
 \right)^2\right],
 \f]
 where
 \f[
 \frac{\nabla^2 \Psi_T({\bf R})}{\Psi_T({\bf R})} =
 \nabla^2 \ln \Psi_T({\bf R}) +
 ({\bf \nabla} \ln \Psi_T({\bf R}))^2
 \f]
 \f[
 \frac{{\bf \nabla} \Psi_T({\bf R})}{\Psi_T({\bf R})} =
 {\bf \nabla} \ln \Psi_T({\bf R})
 \f]
 it is possible to check the sampling and the evaluation
 of \f$ \Psi_T, \f$ e.g. the gradient and laplacian.
 The expectation value of this estimator should fluctuate
 around zero.
*/
struct ConservedEnergy: public QMCHamiltonianBase
{

  ConservedEnergy() {}
  ~ConservedEnergy() { }

  void resetTargetParticleSet(ParticleSet& P) { }

  Return_t
  evaluate(ParticleSet& P)
  {
    RealType gradsq = Dot(P.G,P.G);
    RealType lap = Sum(P.L);
    Value = lap+2.0*gradsq;
    return 0.0;
  }

  inline Return_t evaluate(ParticleSet& P, std::vector<NonLocalData>& Txy)
  {
    return evaluate(P);
  }

  /** Do nothing */
  bool put(xmlNodePtr cur)
  {
    return true;
  }

  bool get(std::ostream& os) const
  {
    os << "ConservedEnergy";
    return true;
  }

  QMCHamiltonianBase* makeClone(ParticleSet& qp, TrialWaveFunction& psi)
  {
    return new ConservedEnergy;
  }

#ifdef QMC_CUDA
  ////////////////////////////////
  // Vectorized version for GPU //
  ////////////////////////////////
  // Nothing is done on GPU here, just copy into vector
  void addEnergy(MCWalkerConfiguration &W,
                 std::vector<RealType> &LocalEnergy)
  {
    // Value of LocalEnergy is not used in caller because this is auxiliary H.
    std::vector<Walker_t*> &walkers = W.WalkerList;
    for (int iw=0; iw<walkers.size(); iw++)
    {
      Walker_t &w = *(walkers[iw]);
      double flux = 0.0;
      for (int ptcl=0; ptcl<w.G.size(); ptcl++)
        flux +=2.0 * dot(w.G[ptcl],w.G[ptcl]) + w.L[ptcl];
      w.getPropertyBase()[NUMPROPERTIES+myIndex] = flux;
    }
  }
#endif

};
}
#endif

/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$
 ***************************************************************************/

