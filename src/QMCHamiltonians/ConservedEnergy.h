//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//                    Mark Dewing, markdewing@gmail.com, University of Illinois at Urbana-Champaign
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    
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

 For complex wavefunctions, the
 \f[ \nabla \Psi_T \cdot \nabla \Psi_T^* \f] term from the divergence
 theorem should use the complex conjugate.
 The \f[ \nabla \Psi_T \cdot \nabla \Psi_T \f] term from expressing
 \f[\Psi\f] in terms of \f[\ln \Psi\] should use normal complex
 multiplication.
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
#ifdef QMC_COMPLEX
    RealType gradsq_cc = Dot_CC(P.G,P.G);
    Value = lap + gradsq + gradsq_cc;
#else
    Value = lap + 2*gradsq;
#endif
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
      RealType flux;
      RealType gradsq = Dot(w.G,w.G);
      RealType lap = Sum(w.L);
#ifdef QMC_COMPLEX
      RealType gradsq_cc = Dot_CC(w.G,w.G);
      flux = lap + gradsq + gradsq_cc;
#else
      flux = lap + 2*gradsq;
#endif
      w.getPropertyBase()[NUMPROPERTIES+myIndex] = flux;
    }
  }
#endif

};
}
#endif


