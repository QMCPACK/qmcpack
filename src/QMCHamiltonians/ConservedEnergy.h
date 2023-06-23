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
#include "QMCHamiltonians/OperatorBase.h"
#include "ParticleBase/ParticleAttribOps.h"

namespace qmcplusplus
{
using WP = WalkerProperties::Indexes;

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
struct ConservedEnergy : public OperatorBase
{
  ConservedEnergy() {}
  ~ConservedEnergy() override {}

  void resetTargetParticleSet(ParticleSet& P) override {}

  std::string getClassName() const override { return "ConservedEnergy"; }

  Return_t evaluate(ParticleSet& P) override
  {
    RealType gradsq = Dot(P.G, P.G);
    RealType lap    = Sum(P.L);
#ifdef QMC_COMPLEX
    RealType gradsq_cc = Dot_CC(P.G, P.G);
    value_             = lap + gradsq + gradsq_cc;
#else
    value_ = lap + 2 * gradsq;
#endif
    return 0.0;
  }

  /** Do nothing */
  bool put(xmlNodePtr cur) override { return true; }

  bool get(std::ostream& os) const override
  {
    os << "ConservedEnergy";
    return true;
  }

  std::unique_ptr<OperatorBase> makeClone(ParticleSet& qp, TrialWaveFunction& psi) final
  {
    return std::make_unique<ConservedEnergy>();
  }
};
} // namespace qmcplusplus
#endif
