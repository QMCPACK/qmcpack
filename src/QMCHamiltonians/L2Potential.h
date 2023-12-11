//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//////////////////////////////////////////////////////////////////////////////////////


#ifndef QMCPLUSPLUS_L2POTENTIAL_H
#define QMCPLUSPLUS_L2POTENTIAL_H
#include "Particle/ParticleSet.h"
#include "QMCHamiltonians/OperatorBase.h"
#include "Numerics/OneDimGridBase.h"
#include "Numerics/OneDimGridFunctor.h"
#include "Numerics/OneDimLinearSpline.h"
#include "Numerics/OneDimCubicSpline.h"
#include "QMCWaveFunctions/TrialWaveFunction.h"

namespace qmcplusplus
{
struct L2RadialPotential : public QMCTraits
{
  using GridType            = OneDimGridBase<RealType>;
  using RadialPotentialType = OneDimCubicSpline<RealType>;

  std::unique_ptr<RadialPotentialType> vL2;
  RealType rcut;

  RealType evaluate(RealType r) { return vL2->splint(r); }

  RealType evaluate_guard(RealType r)
  {
    if (r < rcut)
      return vL2->splint(r);
    else
      return 0.0;
  }

  L2RadialPotential* makeClone()
  {
    auto c = new L2RadialPotential();
    c->vL2.reset(vL2->makeClone());
    c->rcut = rcut;
    return c;
  }
};


/** @ingroup hamiltonian
 * \brief Evaluate the L2 potentials around each ion.
 */

struct L2Potential : public OperatorBase
{
  ///reference to the ionic configuration
  const ParticleSet& IonConfig;
  ///the number of ions
  int NumIons;
  ///distance table index
  int myTableIndex;
  ///unique set of L2 PP to cleanup
  std::vector<std::unique_ptr<L2RadialPotential>> PPset;
  ///PP[iat] is the L2 potential for the iat-th particle
  std::vector<L2RadialPotential*> PP;
  ///Associated trial wavefunction
  TrialWaveFunction* psi_ref;

  L2Potential(const ParticleSet& ions, ParticleSet& els, TrialWaveFunction& psi);

  bool dependsOnWaveFunction() const override { return true; }
  std::string getClassName() const override { return "L2Potential"; }

  void resetTargetParticleSet(ParticleSet& P) override;

  Return_t evaluate(ParticleSet& P) override;

  void evaluateDK(ParticleSet& P, int iel, TensorType& D, PosType& K);
  void evaluateD(ParticleSet& P, int iel, TensorType& D);

  bool put(xmlNodePtr cur) override { return true; }

  bool get(std::ostream& os) const override
  {
    os << "L2Potential: " << IonConfig.getName();
    return true;
  }

  std::unique_ptr<OperatorBase> makeClone(ParticleSet& qp, TrialWaveFunction& psi) final;

  /** Add a RadialPotentialType of a species
   * @param groupID index of the ion species
   * @param ppot L2 pseudopotential
   */
  void add(int groupID, std::unique_ptr<L2RadialPotential>&& ppot);
};
} // namespace qmcplusplus
#endif
