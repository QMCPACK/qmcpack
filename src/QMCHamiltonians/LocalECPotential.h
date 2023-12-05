//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


#ifndef QMCPLUSPLUS_LOCALECPPOTENTIAL_H
#define QMCPLUSPLUS_LOCALECPPOTENTIAL_H
#include "Particle/ParticleSet.h"
#include "QMCHamiltonians/OperatorBase.h"
#include "Numerics/OneDimGridBase.h"
#include "Numerics/OneDimGridFunctor.h"
#include "Numerics/OneDimLinearSpline.h"
#include "Numerics/OneDimCubicSpline.h"
#include "Particle/DistanceTable.h"

namespace qmcplusplus
{
/** @ingroup hamiltonian
 * \brief Evaluate the local potentials (either pseudo or full core) around each ion.
 */

struct LocalECPotential : public OperatorBase
{
  using GridType            = OneDimGridBase<RealType>;
  using RadialPotentialType = OneDimCubicSpline<RealType>;

  ///reference to the ionic configuration
  const ParticleSet& IonConfig;
  ///the number of ioncs
  int NumIons;
  ///distance table index
  int myTableIndex;
  ///temporary energy per particle for pbyp move
  RealType PPtmp;
  ///unique set of local ECP to cleanup
  std::vector<std::unique_ptr<RadialPotentialType>> PPset;
  ///PP[iat] is the local potential for the iat-th particle
  std::vector<RadialPotentialType*> PP;
  ///effective charge per ion
  std::vector<RealType> Zeff;
  ///effective charge per species
  std::vector<RealType> gZeff;
  ///energy per particle
  Vector<RealType> PPart;
#if !defined(REMOVE_TRACEMANAGER)
  ///single particle trace samples
  Array<TraceReal, 1>* Ve_sample;
  Array<TraceReal, 1>* Vi_sample;
#endif
  const ParticleSet& Peln;
  const ParticleSet& Pion;

  LocalECPotential(const ParticleSet& ions, ParticleSet& els);

  std::string getClassName() const override { return "LocalECPotential"; }
  void resetTargetParticleSet(ParticleSet& P) override;

#if !defined(REMOVE_TRACEMANAGER)
  void contributeParticleQuantities() override;
  void checkoutParticleQuantities(TraceManager& tm) override;
  Return_t evaluate_sp(ParticleSet& P); //collect
  void deleteParticleQuantities() override;
#endif

  Return_t evaluate(ParticleSet& P) override;

  Return_t evaluateWithIonDerivs(ParticleSet& P,
                                 ParticleSet& ions,
                                 TrialWaveFunction& psi,
                                 ParticleSet::ParticlePos& hf_terms,
                                 ParticleSet::ParticlePos& pulay_terms) override;


  Return_t evaluate_orig(ParticleSet& P);

  bool put(xmlNodePtr cur) override { return true; }

  bool get(std::ostream& os) const override
  {
    os << "LocalECPotential: " << IonConfig.getName();
    return true;
  }

  std::unique_ptr<OperatorBase> makeClone(ParticleSet& qp, TrialWaveFunction& psi) override;

  /** Add a RadialPotentialType of a species
   * @param groupID index of the ion species
   * @param ppot local pseudopotential
   * @param z effective charge of groupID particle
   */
  void add(int groupID, std::unique_ptr<RadialPotentialType>&& ppot, RealType z);
};
} // namespace qmcplusplus
#endif
