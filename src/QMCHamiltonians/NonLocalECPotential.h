//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


#ifndef QMCPLUSPLUS_NONLOCAL_ECPOTENTIAL_H
#define QMCPLUSPLUS_NONLOCAL_ECPOTENTIAL_H
#include "QMCHamiltonians/NonLocalTOperator.h"
#include "QMCHamiltonians/ForceBase.h"
#include "Particle/NeighborLists.h"

namespace qmcplusplus
{
class NonLocalECPComponent;
template<typename T>
struct NLPPJob;

/** @ingroup hamiltonian
 * \brief Evaluate the semi local potentials
 */
class NonLocalECPotential : public OperatorBase, public ForceBase
{
public:
  NonLocalECPotential(ParticleSet& ions, ParticleSet& els, TrialWaveFunction& psi, bool computeForces, bool enable_DLA);

  ~NonLocalECPotential();

  void resetTargetParticleSet(ParticleSet& P) override;

#if !defined(REMOVE_TRACEMANAGER)
  virtual void contribute_particle_quantities() override;
  virtual void checkout_particle_quantities(TraceManager& tm) override;
  virtual void delete_particle_quantities() override;
#endif

  Return_t evaluate(ParticleSet& P) override;

  void mw_evaluate(const RefVector<OperatorBase>& O_list, const RefVector<ParticleSet>& P_list) override;

  Return_t evaluateWithToperator(ParticleSet& P) override;

  void mw_evaluateWithToperator(const RefVector<OperatorBase>& O_list, const RefVector<ParticleSet>& P_list) override;

  Return_t evaluateWithIonDerivs(ParticleSet& P,
                                 ParticleSet& ions,
                                 TrialWaveFunction& psi,
                                 ParticleSet::ParticlePos_t& hf_terms,
                                 ParticleSet::ParticlePos_t& pulay_terms) override;

  /** set non local moves options
   * @param cur the xml input
   */
  void setNonLocalMoves(xmlNodePtr cur) { UseTMove = nonLocalOps.put(cur); }

  void setNonLocalMoves(const std::string& non_local_move_option,
                        const double tau,
                        const double alpha,
                        const double gamma)
  {
    UseTMove = nonLocalOps.thingsThatShouldBeInMyConstructor(non_local_move_option, tau, alpha, gamma);
  }
  /** make non local moves with particle-by-particle moves
   * @param P particle set
   * @return the number of accepted moves
   */
  int makeNonLocalMovesPbyP(ParticleSet& P);

  Return_t evaluateValueAndDerivatives(ParticleSet& P,
                                       const opt_variables_type& optvars,
                                       const std::vector<ValueType>& dlogpsi,
                                       std::vector<ValueType>& dhpsioverpsi) override;

  /** Do nothing */
  bool put(xmlNodePtr cur) override { return true; }

  bool get(std::ostream& os) const override
  {
    os << "NonLocalECPotential: " << IonConfig.getName();
    return true;
  }

  OperatorBase* makeClone(ParticleSet& qp, TrialWaveFunction& psi) override;

  void addComponent(int groupID, NonLocalECPComponent* pp);

  /** set the internal RNG pointer as the given pointer
   * @param rng input RNG pointer
   */
  void setRandomGenerator(RandomGenerator_t* rng) override { myRNG = rng; }

  void addObservables(PropertySetType& plist, BufferType& collectables) override;

  void setObservables(PropertySetType& plist) override;

  void setParticlePropertyList(PropertySetType& plist, int offset) override;

  void registerObservables(std::vector<observable_helper*>& h5list, hid_t gid) const override;

protected:
  ///random number generator
  RandomGenerator_t* myRNG;
  ///the set of local-potentials (one for each ion)
  std::vector<NonLocalECPComponent*> PP;
  ///unique NonLocalECPComponent to remove
  std::vector<NonLocalECPComponent*> PPset;
  ///reference to the center ion
  ParticleSet& IonConfig;
  ///target TrialWaveFunction
  TrialWaveFunction& Psi;

private:
  ///number of ions
  int NumIons;
  ///index of distance table for the ion-el pair
  int myTableIndex;
  ///reference to the electrons
  ParticleSet& Peln;
  ///neighborlist of electrons
  NeighborLists ElecNeighborIons;
  ///neighborlist of ions
  NeighborLists IonNeighborElecs;
  ///use T-moves
  int UseTMove;
  ///ture if an electron is affected by other electrons moved by T-moves
  std::vector<bool> elecTMAffected;
  ///non local operator
  NonLocalTOperator nonLocalOps;
  ///true if we should compute forces
  bool ComputeForces;
  ///true, determinant localization approximation(DLA) is enabled
  bool use_DLA;
  ///Pulay force vector
  ParticleSet::ParticlePos_t PulayTerm;
#if !defined(REMOVE_TRACEMANAGER)
  ///single particle trace samples
  Array<TraceReal, 1>* Ve_sample;
  Array<TraceReal, 1>* Vi_sample;
#endif
  ///NLPP job list of ion-electron pairs by spin group
  std::vector<std::vector<NLPPJob<RealType>>> nlpp_jobs;

  /** the actual implementation, used by evaluate and evaluateWithToperator
   * @param P particle set
   * @param Tmove whether Txy for Tmove is updated
   */
  void evaluateImpl(ParticleSet& P, bool Tmove);

  /** the actual implementation for batched walkers, used by mw_evaluate and mw_evaluateWithToperator
   * @param O_list the list of NonLocalECPotential in a walker batch
   * @param P_list the list of ParticleSet in a walker batch
   * @param Tmove whether Txy for Tmove is updated
   */
  void mw_evaluateImpl(const RefVector<OperatorBase>& O_list, const RefVector<ParticleSet>& P_list, bool Tmove);

  /** compute the T move transition probability for a given electron
   * member variable nonLocalOps.Txy is updated
   * @param P particle set
   * @param ref_elec reference electron id
   */
  void computeOneElectronTxy(ParticleSet& P, const int ref_elec);

  /** mark all the electrons affected by Tmoves and update ElecNeighborIons and IonNeighborElecs
   * @param myTable electron ion distance table
   * @param iel reference electron
   * Note this funtion should be called before acceptMove for a Tmove
   */
  void markAffectedElecs(const DistanceTableData& myTable, int iel);
};
} // namespace qmcplusplus
#endif
