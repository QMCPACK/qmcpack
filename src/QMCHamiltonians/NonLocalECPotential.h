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
/** @ingroup hamiltonian
 * \brief Evaluate the semi local potentials
 */
class NonLocalECPotential : public OperatorBase, public ForceBase
{
public:
  NonLocalECPotential(ParticleSet& ions,
                      ParticleSet& els,
                      TrialWaveFunction& psi,
                      bool computeForces = false,
                      bool useVP         = false);

  ~NonLocalECPotential();

  void resetTargetParticleSet(ParticleSet& P);

#if !defined(REMOVE_TRACEMANAGER)
  virtual void contribute_particle_quantities();
  virtual void checkout_particle_quantities(TraceManager& tm);
  virtual void delete_particle_quantities();
#endif

  Return_t evaluate(ParticleSet& P);

  Return_t evaluateWithIonDerivs(ParticleSet& P,
                                 ParticleSet& ions,
                                 TrialWaveFunction& psi,
                                 ParticleSet::ParticlePos_t& hf_terms,
                                 ParticleSet::ParticlePos_t& pulay_terms);

  Return_t evaluateWithToperator(ParticleSet& P);

  /** set non local moves options
   * @param cur the xml input
   */
  void setNonLocalMoves(xmlNodePtr cur) { UseTMove = nonLocalOps.put(cur); }

  void setNonLocalMoves(const std::string& non_local_move_option,
                                        const double tau,
                                        const double alpha,
                                        const double gamma)
  {
    UseTMove = nonLocalOps.thingsThatShouldBeInMyConstructor(non_local_move_option,
                                        tau,
                                        alpha,
                                        gamma);
  }
  /** make non local moves with particle-by-particle moves
   * @param P particle set
   * @return the number of accepted moves
   */
  int makeNonLocalMovesPbyP(ParticleSet& P);

  Return_t evaluateValueAndDerivatives(ParticleSet& P,
                                       const opt_variables_type& optvars,
                                       const std::vector<ValueType>& dlogpsi,
                                       std::vector<ValueType>& dhpsioverpsi);

  /** Do nothing */
  bool put(xmlNodePtr cur) { return true; }

  bool get(std::ostream& os) const
  {
    os << "NonLocalECPotential: " << IonConfig.getName();
    return true;
  }

  OperatorBase* makeClone(ParticleSet& qp, TrialWaveFunction& psi);

  void addComponent(int groupID, NonLocalECPComponent* pp);

  /** set the internal RNG pointer as the given pointer
   * @param rng input RNG pointer
   */
  void setRandomGenerator(RandomGenerator_t* rng) { myRNG = rng; }

  void addObservables(PropertySetType& plist, BufferType& collectables);

  void setObservables(PropertySetType& plist);

  void setParticlePropertyList(PropertySetType& plist, int offset);

  void registerObservables(std::vector<observable_helper*>& h5list, hid_t gid) const;

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
  ///Pulay force vector
  ParticleSet::ParticlePos_t PulayTerm;
#if !defined(REMOVE_TRACEMANAGER)
  ///single particle trace samples
  Array<TraceReal, 1> Ve_samp_tmp;
  Array<TraceReal, 1> Vi_samp_tmp;
  Array<TraceReal, 1>* Ve_sample;
  Array<TraceReal, 1>* Vi_sample;
#endif

  /** the actual implementation, used by evaluate and evaluateWithToperator
   * @param P particle set
   * @param Tmove whether Txy for Tmove is updated
   */
  void evaluate(ParticleSet& P, bool Tmove);


  /** compute the T move transition probability for a given electron
   * member variable nonLocalOps.Txy is updated
   * @param P particle set
   * @param ref_elec reference electron id
   */
  void computeOneElectronTxy(ParticleSet& P, const int ref_elec);

  /** mark all the electrons affected by Tmoves
   * @param myTable electron ion distance table
   * @param iel reference electron
   */
  void markAffectedElecs(const DistanceTableData& myTable, int iel);
};
} // namespace qmcplusplus
#endif
