//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: D. Das, University of Illinois at Urbana-Champaign
//                    John R. Gergely,  University of Illinois at Urbana-Champaign
//                    Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


/**@file OperatorBase.h
 *@brief Declaration of OperatorBase
 */
#ifndef QMCPLUSPLUS_HAMILTONIANBASE_H
#define QMCPLUSPLUS_HAMILTONIANBASE_H

#include "Particle/ParticleSet.h"
#include "OhmmsData/RecordProperty.h"
#include "Utilities/RandomGenerator.h"
#include "QMCHamiltonians/ObservableHelper.h"
#include "Containers/MinimalContainers/RecordArray.hpp"
#if !defined(REMOVE_TRACEMANAGER)
#include "Estimators/TraceManager.h"
#endif
#include "QMCWaveFunctions/OrbitalSetTraits.h"

#include <bitset>
#include <memory> // std::unique_ptr

namespace qmcplusplus
{
/**
 * @defgroup hamiltonian Hamiltonian group
 * @brief QMCHamiltonian and its component, OperatorBase
 */
class MCWalkerConfiguration;
class DistanceTableData;
class TrialWaveFunction;
class QMCHamiltonian;
class ResourceCollection;

// TODO add documentation
struct NonLocalData : public QMCTraits
{
  IndexType pid;
  RealType weight;
  PosType delta;
  NonLocalData();
  NonLocalData(IndexType id, RealType w, const PosType& d);
};

/** @ingroup hamiltonian
 * @brief An abstract class for Local Energy operators
 *
 * Return_t is defined as RealTye.
 * The types should be checked when using complex wave functions.
 */
class OperatorBase : public QMCTraits
{
public:
  /// type of return value of evaluate
  using Return_t = FullPrecRealType;
  /**
   * alias for the serialized buffer
   * PooledData<RealType> is used to serialized an anonymous buffer
   */
  using BufferType = ParticleSet::Buffer_t;
  /// alias for the walker
  using Walker_t = ParticleSet::Walker_t;
  /// alias for the ParticleScalar
  using ParticleScalar_t = ParticleSet::Scalar_t;

  ///enum class to denote energy domain of operators
  enum class energy_domains
  {
    kinetic = 0,
    potential,
    no_energy_domain
  };

  ///enum class to denote quantum domain of operators
  enum class quantum_domains
  {
    no_quantum_domain = 0,
    classical,
    quantum,
    classical_classical,
    quantum_classical,
    quantum_quantum
  };

  /// UpdateModes
  static constexpr int PRIMARY     = 0;
  static constexpr int OPTIMIZABLE = 1;
  static constexpr int RATIOUPDATE = 2;
  static constexpr int PHYSICAL    = 3;
  static constexpr int COLLECTABLE = 4;
  static constexpr int NONLOCAL    = 5;

  ///set the current update mode
  std::bitset<8> updateMode;

  ///current value TODO: add better docs
  Return_t value;

  ///name of this object
  std::string myName;

#if !defined(REMOVE_TRACEMANAGER)
  ///whether traces are being collected
  TraceRequest request;
#endif

  ///constructor
  OperatorBase();

  ///virtual destructor
  virtual ~OperatorBase() = default;

  /** write about the class */
  virtual bool get(std::ostream& os) const = 0;

  /**
   * Reset the data with the target ParticleSet
   * @param P new target ParticleSet
   */
  virtual void resetTargetParticleSet(ParticleSet& P) = 0;

  /**
   * Evaluate the local energy contribution of this component
   * @param P input configuration containing N particles
   * @return the value of the Hamiltonian component
   */
  virtual Return_t evaluate(ParticleSet& P) = 0;

  // TODO add documentation
  // FIXME this should be a protected function, only QMCDrivers/WFOpt/CostFunctionCrowdData calls it publicly
  virtual std::unique_ptr<OperatorBase> makeClone(ParticleSet& qp, TrialWaveFunction& psi) = 0;

  /**
   * Return an average value by collective operation
   */
  virtual Return_t getEnsembleAverage();

  /**
   * update data associated with a particleset
   * @param s
   * Default implementation does nothing. Only A-A interactions for s needs to implement its own method.
   */
  virtual void updateSource(ParticleSet& s);

  /**
   * set the values evaluated by this object to plist
   * @param plist RecordNameProperty
   *
   * Default implementation is to assign Value which is updated
   * by evaluate  function using myIndex.
   */
  // TODO() noexcept
  virtual void setObservables(PropertySetType& plist);

  /**
   * named values to  the property list
   * @param plist RecordNameProperty
   * @param collectables Observables that are accumulated by evaluate
   * Default implementaton uses addValue(plist)
   */
  virtual void addObservables(PropertySetType& plist, BufferType& collectables);

  /**
   * add to observable descriptor for hdf5
   * @param h5desc contains a set of hdf5 descriptors for a scalar observable
   * @param gid hdf5 group to which the observables belong
   *
   * The default implementation is to register a scalar for this->Value
   */
  virtual void registerObservables(std::vector<ObservableHelper>& h5desc, hid_t gid) const;

  /**
   * add to collectables descriptor for hdf5
   * @param h5desc contains a set of hdf5 descriptors for a scalar observable
   * @param gid hdf5 group to which the observables belong
   *
   * The default implementation does nothing. The derived classes which compute
   * big data, e.g. density, should overwrite this function.
   */
  virtual void registerCollectables(std::vector<ObservableHelper>& h5desc, hid_t gid) const;

  // TODO add documentation
  virtual void getRequiredTraces(TraceManager& tm);

  // TODO add documentation
  virtual void setParticlePropertyList(PropertySetType& plist, int offset);

  // TODO add documentation
  virtual void setRandomGenerator(RandomGenerator_t* rng);

  /**
     * Evaluate the local energy contribution of this component, deterministically based on current state.
     * The correct behavior of this routine requires estimators with non-deterministic components
     * in their evaluate() function to override this function.
     * @param P input configuration containing N particles
     * @return the value of the Hamiltonian component
     */
  virtual Return_t evaluateDeterministic(ParticleSet& P);

  /**
     * Evaluate the local energy contribution of this component with Toperators updated if requested
     * @param P input configuration containing N particles
     * @return the value of the Hamiltonian component
     */
  virtual Return_t evaluateWithToperator(ParticleSet& P);

  /**
   * evaluate contribution to local energy  and derivatives w.r.t ionic coordinates from OperatorBase.
   * @param P target particle set (electrons)
   * @param ions source particle set (ions)
   * @param psi Trial wave function
   * @param hf_terms  Adds OperatorBase's contribution to Re [(dH)Psi]/Psi
   * @param pulay_terms Adds OperatorBase's contribution to Re [(H-E_L)dPsi]/Psi
   * @return Contribution of OperatorBase to Local Energy.
   * FIXME only P argument is used
   */
  virtual Return_t evaluateWithIonDerivs(ParticleSet& P,
                                         ParticleSet& ions,
                                         TrialWaveFunction& psi,
                                         ParticleSet::ParticlePos_t& hf_term,
                                         ParticleSet::ParticlePos_t& pulay_term);

  /**
     * Evaluate contribution to local energy  and derivatives w.r.t ionic coordinates from OperatorBase.
     * @param P target particle set (electrons)
     * @param ions source particle set (ions)
     * @param psi Trial wave function
     * @param hf_terms  Adds OperatorBase's contribution to Re [(dH)Psi]/Psi
     * @param pulay_terms Adds OperatorBase's contribution to Re [(H-E_L)dPsi]/Psi
     * @return Contribution of OperatorBase to Local Energy.
     */
  virtual Return_t evaluateWithIonDerivsDeterministic(ParticleSet& P,
                                                      ParticleSet& ions,
                                                      TrialWaveFunction& psi,
                                                      ParticleSet::ParticlePos_t& hf_term,
                                                      ParticleSet::ParticlePos_t& pulay_term);

  /**
     * Evaluate value and derivatives wrt the optimizables
     *
     * Default uses evaluate
     * FIXME only P argument is used
     */
  virtual Return_t evaluateValueAndDerivatives(ParticleSet& P,
                                               const opt_variables_type& optvars,
                                               const std::vector<ValueType>& dlogpsi,
                                               std::vector<ValueType>& dhpsioverpsi);

  /**
     * Evaluate the contribution of this component of multiple walkers
     * Take o_list and p_list update evaluation result variables in o_list?
     * really should reduce vector of local_energies. matching the ordering and size of o list
     * the this can be call for 1 or more QMCHamiltonians
     * @param O_list
     * @param P_list
     * TODO add parameter documentation above
     */
  virtual void mwEvaluate(const RefVectorWithLeader<OperatorBase>& O_list,
                           const RefVectorWithLeader<ParticleSet>& P_list) const;

  // TODO add documentation
  virtual void mwEvaluateWithParameterDerivatives(const RefVectorWithLeader<OperatorBase>& O_list,
                                                   const RefVectorWithLeader<ParticleSet>& P_list,
                                                   const opt_variables_type& optvars,
                                                   RecordArray<ValueType>& dlogpsi,
                                                   RecordArray<ValueType>& dhpsioverpsi) const;

  /**
     * Evaluate the contribution of this component of multiple walkers
     * @param O_list
     * @param P_list
     * TODO add documentation
     */
  virtual void mwEvaluateWithToperator(const RefVectorWithLeader<OperatorBase>& O_list,
                                        const RefVectorWithLeader<ParticleSet>& P_list) const;

  // TODO add documentation
  virtual void setHistories(Walker_t& ThisWalker);

  // TODO add documentation
  virtual Return_t rejectedMove(ParticleSet& P);

  // TODO add documentation
  virtual void add2Hamiltonian(ParticleSet& qp, TrialWaveFunction& psi, QMCHamiltonian& targetH);

  /**
     * initialize a shared resource and hand it to a collection, default does nothing
     * @param collection
     */
  virtual void createResource(ResourceCollection& collection) const;

  /**
     * acquire a shared resource from a collection, default does nothing
     * @param collection
     * @param O_list
     */
  virtual void acquireResource(ResourceCollection& collection, const RefVectorWithLeader<OperatorBase>& O_list) const;

  /**
     * return a shared resource to a collection, default does nothing
     * @param collection
     * @param O_list
     */
  virtual void releaseResource(ResourceCollection& collection, const RefVectorWithLeader<OperatorBase>& O_list) const;

  // TODO add documentation
  bool isClassical() const noexcept;
  bool isQuantum() const noexcept;
  bool isClassicalClassical() const noexcept;
  bool isQuantumClassical() const noexcept;
  bool isQuantumQuantum() const noexcept;
  bool isNonLocal() const noexcept;

  /**
   * Return the mode i
   * @param i index among PRIMARY, OPTIMIZABLE, RATIOUPDATE, PHYSICAL
   */
  bool getMode(int i) const noexcept;

#if !defined(REMOVE_TRACEMANAGER)
  ///make trace quantities available
  void contributeTraceQuantities();

  ///collect scalar trace data
  void collectScalarTraces();
#endif

  ///checkout trace arrays
  void checkoutTraceQuantities(TraceManager& tm);

  ///delete trace arrays
  void deleteTraceQuantities();


protected:
  /// starting index of this object
  int myIndex;

  /// a new value for a proposed move
  Return_t newValue;

  ///reference to the current walker
  Walker_t* tWalker = nullptr;

#if !defined(REMOVE_TRACEMANAGER)
  // TODO add documentation
  bool haveRequiredTraces;

  bool streamingParticles;
#endif

  /**
   * Read the input parameter
   * @param cur xml node for a OperatorBase object
   */
  virtual bool put(xmlNodePtr cur) = 0;

  //TODO add documentation
#if !defined(REMOVE_TRACEMANAGER)
  virtual void contributeScalarQuantities();

  virtual void checkoutScalarQuantities(TraceManager& tm);

  virtual void collectScalarQuantities();

  virtual void deleteScalarQuantities();

  virtual void contributeParticleQuantities();

  virtual void checkoutParticleQuantities(TraceManager& tm);

  virtual void deleteParticleQuantities();
#endif

  virtual void addEnergy(MCWalkerConfiguration& W, std::vector<RealType>& LocalEnergy);

  virtual void addEnergy(MCWalkerConfiguration& W,
                         std::vector<RealType>& LocalEnergy,
                         std::vector<std::vector<NonLocalData>>& Txy);

  //FIXME: only implemented by NonLocalECPotential
  virtual void setComputeForces(bool compute);

  ///set energy domain
  void setEnergyDomain(energy_domains edomain);

  ///set quantum domain
  void setQuantumDomain(quantum_domains qdomain);

  ///set quantum domain for one-body operator
  void oneBodyQuantumDomain(const ParticleSet& P);

  ///set quantum domain for two-body operator
  void twoBodyQuantumDomain(const ParticleSet& P);

  ///set quantum domain for two-body operator
  void twoBodyQuantumDomain(const ParticleSet& P1, const ParticleSet& P2);

  /**
     * Named values to  the property list
     * Previously addObservables but it is renamed and a non-virtial function.
     * @param plist RecordNameProperty
     */
  void addValue(PropertySetType& plist);


private:
  ///quantum_domain of the (particle) operator, default = no_quantum_domain
  quantum_domains quantumDomain;

  ///energy domain of the operator (kinetic/potential), default = no_energy_domain
  energy_domains energyDomain;

#if !defined(REMOVE_TRACEMANAGER)
  //TODO add documentation
  bool streamingScalars;

  std::vector<RealType> valueVector;

  ///array to store sample value
  Array<RealType, 1>* valueSample;
#endif

  ///return whether the energy domain is valid
  bool energyDomainValid(const energy_domains edomain) const noexcept;

  ///return whether the energy domain is valid
  bool energyDomainValid() const noexcept;

  ///return whether the quantum domain is valid
  bool quantumDomainValid(const quantum_domains qdomain) const noexcept;

  ///return whether the quantum domain is valid
  bool quantumDomainValid() const noexcept;
};
} // namespace qmcplusplus
#endif
