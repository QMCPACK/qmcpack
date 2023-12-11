//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2022 QMCPACK developers.
//
// File developed by: D. Das, University of Illinois at Urbana-Champaign
//                    John R. Gergely,  University of Illinois at Urbana-Champaign
//                    Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//                    Peter W. Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


/**@file
 *@brief Declaration of OperatorBase
 */
#ifndef QMCPLUSPLUS_HAMILTONIANBASE_H
#define QMCPLUSPLUS_HAMILTONIANBASE_H

#include "Particle/ParticleSet.h"
#include "OhmmsData/RecordProperty.h"
#include "Utilities/RandomGenerator.h"
#include "QMCHamiltonians/ObservableHelper.h"
#include "Containers/MinimalContainers/RecordArray.hpp"
#include "QMCWaveFunctions/TWFFastDerivWrapper.h"
#if !defined(REMOVE_TRACEMANAGER)
#include "Estimators/TraceManager.h"
#endif
#include "QMCHamiltonians/Listener.hpp"
#include "QMCWaveFunctions/OrbitalSetTraits.h"
#include <bitset>
#include <memory> // std::unique_ptr

namespace qmcplusplus
{
class MCWalkerConfiguration;

/**@defgroup hamiltonian Hamiltonian group
 * @brief QMCHamiltonian and its component, OperatorBase
 *
 */
class TrialWaveFunction;
class QMCHamiltonian;
class ResourceCollection;
struct NonLocalData;

/** @ingroup hamiltonian
 * @brief An abstract class for Local Energy operators
 *
 * Return_t is defined as RealTye.
 * The types should be checked when using complex wave functions.
 */
class OperatorBase : public QMCTraits
{
public:
  /** type of return value of evaluate
   */
  using Return_t = FullPrecRealType;

  /** For fast derivative evaluation
   */
  using ValueMatrix = SPOSet::ValueMatrix;
  using GradMatrix  = SPOSet::GradMatrix;

  /** typedef for the serialized buffer
   *
   * PooledData<RealType> is used to serialized an anonymous buffer
   */
  using BufferType = ParticleSet::Buffer_t;

  ///typedef for the walker
  using Walker_t = ParticleSet::Walker_t;

  ///typedef for the ParticleScalar
  using ParticleScalar = ParticleSet::Scalar_t;

  ///typedef for SPOMap
  using SPOMap = SPOSet::SPOMap;

  ///enum to denote energy domain of operators
  enum EnergyDomains
  {
    KINETIC = 0,
    POTENTIAL,
    NO_ENERGY_DOMAIN
  };

  enum QuantumDomains
  {
    NO_QUANTUM_DOMAIN = 0,
    CLASSICAL,
    QUANTUM,
    CLASSICAL_CLASSICAL,
    QUANTUM_CLASSICAL,
    QUANTUM_QUANTUM
  };

  ///enum for update_mode
  enum
  {
    PRIMARY     = 0,
    OPTIMIZABLE = 1,
    RATIOUPDATE = 2,
    PHYSICAL    = 3,
    COLLECTABLE = 4,
    NONLOCAL    = 5,
  };

  /**
   * @brief Construct a new Operator Base object
   * Default and unique empty constructor. Initializes with default values.
   */
  OperatorBase();

  virtual ~OperatorBase() = default;

  /// return true if this operator depends on a wavefunction
  virtual bool dependsOnWaveFunction() const { return false; }

  //////// GETTER AND SETTER FUNCTIONS ////////////////

  /**
   * @brief get update_mode_ reference
   * 
   * @return std::bitset<8>& reference of get_update_mode_
   */
  std::bitset<8>& getUpdateMode() noexcept;

  /**
   * @brief get a copy of value_
   * 
   * @return Return_t copy of value_
   */
  Return_t getValue() const noexcept;

  /**
   * @brief getter a copy of my_name_, rvalue small string optimization
   * 
   * @return std::string copy of my_name_ member
   */
  std::string getName() const noexcept;

  /// return class name
  virtual std::string getClassName() const = 0;

  /**
   * @brief Set my_name member, uses small string optimization (pass by value)
   * 
   * @param name input
   */
  void setName(const std::string name) noexcept;

#if !defined(REMOVE_TRACEMANAGER)
  /**
   * @brief Get request_ member
   * 
   * @return TraceRequest& reference to request_
   */
  TraceRequest& getRequest() noexcept;
#endif

  //////// PURELY VIRTUAL FUNCTIONS ////////////////
  /** 
   * @brief Reset the data with the target ParticleSet
   * @param P new target ParticleSet
   */
  virtual void resetTargetParticleSet(ParticleSet& P) = 0;

  /** 
   * @brief Evaluate the local energy contribution of this component
   * @param P input configuration containing N particles
   * @return the value of the Hamiltonian component
   */
  virtual Return_t evaluate(ParticleSet& P) = 0;

  /** write about the class */
  virtual bool get(std::ostream& os) const = 0;

  // TODO: add docs
  virtual std::unique_ptr<OperatorBase> makeClone(ParticleSet& qp, TrialWaveFunction& psi) = 0;

  //////// VIRTUAL FUNCTIONS ////////////////

  /** 
   * @brief named values to  the property list
   * Default implementaton uses addValue(plist_)
   * 
   * @param plist RecordNameProperty
   * @param collectables Observables that are accumulated by evaluate
   */
  virtual void addObservables(PropertySetType& plist, BufferType& collectables);

  /** 
   * @brief add to observable descriptor for hdf5
   * The default implementation is to register a scalar for this->value_
   * 
   * @param h5desc contains a set of hdf5 descriptors for a scalar observable
   * @param gid hdf5 group to which the observables belong
   */
  virtual void registerObservables(std::vector<ObservableHelper>& h5desc, hdf_archive& file) const;

  /*** 
   * @brief add to collectables descriptor for hdf5
   * The default implementation does nothing. The derived classes which compute
   * big data, e.g. density, should overwrite this function.
   * 
   * @param h5desc contains a set of hdf5 descriptors for a scalar observable
   * @param gid hdf5 group to which the observables belong
   */
  virtual void registerCollectables(std::vector<ObservableHelper>& h5desc, hdf_archive& file) const;

  /** 
   * @brief Set the values evaluated by this object to plist
   * Default implementation is to assign Value which is updated
   * by evaluate function using my_index_.
   *
   * @param plist RecordNameProperty
   */
  virtual void setObservables(PropertySetType& plist);

  // TODO: add docs
  virtual void setParticlePropertyList(PropertySetType& plist, int offset);

  // TODO: add docs
  virtual void setHistories(Walker_t& ThisWalker);

  /** 
   * @brief Evaluate the local energy contribution of this component, deterministically based on current state.
   * The correct behavior of this routine requires estimators with non-deterministic components
   * in their evaluate() function to override this function.

   * @param P input configuration containing N particles
   * @return the value of the Hamiltonian component
   */
  virtual Return_t evaluateDeterministic(ParticleSet& P);

  /**
   * @brief Evaluate the contribution of this component of multiple walkers.
   * Take o_list and p_list update evaluation result variables in o_list?
   * really should reduce vector of local_energies. matching the ordering and size of o list
   * the this can be call for 1 or more QMCHamiltonians

   * @param o_list 
   * @param wf_list 
   * @param p_list 
   */
  virtual void mw_evaluate(const RefVectorWithLeader<OperatorBase>& o_list,
                           const RefVectorWithLeader<TrialWaveFunction>& wf_list,
                           const RefVectorWithLeader<ParticleSet>& p_list) const;

  /**
   * @brief Evaluate the contribution of this component of multiple walkers per particle and report
   * to registerd listeners from objects in Estimators
   *
   * Base class implementation decays to the mw_evaluate so if not overridden the estimator doesn't
   * hear from this operator.
   *
   * specialized versions of this should take advantage of multiwalker resources
   * to reduce the resource cost of collecting these values. 
   */
  virtual void mw_evaluatePerParticle(const RefVectorWithLeader<OperatorBase>& o_list,
                                      const RefVectorWithLeader<TrialWaveFunction>& wf_list,
                                      const RefVectorWithLeader<ParticleSet>& p_list,
                                      const std::vector<ListenerVector<RealType>>& listeners,
                                      const std::vector<ListenerVector<RealType>>& listeners_ions) const;

  /**
   * @brief TODO: add docs

   * @param o_list 
   * @param p_list 
   * @param optvars 
   * @param dlogpsi 
   * @param dhpsioverpsi 
   */
  virtual void mw_evaluateWithParameterDerivatives(const RefVectorWithLeader<OperatorBase>& o_list,
                                                   const RefVectorWithLeader<ParticleSet>& p_list,
                                                   const opt_variables_type& optvars,
                                                   const RecordArray<ValueType>& dlogpsi,
                                                   RecordArray<ValueType>& dhpsioverpsi) const;

  /**
   * @brief TODO: add docs
   * 
   * @param P 
   * @return Return_t 
   */
  virtual Return_t rejectedMove(ParticleSet& P);

  /** 
   * @brief Evaluate the local energy contribution of this component with Toperators updated if requested

   * @param P input configuration containing N particles
   * @return the value of the Hamiltonian component
   */
  virtual Return_t evaluateWithToperator(ParticleSet& P);

  /**
   * @brief Evaluate the contribution of this component of multiple walkers

   * @param o_list 
   * @param wf_list 
   * @param p_list 
   */
  virtual void mw_evaluateWithToperator(const RefVectorWithLeader<OperatorBase>& o_list,
                                        const RefVectorWithLeader<TrialWaveFunction>& wf_list,
                                        const RefVectorWithLeader<ParticleSet>& p_list) const;

  /**
   * @brief Evaluate the contribution of this component of multiple walkers per particle and report
   * to registerd listeners from objects in Estimators
   *
   * default implementation decays to the mw_evaluatePerParticle.
   *
   * specialized versions of this should take advantage of multiwalker resources
   * to reduce the resource cost of collecting these values. 
   */
  virtual void mw_evaluatePerParticleWithToperator(const RefVectorWithLeader<OperatorBase>& o_list,
                                                   const RefVectorWithLeader<TrialWaveFunction>& wf_list,
                                                   const RefVectorWithLeader<ParticleSet>& p_list,
                                                   const std::vector<ListenerVector<RealType>>& listeners,
                                                   const std::vector<ListenerVector<RealType>>& listeners_ions) const;


  /**
   * @brief Evaluate value and derivatives wrt the optimizables. Default uses evaluate.

   * @param P 
   * @param optvars 
   * @param dlogpsi 
   * @param dhpsioverpsi 
   * @return Return_t 
   */
  virtual Return_t evaluateValueAndDerivatives(ParticleSet& P,
                                               const opt_variables_type& optvars,
                                               const Vector<ValueType>& dlogpsi,
                                               Vector<ValueType>& dhpsioverpsi);

  /** 
   * @brief Evaluate contribution to local energy  and derivatives w.r.t ionic coordinates from OperatorBase.  

   * @param P target particle set (electrons)
   * @param ions source particle set (ions)
   * @param psi Trial wave function
   * @param hf_terms  Adds OperatorBase's contribution to Re [(dH)Psi]/Psi
   * @param pulay_terms Adds OperatorBase's contribution to Re [(H-E_L)dPsi]/Psi 
   * @return Contribution of OperatorBase to Local Energy.
   */
  virtual Return_t evaluateWithIonDerivs(ParticleSet& P,
                                         ParticleSet& ions,
                                         TrialWaveFunction& psi,
                                         ParticleSet::ParticlePos& hf_term,
                                         ParticleSet::ParticlePos& pulay_term);

  /** 
   * @brief Evaluate contribution to local energy  and derivatives w.r.t ionic coordinates from OperatorBase.
   * If there's no stochastic component, defaults to evaluateWithIonDerivs.
   * If not otherwise specified, this defaults to evaluate().

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
                                                      ParticleSet::ParticlePos& hf_term,
                                                      ParticleSet::ParticlePos& pulay_term);

  /** 
   * @brief Evaluate "B" matrix for observable.  Filippi scheme for computing fast derivatives.

   * @param[in] P target particle set (electrons)
   * @param[in] psi, Trial Wavefunction wrapper for fast derivatives.
   * @param[in,out] B.  List of B matrices for each species.  
   * @return Void
   */
  inline virtual void evaluateOneBodyOpMatrix(ParticleSet& P,
                                              const TWFFastDerivWrapper& psi,
                                              std::vector<ValueMatrix>& B)
  {}

  /** 
   * @brief Evaluate "dB/dR" matrices for observable.  Filippi scheme for computing fast derivatives.

   * @param[in] P, target particle set (electrons)
   * @param[in] source, ion particle set 
   * @param[in] psi, Trial Wavefunction wrapper for fast derivatives.
   * @param[in] iat, 
   * @param[in,out] dB/dR. Specifically, [ dB/dx_iat, dB/dy_iat, dB/dz_iat ], B is defined above.
   * @return Void
   */
  inline virtual void evaluateOneBodyOpMatrixForceDeriv(ParticleSet& P,
                                                        ParticleSet& source,
                                                        const TWFFastDerivWrapper& psi,
                                                        const int iat,
                                                        std::vector<std::vector<ValueMatrix>>& Bforce)
  {}
  /** 
   * @brief Update data associated with a particleset.
   * Default implementation does nothing. Only A-A interactions for s needs to implement its own method.

   * @param s source particle set
   */
  virtual void updateSource(ParticleSet& s);

  /** 
   * @brief Return an average value by collective operation
   */
  virtual Return_t getEnsembleAverage();

  /**
   * @brief Initialize a shared resource and hand it to a collection

   * @param collection 
   */
  virtual void createResource(ResourceCollection& collection) const;

  /**
   * @brief Acquire a shared resource from a collection

   * @param collection 
   * @param o_list 
   */
  virtual void acquireResource(ResourceCollection& collection, const RefVectorWithLeader<OperatorBase>& o_list) const;

  /**
   * @brief Return a shared resource to a collection
   * 
   * @param collection 
   * @param o_list 
   */
  virtual void releaseResource(ResourceCollection& collection, const RefVectorWithLeader<OperatorBase>& o_list) const;

  /**
   * @brief Set the Random Generator object
   * TODO: add docs
   * @param rng 
   */
  virtual void setRandomGenerator(RandomBase<FullPrecRealType>* rng);

  /**
   * @brief TODO: add docs
   * 
   * @param qp 
   * @param psi 
   * @param targetH 
   */
  virtual void add2Hamiltonian(ParticleSet& qp, TrialWaveFunction& psi, QMCHamiltonian& targetH);

#if !defined(REMOVE_TRACEMANAGER)
  /**
   * @brief TODO: add docs
   * 
   * @param tm 
   */
  virtual void getRequiredTraces(TraceManager& tm);
#endif

  // TODO: add docs

  virtual void informOfPerParticleListener() { has_listener_ = true; }

  bool isClassical() const noexcept;
  bool isQuantum() const noexcept;
  bool isClassicalClassical() const noexcept;
  bool isQuantumClassical() const noexcept;
  bool isQuantumQuantum() const noexcept;

  /** 
   * @brief Return the mode i
   * @param i index among PRIMARY, OPTIMIZABLE, RATIOUPDATE, PHYSICAL
   */
  bool getMode(const int i) const noexcept;

  /**
   * @brief TODO: add docs
   * 
   * @return true 
   * @return false 
   */
  bool isNonLocal() const noexcept;

  bool hasListener() const noexcept;
#if !defined(REMOVE_TRACEMANAGER)

  /**
   * @brief Make trace quantities available
   */
  void contributeTraceQuantities();

  /**
   * @brief Checkout trace arrays 
   * Derived classes must guard individual checkouts using request info 
   * @param tm 
   */
  void checkoutTraceQuantities(TraceManager& tm);

  /**
   * @brief Collect scalar trace data
   */
  void collectScalarTraces();

  /**
   * @brief delete trace arrays
   */
  void deleteTraceQuantities();
#endif

protected:
  ///set the current update mode
  std::bitset<8> update_mode_;

  ///current value
  Return_t value_;

  ///name of this object
  std::string name_;

#if !defined(REMOVE_TRACEMANAGER)
  ///whether traces are being collected
  TraceRequest request_;
#endif

  ///starting index of this object
  int my_index_;

  ///a new value for a proposed move
  Return_t new_value_;

  ///reference to the current walker
  Walker_t* t_walker_;

#if !defined(REMOVE_TRACEMANAGER)
  bool streaming_particles_;
  bool have_required_traces_;
#endif

  /////PURELY VIRTUAL FUNCTIONS

  /**
   * Read the input parameter
   * @param cur xml node for a OperatorBase object
   */
  virtual bool put(xmlNodePtr cur) = 0;

  //////VIRTUAL FUNCTIONS

#if !defined(REMOVE_TRACEMANAGER)
  virtual void contributeScalarQuantities();

  virtual void checkoutScalarQuantities(TraceManager& tm);

  virtual void collectScalarQuantities();

  virtual void deleteScalarQuantities();

  virtual void contributeParticleQuantities();
  virtual void checkoutParticleQuantities(TraceManager& tm);
  virtual void deleteParticleQuantities();
#endif

  virtual void setComputeForces(bool compute);

  /**
   * @brief Set the Energy Domain
   * 
   * @param edomain 
   */
  void setEnergyDomain(EnergyDomains edomain);

  ///set quantum domain
  void setQuantumDomain(QuantumDomains qdomain);

  ///set quantum domain for one-body operator
  void oneBodyQuantumDomain(const ParticleSet& P);

  ///set quantum domain for two-body operator
  void twoBodyQuantumDomain(const ParticleSet& P);

  ///set quantum domain for two-body operator
  void twoBodyQuantumDomain(const ParticleSet& P1, const ParticleSet& P2);

  /**
   * @brief named values to  the property list
   * @param plist RecordNameProperty
   *
   * Previously addObservables but it is renamed and a non-virtial function.
   */
  void addValue(PropertySetType& plist);

private:
#if !defined(REMOVE_TRACEMANAGER)
  bool streaming_scalars_;

  ///array to store sample value
  Array<RealType, 1>* value_sample_;
#endif

  /** Is there a per particle listener
   *  sadly this is necessary due to state machines
   */
  bool has_listener_ = false;

  ///quantum_domain_ of the (particle) operator, default = no_quantum_domain
  QuantumDomains quantum_domain_;
  ///energy domain of the operator (kinetic/potential), default = no_energy_domain
  EnergyDomains energy_domain_;

  ///return whether the energy domain is valid
  bool energyDomainValid(EnergyDomains edomain) const noexcept;

  ///return whether the quantum domain is valid
  bool quantumDomainValid(QuantumDomains qdomain) const noexcept;
};
} // namespace qmcplusplus
#endif
