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
class MCWalkerConfiguration;

/**@defgroup hamiltonian Hamiltonian group
 * @brief QMCHamiltonian and its component, OperatorBase
 *
 */
class DistanceTableData;
class TrialWaveFunction;
class QMCHamiltonian;
class ResourceCollection;

struct NonLocalData : public QMCTraits
{
  IndexType PID;
  RealType Weight;
  PosType Delta;
  inline NonLocalData() : PID(-1), Weight(1.0) {}
  inline NonLocalData(IndexType id, RealType w, const PosType& d) : PID(id), Weight(w), Delta(d) {}
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
  /** type of return value of evaluate
   */
  using Return_t = FullPrecRealType;

  /** typedef for the serialized buffer
   *
   * PooledData<RealType> is used to serialized an anonymous buffer
   */
  using BufferType = ParticleSet::Buffer_t;

  ///typedef for the walker
  using Walker_t = ParticleSet::Walker_t;

  ///typedef for the ParticleScalar
  using ParticleScalar_t = ParticleSet::Scalar_t;

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

  ///constructor
  OperatorBase();

  ///virtual destructor
  virtual ~OperatorBase() {}

  // getter for update_mode member
  /**
   * @brief get update_mode_ 
   * @return std::bitset<8>& reference of get_update_mode_
   */
  std::bitset<8>& getUpdateMode() noexcept;

  /**
   * @brief get a copy of value_
   * @return Return_t copy of value_
   */
  Return_t getValue() const noexcept;

  /**
   * @brief getter a copy of my_name_, rvalue small string optimization
   * @return std::string copy of my_name_ member
   */
  std::string getName() const noexcept;

  /**
   * @brief Set my_name member, uses small string optimization (pass by value)
   * @param name input
   */
  void setName(const std::string name) noexcept;

#if !defined(REMOVE_TRACEMANAGER)
  /**
   * @brief Get request_ member
   * @return TraceRequest& reference to request_
   */
  TraceRequest& getRequest() noexcept;
#endif

  inline bool isClassical() { return quantum_domain_ == CLASSICAL; }
  inline bool isQuantum() { return quantum_domain_ == QUANTUM; }
  inline bool isClassicalClassical() { return quantum_domain_ == CLASSICAL_CLASSICAL; }
  inline bool isQuantumClassical() { return quantum_domain_ == QUANTUM_CLASSICAL; }
  inline bool isQuantumQuantum() { return quantum_domain_ == QUANTUM_QUANTUM; }

  /** return the mode i
   * @param i index among PRIMARY, OPTIMIZABLE, RATIOUPDATE, PHYSICAL
   */
  inline bool getMode(int i) { return update_mode_[i]; }

  inline bool isNonLocal() const { return update_mode_[NONLOCAL]; }

  /** named values to  the property list
   * @param plist RecordNameProperty
   * @param collectables Observables that are accumulated by evaluate
   *
   * Default implementaton uses addValue(plist)
   */
  virtual void addObservables(PropertySetType& plist, BufferType& collectables) { addValue(plist); }

  /*** add to observable descriptor for hdf5
   * @param h5desc contains a set of hdf5 descriptors for a scalar observable
   * @param gid hdf5 group to which the observables belong
   *
   * The default implementation is to register a scalar for this->Value
   */
  virtual void registerObservables(std::vector<ObservableHelper>& h5desc, hid_t gid) const;

  /*** add to collectables descriptor for hdf5
   * @param h5desc contains a set of hdf5 descriptors for a scalar observable
   * @param gid hdf5 group to which the observables belong
   *
   * The default implementation does nothing. The derived classes which compute
   * big data, e.g. density, should overwrite this function.
   */
  virtual void registerCollectables(std::vector<ObservableHelper>& h5desc, hid_t gid) const {}

  /** set the values evaluated by this object to plist
   * @param plist RecordNameProperty
   *
   * Default implementation is to assign Value which is updated
   * by evaluate  function using myIndex.
   */
  virtual void setObservables(PropertySetType& plist) { plist[my_index_] = value_; }

  virtual void setParticlePropertyList(PropertySetType& plist, int offset) { plist[my_index_ + offset] = value_; }

  //virtual void setHistories(Walker<Return_t, ParticleSet::ParticleGradient_t>& ThisWalker)
  virtual void setHistories(Walker_t& ThisWalker) { t_walker_ = &(ThisWalker); }

  /** reset the data with the target ParticleSet
   * @param P new target ParticleSet
   */
  virtual void resetTargetParticleSet(ParticleSet& P) = 0;

  /** Evaluate the local energy contribution of this component
   *@param P input configuration containing N particles
   *@return the value of the Hamiltonian component
   */
  virtual Return_t evaluate(ParticleSet& P) = 0;
  /** Evaluate the local energy contribution of this component, deterministically based on current state.
   *@param P input configuration containing N particles
   *@return the value of the Hamiltonian component
   */
  virtual Return_t evaluateDeterministic(ParticleSet& P);
  /** Evaluate the contribution of this component of multiple walkers */
  virtual void mw_evaluate(const RefVectorWithLeader<OperatorBase>& o_list,
                           const RefVectorWithLeader<TrialWaveFunction>& wf_list,
                           const RefVectorWithLeader<ParticleSet>& p_list) const;

  virtual void mw_evaluateWithParameterDerivatives(const RefVectorWithLeader<OperatorBase>& o_list,
                                                   const RefVectorWithLeader<ParticleSet>& p_list,
                                                   const opt_variables_type& optvars,
                                                   RecordArray<ValueType>& dlogpsi,
                                                   RecordArray<ValueType>& dhpsioverpsi) const;


  virtual Return_t rejectedMove(ParticleSet& P) { return 0; }
  /** Evaluate the local energy contribution of this component with Toperators updated if requested
   *@param P input configuration containing N particles
   *@return the value of the Hamiltonian component
   */
  virtual Return_t evaluateWithToperator(ParticleSet& P) { return evaluate(P); }

  /** Evaluate the contribution of this component of multiple walkers */
  virtual void mw_evaluateWithToperator(const RefVectorWithLeader<OperatorBase>& o_list,
                                        const RefVectorWithLeader<TrialWaveFunction>& wf_list,
                                        const RefVectorWithLeader<ParticleSet>& p_list) const
  {
    mw_evaluate(o_list, wf_list, p_list);
  }

  /** evaluate value and derivatives wrt the optimizables
   *
   * Default uses evaluate
   */
  virtual Return_t evaluateValueAndDerivatives(ParticleSet& P,
                                               const opt_variables_type& optvars,
                                               const std::vector<ValueType>& dlogpsi,
                                               std::vector<ValueType>& dhpsioverpsi)
  {
    return evaluate(P);
  }

  /** evaluate contribution to local energy  and derivatives w.r.t ionic coordinates from OperatorBase.  
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
                                         ParticleSet::ParticlePos_t& hf_term,
                                         ParticleSet::ParticlePos_t& pulay_term)
  {
    return evaluate(P);
  }

  /** evaluate contribution to local energy  and derivatives w.r.t ionic coordinates from OperatorBase.  
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
                                                      ParticleSet::ParticlePos_t& pulay_term)
  {
    //If there's no stochastic component, defaults to above defined evaluateWithIonDerivs.
    //If not otherwise specified, this defaults to evaluate().
    return evaluateWithIonDerivs(P, ions, psi, hf_term, pulay_term);
  }
  /** update data associated with a particleset
   * @param s source particle set
   *
   * Default implementation does nothing. Only A-A interactions for s needs to implement its own method.
   */
  virtual void updateSource(ParticleSet& s) {}

  /** return an average value by collective operation
   */
  virtual Return_t getEnsembleAverage() { return 0.0; }

  /** write about the class */
  virtual bool get(std::ostream& os) const = 0;

  /** initialize a shared resource and hand it to a collection
   */
  virtual void createResource(ResourceCollection& collection) const {}

  /** acquire a shared resource from a collection
   */
  virtual void acquireResource(ResourceCollection& collection, const RefVectorWithLeader<OperatorBase>& o_list) const {}

  /** return a shared resource to a collection
   */
  virtual void releaseResource(ResourceCollection& collection, const RefVectorWithLeader<OperatorBase>& o_list) const {}

  virtual std::unique_ptr<OperatorBase> makeClone(ParticleSet& qp, TrialWaveFunction& psi) = 0;

  //virtual std::unique_ptr<OperatorBase> makeClone(ParticleSet& qp, TrialWaveFunction& psi, QMCHamiltonian& H);

  virtual void setRandomGenerator(RandomGenerator_t* rng)
  {
    //empty
  }

  virtual void add2Hamiltonian(ParticleSet& qp, TrialWaveFunction& psi, QMCHamiltonian& targetH);

#if !defined(REMOVE_TRACEMANAGER)
  ///make trace quantities available
  inline void contributeTraceQuantities()
  {
    contributeScalarQuantities();
    contributeParticleQuantities();
  }

  ///checkout trace arrays
  inline void checkoutTraceQuantities(TraceManager& tm)
  {
    //derived classes must guard individual checkouts using request info
    checkoutScalarQuantities(tm);
    checkoutParticleQuantities(tm);
  }

  ///collect scalar trace data
  inline void collectScalarTraces()
  {
    //app_log()<<"OperatorBase::collectScalarTraces"<< std::endl;
    collectScalarQuantities();
  }

  ///delete trace arrays
  inline void deleteTraceQuantities()
  {
    deleteScalarQuantities();
    deleteParticleQuantities();
    streaming_scalars_    = false;
    streaming_particles_  = false;
    have_required_traces_ = false;
    request_.reset();
  }

  virtual void getRequiredTraces(TraceManager& tm){};
#endif

  virtual void addEnergy(MCWalkerConfiguration& W, std::vector<RealType>& LocalEnergy);

  virtual void addEnergy(MCWalkerConfiguration& W,
                         std::vector<RealType>& LocalEnergy,
                         std::vector<std::vector<NonLocalData>>& Txy)
  {
    addEnergy(W, LocalEnergy);
  }

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

  /**
   * Read the input parameter
   * @param cur xml node for a OperatorBase object
   */
  virtual bool put(xmlNodePtr cur) = 0;

#if !defined(REMOVE_TRACEMANAGER)
  virtual void contributeScalarQuantities() { request_.contribute_scalar(name_); }

  virtual void checkoutScalarQuantities(TraceManager& tm)
  {
    streaming_scalars_ = request_.streaming_scalar(name_);
    if (streaming_scalars_)
      value_sample_ = tm.checkout_real<1>(name_);
  }

  virtual void collectScalarQuantities()
  {
    if (streaming_scalars_)
      (*value_sample_)(0) = value_;
  }

  virtual void deleteScalarQuantities()
  {
    if (streaming_scalars_)
      delete value_sample_;
  }

  virtual void contributeParticleQuantities(){};
  virtual void checkoutParticleQuantities(TraceManager& tm){};
  virtual void deleteParticleQuantities(){};
#endif

  virtual void setComputeForces(bool compute)
  {
    // empty
  }

  ///set energy domain
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
   * named values to  the property list
   * @param plist RecordNameProperty
   *
   * Previously addObservables but it is renamed and a non-virtial function.
   */
  inline void addValue(PropertySetType& plist)
  {
    if (!update_mode_[COLLECTABLE])
      my_index_ = plist.add(name_.c_str());
  }

private:
  ///quantum_domain_ of the (particle) operator, default = no_quantum_domain
  QuantumDomains quantum_domain_;
  ///energy domain of the operator (kinetic/potential), default = no_energy_domain
  EnergyDomains energy_domain_;

#if !defined(REMOVE_TRACEMANAGER)
  bool streaming_scalars_;

  ///array to store sample value
  Array<RealType, 1>* value_sample_;
#endif

  ///return whether the energy domain is valid
  inline bool energyDomainValid(EnergyDomains edomain) const { return edomain != NO_ENERGY_DOMAIN; }

  ///return whether the energy domain is valid
  inline bool energyDomainValid() const { return energyDomainValid(energy_domain_); }

  ///return whether the quantum domain is valid
  bool quantumDomainValid(QuantumDomains qdomain);

  ///return whether the quantum domain is valid
  inline bool quantumDomainValid() { return quantumDomainValid(quantum_domain_); }
};
} // namespace qmcplusplus
#endif
