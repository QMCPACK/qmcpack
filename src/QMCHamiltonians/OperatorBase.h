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

  ///set the current update mode
  std::bitset<8> UpdateMode;

  ///number of dependents: to be removed
  int Dependants;
  ///current value
  Return_t Value;

  /// This is used to store the value for force on the source
  /// ParticleSet.  It is accumulated if setComputeForces(true).
  ParticleSet::ParticlePos_t IonForce;

  //Walker<Return_t, ParticleSet::ParticleGradient_t>* tWalker;
  ///name of this object
  std::string myName;
  ///name of dependent object: to be removed
  std::string depName;

#if !defined(REMOVE_TRACEMANAGER)
  ///whether traces are being collected
  TraceRequest request;

  std::vector<RealType> ValueVector;

#endif

  ///constructor
  OperatorBase();

  ///virtual destructor
  virtual ~OperatorBase() {}

  inline bool is_classical() { return quantum_domain == CLASSICAL; }
  inline bool is_quantum() { return quantum_domain == QUANTUM; }
  inline bool is_classical_classical() { return quantum_domain == CLASSICAL_CLASSICAL; }
  inline bool is_quantum_classical() { return quantum_domain == QUANTUM_CLASSICAL; }
  inline bool is_quantum_quantum() { return quantum_domain == QUANTUM_QUANTUM; }

  /** return the mode i
   * @param i index among PRIMARY, OPTIMIZABLE, RATIOUPDATE, PHYSICAL
   */
  inline bool getMode(int i) { return UpdateMode[i]; }

  inline bool isNonLocal() const { return UpdateMode[NONLOCAL]; }

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
  virtual void setObservables(PropertySetType& plist) { plist[myIndex] = Value; }

  virtual void setParticlePropertyList(PropertySetType& plist, int offset) { plist[myIndex + offset] = Value; }

  //virtual void setHistories(Walker<Return_t, ParticleSet::ParticleGradient_t>& ThisWalker)
  virtual void setHistories(Walker_t& ThisWalker) { tWalker = &(ThisWalker); }

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
  virtual void update_source(ParticleSet& s) {}

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
  inline void contribute_trace_quantities()
  {
    contribute_scalar_quantities();
    contribute_particle_quantities();
  }

  ///checkout trace arrays
  inline void checkout_trace_quantities(TraceManager& tm)
  {
    //derived classes must guard individual checkouts using request info
    checkout_scalar_quantities(tm);
    checkout_particle_quantities(tm);
  }

  ///collect scalar trace data
  inline void collect_scalar_traces()
  {
    //app_log()<<"OperatorBase::collect_scalar_traces"<< std::endl;
    collect_scalar_quantities();
  }

  ///delete trace arrays
  inline void delete_trace_quantities()
  {
    delete_scalar_quantities();
    delete_particle_quantities();
    streaming_scalars    = false;
    streaming_particles  = false;
    have_required_traces = false;
    request.reset();
  }

  virtual void get_required_traces(TraceManager& tm){};
#endif

  virtual void addEnergy(MCWalkerConfiguration& W, std::vector<RealType>& LocalEnergy);

  virtual void addEnergy(MCWalkerConfiguration& W,
                         std::vector<RealType>& LocalEnergy,
                         std::vector<std::vector<NonLocalData>>& Txy)
  {
    addEnergy(W, LocalEnergy);
  }

protected:
  ///starting index of this object
  int myIndex;

  ///a new value for a proposed move
  Return_t NewValue;

  ///reference to the current walker
  Walker_t* tWalker;

#if !defined(REMOVE_TRACEMANAGER)
  bool streaming_particles;
  bool have_required_traces;
#endif

  /**
   * Read the input parameter
   * @param cur xml node for a OperatorBase object
   */
  virtual bool put(xmlNodePtr cur) = 0;

#if !defined(REMOVE_TRACEMANAGER)
  virtual void contribute_scalar_quantities() { request.contribute_scalar(myName); }

  virtual void checkout_scalar_quantities(TraceManager& tm)
  {
    streaming_scalars = request.streaming_scalar(myName);
    if (streaming_scalars)
      value_sample = tm.checkout_real<1>(myName);
  }

  virtual void collect_scalar_quantities()
  {
    if (streaming_scalars)
      (*value_sample)(0) = Value;
  }

  virtual void delete_scalar_quantities()
  {
    if (streaming_scalars)
      delete value_sample;
  }

  virtual void contribute_particle_quantities(){};
  virtual void checkout_particle_quantities(TraceManager& tm){};
  virtual void delete_particle_quantities(){};
#endif

  virtual void setComputeForces(bool compute)
  {
    // empty
  }

  ///set energy domain
  void set_energy_domain(EnergyDomains edomain);

  ///set quantum domain
  void set_quantum_domain(QuantumDomains qdomain);

  ///set quantum domain for one-body operator
  void one_body_quantum_domain(const ParticleSet& P);

  ///set quantum domain for two-body operator
  void two_body_quantum_domain(const ParticleSet& P);

  ///set quantum domain for two-body operator
  void two_body_quantum_domain(const ParticleSet& P1, const ParticleSet& P2);

  /**
   * named values to  the property list
   * @param plist RecordNameProperty
   *
   * Previously addObservables but it is renamed and a non-virtial function.
   */
  inline void addValue(PropertySetType& plist)
  {
    if (!UpdateMode[COLLECTABLE])
      myIndex = plist.add(myName.c_str());
  }

private:
  ///quantum_domain of the (particle) operator, default = no_quantum_domain
  QuantumDomains quantum_domain;
  ///energy domain of the operator (kinetic/potential), default = no_energy_domain
  EnergyDomains energy_domain;

#if !defined(REMOVE_TRACEMANAGER)
  bool streaming_scalars;

  ///array to store sample value
  Array<RealType, 1>* value_sample;
#endif

  ///return whether the energy domain is valid
  inline bool energy_domain_valid(EnergyDomains edomain) const { return edomain != NO_ENERGY_DOMAIN; }

  ///return whether the energy domain is valid
  inline bool energy_domain_valid() const { return energy_domain_valid(energy_domain); }

  ///return whether the quantum domain is valid
  bool quantum_domain_valid(QuantumDomains qdomain);

  ///return whether the quantum domain is valid
  inline bool quantum_domain_valid() { return quantum_domain_valid(quantum_domain); }
};
} // namespace qmcplusplus
#endif
