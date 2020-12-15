//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2020 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//
// File refactored from: OperatorBase.h
//////////////////////////////////////////////////////////////////////////////////////


/**@file
 */
#ifndef QMCPLUSPLUS_OPERATORESTBASE_H
#define QMCPLUSPLUS_OPERATORESTBASE_H

#include "Particle/ParticleSet.h"
#include "OhmmsData/RecordProperty.h"
#include "Utilities/RandomGenerator.h"
#include "QMCHamiltonians/observable_helper.h"
#include "QMCWaveFunctions/OrbitalSetTraits.h"
#include <bitset>

namespace qmcplusplus
{
class DistanceTableData;
class TrialWaveFunction;

/** @ingroup Estimators
 * @brief An abstract class for gridded estimators
 *
 */
class OperatorEstBase
{
public:
  using QMCT = QMCTraits;
  typedef ParticleSet::Buffer_t BufferType;
  ///typedef for the walker
  typedef ParticleSet::Walker_t Walker_t;
  ///typedef for the ParticleScalar
  typedef ParticleSet::Scalar_t ParticleScalar_t;

  ///enum to denote energy domain of operators
  enum energy_domains
  {
    kinetic = 0,
    potential,
    no_energy_domain
  };

  enum quantum_domains
  {
    no_quantum_domain = 0,
    classical,
    quantum,
    classical_classical,
    quantum_classical,
    quantum_quantum
  };

  enum data_locality
  {
    main_block = 0,
    separate
  };

  ///quantum_domain of the (particle) operator, default = no_quantum_domain
  quantum_domains quantum_domain;
  ///energy domain of the operator (kinetic/potential), default = no_energy_domain
  energy_domains energy_domain;

  ///data locality with respect to walker buffer

  ///enum for UpdateMode
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
  ///starting index of this object
  int myIndex;
  ///number of dependents: to be removed
  int Dependants;
  ///current value
  QMCT::FullPrecRealType Value;
  ///a new value for a proposed move
  QMCT::FullPrecRealType NewValue;
  /// This is used to store the value for force on the source
  /// ParticleSet.  It is accumulated if setComputeForces(true).
  ParticleSet::ParticlePos_t IonForce;
  ///reference to the current walker
  Walker_t* tWalker;
  //Walker<Return_t, ParticleSet::ParticleGradient_t>* tWalker;
  ///name of this object
  std::string myName;
  ///name of dependent object: to be removed
  std::string depName;

  ///constructor
  OperatorEstBase();

  ///virtual destructor
  virtual ~OperatorEstBase() {}

  ///set energy domain
  void set_energy_domain(energy_domains edomain);

  ///return whether the energy domain is valid
  inline bool energy_domain_valid(energy_domains edomain) const { return edomain != no_energy_domain; }

  ///return whether the energy domain is valid
  inline bool energy_domain_valid() const { return energy_domain_valid(energy_domain); }

  ///set quantum domain
  void set_quantum_domain(quantum_domains qdomain);

  ///set quantum domain for one-body operator
  void one_body_quantum_domain(const ParticleSet& P);

  ///set quantum domain for two-body operator
  void two_body_quantum_domain(const ParticleSet& P);

  ///set quantum domain for two-body operator
  void two_body_quantum_domain(const ParticleSet& P1, const ParticleSet& P2);

  ///return whether the quantum domain is valid
  bool quantum_domain_valid(quantum_domains qdomain);

  ///return whether the quantum domain is valid
  inline bool quantum_domain_valid() { return quantum_domain_valid(quantum_domain); }

  inline bool is_classical() { return quantum_domain == classical; }
  inline bool is_quantum() { return quantum_domain == quantum; }
  inline bool is_classical_classical() { return quantum_domain == classical_classical; }
  inline bool is_quantum_classical() { return quantum_domain == quantum_classical; }
  inline bool is_quantum_quantum() { return quantum_domain == quantum_quantum; }

  /** return the mode i
   * @param i index among PRIMARY, OPTIMIZABLE, RATIOUPDATE, PHYSICAL
   */
  inline bool getMode(int i) { return UpdateMode[i]; }

  inline bool isNonLocal() const { return UpdateMode[NONLOCAL]; }

  /** named values to  the property list
   * @param plist RecordNameProperty
   *
   * Previously addObservables but it is renamed and a non-virtial function.
   */
  inline void addValue(QMCT::PropertySetType& plist)
  {
    if (!UpdateMode[COLLECTABLE])
      myIndex = plist.add(myName.c_str());
  }

  /** Evaluate the local energy contribution of this component
   *@param P input configuration containing N particles
   *@return the value of the Hamiltonian component
   */
  virtual QMCT::FullPrecRealType evaluate(ParticleSet& P) = 0;
  /** Evaluate the contribution of this component of multiple walkers */
  virtual void mw_evaluate(const RefVector<OperatorEstBase>& O_list, const RefVector<ParticleSet>& P_list);

  virtual QMCT::FullPrecRealType rejectedMove(ParticleSet& P) { return 0; }
  /** Evaluate the local energy contribution of this component with Toperators updated if requested
   *@param P input configuration containing N particles
   *@return the value of the Hamiltonian component
   */
  virtual QMCT::FullPrecRealType evaluateWithToperator(ParticleSet& P) { return evaluate(P); }

  /** Evaluate the contribution of this component of multiple walkers */
  virtual void mw_evaluateWithToperator(const RefVector<OperatorEstBase>& O_list, const RefVector<ParticleSet>& P_list)
  {
    mw_evaluate(O_list, P_list);
  }

  /** evaluate value and derivatives wrt the optimizables
   *
   * Default uses evaluate
   */
  virtual QMCT::FullPrecRealType evaluateValueAndDerivatives(ParticleSet& P,
                                                             const opt_variables_type& optvars,
                                                             const std::vector<QMCT::ValueType>& dlogpsi,
                                                             std::vector<QMCT::ValueType>& dhpsioverpsi)
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
  virtual QMCT::FullPrecRealType evaluateWithIonDerivs(ParticleSet& P,
                                                       ParticleSet& ions,
                                                       TrialWaveFunction& psi,
                                                       ParticleSet::ParticlePos_t& hf_term,
                                                       ParticleSet::ParticlePos_t& pulay_term)
  {
    return evaluate(P);
  }
  /** update data associated with a particleset
   * @param s source particle set
   *
   * Default implementation does nothing. Only A-A interactions for s needs to implement its own method.
   */
  virtual void update_source(ParticleSet& s) {}

  /** return an average value by collective operation
   */
  virtual QMCT::FullPrecRealType getEnsembleAverage() { return 0.0; }

  /** write about the class */
  virtual bool get(std::ostream& os) const = 0;

  virtual OperatorEstBase* makeClone(ParticleSet& qp, TrialWaveFunction& psi) = 0;
};
} // namespace qmcplusplus
#endif
