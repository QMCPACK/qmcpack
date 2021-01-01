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

#include <variant>

#include "Particle/ParticleSet.h"
#include "OhmmsData/RecordProperty.h"
#include "Utilities/RandomGenerator.h"
#include "QMCHamiltonians/observable_helper.h"
#include "QMCWaveFunctions/OrbitalSetTraits.h"
#include "type_traits/DataLocality.h"
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
  using QMCT      = QMCTraits;
  using MCPWalker = Walker<QMCTraits, PtclOnLatticeTraits>;

  /** the type in this variant changes based on data locality
   */
  using Data = std::variant<std::unique_ptr<std::vector<QMCT::RealType>>, std::shared_ptr<std::vector<QMCT::RealType>>>;

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

  /// locality for accumulation data
  DataLocality data_locality_;

  /// quantum_domain of the (particle) operator, default = no_quantum_domain
  quantum_domains quantum_domain;
  /// energy domain of the operator (kinetic/potential), default = no_energy_domain
  energy_domains energy_domain;

  /// enum for UpdateMode
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
  ///current value
  QMCT::FullPrecRealType Value;
  ///a new value for a proposed move
  QMCT::FullPrecRealType NewValue;
  /// This is used to store the value for force on the source
  /// ParticleSet.  It is accumulated if setComputeForces(true).
  ParticleSet::ParticlePos_t IonForce;
  //Walker<Return_t, ParticleSet::ParticleGradient_t>* tWalker;
  ///name of this object
  std::string myName;

  QMCT::FullPrecRealType walkers_weight_;

  QMCT::FullPrecRealType get_walkers_weight() const { return walkers_weight_; }
  ///constructor
  OperatorEstBase();
  OperatorEstBase(const OperatorEstBase& oth);
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

  /** Accumulate whatever it is you are accumulating with respect to walkers
   */
  virtual void accumulate(RefVector<MCPWalker>& walkers, RefVector<ParticleSet>& psets) = 0;

  virtual void collect(const OperatorEstBase& oeb) = 0;

  virtual void write() = 0;

  Data& get_data() { return data_; };
  /*** add to OperatorEstimator descriptor for hdf5
   * @param h5desc contains a set of hdf5 descriptors for a scalar observable
   * @param gid hdf5 group to which the observables belong
   *
   * The default implementation does nothing. The derived classes which compute
   * big data, e.g. density, should overwrite this function.
   */
  virtual void registerOperatorEstimator(std::vector<observable_helper*>& h5desc, hid_t gid) const {}

  virtual OperatorEstBase* clone() = 0;


protected:
  Data data_;
};
} // namespace qmcplusplus
#endif
