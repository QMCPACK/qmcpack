//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2020 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//
// File refactored from: QMCDriver.cpp
//////////////////////////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_QMCDRIVERINPUT_H
#define QMCPLUSPLUS_QMCDRIVERINPUT_H

#include <optional>

#include "Configuration.h"
#include "OhmmsData/ParameterSet.h"
#include "InputTypes.hpp"
#include "DriverDebugChecks.h"
#include "EstimatorManagerInput.h"
#include "type_traits/template_types.hpp"

namespace qmcplusplus
{
/** Input representation for Driver base class runtime parameters
 */
class QMCDriverInput
{
public:
  using IndexType             = QMCTraits::IndexType;
  using RealType              = QMCTraits::RealType;
  using FullPrecisionRealType = QMCTraits::FullPrecRealType;

  void readXML(xmlNodePtr cur);

  // To allow compile check if move constructor is still implicit
  QMCDriverInput()                      = default;
  QMCDriverInput(const QMCDriverInput&) = default;
  QMCDriverInput& operator=(const QMCDriverInput&) = default;
  QMCDriverInput(QMCDriverInput&&) noexcept;
  QMCDriverInput& operator=(QMCDriverInput&&) noexcept;

protected:
  bool scoped_profiling_ = false;
  /// determine additional checks for debugging purpose
  DriverDebugChecks debug_checks_ = DriverDebugChecks::ALL_OFF;
  /// measure load imbalance (add a barrier) before data aggregation (obvious synchronization)
  bool measure_imbalance_ = false;

  /** @ingroup Input Parameters for QMCDriver base class
   *  @{
   *  All input determined variables should be here
   *  They are read only for Drivers
   *  Do not write out blocks of gets for variables like this
   *  there will be code generation snippets soon in utils
   */

  /// if true, batched operations are serialized over walkers
  bool crowd_serialize_walkers_ = false;
  /// period to recalculate the walker properties from scratch.
  int recalculate_properties_period_ = 100;
  /// period of recording walker positions and IDs for forward walking afterwards
  input::PeriodStride config_dump_period_;
  IndexType starting_step_ = 0;
  IndexType num_crowds_    = 0;
  // This is the global walkers it is a hard limit for VMC and the target for DMC
  IndexType total_walkers_     = 0;
  IndexType walkers_per_rank_  = 0;
  IndexType requested_samples_ = 0;
  IndexType sub_steps_         = 1;
  // max unnecessary in this context
  IndexType max_blocks_            = 1;
  IndexType max_steps_             = 1;
  IndexType warmup_steps_          = 0;
  IndexType steps_between_samples_ = 1;
  IndexType samples_per_thread_    = 0;
  RealType tau_                    = 0.1;
  RealType spin_mass_              = 1.0;
  // call recompute at the end of each block in the full/mixed precision case.
  IndexType blocks_between_recompute_ = std::is_same<RealType, FullPrecisionRealType>::value ? 0 : 1;
  bool append_run_                    = false;

  // from QMCDriverFactory
  std::string qmc_method_{"invalid"};
  std::string update_mode_{"pbyp"};

  /** The EstimatorManagerInput for batched version input
   */
  std::optional<EstimatorManagerInput> estimator_manager_input_;

  // from putQMCInfo
  input::PeriodStride walker_dump_period_{0, 0};
  input::PeriodStride check_point_period_{0, 0};
  bool dump_config_  = false;
  IndexType k_delay_ = 0;
  bool reset_random_ = false;

  // from QMCUpdateBase
  RealType max_disp_sq_ = -1.0;

  // for drift modifer
  std::string drift_modifier_{"UNR"};
  RealType drift_modifier_unr_a_ = 1.0;

  /** @}
   */

public:
  int get_recalculate_properties_period() const { return recalculate_properties_period_; }
  input::PeriodStride get_config_dump_period() const { return config_dump_period_; }
  IndexType get_starting_step() const { return starting_step_; }
  IndexType get_num_crowds() const { return num_crowds_; }
  IndexType get_walkers_per_rank() const { return walkers_per_rank_; }
  IndexType get_total_walkers() const { return total_walkers_; }
  IndexType get_requested_samples() const { return requested_samples_; }
  IndexType get_sub_steps() const { return sub_steps_; }
  RealType get_max_disp_sq() const { return max_disp_sq_; }
  IndexType get_max_blocks() const { return max_blocks_; }
  IndexType get_max_steps() const { return max_steps_; }
  IndexType get_warmup_steps() const { return warmup_steps_; }
  IndexType get_steps_between_samples() const { return steps_between_samples_; }
  IndexType get_samples_per_thread() const { return samples_per_thread_; }
  RealType get_tau() const { return tau_; }
  RealType get_spin_mass() const { return spin_mass_; }
  IndexType get_blocks_between_recompute() const { return blocks_between_recompute_; }
  bool get_append_run() const { return append_run_; }
  input::PeriodStride get_walker_dump_period() const { return walker_dump_period_; }
  input::PeriodStride get_check_point_period() const { return check_point_period_; }
  IndexType get_k_delay() const { return k_delay_; }
  bool get_reset_random() const { return reset_random_; }
  bool get_dump_config() const { return dump_config_; }

  const std::string& get_qmc_method() const { return qmc_method_; }
  const std::string& get_update_mode() const { return update_mode_; }
  DriverDebugChecks get_debug_checks() const { return debug_checks_; }
  bool get_scoped_profiling() const { return scoped_profiling_; }
  bool areWalkersSerialized() const { return crowd_serialize_walkers_; }
  bool get_measure_imbalance() const { return measure_imbalance_; }

  const std::string get_drift_modifier() const { return drift_modifier_; }
  RealType get_drift_modifier_unr_a() const { return drift_modifier_unr_a_; }

  const std::optional<EstimatorManagerInput>& get_estimator_manager_input() const { return estimator_manager_input_; }
};

// These will cause a compiler error if the implicit move constructor has been broken
inline QMCDriverInput::QMCDriverInput(QMCDriverInput&&) noexcept = default;
inline QMCDriverInput& QMCDriverInput::operator=(QMCDriverInput&&) noexcept = default;

} // namespace qmcplusplus

#endif
