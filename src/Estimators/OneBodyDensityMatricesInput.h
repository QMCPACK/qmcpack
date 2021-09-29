//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2021 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//
// Some code refactored from: DensityMatrices1b.h
//////////////////////////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_ONE_BODY_DENSITY_MATRICES_INPUT_H
#define QMCPLUSPLUS_ONE_BODY_DENSITY_MATRICES_INPUT_H

#include "Configuration.h"
#include "InputSection.h"

namespace qmcplusplus
{
/** Native representation for DensityMatrices1B Estimator's inputs
 */
class OneBodyDensityMatricesInput
{
public:
  enum class Integrators
  {
    UNIFORM_GRID,
    UNIFORM,
    DENSITY,
    NO_INTEGRATOR
  };

  enum class Evaluators
  {
    LOOP,
    MATRIX,
    NO_EVALUATOR
  };

  enum class Samplings
  {
    VOLUME_BASED,
    METROPOLIS,
    NO_SAMPLING
  };

  class OneBodyDensityMatrixInputSection : public InputSection
  {
  public:
    /** parse time definition of input parameters */
    OneBodyDensityMatrixInputSection()
    {
      section_name = "OneBodyDensityMatrix";
      attributes   = {"name", "type"};
      parameters   = {"basis",      "energy_matrix", "integrator",        "evaluator",        "scale",
                    "center",     "points",        "samples",           "warmup",           "timestep",
                    "use_drift",  "check_overlap", "check_derivatives", "acceptance_ratio", "rstats",
                    "normalized", "volumed_normed"};
      bools        = {"energy_matrix", "use_drift",         "normalized", "volume_normed",
               "check_overlap", "check_derivatives", "rstats",     "acceptance_ratio"};
      strings      = {"name", "type", "basis", "integrator", "evaluator"};
      integers     = {"points", "samples"};
      reals        = {"scale", "timestep"};
      required     = {"name", "basis"};
      // I'd much rather see the default defined in simple native c++ as below
    }

    /** do parse time checks of input */
    void checkParticularValidity() override;
  };

  using Position = QMCTraits::PosType;
  using Real     = QMCTraits::FullPrecRealType;

  OneBodyDensityMatricesInput() = default;
  OneBodyDensityMatricesInput(xmlNodePtr cur);

private:
  OneBodyDensityMatrixInputSection input_section_;

  // Default parameters for OneBodyDensityMatrices
  bool energy_matrix_     = false;
  bool use_drift_         = false;
  bool normalized_        = true;
  bool volume_normalized_ = true;
  bool check_overlap_     = false;
  bool check_derivatives_ = false;
  bool rstats_            = false;
  bool acceptance_ratio_  = false;
  Integrators integrator_ = Integrators::UNIFORM_GRID;
  Samplings sampling_     = Samplings::VOLUME_BASED;
  Evaluators evaluator_   = Evaluators::LOOP;
  Real scale_             = 1.0;
  Position center_        = 0.0;
  Real timestep_          = 0.5;
  int points_             = 10;
  int samples_            = 10;
  int warmup_samples_     = 30;
public:
  bool get_energy_matrix() const { return energy_matrix_; }
  bool get_use_drift() const { return use_drift_; }
  bool get_normalized() const { return normalized_; }
  bool get_volume_normalized() const { return volume_normalized_; }
  bool get_check_overlap() const { return check_overlap_; }
  bool get_check_derivatives() const { return check_derivatives_; }
  bool get_rstats() const { return rstats_; }
  bool get_acceptance_ratio() const { return acceptance_ratio_; }
  Integrators get_integrator() const { return integrator_; }
  Samplings get_sampling() const { return sampling_; }
  Evaluators get_evaluator() const { return evaluator_; }
  Real get_scale() const { return scale_; }
  Position get_center() const { return center_; }
  Real get_timestep() const { return timestep_; }
  int get_points() const { return points_; }
  int get_samples() const { return samples_; }
  int get_warmup_samples() const { return warmup_samples_; }
};

} // namespace qmcplusplus

#endif
