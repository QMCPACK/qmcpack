//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2026 QMCPACK developers.
//
//////////////////////////////////////////////////////////////////////////////////////

#include "QMCDrivers/WaveFunctionTesterBatched.h"

#include "Concurrency/ParallelExecutor.hpp"
#include "EstimatorInputDelegates.h"
#include "Message/CommOperators.h"
#include "Message/UniformCommunicateError.h"

#include <array>
#include <fstream>
#include <iomanip>

namespace qmcplusplus
{
WaveFunctionTesterBatched::WaveFunctionTesterBatched(const ProjectData& project_data,
                                                     QMCDriverInput&& qmcdriver_input,
                                                     UPtr<EstimatorManagerNew>&& estimator_manager,
                                                     WaveFunctionTesterBatchedInput&& input,
                                                     WalkerConfigurations& wc,
                                                     MCPopulation&& population,
                                                     const RefVector<RandomBase<FullPrecRealType>>& rng_refs,
                                                     Communicate* comm)
    : QMCDriverNew(project_data,
                   std::move(qmcdriver_input),
                   std::move(estimator_manager),
                   wc,
                   std::move(population),
                   rng_refs,
                   "WaveFunctionTesterBatched::",
                   comm,
                   "WaveFunctionTesterBatched"),
      input_(std::move(input))
{}

void WaveFunctionTesterBatched::process(xmlNodePtr)
{
  const int num_crowds = determineNumCrowds(qmcdriver_input_.get_num_crowds(), rngs_.size());
  const AdjustedWalkerCounts walker_counts =
      adjustGlobalWalkerCount(*myComm, walker_configs_ref_.getActiveWalkers(), qmcdriver_input_.get_total_walkers(),
                              qmcdriver_input_.get_walkers_per_rank(), 1.0, num_crowds);

  steps_per_block_ = 1;
  initPopulationAndCrowds(walker_counts);

  test_contexts_.resize(crowds_.size());
  for (int crowd_index = 0; crowd_index < crowds_.size(); ++crowd_index)
    test_contexts_[crowd_index] = std::make_unique<TestContext>(rngs_[crowd_index]);
}

bool WaveFunctionTesterBatched::run()
{
  RefVector<ContextForSteps> context_refs;
  context_refs.reserve(test_contexts_.size());
  for (auto& context : test_contexts_)
    context_refs.push_back(*context);

  ParallelExecutor<> initialize_task;
  initialize_task(crowds_.size(), initialLogEvaluation, crowds_, context_refs, qmcdriver_input_.areWalkersSerialized());

  std::array<char, 16> filename;
  if (std::snprintf(filename.data(), filename.size(), "wftest.%03d", myComm->rank()) < 0)
    throw UniformCommunicateError("Failed to generate batched wavefunction tester output filename.");

  std::ofstream output(filename.data());
  output.precision(15);
  output << "Numerical gradient and Laplacian test\n";

  bool all_okay                         = true;
  const RealType delta                  = input_.get_delta() > 0.0 ? input_.get_delta() : 1.0e-4;
  const RealType tolerance              = input_.get_tolerance() > 0.0 ? input_.get_tolerance() : 1.0e-3;
  constexpr RealType absolute_tolerance = 1.0e-7;

  for (int crowd_index = 0; crowd_index < crowds_.size(); ++crowd_index)
  {
    Crowd& crowd = *crowds_[crowd_index];
    for (int walker_index = 0; walker_index < crowd.size(); ++walker_index)
    {
      ParticleSet& electrons          = crowd.get_walker_elecs()[walker_index];
      TrialWaveFunction& wavefunction = crowd.get_walker_twfs()[walker_index];
      const int particle_count        = electrons.getTotalNum();
      ParticleSet::ParticleGradient analytic_gradient(electrons.G);
      ParticleSet::ParticleLaplacian analytic_laplacian(electrons.L);

      output << "Walker # " << walker_index << '\n';
      for (int particle_index = 0; particle_index < particle_count; ++particle_index)
      {
        ParticleSet::GradType finite_gradient;
        ParticleSet::ValueType finite_laplacian      = 0.0;
        const ParticleSet::PosType original_position = electrons.R[particle_index];
        const RealType log_center                    = wavefunction.evaluateLog(electrons);

        for (int dimension = 0; dimension < OHMMS_DIM; ++dimension)
        {
          ParticleSet::PosType minus_position(original_position);
          minus_position[dimension] -= delta;
          electrons.R[particle_index] = minus_position;
          electrons.update();
          const RealType log_minus = wavefunction.evaluateLog(electrons);

          ParticleSet::PosType plus_position(original_position);
          plus_position[dimension] += delta;
          electrons.R[particle_index] = plus_position;
          electrons.update();
          const RealType log_plus = wavefunction.evaluateLog(electrons);

          finite_gradient[dimension] = (log_plus - log_minus) / (2.0 * delta);
          finite_laplacian += (log_plus + log_minus - 2.0 * log_center) / (delta * delta);
        }

        electrons.R[particle_index] = original_position;
        electrons.update();
        wavefunction.evaluateLog(electrons);

        const auto analytic_lap       = analytic_laplacian[particle_index];
        const auto lap_error          = std::abs(analytic_lap - finite_laplacian);
        const auto lap_denom          = std::max(std::abs(analytic_lap), std::abs(finite_laplacian));
        const auto lap_relative_error = lap_denom > 0.0 ? lap_error / lap_denom : 0.0;

        std::array<RealType, OHMMS_DIM> gradient_relative_errors;
        for (int dimension = 0; dimension < OHMMS_DIM; ++dimension)
        {
          const auto gradient_error =
              std::abs(analytic_gradient[particle_index][dimension] - finite_gradient[dimension]);
          const auto gradient_denom =
              std::max(std::abs(analytic_gradient[particle_index][dimension]), std::abs(finite_gradient[dimension]));
          gradient_relative_errors[dimension] = gradient_denom > 0.0 ? gradient_error / gradient_denom : 0.0;
          all_okay =
              all_okay && !(gradient_relative_errors[dimension] > tolerance && gradient_error > absolute_tolerance);
        }
        all_okay = all_okay && !(lap_relative_error > tolerance && lap_error > absolute_tolerance);

        output << "For particle #" << particle_index << " at " << original_position << '\n';
        output << "Gradient      = " << analytic_gradient[particle_index] << '\n';
        output << "  Finite diff = " << finite_gradient << '\n';
        output << "  Error       = " << analytic_gradient[particle_index] - finite_gradient << '\n';
        output << "  Relative Error = ";
        for (const RealType relative_error : gradient_relative_errors)
          output << relative_error << ' ';
        output << "\n\n";
        output << "Laplacian     = " << analytic_lap << '\n';
        output << "  Finite diff = " << finite_laplacian << '\n';
        output << "  Error       = " << analytic_lap - finite_laplacian << "  Relative Error = " << lap_relative_error
               << "\n\n";
      }

      output << "Particle       Ratio of Ratios     Computed Ratio   Internal Ratio\n";
      for (int particle_index = 0; particle_index < particle_count; ++particle_index)
      {
        ParticleSet::PosType displacement;
        displacement               = 1.0e-3;
        const RealType initial_log = wavefunction.evaluateLog(electrons);
        electrons.makeMove(particle_index, displacement);
        const TrialWaveFunction::PsiValue internal_ratio = wavefunction.calcRatio(electrons, particle_index);
        electrons.rejectMove(particle_index);
        wavefunction.rejectMove(particle_index);

        const ParticleSet::PosType original_position = electrons.R[particle_index];
        electrons.R[particle_index] += displacement;
        electrons.update();
        const RealType displaced_log = wavefunction.evaluateLog(electrons);
        electrons.R[particle_index]  = original_position;
        electrons.update();
        wavefunction.evaluateLog(electrons);

        const RealType computed_ratio  = std::exp(displaced_log - initial_log);
        const RealType ratio_of_ratios = std::abs(internal_ratio) / computed_ratio;
        output << particle_index << ' ' << ratio_of_ratios << ' ' << computed_ratio << ' ' << internal_ratio << '\n';
        all_okay = all_okay && std::abs(ratio_of_ratios - 1.0) <= 1.0e-6;
      }
    }
  }

  int global_pass = all_okay ? 1 : 0;
  myComm->allreduce(global_pass);
  const bool passed = global_pass == myComm->size();
  app_log() << "Finite difference test: " << (passed ? "PASS" : "FAIL") << std::endl;
  app_log() << "Ratio test: " << (passed ? "PASS" : "FAIL") << std::endl;
  return passed;
}
} // namespace qmcplusplus