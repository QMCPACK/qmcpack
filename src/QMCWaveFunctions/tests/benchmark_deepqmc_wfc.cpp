//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2026 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//////////////////////////////////////////////////////////////////////////////////////

#include <catch2/catch_test_macros.hpp>
#include <catch2/benchmark/catch_benchmark.hpp>

#include <cstdlib>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

#include "Configuration.h"
#include "Particle/ParticleSet.h"
#include "QMCWaveFunctions/DeepQMC/DeepQMCBridge.h"
#include "QMCWaveFunctions/DeepQMC/DeepQMCWF.h"
#include "QMCWaveFunctions/TrialWaveFunction.h"
#include "Utilities/RuntimeOptions.h"

namespace qmcplusplus
{
namespace
{
const char* getRequiredEnv(const char* name)
{
  const char* value = std::getenv(name);
  INFO(std::string("Set ") + name + " to run the DeepQMC benchmark");
  REQUIRE(value != nullptr);
  REQUIRE(value[0] != '\0');
  return value;
}

ParticleSet makeIons(const SimulationCell& simulation_cell)
{
  ParticleSet ions(simulation_cell);
  ions.setName("ion0");
  ions.create({1});
  ions.R[0] = {0.0, 0.0, 0.0};
  ions.update();
  return ions;
}

std::unique_ptr<ParticleSet> makeElectrons(const SimulationCell& simulation_cell, int walker_index)
{
  auto electrons = std::make_unique<ParticleSet>(simulation_cell);
  electrons->setName("e");
  electrons->create({1, 1});

  const DeepQMCBridge::RealType shift = 0.01 * static_cast<DeepQMCBridge::RealType>(walker_index % 17);
  electrons->R[0]                     = {0.25 + shift, -0.10 + 0.5 * shift, 0.05 - 0.25 * shift};
  electrons->R[1]                     = {-0.20 + 0.25 * shift, 0.15 - shift, -0.05 + 0.5 * shift};
  electrons->update();
  return electrons;
}

std::unique_ptr<TrialWaveFunction> makeTrialWaveFunctionWithDeepQMC(const RuntimeOptions& runtime_options,
                                                                    const ParticleSet& ions,
                                                                    std::unique_ptr<const DeepQMCBridge> bridge)
{
  auto twf = std::make_unique<TrialWaveFunction>(runtime_options, "deepqmc_benchmark");
  twf->addComponent(std::make_unique<DeepQMCWF>("DNN", ions, std::move(bridge), 0));
  return twf;
}

struct DeepQMCBenchmarkBatch
{
  DeepQMCBenchmarkBatch(int batch_size,
                        const RuntimeOptions& runtime_options,
                        const SimulationCell& simulation_cell,
                        const ParticleSet& ions,
                        const std::string& checkpoint_path)
  {
    electrons.reserve(batch_size);
    wavefunctions.reserve(batch_size);
    for (int iw = 0; iw < batch_size; ++iw)
    {
      electrons.push_back(makeElectrons(simulation_cell, iw));
      auto bridge = iw == 0
          ? makePythonDeepQMCBridge(checkpoint_path, "")
          : makeUnavailableDeepQMCBridge("Only the leader bridge is used by batched DeepQMCWF evaluation");
      wavefunctions.push_back(makeTrialWaveFunctionWithDeepQMC(runtime_options, ions, std::move(bridge)));
    }

    wf_refs = std::make_unique<RefVectorWithLeader<TrialWaveFunction>>(*wavefunctions.front());
    p_refs  = std::make_unique<RefVectorWithLeader<ParticleSet>>(*electrons.front());
    for (int iw = 0; iw < batch_size; ++iw)
    {
      wf_refs->push_back(*wavefunctions[iw]);
      p_refs->push_back(*electrons[iw]);
    }
  }

  void evaluate() { TrialWaveFunction::mw_evaluateLog(*wf_refs, *p_refs); }

  std::vector<std::unique_ptr<ParticleSet>> electrons;
  std::vector<std::unique_ptr<TrialWaveFunction>> wavefunctions;
  std::unique_ptr<RefVectorWithLeader<TrialWaveFunction>> wf_refs;
  std::unique_ptr<RefVectorWithLeader<ParticleSet>> p_refs;
};

std::string benchmarkName(int batch_size)
{
  std::ostringstream name;
  name << "TrialWaveFunction::mw_evaluateLog DeepQMC He batch_size=" << batch_size;
  return name.str();
}
} // namespace

TEST_CASE("DeepQMC TrialWaveFunction mw_evaluateLog batching benchmark", "[wavefunction][deepqmc][.benchmark]")
{
  const std::string checkpoint_path(getRequiredEnv("DEEPQMC_HE_CHECKPOINT"));

  // The benchmark executable is intentionally constructed directly from the DeepQMC
  // WaveFunctionComponent rather than through WaveFunctionFactory/XML. The external
  // Python environment still needs to be supplied by the launcher, e.g. PYTHONPATH
  // containing DeepQMC and its dependencies.
  const SimulationCell simulation_cell;
  ParticleSet ions = makeIons(simulation_cell);
  RuntimeOptions runtime_options;

  const std::vector<int> batch_sizes{1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048};
  for (int batch_size : batch_sizes)
  {
    DYNAMIC_SECTION("batch size " << batch_size)
    {
      DeepQMCBenchmarkBatch batch(batch_size, runtime_options, simulation_cell, ions, checkpoint_path);

      // Trigger Python/JAX setup and shape-specialized compilation outside the timed region.
      batch.evaluate();

      BENCHMARK_ADVANCED(benchmarkName(batch_size).c_str())(Catch::Benchmark::Chronometer meter)
      {
        meter.measure([&batch]() { batch.evaluate(); });
      };
    }
  }
}
} // namespace qmcplusplus
