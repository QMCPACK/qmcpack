//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//////////////////////////////////////////////////////////////////////////////////////
#include "QMCDrivers/tests/CountingWaveFunctionComponent.h"

#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <memory>

namespace qmcplusplus
{
namespace
{
struct CountingCrowd
{
  explicit CountingCrowd(std::size_t size)
  {
    for (std::size_t i = 0; i < size; ++i)
    {
      particles.push_back(std::make_unique<ParticleSet>(simulation_cell));
      components.push_back(std::make_unique<testing::CountingWaveFunctionComponent>());
    }
    for (auto& component : components)
      wfc_list.push_back(*component);
    for (auto& particle : particles)
    {
      p_list.push_back(*particle);
      G_list.push_back(particle->G);
      L_list.push_back(particle->L);
    }
    component_list = std::make_unique<RefVectorWithLeader<WaveFunctionComponent>>(*components.front(), wfc_list);
    particle_list  = std::make_unique<RefVectorWithLeader<ParticleSet>>(*particles.front(), p_list);
  }

  const SimulationCell simulation_cell;
  std::vector<std::unique_ptr<ParticleSet>> particles;
  std::vector<std::unique_ptr<testing::CountingWaveFunctionComponent>> components;
  RefVectorWithLeader<WaveFunctionComponent>::BaseVec wfc_list;
  RefVectorWithLeader<ParticleSet>::BaseVec p_list;
  RefVector<ParticleSet::ParticleGradient> G_list;
  RefVector<ParticleSet::ParticleLaplacian> L_list;
  std::unique_ptr<RefVectorWithLeader<WaveFunctionComponent>> component_list;
  std::unique_ptr<RefVectorWithLeader<ParticleSet>> particle_list;
};

TEST_CASE("CountingWaveFunctionComponent clones have independent state", "[drivers]")
{
  testing::CountingWaveFunctionComponent original;
  const SimulationCell simulation_cell;
  ParticleSet particles(simulation_cell);
  ParticleSet::ParticleGradient gradient;
  ParticleSet::ParticleLaplacian laplacian;
  original.evaluateLog(particles, gradient, laplacian);
  auto clone          = original.makeClone(particles);
  auto& counted_clone = static_cast<testing::CountingWaveFunctionComponent&>(*clone);
  REQUIRE(original.getEvaluationCount() == 1);
  REQUIRE(counted_clone.getEvaluationCount() == 0);
  counted_clone.evaluateLog(particles, gradient, laplacian);
  REQUIRE(original.getEvaluationCount() == 1);
  REQUIRE(counted_clone.getEvaluationCount() == 1);
  REQUIRE(original.getProposals().empty());
  REQUIRE(original.getResolutions().empty());
  REQUIRE(counted_clone.getProposals().empty());
  REQUIRE(counted_clone.getResolutions().empty());
}

TEST_CASE("CountingWaveFunctionComponent records all-particle masks and resolutions", "[drivers]")
{
  for (const std::size_t size : {1, 2, 3, 5})
  {
    std::vector<std::vector<bool>> masks(4, std::vector<bool>(size, false));
    for (std::size_t i = 0; i < size; ++i)
    {
      masks[0][i] = true;
      masks[2][i] = i % 2 == 0;
      masks[3][i] = i % 2 != 0;
    }
    CAPTURE(size);
    for (const auto& mask : masks)
    {
      CountingCrowd crowd(size);
      const auto& leader = *crowd.components.front();
      leader.mw_evaluateLogAllParticles(*crowd.component_list, *crowd.particle_list, crowd.G_list, crowd.L_list, mask);
      for (std::size_t i = 0; i < size; ++i)
      {
        const auto& proposals = crowd.components[i]->getProposals();
        REQUIRE(proposals.size() == 1);
        CHECK(proposals[0].evaluate == mask[i]);
        CHECK(proposals[0].evaluations == static_cast<int>(mask[i]));
        CHECK(crowd.components[i]->getEvaluationCount() == static_cast<int>(mask[i]));
      }
    }

    for (const auto& accepted : masks)
    {
      const bool mixed =
          std::any_of(accepted.begin(), accepted.end(), [&](bool value) { return value != accepted.front(); });
      CountingCrowd crowd(size);
      const auto& leader = *crowd.components.front();
      leader.mw_evaluateLogAllParticles(*crowd.component_list, *crowd.particle_list, crowd.G_list, crowd.L_list,
                                        std::vector<bool>(size, true));
      leader.mw_accept_rejectMoveAllParticles(*crowd.component_list, *crowd.particle_list, crowd.G_list, crowd.L_list,
                                              accepted);
      for (std::size_t i = 0; i < size; ++i)
      {
        const auto& resolutions = crowd.components[i]->getResolutions();
        REQUIRE(resolutions.size() == 1);
        CHECK(resolutions[0].accepted == accepted[i]);
        CHECK(resolutions[0].mixed == mixed);
        CHECK(resolutions[0].evaluations == static_cast<int>(!accepted[i]));
        CHECK(crowd.components[i]->getEvaluationCount() == 1 + static_cast<int>(!accepted[i]));
      }
    }
  }
}
} // namespace
} // namespace qmcplusplus
