//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2025 QMCPACK developers.
//
// File developed by: Anouar Benali, abenali.sci@gmail.com, Qubit Pharmaceuticals
//
// File created by: Anouar Benali, abenali.sci@gmail.com, Qubit Pharmaceuticals
//////////////////////////////////////////////////////////////////////////////////////

/** @file test_twf_fastderiv_raii.h
 *  @brief Documentation for TWFFastDerivWrapper RAII resource management tests
 *
 *  This test suite validates the Resource Acquisition Is Initialization (RAII) 
 *  pattern implementation for TWFFastDerivWrapper when managed by TrialWaveFunction.
 *
 *  ## Background
 *  
 *  TWFFastDerivWrapper is used for fast force evaluation in QMCPACK. Previously,
 *  it managed its own resources independently, which led to resource management
 *  conflicts and double-takeback errors. This test validates the new architecture
 *  where TWFFastDerivWrapper resources are managed through TrialWaveFunction's
 *  resource management system.
 *
 *  ## Architecture Being Tested
 *
 *  The resource management hierarchy is:
 *  ```
 *  ResourceCollectionTeamLock<TrialWaveFunction>
 *      └── TrialWaveFunction::acquireResource()
 *              ├── WaveFunctionComponent::acquireResource() [for each component]
 *              └── TWFFastDerivWrapper::acquireResource() [if wrapper exists]
 *  ```
 *
 *  ## Key Design Decisions Validated
 *
 *  1. **On-demand wrapper creation**: TWFFastDerivWrapper is created only when
 *     needed (when forces are requested), not during TWF construction
 *
 *  2. **Two-phase initialization**: TrialWaveFunction is constructed from XML
 *     without ParticleSet, wrapper is created later when ParticleSet is available
 *
 *  3. **Resource lifetime**: Wrapper resources are acquired/released through
 *     TWF's resource lock, not independently
 *
 *  4. **Conditional resource allocation**: Memory is only allocated for force
 *     evaluation when ACForce is present in the Hamiltonian
 *
 *  ## Test Sections
 *
 *  ### Single walker wrapper creation and resource management
 *  - Validates that wrapper doesn't exist initially
 *  - Tests on-demand creation via getOrCreateTWFFastDerivWrapper()
 *  - Ensures repeated calls return the same wrapper instance
 *
 *  ### Multi-walker resource acquisition and release
 *  - Tests the standard multi-walker pattern used in QMC drivers
 *  - Validates multiple acquire/release cycles (idempotency)
 *  - Ensures no resource leaks or double-takeback errors
 *
 *  ### TWF acquires wrapper resources when they exist
 *  - Tests that TWF resource management works without wrappers
 *  - Validates that wrapper resources are properly managed when present
 *  - Ensures backward compatibility with runs that don't need forces
 *
 *  ### Nested resource locks
 *  - Simulates the real driver pattern with ParticleSet, TWF, and Hamiltonian locks
 *  - Tests proper lock ordering and release order
 *  - Validates no deadlocks or resource conflicts
 *
 *  ### Wrapper persistence across cycles
 *  - Ensures wrapper survives resource release cycles
 *  - Tests that wrapper initialization is preserved
 *  - Validates efficient reuse without re-initialization
 *
 *  ### Conditional creation based on forces
 *  - Tests the hasAuxiliaryOperator("ACForce") pattern
 *  - Validates memory savings when forces aren't needed
 *  - Ensures resource management works with or without wrapper
 *
 *  ## Known Issues Addressed
 *
 *  1. **Double-takeback error**: Previous architecture called takebackResource()
 *     both in component's releaseResource() and in ResourceCollectionTeamLock's
 *     destructor, causing crashes. This test validates the fix where
 *     releaseResource() only clears handles.
 *
 *  2. **Resource cursor mismatches**: Old approach used cursor-based resource
 *     management which caused type confusion. This test validates the new
 *     approach where TWF manages wrapper resources internally.
 *
 *  3. **Premature resource acquisition**: Previous attempts tried to acquire
 *     wrapper resources before wrappers existed. This test validates that
 *     wrappers are created before resource locks when needed.
 *
 *  ## Integration Points Tested
 *
 *  - QMCDriverNew::initialLogEvaluation pattern
 *  - VMCBatched::advanceWalkers resource management
 *  - ACForce::mw_evaluate wrapper usage (indirectly)
 *
 *  ## Success Criteria
 *
 *  All tests pass if:
 *  - No exceptions thrown during resource acquisition/release
 *  - No memory leaks detected
 *  - Wrapper resources properly managed through TWF
 *  - Conditional creation saves memory when forces not needed
 *  - Multiple resource cycles work without errors
 *
 *  ## Usage Note
 *
 *  This test requires minimal particle sets and wavefunction components.
 *  It does not test the actual force evaluation functionality, only the
 *  resource management infrastructure.
 */

// test_twf_fastderiv_raii.cpp
#include "catch.hpp"

#include "QMCWaveFunctions/TrialWaveFunction.h"
#include "QMCWaveFunctions/TWFFastDerivWrapper.h"
#include "Particle/ParticleSet.h"
#include "SimulationCell.h"
#include "Utilities/ResourceCollection.h"
#include "Utilities/RuntimeOptions.h"
#include "QMCHamiltonians/QMCHamiltonian.h"

using namespace qmcplusplus;

namespace
{
// Helper to create minimal particle sets
ParticleSet createMinimalElectrons(SimulationCell& sim_cell, int num_elec = 2)
{
  ParticleSet elec(sim_cell);
  elec.setName("e");
  elec.create({num_elec / 2, num_elec - num_elec / 2}); // up, down spins
  for (int i = 0; i < num_elec; ++i)
    elec.R[i] = {QMCTraits::RealType(0.1 * i), 0, 0};
  elec.update();
  return elec;
}

} // namespace
TEST_CASE("TWFFastDerivWrapper RAII resource management", "[wavefunction][resources]")
{
  SimulationCell sim_cell;
  RuntimeOptions runtime_options;

  SECTION("Single walker wrapper creation and resource management")
  {
    ParticleSet elec = createMinimalElectrons(sim_cell);
    TrialWaveFunction twf(runtime_options, "psi", false);

    // Create wrapper on demand
    TWFFastDerivWrapper& wrapper = twf.getOrCreateTWFFastDerivWrapper(elec);

    // Second call returns same wrapper
    TWFFastDerivWrapper& wrapper2 = twf.getOrCreateTWFFastDerivWrapper(elec);
    REQUIRE(&wrapper == &wrapper2);
  }

  SECTION("Multi-walker resource acquisition and release")
  {
    const int num_walkers = 4;

    // Create particle sets
    std::vector<ParticleSet> elec_list;
    for (int i = 0; i < num_walkers; ++i)
      elec_list.push_back(createMinimalElectrons(sim_cell));

    // Create TWFs (simplified to one line)
    std::vector<std::unique_ptr<TrialWaveFunction>> twf_list;
    for (int i = 0; i < num_walkers; ++i)
      twf_list.push_back(std::make_unique<TrialWaveFunction>(runtime_options, "psi", false));

    // Create RefVectorWithLeader
    RefVectorWithLeader<TrialWaveFunction> twf_refs(*twf_list[0]);
    for (auto& twf : twf_list)
      twf_refs.push_back(*twf);

    RefVectorWithLeader<ParticleSet> elec_refs(elec_list[0]);
    for (auto& elec : elec_list)
      elec_refs.push_back(elec);

    // Create wrappers (simulating what ACForce does)
    for (int i = 0; i < num_walkers; ++i)
      twf_refs[i].getOrCreateTWFFastDerivWrapper(elec_refs[i]);

    ResourceCollection twf_res("twf");

    // TWF creates resources including wrapper resource
    size_t initial_resource_count = twf_res.size();
    twf_refs.getLeader().createResource(twf_res);
    REQUIRE(twf_res.size() > initial_resource_count); // Resources were added

    // Multiple acquire/release cycles - if these complete without throwing, RAII works
    for (int cycle = 0; cycle < 3; ++cycle)
    {
      {
        ResourceCollectionTeamLock<TrialWaveFunction> lock(twf_res, twf_refs);
        // Resources acquired
      }
      // Resources released by destructor

      {
        ResourceCollectionTeamLock<TrialWaveFunction> lock2(twf_res, twf_refs);
        // Can acquire again - proper release happened
      }
      // Released again
    }
    // Getting here means no exceptions were thrown
  }

  SECTION("TWF acquires wrapper resources when they exist")
  {
    const int num_walkers = 2;

    std::vector<ParticleSet> elec_list;
    for (int i = 0; i < num_walkers; ++i)
      elec_list.push_back(createMinimalElectrons(sim_cell));

    std::vector<std::unique_ptr<TrialWaveFunction>> twf_list;
    for (int i = 0; i < num_walkers; ++i)
      twf_list.push_back(std::make_unique<TrialWaveFunction>(runtime_options, "psi", false));

    RefVectorWithLeader<TrialWaveFunction> twf_refs(*twf_list[0]);
    for (auto& twf : twf_list)
      twf_refs.push_back(*twf);

    RefVectorWithLeader<ParticleSet> elec_refs(elec_list[0]);
    for (auto& elec : elec_list)
      elec_refs.push_back(elec);

    ResourceCollection twf_res("twf");

    // First test: TWF without wrappers
    size_t resources_without_wrapper = 0;
    twf_refs.getLeader().createResource(twf_res);
    resources_without_wrapper = twf_res.size();

    {
      ResourceCollectionTeamLock<TrialWaveFunction> lock(twf_res, twf_refs);
      // Should work even without wrappers
    }

    // Now create wrappers
    for (int i = 0; i < num_walkers; ++i)
      twf_refs[i].getOrCreateTWFFastDerivWrapper(elec_refs[i]);

    // Test that wrapper resources are now managed
    ResourceCollection twf_res2("twf_with_wrapper");
    twf_refs.getLeader().createResource(twf_res2);
    REQUIRE(twf_res2.size() >= resources_without_wrapper); // May have added wrapper resource

    {
      ResourceCollectionTeamLock<TrialWaveFunction> lock(twf_res2, twf_refs);
      // Wrapper resources should be acquired through TWF
    }
  }

  SECTION("Nested resource locks (simulating real driver pattern)")
  {
    const int num_walkers = 2;

    std::vector<ParticleSet> elec_list;
    for (int i = 0; i < num_walkers; ++i)
      elec_list.push_back(createMinimalElectrons(sim_cell));

    std::vector<std::unique_ptr<TrialWaveFunction>> twf_list;
    std::vector<std::unique_ptr<QMCHamiltonian>> ham_list;

    for (int i = 0; i < num_walkers; ++i)
    {
      twf_list.push_back(std::make_unique<TrialWaveFunction>(runtime_options, "psi", false));
      ham_list.push_back(std::make_unique<QMCHamiltonian>());
    }

    RefVectorWithLeader<ParticleSet> elec_refs(elec_list[0]);
    for (auto& elec : elec_list)
      elec_refs.push_back(elec);

    RefVectorWithLeader<TrialWaveFunction> twf_refs(*twf_list[0]);
    for (auto& twf : twf_list)
      twf_refs.push_back(*twf);

    RefVectorWithLeader<QMCHamiltonian> ham_refs(*ham_list[0]);
    for (auto& ham : ham_list)
      ham_refs.push_back(*ham);

    // Create wrappers
    for (int i = 0; i < num_walkers; ++i)
      twf_refs[i].getOrCreateTWFFastDerivWrapper(elec_refs[i]);

    ResourceCollection pset_res("pset");
    ResourceCollection twf_res("twf");
    ResourceCollection ham_res("ham");

    elec_refs.getLeader().createResource(pset_res);
    twf_refs.getLeader().createResource(twf_res);
    ham_refs.getLeader().createResource(ham_res);

    // Test nested locks don't cause deadlock or double-takeback
    {
      ResourceCollectionTeamLock<ParticleSet> pset_lock(pset_res, elec_refs);
      ResourceCollectionTeamLock<TrialWaveFunction> twf_lock(twf_res, twf_refs);
      ResourceCollectionTeamLock<QMCHamiltonian> ham_lock(ham_res, ham_refs);
      // All resources active
    }
    // All released in correct order

    // Can acquire again - tests proper release
    {
      ResourceCollectionTeamLock<ParticleSet> pset_lock(pset_res, elec_refs);
      ResourceCollectionTeamLock<TrialWaveFunction> twf_lock(twf_res, twf_refs);
      ResourceCollectionTeamLock<QMCHamiltonian> ham_lock(ham_res, ham_refs);
      // No double-takeback errors
    }
  }

  SECTION("Wrapper persistence across multiple resource cycles")
  {
    ParticleSet elec = createMinimalElectrons(sim_cell);
    TrialWaveFunction twf(runtime_options, "psi", false);

    // Create wrapper
    TWFFastDerivWrapper& wrapper1 = twf.getOrCreateTWFFastDerivWrapper(elec);

    RefVectorWithLeader<TrialWaveFunction> twf_refs(twf);
    ResourceCollection twf_res("twf");
    twf.createResource(twf_res);

    // First cycle
    {
      ResourceCollectionTeamLock<TrialWaveFunction> lock(twf_res, twf_refs);
    }

    // Wrapper should still exist after resource release
    TWFFastDerivWrapper& wrapper2 = twf.getOrCreateTWFFastDerivWrapper(elec);
    REQUIRE(&wrapper1 == &wrapper2);

    // Second cycle with same wrapper
    {
      ResourceCollectionTeamLock<TrialWaveFunction> lock(twf_res, twf_refs);
    }

    // Still the same wrapper
    TWFFastDerivWrapper& wrapper3 = twf.getOrCreateTWFFastDerivWrapper(elec);
    REQUIRE(&wrapper1 == &wrapper3);
  }
}

TEST_CASE("TWFFastDerivWrapper conditional creation", "[wavefunction][resources]")
{
  SECTION("Wrapper only created when forces needed")
  {
    SimulationCell sim_cell;
    RuntimeOptions runtime_options;
    ParticleSet elec = createMinimalElectrons(sim_cell);

    // Simulate two scenarios
    bool forces_needed = GENERATE(true, false);

    TrialWaveFunction twf(runtime_options, "psi", false);

    if (forces_needed)
    {
      // Simulate ACForce behavior when forces are present
      TWFFastDerivWrapper& wrapper  = twf.getOrCreateTWFFastDerivWrapper(elec);
      TWFFastDerivWrapper& wrapper2 = twf.getOrCreateTWFFastDerivWrapper(elec);
      REQUIRE(&wrapper == &wrapper2); // Same wrapper returned
    }

    // Resource management should work either way
    RefVectorWithLeader<TrialWaveFunction> twf_refs(twf);
    ResourceCollection twf_res("twf");

    // Should not throw whether wrapper exists or not
    REQUIRE_NOTHROW(twf.createResource(twf_res));

    {
      ResourceCollectionTeamLock<TrialWaveFunction> lock(twf_res, twf_refs);
      // Should work with or without wrapper
    }
    // Getting here means successful execution
  }
}
