//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2022 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//////////////////////////////////////////////////////////////////////////////////////


#include "catch.hpp"

#include "type_traits/template_types.hpp"
#include "QMCHamiltonians/QMCHamiltonian.h"
#include "Particle/tests/MinimalParticlePool.h"
#include "QMCWaveFunctions/tests/MinimalWaveFunctionPool.h"
#include "QMCHamiltonians/tests/MinimalHamiltonianPool.h"
#include "TestListenerFunction.h"
#include "Utilities/ResourceCollection.h"
#include "Utilities/StlPrettyPrint.hpp"
#include "Utilities/RuntimeOptions.h"
#include "Utilities/for_testing/NativeInitializerPrint.hpp"

namespace qmcplusplus
{
using QMCT = QMCTraits;
using Real = QMCT::RealType;

constexpr bool generate_test_data = false;

TEST_CASE("QMCHamiltonian::flex_evaluate", "[hamiltonian]")
{
  RuntimeOptions runtime_options;
  Communicate* comm;
  comm = OHMMS::Controller;

  auto particle_pool     = MinimalParticlePool::make_diamondC_1x1x1(comm);
  auto wavefunction_pool = MinimalWaveFunctionPool::make_diamondC_1x1x1(runtime_options, comm, particle_pool);
  auto hamiltonian_pool  = MinimalHamiltonianPool::make_hamWithEE(comm, particle_pool, wavefunction_pool);

  TrialWaveFunction twf(runtime_options);

  std::vector<QMCHamiltonian> hamiltonians;
  //hamiltonians.emplace_back(*(hamiltonian_pool.getPrimary()));
  //hamiltonians.emplace_back(*(hamiltonian_pool.getPrimary()));

  std::vector<ParticleSet> elecs;
  elecs.emplace_back(*(particle_pool.getParticleSet("e")));
  elecs.emplace_back(*(particle_pool.getParticleSet("e")));

  // TODO: finish initializing the elecs.
  //std::vector<QMCHamiltonian::RealType> local_energies(QMCHamiltonian::flex_evaluate(makeRefVector<QMCHamiltonian>(hamiltonians), makeRefVector<ParticleSet>(elecs)));

  //TODO: Would be nice to check some values but I think the system needs a little more setup
}

/** QMCHamiltonian + Hamiltonians with listeners integration test
 */
TEST_CASE("integrateListeners", "[hamiltonian]")
{
  RuntimeOptions runtime_options;
  Communicate* comm = OHMMS::Controller;

  auto particle_pool     = MinimalParticlePool::make_diamondC_1x1x1(comm);
  auto wavefunction_pool = MinimalWaveFunctionPool::make_diamondC_1x1x1(runtime_options, comm, particle_pool);
  auto hamiltonian_pool  = MinimalHamiltonianPool::makeHamWithEEEI(comm, particle_pool, wavefunction_pool);
  auto& pset_target      = *(particle_pool.getParticleSet("e"));
  //auto& species_set        = pset_target.getSpeciesSet();
  //auto& spo_map            = wavefunction_pool.getWaveFunction("wavefunction")->getSPOMap();
  auto& trial_wavefunction = *(wavefunction_pool.getPrimary());

  UPtrVector<QMCHamiltonian> hams;
  UPtrVector<TrialWaveFunction> twfs;
  std::vector<ParticleSet> psets;

  // This must be done before clones otherwise the clone particle sets do not have the correct state.
  hamiltonian_pool.getPrimary()->informOperatorsOfListener();


  int num_walkers   = 4;
  int num_electrons = particle_pool.getParticleSet("e")->getTotalNum();
  int num_ions      = particle_pool.getParticleSet("ion")->getTotalNum();

  for (int iw = 0; iw < num_walkers; ++iw)
  {
    psets.emplace_back(pset_target);
    psets.back().randomizeFromSource(*particle_pool.getParticleSet("ion"));
    twfs.emplace_back(trial_wavefunction.makeClone(psets.back()));
    hams.emplace_back(hamiltonian_pool.getPrimary()->makeClone(psets.back(), *twfs.back()));
  }

  RefVector<QMCHamiltonian> ham_refs = convertUPtrToRefVector(hams);
  RefVectorWithLeader<QMCHamiltonian> ham_list{ham_refs[0], ham_refs};

  ResourceCollection ham_res("test_ham_res");
  ham_list.getLeader().createResource(ham_res);
  ResourceCollectionTeamLock<QMCHamiltonian> ham_lock(ham_res, ham_list);

  Matrix<Real> kinetic(num_walkers, num_electrons);
  Matrix<Real> local_pots(num_walkers, num_electrons);
  Matrix<Real> local_nrg(num_walkers, num_electrons);
  Matrix<Real> ion_pots(num_walkers, num_ions);

  using testing::getParticularListener;
  using testing::getSummingListener;

  ListenerVector<Real> listener_kinetic("kinetic", getSummingListener(kinetic));
  QMCHamiltonian::mw_registerKineticListener(ham_list.getLeader(), listener_kinetic);
  ListenerVector<Real> listener_potential("potential", getSummingListener(local_pots));
  QMCHamiltonian::mw_registerLocalPotentialListener(ham_list.getLeader(), listener_potential);
  ListenerVector<Real> listener_energy("local_energy", getSummingListener(local_nrg));
  QMCHamiltonian::mw_registerLocalEnergyListener(ham_list.getLeader(), listener_energy);
  ListenerVector<Real> listener_ion_potential("ion_potential", getSummingListener(ion_pots));
  QMCHamiltonian::mw_registerLocalIonPotentialListener(ham_list.getLeader(), listener_ion_potential);

  auto p_refs = makeRefVector<ParticleSet>(psets);
  RefVectorWithLeader<ParticleSet> p_list{p_refs[0], p_refs};

  ResourceCollection pset_res("test_pset_res");
  p_list.getLeader().createResource(pset_res);
  ResourceCollectionTeamLock<ParticleSet> pset_lock(pset_res, p_list);

  auto twf_refs = convertUPtrToRefVector(twfs);
  RefVectorWithLeader<TrialWaveFunction> twf_list{twf_refs[0], twf_refs};

  ResourceCollection wfc_res("test_wfc_res");
  twf_list.getLeader().createResource(wfc_res);
  ResourceCollectionTeamLock<TrialWaveFunction> mw_wfc_lock(wfc_res, twf_list);

  // Otherwise we'll get no kinetic energy values
  p_refs[0].get().L[0] = 1.0;
  p_refs[1].get().L[1] = 1.0;
  p_refs[2].get().L[2] = 1.0;
  p_refs[3].get().L[3] = 1.0;

  if constexpr (generate_test_data)
  {
    std::cout << "QMCHamiltonian-registerListeners initialize psets with:\n{";

    for (int iw = 0; iw < num_walkers; ++iw)
    {
      //psets.emplace_back(pset_target);
      std::cout << "{";
      for (auto r : p_refs[iw].get().R)
        std::cout << NativePrint(r) << ",";
      std::cout << "},\n";
    }
    std::cout << "}\n";
  }
  else
  {
    std::vector<ParticleSet::ParticlePos> deterministic_rs = {
        {
            {0.515677886, 0.9486072745, -1.17423246},
            {-0.3166678423, 1.376550506, 1.445290031},
            {1.96071365, 2.47265689, 1.051449486},
            {0.745853269, 0.5551359072, 4.080774681},
            {-0.3515016103, -0.5192222523, 0.9941510909},
            {-0.8354426872, 0.7071638258, -0.3409843552},
            {0.4386044751, 1.237378731, 2.331874152},
            {2.125850717, 0.3221067321, 0.5825731561},
        },
        {
            {-0.4633736785, 0.06318772224, -0.8580153742},
            {-1.174926354, -0.6276503679, 0.07458759314},
            {1.327618206, 2.085829379, 1.415749862},
            {0.9114727103, 0.1789183931, -0.08135540251},
            {-2.267908723, 0.802928773, 0.9522812957},
            {1.502715257, -1.84493529, 0.2805620469},
            {3.168934617, 0.1348337978, 1.371092768},
            {0.8310229518, 1.070827168, 1.18016733},
        },
        {
            {-0.04847732172, -1.201739871, -1.700527771},
            {0.1589259538, -0.3096047065, -2.066626415},
            {2.255976232, 1.629132391, -0.8024446773},
            {2.534792993, 3.121092901, 1.609780703},
            {-0.2892376071, -0.152022511, -2.727613712},
            {0.2477154804, 0.5039232765, 2.995702733},
            {3.679345099, 3.037770313, 2.808899306},
            {0.6418578532, 1.935944544, 1.515637954},
        },
        {
            {0.91126951, 0.0234699242, 1.442297821},
            {-0.9240061217, -0.1014997844, 0.9081020061},
            {1.887794866, 2.210192703, 2.209118551},
            {2.758945014, -1.21140421, 1.3337907},
            {0.376540703, 0.3485486555, 0.9056881595},
            {-0.3512770187, -0.4056820917, -2.068499576},
            {0.5358460986, 2.720153363, 1.41999706},
            {2.284020089, 1.173071915, 1.044597715},
        },
    };
    for (int iw = 0; iw < num_walkers; ++iw)
    {
      p_refs[iw].get().R = deterministic_rs[iw];
    }
  }

  ParticleSet::mw_update(p_list);

  QMCHamiltonian::mw_evaluate(ham_list, twf_list, p_list);


  if constexpr (generate_test_data)
  {
    std::vector<Real> vector_kinetic;
    std::copy(kinetic.begin(), kinetic.end(), std::back_inserter(vector_kinetic));
    std::cout << " size kinetic: " << vector_kinetic.size() << '\n';
    std::cout << " std::vector<Real> kinetic_ref_vector = " << NativePrint(vector_kinetic) << ";\n";

    std::vector<Real> vector_pots;
    std::copy(local_pots.begin(), local_pots.end(), std::back_inserter(vector_pots));
    std::cout << " size potentials: " << vector_pots.size() << '\n';
    std::cout << " std::vector<Real> potential_ref_vector = " << NativePrint(vector_pots) << ";\n";

    std::vector<Real> vector_ions;
    std::copy(ion_pots.begin(), ion_pots.end(), std::back_inserter(vector_ions));
    std::cout << " size ion potentials: " << vector_ions.size() << '\n';
    std::cout << " std::vector<Real> ion_potential_ref_vector = " << NativePrint(vector_ions) << ";\n";
  }
  else
  {
    std::vector<Real> kinetic_ref_vector = {
        -0.5, -0, -0,   -0, -0, -0, -0, -0, -0, -0.5, -0, -0,   -0, -0, -0, -0,
        -0,   -0, -0.5, -0, -0, -0, -0, -0, -0, -0,   -0, -0.5, -0, -0, -0, -0,
    };
    std::vector<Real> potential_ref_vector = {
        -0.4304472804, -0.2378402054, -0.9860844016, -0.08201235533, -0.8099648952, -0.9536881447, -0.4498806596,
        -0.4415832162, -0.9407452345, -0.6449454427, -2.240605593,   -1.194514155,  0.04097786546, -0.1860728562,
        -0.1374063194, -0.8808210492, 0.3020428121,  0.3772183955,   -0.112801373,  -0.5531326532, 0.4633262753,
        0.3185032904,  0.3001851439,  -0.4555109739, -0.2190704495,  -0.7140043378, -1.641614318,  0.1718038917,
        -0.7621642947, 0.2606962323,  -0.1486036181, -1.001747012,
    };
    std::vector<Real> ion_potential_ref_vector = {
        0.04332125187, 0.173753351, -0.9332901239, -1.323815107, 1.695662975, 0.5990064144, 0.1634206772, -1.207304001,
    };

    std::size_t num_data = num_electrons * num_walkers;
    for (std::size_t id = 0; id < num_data; ++id)
    {
      INFO("id : " << id);
      CHECK(kinetic_ref_vector[id] == Approx(kinetic(id)));
    }

    for (std::size_t id = 0; id < num_data; ++id)
    {
      INFO("id : " << id);
      CHECK(potential_ref_vector[id] == Approx(local_pots(id)));
    }

    std::size_t num_data_ions = num_ions * num_walkers;
    for (std::size_t id = 0; id < num_data_ions; ++id)
    {
      INFO("id : " << id);
      CHECK(ion_potential_ref_vector[id] == Approx(ion_pots(id)));
    }

    // When only EE and EI coulomb potentials the local energy is just the sum of
    // the local_pots and kinetic
    auto sum_local_pots = std::accumulate(local_pots.begin(), local_pots.end(), 0.0);
    auto sum_kinetic    = std::accumulate(kinetic.begin(), kinetic.end(), 0.0);
    auto sum_local_nrg  = std::accumulate(local_nrg.begin(), local_nrg.end(), 0.0);
    CHECK(sum_local_nrg == Approx(sum_local_pots + sum_kinetic));
  }
}

} // namespace qmcplusplus
