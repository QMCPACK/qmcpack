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
#include "Utilities/StlPrettyPrint.hpp"
#include "Utilities/for_testing/NativeInitializerPrint.hpp"

namespace qmcplusplus
{
using QMCT = QMCTraits;
using Real = QMCT::RealType;

constexpr bool generate_test_data = false;

TEST_CASE("QMCHamiltonian::flex_evaluate", "[hamiltonian]")
{
  Communicate* comm;
  comm = OHMMS::Controller;

  outputManager.pause();

  auto particle_pool     = MinimalParticlePool::make_diamondC_1x1x1(comm);
  auto wavefunction_pool = MinimalWaveFunctionPool::make_diamondC_1x1x1(comm, particle_pool);
  auto hamiltonian_pool  = MinimalHamiltonianPool::make_hamWithEE(comm, particle_pool, wavefunction_pool);

  TrialWaveFunction twf;

  std::vector<QMCHamiltonian> hamiltonians;
  //hamiltonians.emplace_back(*(hamiltonian_pool.getPrimary()));
  //hamiltonians.emplace_back(*(hamiltonian_pool.getPrimary()));

  std::vector<ParticleSet> elecs;
  elecs.emplace_back(*(particle_pool.getParticleSet("e")));
  elecs.emplace_back(*(particle_pool.getParticleSet("e")));

  // TODO: finish initializing the elecs.
  //std::vector<QMCHamiltonian::RealType> local_energies(QMCHamiltonian::flex_evaluate(makeRefVector<QMCHamiltonian>(hamiltonians), makeRefVector<ParticleSet>(elecs)));

  //TODO: Would be nice to check some values but I think the system needs a little more setup
  outputManager.resume();
}

TEST_CASE("QMCHamiltonian-registerListeners", "[hamiltonian]")
{
  using testing::getParticularListener;
  Communicate* comm;
  comm = OHMMS::Controller;

  outputManager.pause();

  auto particle_pool       = MinimalParticlePool::make_diamondC_1x1x1(comm);
  auto wavefunction_pool   = MinimalWaveFunctionPool::make_diamondC_1x1x1(comm, particle_pool);
  auto hamiltonian_pool    = MinimalHamiltonianPool::makeHamWithEEEI(comm, particle_pool, wavefunction_pool);
  auto& pset_target        = *(particle_pool.getParticleSet("e"));
  //auto& species_set        = pset_target.getSpeciesSet();
  //auto& spo_map            = wavefunction_pool.getWaveFunction("wavefunction")->getSPOMap();
  auto& trial_wavefunction = *(wavefunction_pool.getPrimary());

  UPtrVector<QMCHamiltonian> hams;
  UPtrVector<TrialWaveFunction> twfs;
  std::vector<ParticleSet> psets;

  // This must be done before clones otherwise the clone particle sets do not have the correct state.
  hamiltonian_pool.getPrimary()->informOperatorsOfListener();


  int num_walkers   = 4;
  int num_electrons = 8;
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

  ListenerVector<Real> listener_kinetic("kinetic", getParticularListener(kinetic));
  QMCHamiltonian::mw_registerKineticListener(ham_list.getLeader(), listener_kinetic);
  ListenerVector<Real> listener_potential("kinetic", getParticularListener(local_pots));
  QMCHamiltonian::mw_registerLocalPotentialListener(ham_list.getLeader(), listener_potential);

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
  }
  else
  {
    std::vector<Real> kinetic_ref_vector = {
        -0.5, -0, -0,   -0, -0, -0, -0, -0, -0, -0.5, -0, -0,   -0, -0, -0, -0,
        -0,   -0, -0.5, -0, -0, -0, -0, -0, -0, -0,   -0, -0.5, -0, -0, -0, -0,
    };
    std::vector<Real> potential_ref_vector = {
        0.1902921988,  0.3051115892,  -0.4795567016, 0.4963486852,  -0.3115299199, -0.3619893752, 0.09889832988,
        0.2794993746,  -0.5591061164, -0.1018774556, -1.758163533,  -0.6303051358, 0.4878929502,  0.3358202323,
        0.3278906128,  -0.3592569681, 0.2308181408,  0.3380463886,  0.4502142541,  0.09129681648, 0.5120332571,
        0.4968745532,  0.52070839,    -0.3453223967, 0.1190383566,  -0.1566598479, -1.061808709,  0.4489104638,
        -0.4774404557, 0.3925164499,  0.1283951706,  -0.4368348731,
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
  }
  outputManager.resume();
}
} // namespace qmcplusplus
