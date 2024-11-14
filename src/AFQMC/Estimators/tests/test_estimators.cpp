//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Miguel Morales, moralessilva2@llnl.gov, Lawrence Livermore National Laboratory
//
// File created by: Fionn Malone, malone14@llnl.gov, Lawrence Livermore National Laboratory
//////////////////////////////////////////////////////////////////////////////////////

#include "catch.hpp"

#include "Configuration.h"

#include "OhmmsData/Libxml2Doc.h"
#include "ProjectData.h"
#include "Utilities/TimerManager.h"
#include "hdf/hdf_archive.h"

#undef APP_ABORT
#define APP_ABORT(x)             \
  {                              \
    std::cout << x << std::endl; \
    throw;                       \
  }

#include <stdio.h>
#include <string>

#include "AFQMC/config.h"
#include "AFQMC/Matrix/tests/matrix_helpers.h"
#include "AFQMC/Hamiltonians/HamiltonianFactory.h"
#include "AFQMC/Hamiltonians/Hamiltonian.hpp"
#include "AFQMC/Wavefunctions/WavefunctionFactory.h"
#include "AFQMC/Wavefunctions/Wavefunction.hpp"
#include "AFQMC/Walkers/WalkerSet.hpp"
#include "AFQMC/Estimators/EstimatorBase.h"
#include "AFQMC/Propagators/PropagatorFactory.h"
#include "AFQMC/Propagators/Propagator.hpp"
#include "AFQMC/Estimators/BackPropagatedEstimator.hpp"
#include "AFQMC/Utilities/test_utils.hpp"
#include "AFQMC/Memory/buffer_managers.h"

using std::cerr;
using std::complex;
using std::cout;
using std::endl;
using std::ifstream;
using std::setprecision;
using std::string;

extern std::string UTEST_HAMIL, UTEST_WFN;

namespace qmcplusplus
{
using namespace afqmc;

template<class Allocator>
void reduced_density_matrix(boost::mpi3::communicator& world)
{
  using pointer = typename Allocator::pointer;

  if (check_hamil_wfn_for_utest("reduced_density_matrix", UTEST_WFN, UTEST_HAMIL))
  {
    getGlobalTimerManager().set_timer_threshold(timer_level_coarse);

    // Global Task Group
    afqmc::GlobalTaskGroup gTG(world);

    int NMO, NAEA, NAEB;
    std::tie(NMO, NAEA, NAEB) = read_info_from_hdf(UTEST_HAMIL);

    std::map<std::string, AFQMCInfo> InfoMap;
    InfoMap.insert(std::pair<std::string, AFQMCInfo>("info0", AFQMCInfo{"info0", NMO, NAEA, NAEB}));
    HamiltonianFactory HamFac(InfoMap);
    std::string hamil_xml = R"(<Hamiltonian name="ham0" info="info0">
      <parameter name="filetype">hdf5</parameter>
      <parameter name="filename">)" +
        UTEST_HAMIL + R"(</parameter>
      <parameter name="cutoff_decomposition">1e-5</parameter>
    </Hamiltonian>
    )";
    const char* ham_xml_block = hamil_xml.c_str();
    Libxml2Document doc;
    bool okay = doc.parseFromString(ham_xml_block);
    REQUIRE(okay);
    std::string ham_name("ham0");
    HamFac.push(ham_name, doc.getRoot());
    Hamiltonian& ham = HamFac.getHamiltonian(gTG, ham_name);

    WALKER_TYPES type                = afqmc::getWalkerType(UTEST_WFN);
    const char* wlk_xml_block_closed = R"(<WalkerSet name="wset0">
      <parameter name="walker_type">closed</parameter>
      <parameter name="back_propagation_steps">10</parameter>
    </WalkerSet>
    )";
    const char* wlk_xml_block_coll   = R"(<WalkerSet name="wset0">
      <parameter name="walker_type">collinear</parameter>
      <parameter name="back_propagation_steps">10</parameter>
    </WalkerSet>
    )";
    const char* wlk_xml_block_noncol = R"(<WalkerSet name="wset0">
      <parameter name="walker_type">noncollinear</parameter>
      <parameter name="back_propagation_steps">10</parameter>
    </WalkerSet>
    )";

    const char* wlk_xml_block =
        ((type == CLOSED) ? (wlk_xml_block_closed) : (type == COLLINEAR ? wlk_xml_block_coll : wlk_xml_block_noncol));
    Libxml2Document doc3;
    okay = doc3.parseFromString(wlk_xml_block);
    REQUIRE(okay);

    std::string wfn_xml = R"(<Wavefunction name="wfn0" info="info0">
      <parameter name="filetype">ascii</parameter>
      <parameter name="filename">)" +
        UTEST_WFN + R"(</parameter>
      <parameter name="cutoff">1e-6</parameter>
    </Wavefunction>
    )";

    const char* wfn_xml_block = wfn_xml.c_str();
    auto TG                   = TaskGroup_(gTG, std::string("WfnTG"), 1, gTG.getTotalCores());
    Allocator alloc_(make_localTG_allocator<ComplexType>(TG));
    int nwalk = 1; // choose prime number to force non-trivial splits in shared routines
    RandomGenerator rng;
    Libxml2Document doc2;
    okay = doc2.parseFromString(wfn_xml_block);
    REQUIRE(okay);
    std::string wfn_name("wfn0");
    WavefunctionFactory WfnFac(InfoMap);
    WfnFac.push(wfn_name, doc2.getRoot());
    Wavefunction& wfn = WfnFac.getWavefunction(TG, TG, wfn_name, type, &ham, 1e-6, nwalk);

    const char* propg_xml_block = R"(<Propagator name="prop0"></Propagator>)";
    Libxml2Document doc5;
    okay = doc5.parseFromString(propg_xml_block);
    REQUIRE(okay);
    std::string prop_name("prop0");
    PropagatorFactory PropgFac(InfoMap);
    PropgFac.push(prop_name, doc5.getRoot());
    Propagator& prop = PropgFac.getPropagator(TG, prop_name, wfn, rng);

    WalkerSet wset(TG, doc3.getRoot(), InfoMap["info0"], rng);
    auto initial_guess = WfnFac.getInitialGuess(wfn_name);
    REQUIRE(std::get<0>(initial_guess.sizes()) == 2);
    REQUIRE(std::get<1>(initial_guess.sizes()) == NMO);
    REQUIRE(std::get<2>(initial_guess.sizes()) == NAEA);
    wset.resize(nwalk, initial_guess[0], initial_guess[0]);
    using EstimPtr = std::shared_ptr<EstimatorBase>;
    std::vector<EstimPtr> estimators;
    const char* est_xml_block = R"(<Estimator name="back_propagation">
      <parameter name="nsteps">1</parameter>
      <parameter name="block_size">2</parameter>
      <OneRDM> </OneRDM>
      <parameter name="path_restoration">false</parameter>
    </Estimator>
    )";
    Libxml2Document doc4;
    okay = doc4.parseFromString(est_xml_block);
    REQUIRE(okay);
    bool impsamp = true;
    estimators.push_back(std::make_shared<BackPropagatedEstimator>(TG, InfoMap["info0"], "none", doc4.getRoot(), type,
                                                                   wset, wfn, prop, impsamp));

    // generate P1 with dt=0
    prop.generateP1(0.0, wset.getWalkerType());

    std::string file = create_test_hdf(UTEST_WFN, UTEST_HAMIL);
    hdf_archive dump;
    std::ofstream out;
    dump.create(file);
    dump.open(file);
    for (int iblock = 0; iblock < 10; iblock++)
    {
      wset.advanceBPPos();
      estimators[0]->accumulate_block(wset);
      estimators[0]->print(out, dump, wset);
    }
    dump.close();
    boost::multi::array<ComplexType, 1> read_data(boost::multi::iextensions<1u>{2 * NMO * NMO});

    ComplexType denom;
    hdf_archive reader;
    // Read from a particular block.
    if (!reader.open(file, H5F_ACC_RDONLY))
    {
      app_error() << " Error opening estimates.h5. \n";
      APP_ABORT("");
    }
    reader.read(read_data, "Observables/BackPropagated/FullOneRDM/Average_0/one_rdm_000000004");
    reader.read(denom, "Observables/BackPropagated/FullOneRDM/Average_0/denominator_000000004");
    // Test EstimatorHandler eventually.
    //int NAEA_READ, NAEB_READ, NMO_READ, WALKER_TYPE_READ;
    //reader.read(NAEA_READ, "Metadata/NAEA");
    //REQUIRE(NAEA_READ==NAEA);
    //reader.read(NAEB_READ, "Metadata/NAEB");
    //REQUIRE(NAEB_READ==NAEB);
    //reader.read(NMO_READ, "Metadata/NMO");
    //REQUIRE(NMO_READ==NMO);
    //reader.read(WALKER_TYPE_READ, "Metadata/WALKER_TYPE");
    //REQUIRE(WALKER_TYPE_READ==type);
    reader.close();
    // Test the RDM. Since no back propagation has been performed the RDM should be
    // identical to the mixed estimate.
    if (type == CLOSED)
    {
      REQUIRE(read_data.num_elements() >= NMO * NMO);
      boost::multi::array_ref<ComplexType, 2> BPRDM(read_data.origin(), {NMO, NMO});
      ma::scal(1.0 / denom, BPRDM);
      ComplexType trace = ComplexType(0.0);
      for (int i = 0; i < NMO; i++)
        trace += BPRDM[i][i];
      CHECK(trace.real() == Approx(NAEA));
      boost::multi::array<ComplexType, 2, Allocator> Gw({1, NMO * NMO}, alloc_);
      wfn.MixedDensityMatrix(wset, Gw, false, true);
      boost::multi::array_ref<ComplexType, 2, pointer> G(Gw.origin(), {NMO, NMO});
      verify_approx(G, BPRDM);
    }
    else if (type == COLLINEAR)
    {
      REQUIRE(read_data.num_elements() >= 2 * NMO * NMO);
      boost::multi::array_ref<ComplexType, 3> BPRDM(read_data.origin(), {2, NMO, NMO});
      ma::scal(1.0 / denom, BPRDM[0]);
      ma::scal(1.0 / denom, BPRDM[1]);
      ComplexType trace = ComplexType(0.0);
      for (int i = 0; i < NMO; i++)
        trace += BPRDM[0][i][i] + BPRDM[1][i][i];
      CHECK(trace.real() == Approx(NAEA + NAEB));
      boost::multi::array<ComplexType, 2, Allocator> Gw({1, 2 * NMO * NMO}, alloc_);
      wfn.MixedDensityMatrix(wset, Gw, false, true);
      boost::multi::array_ref<ComplexType, 3, pointer> G(Gw.origin(), {2, NMO, NMO});
      verify_approx(G, BPRDM);
    }
    else
    {
      APP_ABORT(" NONCOLLINEAR Wavefunction found.\n");
    }
  }
}

TEST_CASE("reduced_density_matrix", "[estimators]")
{
  auto world = boost::mpi3::environment::get_world_instance();
  if (not world.root())
    infoLog.pause();
  auto node = world.split_shared(world.rank());

#if defined(ENABLE_CUDA) || defined(ENABLE_HIP)
  arch::INIT(node);
  using Alloc = device::device_allocator<ComplexType>;
#else
  using Alloc = shared_allocator<ComplexType>;
#endif
  setup_memory_managers(node, 10uL * 1024uL * 1024uL);
  reduced_density_matrix<Alloc>(world);
  release_memory_managers();
}

} // namespace qmcplusplus
