//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Miguel Morales, moralessilva2@llnl.gov, Lawrence Livermore National Laboratory
//
// File created by: Miguel Morales, moralessilva2@llnl.gov, Lawrence Livermore National Laboratory
//////////////////////////////////////////////////////////////////////////////////////

//#undef NDEBUG

#include "catch.hpp"

#include "Configuration.h"

#include "OhmmsData/Libxml2Doc.h"
#include "ProjectData.h"
#include "hdf/hdf_archive.h"
#include "Utilities/RandomGenerator.h"
#include "Utilities/Timer.h"
#include "Platforms/Host/OutputManager.h"

#undef APP_ABORT
#define APP_ABORT(x)             \
  {                              \
    std::cout << x << std::endl; \
    throw;                       \
  }

#include <string>
#include <vector>
#include <complex>
#include <iomanip>
#include <random>
#include <algorithm>

#include "AFQMC/Wavefunctions/Excitations.hpp"
#include "AFQMC/Wavefunctions/WavefunctionFactory.h"
#include "AFQMC/Hamiltonians/HamiltonianFactory.h"
#include "AFQMC/Hamiltonians/Hamiltonian.hpp"
#include "AFQMC/Walkers/WalkerSet.hpp"
#include "AFQMC/SlaterDeterminantOperations/SlaterDetOperations.hpp"
#include "AFQMC/Utilities/test_utils.hpp"
#include "AFQMC/Utilities/Utils.hpp"
#include "AFQMC/Utilities/taskgroup.h"
#include "AFQMC/Utilities/readWfn.h"
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
void test_read_phmsd(boost::mpi3::communicator& world)
{
  using pointer = device_ptr<ComplexType>;

  if (check_hamil_wfn_for_utest("test_read_phmsd", UTEST_WFN, UTEST_HAMIL))
  {
    // Global Task Group
    GlobalTaskGroup gTG(world);
    auto TG    = TaskGroup_(gTG, std::string("WfnTG"), 1, gTG.getTotalCores());
    auto TGwfn = TaskGroup_(gTG, std::string("WfnTG"), 1, gTG.getTotalCores());
    Allocator alloc_(make_localTG_allocator<ComplexType>(TG));

    int NMO;
    int NAEA;
    int NAEB;
    std::tie(NMO, NAEA, NAEB) = read_info_from_wfn(UTEST_WFN, "PHMSD");
    hdf_archive dump;
    if (!dump.open(UTEST_WFN, H5F_ACC_RDONLY))
    {
      app_error() << "Error reading wavefunction file.\n";
    }
    dump.push("Wavefunction", false);
    dump.push("PHMSD", false);

    int ndets_to_read = -1;
    std::string wfn_type;
    WALKER_TYPES walker_type = COLLINEAR;
    std::vector<PsiT_Matrix> PsiT_MO; // read_ph_wavefunction returns the full MO matrix
    std::vector<int> occbuff;
    std::vector<ComplexType> coeffs;
    //ph_excitations<int,ComplexType> abij = read_ph_wavefunction_hdf(dump,ndets_to_read,walker_type,
    //gTG.Node(),NMO,NAEA,NAEB,PsiT_MO,
    //wfn_type);
    read_ph_wavefunction_hdf(dump, coeffs, occbuff, ndets_to_read, walker_type, TGwfn.Node(), NMO, NAEA, NAEB, PsiT_MO,
                             wfn_type);
    boost::multi::array_ref<int, 2> occs(to_address(occbuff.data()), {ndets_to_read, NAEA + NAEB});
    ph_excitations<int, ComplexType> abij = build_ph_struct(coeffs, occs, ndets_to_read, TGwfn.Node(), NMO, NAEA, NAEB);
    //std::vector<int> buff(ndets_to_read*(NAEA+NAEB));
    //boost::multi::array_ref<int,2> occs(to_address(buff.data()), {ndets_to_read, NAEA+NAEB});
    //if(!dump.readEntry(buff, "occs"))
    //APP_ABORT("Error reading occs array.\n");
    //std::vector<ComplexType> ci_coeffs(ndets_to_read);
    //if(!dump.readEntry(ci_coeffs, "ci_coeffs"))
    //APP_ABORT("Error reading occs array.\n");
    using std::get;
    auto cit = abij.configurations_begin();
    std::vector<int> configa(NAEA), configb(NAEB);
    // Is it fortuitous that the order of determinants is the same?
    for (int nd = 0; nd < ndets_to_read; nd++, ++cit)
    {
      int alpha_ix = get<0>(*cit);
      int beta_ix  = get<1>(*cit);
      auto ci      = get<2>(*cit);
      abij.get_configuration(0, alpha_ix, configa);
      abij.get_configuration(1, beta_ix, configb);
      std::sort(configa.begin(), configa.end());
      std::sort(configb.begin(), configb.end());
      for (int i = 0; i < NAEA; i++)
      {
        REQUIRE(configa[i] == occs[nd][i]);
      }
      for (int i = 0; i < NAEB; i++)
      {
        REQUIRE(configb[i] == occs[nd][i + NAEA]);
      }
      REQUIRE(std::abs(coeffs[nd]) == std::abs(ci));
    }
    // Check sign of permutation.
    REQUIRE(abij.number_of_configurations() == ndets_to_read);
  }
}

void getBasicWavefunction(std::vector<int>& occs, std::vector<ComplexType>& coeffs, int NEL)
{
  hdf_archive dump;
  if (!dump.open(UTEST_WFN, H5F_ACC_RDONLY))
  {
    app_error() << "Error reading wavefunction file.\n";
  }
  dump.push("Wavefunction", false);
  dump.push("PHMSD", false);

  std::vector<int> Idata(5);
  if (!dump.readEntry(Idata, "dims"))
    APP_ABORT("Errro reading dims array\n");
  int ndets = Idata[4];
  occs.resize(ndets * NEL);
  if (!dump.readEntry(occs, "occs"))
    APP_ABORT("Error reading occs array.\n");
  std::vector<ComplexType> ci_coeffs(ndets);
  if (!dump.readEntry(coeffs, "ci_coeffs"))
    APP_ABORT("Error reading occs array.\n");
}

// Construct PsiT^{dagger}
void getSlaterMatrix(boost::multi::array<ComplexType, 2>& SM, boost::multi::array_ref<int, 1>& occs, int NEL)
{
  using std::fill_n;
  fill_n(SM.origin(), SM.num_elements(), ComplexType(0.0));
  for (int i = 0; i < NEL; i++)
    SM[i][occs[i]] = ComplexType(1.0);
}

template<class Allocator>
void test_phmsd(boost::mpi3::communicator& world)
{
  using pointer = device_ptr<ComplexType>;

  if (check_hamil_wfn_for_utest("read_phmsd", UTEST_WFN, UTEST_HAMIL))
  {
    // Global Task Group
    GlobalTaskGroup gTG(world);
    auto TG    = TaskGroup_(gTG, std::string("WfnTG"), 1, gTG.getTotalCores());
    auto TGwfn = TaskGroup_(gTG, std::string("WfnTG"), 1, gTG.getTotalCores());
    Allocator alloc_(make_localTG_allocator<ComplexType>(TG));

    int NMO;
    int NAEA;
    int NAEB;
    std::tie(NMO, NAEA, NAEB) = read_info_from_wfn(UTEST_WFN, "PHMSD");
    // Test overlap.
    //wfn.Overlap(wset);
    WALKER_TYPES type = afqmc::getWalkerTypeHDF5(UTEST_WFN, "PHMSD");
    std::map<std::string, AFQMCInfo> InfoMap;
    InfoMap.insert(std::make_pair("info0", AFQMCInfo{"info0", NMO, NAEA, NAEB}));
    HamiltonianFactory HamFac(InfoMap);
    std::string hamil_xml = R"(<Hamiltonian name="ham0" info="info0">
      <parameter name="filetype">hdf5</parameter>
      <parameter name="filename">)" +
        UTEST_HAMIL + R"(</parameter>
      <parameter name="cutoff_decomposition">1e-12</parameter>
      <parameter name="cutoff_1bar">1e-12</parameter>
    </Hamiltonian>
    )";
    const char* ham_xml_block = hamil_xml.c_str();
    Libxml2Document doc;
    bool okay = doc.parseFromString(ham_xml_block);
    REQUIRE(okay);
    std::string ham_name("ham0");
    HamFac.push(ham_name, doc.getRoot());
    Hamiltonian& ham = HamFac.getHamiltonian(gTG, ham_name);

    std::string wfn_xml = R"(<Wavefunction name="wfn0" info="info0" type="phmsd">
      <parameter name="filetype">hdf5</parameter>
      <parameter name="filename">)" +
        UTEST_WFN + R"(</parameter>
      <parameter name="rediag">true</parameter>
      <parameter name="cutoff">1e-6</parameter>
    </Wavefunction>
    )";
    const char* wfn_xml_block = wfn_xml.c_str();
    Libxml2Document doc2;
    okay = doc2.parseFromString(wfn_xml_block);
    REQUIRE(okay);
    std::string wfn_name("wfn0");
    WavefunctionFactory WfnFac(InfoMap);
    WfnFac.push(wfn_name, doc2.getRoot());
    int nwalk                 = 1;
    Wavefunction& wfn         = WfnFac.getWavefunction(TGwfn, TGwfn, wfn_name, type, &ham, 1e-6, nwalk);
    const char* wlk_xml_block = R"(<WalkerSet name="wset0">
      <parameter name="walker_type">collinear</parameter>
    </WalkerSet>
    )";
    Libxml2Document doc3;
    okay = doc3.parseFromString(wlk_xml_block);
    REQUIRE(okay);
    RandomGenerator rng;
    WalkerSet wset(TG, doc3.getRoot(), InfoMap["info0"], rng);
    auto initial_guess = WfnFac.getInitialGuess(wfn_name);
    REQUIRE(std::get<0>(initial_guess.sizes()) == 2);
    REQUIRE(std::get<1>(initial_guess.sizes()) == NMO);
    REQUIRE(std::get<2>(initial_guess.sizes()) == NAEA);

    wset.resize(nwalk, initial_guess[0], initial_guess[1](initial_guess.extension(1), {0, NAEB}));
    // 1. Test Overlap Explicitly
    // 1.a Get raw occupancies and coefficients from file.
    std::vector<ComplexType> coeffs;
    std::vector<int> buff;
    getBasicWavefunction(buff, coeffs, NAEA + NAEB);
    int ndets = coeffs.size();
    boost::multi::array_ref<int, 2> occs(buff.data(), {ndets, NAEA + NAEB});
    // 1.b Compute overlap of trial wavefunction compotents.
    boost::multi::array<ComplexType, 2> Orbs({NMO, NMO});
    for (int i = 0; i < NMO; i++)
      Orbs[i][i] = ComplexType(1.0);
    boost::multi::array<ComplexType, 2> TrialA({NAEB, NMO}), TrialB({NAEB, NMO});
    auto sdet = wfn.getSlaterDetOperations();
    ComplexType ovlp;
    ComplexType ovlp_sum = ComplexType(0.0);
    ComplexType logovlp;
    //boost::multi::array<ComplexType,2> GBuff;
    //GBuff.reextent({NAEA+NAEB,NMO});
    for (int idet = 0; idet < coeffs.size(); idet++)
    {
      // Construct slater matrix from given set of occupied orbitals.
      boost::multi::array_ref<int, 1> oa(occs[idet].origin(), {NAEA});
      getSlaterMatrix(TrialA, oa, NAEA);
      boost::multi::array_ref<int, 1> ob(occs[idet].origin() + NAEA, {NAEB});
      for (int i = 0; i < NAEB; i++)
        ob[i] -= NMO;
      getSlaterMatrix(TrialB, ob, NAEB);
      ComplexType ovlpa = sdet->Overlap(TrialA, *wset[0].SlaterMatrix(Alpha), logovlp);
      ComplexType ovlpb = sdet->Overlap(TrialB, *wset[0].SlaterMatrix(Beta), logovlp);
      ovlp_sum += ma::conj(coeffs[idet]) * ovlpa * ovlpb;
      //boost::multi::array_ref<ComplexType,2> GB(to_address(GBuff.origin()), {NAEA,NMO});
      //sdet->MixedDensityMatrix(TrialB, wset[0].SlaterMatrix(Alpha), GA, logovlp, true);
      //boost::multi::array_ref<ComplexType,2> GA(to_address(GBuff.origin()+NAEA*NMO), {NAEA,NMO});
      //sdet->MixedDensityMatrix(TrialA, wset[0].SlaterMatrix(Alpha), GB, logovlp, true);
    }
    wfn.Overlap(wset);
    // TODO: Remove abs here. I want to test rediag but currently no way to get updated ci
    // coefficients as using factory setup.
    for (auto it = wset.begin(); it != wset.end(); ++it)
    {
      CHECK(std::abs(real(*it->overlap())) == Approx(std::abs(real(ovlp_sum))));
      CHECK(std::abs(imag(*it->overlap())) == Approx(std::abs(imag(ovlp_sum))));
    }
    // It's not straightforward to calculate energy directly in unit test due to half
    // rotation.
    //wfn.Energy(wset);
    //for(auto it = wset.begin(); it!=wset.end(); ++it) {
    //CHECK(real(*it->energy()) == Approx(real(energy)));
    //CHECK(imag(*it->energy()) == Approx(imag(energy)));
    //}
    //auto nCV = wfn.local_number_of_cholesky_vectors();
    //boost::multi::array<ComplexType,1> vMF(iextensions<1u>{nCV});
    //std::cout << "NCHOL : " << nCV << " " << NMO*NMO << std::endl;
    //wfn.vMF(vMF);
    //computeVariationalEnergy(wfn, occs, ham, NAEA, NAEB);
    //std::vector<ComplexType> vMF_sc = computeMeanFieldShift(wfn, occs, coeffs, NAEA, NAEB);
    //for(int i=0; i < vMF.size(); i++) {
    //std::cout << vMF[i] << std::endl;
    //}
  }
}

TEST_CASE("test_read_phmsd", "[test_read_phmsd]")
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

  test_read_phmsd<Alloc>(world);

  release_memory_managers();
}

TEST_CASE("test_phmsd", "[read_phmsd]")
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

  //test_phmsd<Alloc,SlaterDetOperations_serial<Alloc>>(world);
  test_phmsd<Alloc>(world);

  release_memory_managers();
}

} // namespace qmcplusplus
