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

#include "AFQMC/Hamiltonians/HamiltonianFactory.h"
#include "AFQMC/Hamiltonians/Hamiltonian.hpp"
#include "AFQMC/Utilities/myTimer.h"
#include "AFQMC/Utilities/readWfn.h"
#include "AFQMC/SlaterDeterminantOperations/SlaterDetOperations.hpp"
#include "AFQMC/Utilities/test_utils.hpp"

#include "AFQMC/Matrix/csr_matrix_construct.hpp"
#include "AFQMC/Numerics/ma_blas.hpp"

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
void ham_factory(boost::mpi3::communicator& world)
{
  using pointer = device_ptr<ComplexType>;
  if (not file_exists(UTEST_HAMIL))
  {
    app_log() << " Skipping ham_ops_basic_serial. Hamiltonian file not found. \n";
    app_log() << " Run unit test with --hamil /path/to/hamil.h5.\n";
  }
  else
  {
    // Global Task Group
    GlobalTaskGroup gTG(world);

    int NMO, NAEA, NAEB;
    std::tie(NMO, NAEA, NAEB) = read_info_from_hdf(UTEST_HAMIL);
    REQUIRE(NAEA == NAEB);

    std::map<std::string, AFQMCInfo> InfoMap;
    InfoMap.insert(std::make_pair("info0", AFQMCInfo{"info0", NMO, NAEA, NAEB}));
    HamiltonianFactory HamFac(InfoMap);
    std::string hamil_xml = R"(<Hamiltonian name="ham0" info="info0">
      <parameter name="filetype">hdf5</parameter>
      <parameter name="filename">)" +
        UTEST_HAMIL + R"(</parameter>
      <parameter name="cutoff_decomposition">1e-5</parameter>
    </Hamiltonian>
    )";
    const char* xml_block = hamil_xml.c_str();
    Libxml2Document doc;
    bool okay = doc.parseFromString(xml_block);
    REQUIRE(okay);
    std::string ham_name("ham0");
    HamFac.push(ham_name, doc.getRoot());

    Hamiltonian& ham_ = HamFac.getHamiltonian(gTG, ham_name);
  }
}

template<class Alloc>
void ham_generation_timing(boost::mpi3::communicator& world)
{
  if (not file_exists("./afqmc_timing.h5"))
  {
    app_log() << " Skipping ham_fac_timing text. afqmc_timing.h5 file not found. \n";
  }
  else
  {
    // Global Task Group
    afqmc::GlobalTaskGroup gTG(world);

    auto [NMO, NAEA, NAEB] = read_info_from_hdf("./afqmc_timing.h5");

    std::map<std::string, AFQMCInfo> InfoMap;
    InfoMap.insert(std::make_pair("info0", AFQMCInfo{"info0", NMO, NAEA, NAEB}));
    HamiltonianFactory HamFac(InfoMap);

    const char* xml_block = R"(<Hamiltonian name="ham0" type="SparseGeneral" info="info0">
      <parameter name="filetype">hdf5</parameter>
      <parameter name="version">new</parameter>
      <parameter name="filename">./afqmc_timing.h5</parameter>
      <parameter name="cutoff_decomposition">1e-5</parameter>
    </Hamiltonian>
    )";
    Libxml2Document doc;
    bool okay = doc.parseFromString(xml_block);
    REQUIRE(okay);
    std::string ham_name("ham0");
    HamFac.push(ham_name, doc.getRoot());

    myTimer Timer_;
    Timer_.start("GenTest0");
    Hamiltonian& ham = HamFac.getHamiltonian(gTG, ham_name);
    Timer_.stop("GenTest0");
    app_log() << "\n*********************************************************************\n"
              << " Time to create hamiltonian in ham_generation_timing_hdf: " << Timer_.total("GenTest0")
              << "\n*********************************************************************\n\n";
  }
}
TEST_CASE("ham_factory", "[hamiltonian_factory]")
{
  auto world = boost::mpi3::environment::get_world_instance();
  if (not world.root())
    infoLog.pause();

#if defined(ENABLE_CUDA) || defined(ENABLE_HIP)
  auto node = world.split_shared(world.rank());

  arch::INIT(node);
  using Alloc = device::device_allocator<ComplexType>;
#else
  auto node   = world.split_shared(world.rank());
  using Alloc = shared_allocator<ComplexType>;
#endif

  ham_factory<Alloc>(world);
}

TEST_CASE("ham_generation_timing_hdf", "[hamiltonian_factory]")
{
  auto world = boost::mpi3::environment::get_world_instance();
  if (not world.root())
    infoLog.pause();

#if defined(ENABLE_CUDA) || defined(ENABLE_HIP)
  auto node = world.split_shared(world.rank());

  arch::INIT(node);
  using Alloc = device::device_allocator<ComplexType>;
#else
  auto node   = world.split_shared(world.rank());
  using Alloc = shared_allocator<ComplexType>;
#endif

  ham_generation_timing<Alloc>(world);
}

} // namespace qmcplusplus
