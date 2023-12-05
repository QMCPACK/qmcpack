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
#include "hdf/hdf_multi.h"

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

#include "AFQMC/config.h"
#include "AFQMC/Hamiltonians/HamiltonianFactory.h"
#include "AFQMC/Hamiltonians/Hamiltonian.hpp"
#include "AFQMC/Hamiltonians/THCHamiltonian.h"
#include "AFQMC/Matrix/csr_hdf5_readers.hpp"
#include "AFQMC/Utilities/readWfn.h"
#include "AFQMC/SlaterDeterminantOperations/SlaterDetOperations.hpp"
#include "AFQMC/Utilities/test_utils.hpp"
#include "AFQMC/Memory/buffer_managers.h"

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

template<class Alloc>
void ham_ops_basic_serial(boost::mpi3::communicator& world)
{
  using pointer = device_ptr<ComplexType>;

  if (check_hamil_wfn_for_utest("ham_ops_basic_serial", UTEST_WFN, UTEST_HAMIL))
  {
    // Global Task Group
    afqmc::GlobalTaskGroup gTG(world);

    // Determine wavefunction type for test results from wavefunction file name which is
    // has the naming convention wfn_(wfn_type).dat.
    // First strip path of filename.
    std::string base_name = UTEST_WFN.substr(UTEST_WFN.find_last_of("\\/") + 1);
    // Remove file extension.
    std::string test_wfn = base_name.substr(0, base_name.find_last_of("."));
    auto file_data       = read_test_results_from_hdf<ValueType>(UTEST_HAMIL, test_wfn);
    int NMO              = file_data.NMO;
    int NAEA             = file_data.NAEA;
    int NAEB             = file_data.NAEB;

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
    const char* xml_block = hamil_xml.c_str();
    Libxml2Document doc;
    bool okay = doc.parseFromString(xml_block);
    REQUIRE(okay);
    std::string ham_name("ham0");
    HamFac.push(ham_name, doc.getRoot());

    Hamiltonian& ham = HamFac.getHamiltonian(gTG, ham_name);

    using CMatrix = ComplexMatrix<Alloc>;
    hdf_archive dump;
    if (!dump.open(UTEST_WFN, H5F_ACC_RDONLY))
    {
      app_error() << " Error opening HDF wavefunction file.\n";
    }
    dump.push("Wavefunction", false);
    dump.push("NOMSD", false);
    std::vector<int> dims(5);
    if (!dump.readEntry(dims, "dims"))
    {
      app_error() << " Error in getCommonInput(): Problems reading dims. \n";
      APP_ABORT("");
    }
    WALKER_TYPES WTYPE(initWALKER_TYPES(dims[3]));
    //    int walker_type = dims[3];
    int NEL  = (WTYPE == CLOSED) ? NAEA : (NAEA + NAEB);
    int NPOL = (WTYPE == NONCOLLINEAR) ? 2 : 1;
    //    WALKER_TYPES WTYPE = CLOSED;
    //    if(walker_type==1) WTYPE = COLLINEAR;
    //    if(walker_type==2) WTYPE = NONCOLLINEAR;

    auto TG = TaskGroup_(gTG, std::string("DummyTG"), 1, gTG.getTotalCores());
    Alloc alloc_(make_localTG_allocator<ComplexType>(TG));
    std::vector<PsiT_Matrix> PsiT;
    PsiT.reserve(2);
    dump.push(std::string("PsiT_0"));
    PsiT.emplace_back(csr_hdf5::HDF2CSR<PsiT_Matrix, shared_allocator<ComplexType>>(dump, gTG.Node()));
    if (WTYPE == COLLINEAR)
    {
      dump.pop();
      dump.push(std::string("PsiT_1"));
      PsiT.emplace_back(csr_hdf5::HDF2CSR<PsiT_Matrix, shared_allocator<ComplexType>>(dump, gTG.Node()));
    }

    dump.pop();
    boost::multi::array<ComplexType, 3> OrbMat({2, NPOL * NMO, NAEA});
    {
      boost::multi::array<ComplexType, 2> Psi0A({NPOL * NMO, NAEA});
      dump.readEntry(Psi0A, "Psi0_alpha");
      for (int i = 0; i < NPOL * NMO; i++)
      {
        for (int j = 0; j < NAEA; j++)
        {
          OrbMat[0][i][j] = Psi0A[i][j];
        }
      }
      if (WTYPE == COLLINEAR)
      {
        boost::multi::array<ComplexType, 2> Psi0B({NMO, NAEA});
        dump.readEntry(Psi0B, "Psi0_beta");
        for (int i = 0; i < NMO; i++)
        {
          for (int j = 0; j < NAEB; j++)
          {
            OrbMat[1][i][j] = Psi0B[i][j];
          }
        }
      }
    }
    dump.close();
    hdf_archive dummy;
    auto HOps(ham.getHamiltonianOperations(false, false, WTYPE, PsiT, 1e-6, 1e-6, TG, TG, dummy));

    // Calculates Overlap, G
// NOTE: Make small factory routine!
#if defined(ENABLE_CUDA) || defined(ENABLE_HIP)
    auto SDet(SlaterDetOperations_serial<ComplexType, DeviceBufferManager>(NPOL * NMO, NAEA, DeviceBufferManager{}));
#else
    auto SDet(SlaterDetOperations_shared<ComplexType>(NPOL * NMO, NAEA));
#endif

    boost::multi::array<ComplexType, 3, Alloc> devOrbMat(OrbMat, alloc_);
    std::vector<devcsr_Matrix> devPsiT(move_vector<devcsr_Matrix>(std::move(PsiT)));

    CMatrix G({NEL, NPOL * NMO}, alloc_);
    ComplexType Ovlp = SDet.MixedDensityMatrix(devPsiT[0], devOrbMat[0], G.sliced(0, NAEA), 0.0, true);
    if (WTYPE == COLLINEAR)
    {
      Ovlp *= SDet.MixedDensityMatrix(devPsiT[1], devOrbMat[1](devOrbMat.extension(1), {0, NAEB}),
                                      G.sliced(NAEA, NAEA + NAEB), 0.0, true);
    }
    CHECK(real(Ovlp) == Approx(1.0));
    CHECK(imag(Ovlp) == Approx(0.0));

    boost::multi::array<ComplexType, 2, Alloc> Eloc({1, 3}, alloc_);
    {
      int nc = 1, nr = NEL * NPOL * NMO;
      if (HOps.transposed_G_for_E())
      {
        nr = 1;
        nc = NEL * NPOL * NMO;
      }
      boost::multi::array_ref<ComplexType, 2, pointer> Gw(make_device_ptr(G.origin()), {nr, nc});
      HOps.energy(Eloc, Gw, 0, TG.getCoreID() == 0);
    }
    Eloc[0][0] = (TG.Node() += ComplexType(Eloc[0][0]));
    Eloc[0][1] = (TG.Node() += ComplexType(Eloc[0][1]));
    Eloc[0][2] = (TG.Node() += ComplexType(Eloc[0][2]));
    if (std::abs(file_data.E0 + file_data.E1) > 1e-8)
    {
      CHECK(real(Eloc[0][0]) == Approx(real(file_data.E0 + file_data.E1)));
      CHECK(imag(Eloc[0][0]) == Approx(imag(file_data.E0 + file_data.E1)));
    }
    else
    {
      app_log() << " E1: " << setprecision(12) << Eloc[0][0] << std::endl;
    }
    if (std::abs(file_data.E2) > 1e-8)
    {
      CHECK(real(Eloc[0][1] + Eloc[0][2]) == Approx(real(file_data.E2)));
      CHECK(imag(Eloc[0][1] + Eloc[0][2]) == Approx(imag(file_data.E2)));
    }
    else
    {
      app_log() << " EJ: " << setprecision(12) << Eloc[0][2] << std::endl;
      app_log() << " EXX: " << setprecision(12) << Eloc[0][1] << std::endl;
      app_log() << " ETotal: " << setprecision(12) << Eloc[0][0] + Eloc[0][1] + Eloc[0][2] << std::endl;
    }

    double sqrtdt = std::sqrt(0.01);
    auto nCV      = HOps.local_number_of_cholesky_vectors();

    CMatrix X({nCV, 1}, alloc_);
    {
      int nc = 1, nr = NEL * NPOL * NMO;
      if (HOps.transposed_G_for_vbias())
      {
        nr = 1;
        nc = NEL * NPOL * NMO;
      }
      boost::multi::array_ref<ComplexType, 2, pointer> Gw(make_device_ptr(G.origin()), {nr, nc});
      HOps.vbias(Gw, X, sqrtdt);
    }
    TG.local_barrier();
    ComplexType Xsum = 0, Xsum2 = 0;
    for (int i = 0; i < X.size(); i++)
    {
      Xsum += X[i][0];
      Xsum2 += ComplexType(0.5) * X[i][0] * X[i][0];
    }
    if (std::abs(file_data.Xsum) > 1e-8)
    {
      CHECK(real(Xsum) == Approx(real(file_data.Xsum)));
      CHECK(imag(Xsum) == Approx(imag(file_data.Xsum)));
    }
    else
    {
      app_log() << " Xsum: " << setprecision(12) << Xsum << std::endl;
      app_log() << " Xsum2 (EJ): " << setprecision(12) << Xsum2 / sqrtdt / sqrtdt << std::endl;
    }

    int vdim1 = (HOps.transposed_vHS() ? 1 : NMO * NMO);
    int vdim2 = (HOps.transposed_vHS() ? NMO * NMO : 1);
    CMatrix vHS({vdim1, vdim2}, alloc_);
    TG.local_barrier();
    HOps.vHS(X, vHS, sqrtdt);
    TG.local_barrier();
    ComplexType Vsum = 0;
    if (HOps.transposed_vHS())
    {
      for (int i = 0; i < std::get<1>(vHS.sizes()); i++)
        Vsum += vHS[0][i];
    }
    else
    {
      for (int i = 0; i < std::get<0>(vHS.sizes()); i++)
        Vsum += vHS[i][0];
    }
    if (std::abs(file_data.Vsum) > 1e-8)
    {
      CHECK(real(Vsum) == Approx(real(file_data.Vsum)));
      CHECK(imag(Vsum) == Approx(imag(file_data.Vsum)));
    }
    else
    {
      app_log() << " Vsum: " << setprecision(12) << Vsum << std::endl;
    }
    // Test Generalised Fock matrix.
    int dm_size;
    CMatrix G2({1, 1}, alloc_);
    dm_size = 2 * NMO * NMO;
    G2.reextent({2 * NMO, NMO});
    boost::multi::array<double, 2> Mat({NMO, NMO});
    hdf_archive ref;
    if (!ref.open("G.h5", H5F_ACC_RDONLY))
    {
      app_error() << " Error opening HDF wavefunction file.\n";
    }
    ref.readEntry(Mat, "Ga");
    for (int i = 0; i < NMO; i++)
    {
      for (int j = 0; j < NMO; j++)
      {
        G2[i][j] = Mat[i][j];
      }
    }
    ref.readEntry(Mat, "Gb");
    for (int i = 0; i < NMO; i++)
    {
      for (int j = 0; j < NMO; j++)
      {
        G2[NMO + i][j] = Mat[i][j];
      }
    }

    //if(WTYPE==COLLINEAR) {
    //dm_size = 2*NMO*NMO;
    //G2.reextent({2*NMO,NMO});
    //} else if(WTYPE==CLOSED) {
    //dm_size = NMO*NMO;
    //G2.reextent({NMO,NMO});
    //} else {
    //APP_ABORT("NON COLLINEAR Wavefunction not implemented.");
    //}
    //Ovlp = SDet.MixedDensityMatrix(devPsiT[0],devOrbMat[0], G2.sliced(0,NMO),0.0,false);
    //if(WTYPE==COLLINEAR) {
    //Ovlp *= SDet.MixedDensityMatrix(devPsiT[1],devOrbMat[1](devOrbMat.extension(1),{0,NAEB}), G.sliced(NMO,2*NMO),0.0,false);
    //}
    int nwalk = 1;
    CMatrix Gw2({nwalk, dm_size}, alloc_);
    for (int nw = 0; nw < nwalk; nw++)
    {
      for (int i = 0; i < NMO; i++)
      {
        for (int j = 0; j < NMO; j++)
        {
          Gw2[nw][i * NMO + j]             = G2[i][j];
          Gw2[nw][NMO * NMO + i * NMO + j] = G2[i + NMO][j];
        }
      }
    }
    boost::multi::array_ref<ComplexType, 2, pointer> Gw2_(make_device_ptr(Gw2.origin()), {nwalk, dm_size});
    boost::multi::array<ComplexType, 3, Alloc> GFock({2, nwalk, dm_size}, alloc_);
    fill_n(GFock.origin(), GFock.num_elements(), ComplexType(0.0));
    //for(int i = 0; i < nwalk; i++) {
    //std::cout << Gw2_[i][0] << std::endl;
    //}
    //std::cout << "INIT: " << Gw2_[0][0] << " " << Gw2_[0][NMO*NMO] << std::endl;
    HOps.generalizedFockMatrix(Gw2_, GFock[0], GFock[1]);
    //boost::multi::array_ref<ComplexType,2,pointer> GR(make_device_ptr(GFock[0][0].origin()), {NMO,NMO});
    //for(int i = 0; i < nwalk; i++) {
    //std::cout << GFock[0][i][0] << std::endl;
    //}
    //std::cout << "GFOCK: " << GFock[0][0][0] << " " << GFock[1][0][0] << std::endl;
    //std::cout << "Fm: " << std::endl;
    std::fill_n(Mat.origin(), Mat.num_elements(), 0.0);
    ref.readEntry(Mat, "Fha");
    for (int i = 0; i < NMO; i++)
    {
      for (int j = 0; j < NMO; j++)
      {
        if (std::abs(Mat[i][j] - real(GFock[1][0][i * NMO + j])) > 1e-5)
        {
          std::cout << "DELTAA: " << i << " " << j << " " << Mat[i][j] << " " << real(GFock[1][0][i * NMO + j])
                    << std::endl;
        }
        //if(std::abs(real(GFock[1][0][i*NMO+j]))>1e-6)
        //std::cout << i << " " << j << " " << real(GFock[1][0][i*NMO+j]) << " " << std::endl;
      }
    }
    //std::cout << "Fp: " << std::endl;
    std::fill_n(Mat.origin(), Mat.num_elements(), 0.0);
    ref.readEntry(Mat, "Fpa");
    for (int i = 0; i < NMO; i++)
    {
      for (int j = 0; j < NMO; j++)
      {
        //std::cout << Mat[i][j] << std::endl;
        //std::cout << Mat[i][j]-real(GFock[0][0][i*NMO+j]) << std::endl;
        if (std::abs(Mat[i][j] - real(GFock[0][0][i * NMO + j])) > 1e-5)
        {
          std::cout << "DELTAB: " << i << " " << j << " " << Mat[i][j] << " " << real(GFock[0][0][i * NMO + j])
                    << std::endl;
        }
        //if(std::abs(real(GFock[0][0][i*NMO+j]))>1e-6)
        //std::cout << i << " " << j << " " << real(GFock[0][0][i*NMO+j]) << " " << real(GFock[0][1][i*NMO+j]) << " " << real(GFock[0][2][i*NMO+j]) << std::endl;
      }
    }
  }
}

TEST_CASE("ham_ops_basic_serial", "[hamiltonian_operations]")
{
  auto world = boost::mpi3::environment::get_world_instance();
  auto node  = world.split_shared(world.rank());

#if defined(ENABLE_CUDA) || defined(ENABLE_HIP)

  arch::INIT(node);
  using Alloc = device::device_allocator<ComplexType>;
#else
  using Alloc = shared_allocator<ComplexType>;
#endif
  setup_memory_managers(node, 10uL * 1024uL * 1024uL);
  ham_ops_basic_serial<Alloc>(world);
  release_memory_managers();
}

} // namespace qmcplusplus
