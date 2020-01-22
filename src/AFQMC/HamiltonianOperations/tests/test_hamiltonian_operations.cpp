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

#include "Message/catch_mpi_main.hpp"

#include "Configuration.h"

#include "OhmmsData/Libxml2Doc.h"
#include "OhmmsApp/ProjectData.h"
#include "io/hdf_archive.h"

#undef APP_ABORT
#define APP_ABORT(x) {std::cout << x <<std::endl; exit(0);}

#include <sys/stat.h>
#include <unistd.h>
#include <stdio.h>
#include <string>
#include <vector>
#include <complex>
#include <iomanip>

#include "AFQMC/config.h"
#include "AFQMC/Hamiltonians/HamiltonianFactory.h"
#include "AFQMC/Hamiltonians/Hamiltonian.hpp"
#include "AFQMC/Hamiltonians/THCHamiltonian.h"
#include "AFQMC/Utilities/readWfn.h"
#include "AFQMC/SlaterDeterminantOperations/SlaterDetOperations.hpp"
#include "AFQMC/Utilities/test_utils.hpp"

#include "AFQMC/Matrix/csr_matrix_construct.hpp"
#include "AFQMC/Numerics/ma_blas.hpp"

using std::string;
using std::complex;
using std::cout;
using std::cerr;
using std::endl;
using std::ifstream;
using std::setprecision;

namespace qmcplusplus
{

using namespace afqmc;


template<class Alloc>
void ham_ops_basic_serial(boost::mpi3::communicator & world)
{

  using pointer = device_ptr<ComplexType>;

  if(not file_exists(UTEST_HAMIL) ||
     not file_exists(UTEST_WFN) ) {
    app_log()<<" Skipping ham_ops_basic_serial. Hamiltonian or wavefunction file not found. \n";
    app_log()<<" Run unit test with --hamil /path/to/hamil.h5 and --wfn /path/to/wfn.dat.\n";
  } else {

    // Global Task Group
    afqmc::GlobalTaskGroup gTG(world);

    // Determine wavefunction type for test results from wavefunction file name which is
    // has the naming convention wfn_(wfn_type).dat.
    // First strip path of filename.
    std::string base_name = UTEST_WFN.substr(UTEST_WFN.find_last_of("\\/")+1);
    // Remove file extension.
    std::string test_wfn = base_name.substr(0, base_name.find_last_of("."));
    auto file_data = read_test_results_from_hdf<ValueType>(UTEST_HAMIL, test_wfn);
    int NMO=file_data.NMO;
    int NAEA=file_data.NAEA;
    int NAEB=file_data.NAEB;

    std::map<std::string,AFQMCInfo> InfoMap;
    InfoMap.insert ( std::pair<std::string,AFQMCInfo>("info0",AFQMCInfo{"info0",NMO,NAEA,NAEB}) );
    HamiltonianFactory HamFac(InfoMap);
    std::string hamil_xml =
"<Hamiltonian name=\"ham0\" info=\"info0\"> \
    <parameter name=\"filetype\">hdf5</parameter> \
    <parameter name=\"filename\">"+UTEST_HAMIL+"</parameter> \
    <parameter name=\"cutoff_decomposition\">1e-5</parameter> \
  </Hamiltonian> \
";
    const char *xml_block = hamil_xml.c_str();
    Libxml2Document doc;
    bool okay = doc.parseFromString(xml_block);
    REQUIRE(okay);
    std::string ham_name("ham0");
    HamFac.push(ham_name,doc.getRoot());

    Hamiltonian& ham = HamFac.getHamiltonian(gTG,ham_name);

    using CMatrix = ComplexMatrix<Alloc>;
    boost::multi::array<ComplexType,3> OrbMat;
    int walker_type = readWfn(std::string(UTEST_WFN),OrbMat,NMO,NAEA,NAEB);
    int NEL = (walker_type==0)?NAEA:(NAEA+NAEB);
    WALKER_TYPES WTYPE = CLOSED;
    if(walker_type==1) WTYPE = COLLINEAR;
    if(walker_type==2) WTYPE = NONCOLLINEAR;

    std::vector<PsiT_Matrix> PsiT;
    PsiT.reserve(2);
    PsiT.emplace_back(csr::shm::construct_csr_matrix_single_input<PsiT_Matrix>(
                                        OrbMat[0],1e-8,'H',gTG.Node()));
    if(walker_type>0)
      PsiT.emplace_back(csr::shm::construct_csr_matrix_single_input<PsiT_Matrix>(
                                        OrbMat[1](OrbMat.extension(1),{0,NAEB}),
                                        1e-8,'H',gTG.Node()));

    hdf_archive dummy;
    auto TG = TaskGroup_(gTG,std::string("DummyTG"),1,gTG.getTotalCores());
    auto HOps(ham.getHamiltonianOperations(false,false,WTYPE,PsiT,1e-6,1e-6,TG,TG,dummy));

    // Calculates Overlap, G
// NOTE: Make small factory routine!
    //SlaterDetOperations SDet( SlaterDetOperations_serial<Alloc>(NMO,NAEA) );
#ifdef ENABLE_CUDA
    auto SDet( SlaterDetOperations_serial<Alloc>(NMO,NAEA) );
#else
    auto SDet( SlaterDetOperations_shared<ComplexType>(NMO,NAEA) );
#endif

    Alloc alloc_(make_localTG_allocator<ComplexType>(TG));
    boost::multi::array<ComplexType,3,Alloc> devOrbMat(OrbMat, alloc_);
    std::vector<devcsr_Matrix> devPsiT(move_vector<devcsr_Matrix>(std::move(PsiT)));

    CMatrix G({NEL,NMO},alloc_);
    ComplexType Ovlp = SDet.MixedDensityMatrix(devPsiT[0],devOrbMat[0],
        G.sliced(0,NAEA),0.0,true);
    if(WTYPE==COLLINEAR) {
      Ovlp *= SDet.MixedDensityMatrix(devPsiT[1],devOrbMat[1](devOrbMat.extension(1),{0,NAEB}),
        G.sliced(NAEA,NAEA+NAEB),0.0,true);
    }
    REQUIRE( real(Ovlp) == Approx(1.0) );
    REQUIRE( imag(Ovlp) == Approx(0.0) );

    boost::multi::array<ComplexType,2,Alloc> Eloc({1,3},alloc_);
    boost::multi::array_ref<ComplexType,2,pointer> Gw(make_device_ptr(G.origin()),{NEL*NMO,1});
    HOps.energy(Eloc,Gw,0,TG.getCoreID()==0);
    Eloc[0][0] = ( TG.Node() += ComplexType(Eloc[0][0]) );
    Eloc[0][1] = ( TG.Node() += ComplexType(Eloc[0][1]) );
    Eloc[0][2] = ( TG.Node() += ComplexType(Eloc[0][2]) );
    if(std::abs(file_data.E0+file_data.E1)>1e-8) {
      REQUIRE( real(Eloc[0][0]) == Approx(real(file_data.E0+file_data.E1)) );
      REQUIRE( imag(Eloc[0][0]) == Approx(imag(file_data.E0+file_data.E1)) );
    } else {
      app_log()<<" E1: " <<setprecision(12) <<Eloc[0][0] <<std::endl;
    }
    if(std::abs(file_data.E2)>1e-8) {
      REQUIRE( real(Eloc[0][1]+Eloc[0][2]) == Approx(real(file_data.E2)));
      REQUIRE( imag(Eloc[0][1]+Eloc[0][2]) == Approx(imag(file_data.E2)));
    } else {
      app_log()<<" EJ: " <<setprecision(12) <<Eloc[0][2] <<std::endl;
      app_log()<<" EXX: " <<setprecision(12) <<Eloc[0][1] <<std::endl;
    }

    double sqrtdt = std::sqrt(0.01);
    auto nCV = HOps.local_number_of_cholesky_vectors();

    CMatrix X({nCV,1},alloc_);
    HOps.vbias(Gw,X,sqrtdt);
    TG.local_barrier();
    ComplexType Xsum=0;
    for(int i=0; i<X.size(); i++)
        Xsum += X[i][0];
    if(std::abs(file_data.Xsum)>1e-8) {
      REQUIRE( real(Xsum) == Approx(real(file_data.Xsum)) );
      REQUIRE( imag(Xsum) == Approx(imag(file_data.Xsum)) );
    } else {
      app_log()<<" Xsum: " <<setprecision(12) <<Xsum <<std::endl;
    }

    int vdim1 = (HOps.transposed_vHS()?1:NMO*NMO);
    int vdim2 = (HOps.transposed_vHS()?NMO*NMO:1);
    CMatrix vHS({vdim1,vdim2},alloc_);
    TG.local_barrier();
    HOps.vHS(X,vHS,sqrtdt);
    TG.local_barrier();
    ComplexType Vsum=0;
    if(HOps.transposed_vHS()) {
      for(int i=0; i<vHS.size(1); i++)
        Vsum += vHS[0][i];
    } else {
      for(int i=0; i<vHS.size(0); i++)
        Vsum += vHS[i][0];
    }
    if(std::abs(file_data.Vsum)>1e-8) {
      REQUIRE( real(Vsum) == Approx(real(file_data.Vsum)) );
      REQUIRE( imag(Vsum) == Approx(imag(file_data.Vsum)) );
    } else {
      app_log()<<" Vsum: " <<setprecision(12) <<Vsum <<std::endl;
    }
    // Test Generalised Fock matrix.
    int dm_size;
    if(WTYPE==COLLINEAR) {
      dm_size = 2*NMO*NMO;
      CMatrix G2({2*NMO,NMO},alloc_);
    } else if(WTYPE==CLOSED) {
      dm_size = NMO*NMO;
      CMatrix G2({NMO,NMO},alloc_);
    } else {
      APP_ABORT("NON COLLINEAR Wavefunction not implemented.");
    }
    Ovlp = SDet.MixedDensityMatrix(devPsiT[0],devOrbMat[0], G2,0.0,false);
    if(WTYPE==COLLINEAR) {
      Ovlp *= SDet.MixedDensityMatrix(devPsiT[1],devOrbMat[1](devOrbMat.extension(1),{0,NAEB}), G.sliced(),0.0,false);
    }
    boost::multi::array_ref<ComplexType,2,pointer> Gw2(make_device_ptr(G2.origin()),{1,dm_size});
    CMatrix GFock({2,1,dm_size},alloc_);
    HOps.generalizedFockMatrix(Gw2,GFock[0],GFock[1]);
  }
}

TEST_CASE("ham_ops_basic_serial", "[hamiltonian_operations]")
{
  OHMMS::Controller->initialize(0, NULL);
  auto world = boost::mpi3::environment::get_world_instance();

#ifdef ENABLE_CUDA
  auto node = world.split_shared(world.rank());

  qmc_cuda::CUDA_INIT(node);
  using Alloc = qmc_cuda::cuda_gpu_allocator<ComplexType>;
#else
  using Alloc = shared_allocator<ComplexType>;
#endif

  ham_ops_basic_serial<Alloc>(world);
}

}
