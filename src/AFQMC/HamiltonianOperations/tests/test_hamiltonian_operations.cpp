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
#include "AFQMC/Utilities/myTimer.h"
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

TEST_CASE("ham_ops_basic_serial", "[hamiltonian_operations]")
{
  OHMMS::Controller->initialize(0, NULL);
  auto world = boost::mpi3::environment::get_world_instance();

  if(not file_exists("./afqmc.h5") ||
     not file_exists("./wfn.dat") ) {
    app_log()<<" Skipping ham_ops_basic_serial. afqmc.h5 and ./wfn.dat files not found. \n";
  } else {

    // Global Task Group
    afqmc::GlobalTaskGroup gTG(world);

    auto file_data = read_test_results_from_hdf<ValueType>("./afqmc.h5");
    int NMO=file_data.NMO;
    int NAEA=file_data.NAEA;
    int NAEB=file_data.NAEB;

    std::map<std::string,AFQMCInfo> InfoMap;
    InfoMap.insert ( std::pair<std::string,AFQMCInfo>("info0",AFQMCInfo{"info0",NMO,NAEA,NAEB}) );
    HamiltonianFactory HamFac(InfoMap);
    const char *xml_block =
"<Hamiltonian name=\"ham0\" type=\"SparseGeneral\" info=\"info0\"> \
    <parameter name=\"filetype\">hdf5</parameter> \
    <parameter name=\"version\">new</parameter> \
    <parameter name=\"filename\">./afqmc.h5</parameter> \
    <parameter name=\"cutoff_decomposition\">1e-5</parameter> \
  </Hamiltonian> \
";
    Libxml2Document doc;
    bool okay = doc.parseFromString(xml_block);
    REQUIRE(okay);
    std::string ham_name("ham0");
    HamFac.push(ham_name,doc.getRoot());

    Hamiltonian& ham = HamFac.getHamiltonian(gTG,ham_name);

    using shm_Alloc = boost::mpi3::intranode::allocator<ComplexType>;
    using shmCMatrix = ComplexMatrix<shared_allocator<ComplexType>>;
    boost::multi::array<ComplexType,3> OrbMat;
    int walker_type = readWfn(std::string("./wfn.dat"),OrbMat,NMO,NAEA,NAEB);
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
    auto HOps(ham.getHamiltonianOperations(false,true,WTYPE,PsiT,1e-6,1e-6,TG,TG,dummy));

    // Calculates Overlap, G
    //SlaterDetOperations SDet( SlaterDetOperations_shared<ComplexType>(NMO,NAEA) );
    SlaterDetOperations_shared<ComplexType> SDet(NMO,NAEA);

    shmCMatrix G({NEL,NMO},shared_allocator<ComplexType>{TG.TG_local()});
    ComplexType Ovlp;
    SDet.MixedDensityMatrix(PsiT[0],OrbMat[0],
        G.sliced(0,NAEA),std::addressof(Ovlp),true);
    if(WTYPE==COLLINEAR) {
      ComplexType Ovlp_;
      SDet.MixedDensityMatrix(PsiT[1],OrbMat[1](OrbMat.extension(1),{0,NAEB}),
        G.sliced(NAEA,NAEA+NAEB),std::addressof(Ovlp_),true);
      Ovlp *= Ovlp_; 
    }
    REQUIRE( real(Ovlp) == Approx(1.0) );
    REQUIRE( imag(Ovlp) == Approx(0.0) );

    boost::multi::array<ComplexType,2> Eloc({1,3});
    boost::multi::array_ref<ComplexType,2> Gw(std::addressof(*G.origin()),{NEL*NMO,1});
    HOps.energy(Eloc,Gw,0,TG.getCoreID()==0);
    Eloc[0][0] = ( TG.Node() += Eloc[0][0] );
    Eloc[0][1] = ( TG.Node() += Eloc[0][1] );
    Eloc[0][2] = ( TG.Node() += Eloc[0][2] );
    if(std::abs(file_data.E0+file_data.E1)>1e-8) {
      REQUIRE( real(Eloc[0][0]) == Approx(real(file_data.E0+file_data.E1)) );
      REQUIRE( imag(Eloc[0][0]) == Approx(imag(file_data.E0+file_data.E1)) );
    } else {
      app_log()<<" E1: " <<setprecision(12) <<Eloc[0][0] <<std::endl;
    }
    if(std::abs(file_data.E2)>1e-8) {
      REQUIRE( real(Eloc[0][1]) == Approx(real(file_data.E2)));
      REQUIRE( imag(Eloc[0][1]) == Approx(imag(file_data.E2)));
    } else {
      app_log()<<" EJ: " <<setprecision(12) <<Eloc[0][2] <<std::endl;
      app_log()<<" EXX: " <<setprecision(12) <<Eloc[0][1] <<std::endl;
    }

    double sqrtdt = std::sqrt(0.01);
    auto nCV = HOps.local_number_of_cholesky_vectors();

    shmCMatrix X({nCV,1},shared_allocator<ComplexType>{TG.TG_local()});
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

    shmCMatrix vHS({NMO*NMO,1},shared_allocator<ComplexType>{TG.TG_local()});
    TG.local_barrier();
    Timer.reset("Generic");
    Timer.start("Generic");
    HOps.vHS(X,vHS,sqrtdt);
    TG.local_barrier();
    Timer.stop("Generic");
    app_log()<<" Time in vHS: " <<Timer.total("Generic") <<std::endl;
    ComplexType Vsum=0;
    for(int i=0; i<vHS.size(); i++)
        Vsum += vHS[i][0];
    if(std::abs(file_data.Vsum)>1e-8) {
      REQUIRE( real(Vsum) == Approx(real(file_data.Vsum)) );
      REQUIRE( imag(Vsum) == Approx(imag(file_data.Vsum)) );
    } else {
      app_log()<<" Vsum: " <<setprecision(12) <<Vsum <<std::endl;
    }
  }
}

TEST_CASE("ham_ops_collinear_distributed", "[hamiltonian_operations]")
{
  OHMMS::Controller->initialize(0, NULL);
  auto world = boost::mpi3::environment::get_world_instance();

  if(not file_exists("./afqmc_collinear.h5") ||
     not file_exists("./wfn_collinear.dat") ) {
    app_log()<<" Skipping ham_ops_collinear_sdet text. afqmc.h5 and ./wfn.dat files not found. \n";
  } else {

    // Global Task Group
    afqmc::GlobalTaskGroup gTG(world);

    auto file_data = read_test_results_from_hdf<ValueType>("./afqmc.h5");
    int NMO=file_data.NMO;
    int NAEA=file_data.NAEA;
    int NAEB=file_data.NAEB;

    std::map<std::string,AFQMCInfo> InfoMap;
    InfoMap.insert ( std::pair<std::string,AFQMCInfo>("info0",AFQMCInfo{"info0",NMO,NAEA,NAEB}) );
    HamiltonianFactory HamFac(InfoMap);
    const char *xml_block =
"<Hamiltonian name=\"ham0\" type=\"SparseGeneral\" info=\"info0\"> \
    <parameter name=\"filetype\">hdf5</parameter> \
    <parameter name=\"version\">new</parameter> \
    <parameter name=\"filename\">./afqmc.h5</parameter> \
    <parameter name=\"cutoff_decomposition\">1e-5</parameter> \
  </Hamiltonian> \
";
    Libxml2Document doc;
    bool okay = doc.parseFromString(xml_block);
    REQUIRE(okay);
    std::string ham_name("ham0");
    HamFac.push(ham_name,doc.getRoot());

    Hamiltonian& ham = HamFac.getHamiltonian(gTG,ham_name);

    using shm_Alloc = boost::mpi3::intranode::allocator<ComplexType>;
    using shmCMatrix = ComplexMatrix<shared_allocator<ComplexType>>;
    boost::multi::array<ComplexType,3> OrbMat;
    int walker_type = readWfn(std::string("./wfn.dat"),OrbMat,NMO,NAEA,NAEB);
    int NEL = (walker_type==0)?NAEA:(NAEA+NAEB);
    WALKER_TYPES WTYPE = CLOSED;
    if(walker_type==1) WTYPE = COLLINEAR;
    if(walker_type==2) WTYPE = NONCOLLINEAR;

    std::vector<PsiT_Matrix> PsiT;
    PsiT.reserve(2);
    PsiT.emplace_back(csr::shm::construct_csr_matrix_single_input<PsiT_Matrix>(
                                        OrbMat[0],1e-8,'H',gTG.Node()));
    if(WTYPE==COLLINEAR)
      PsiT.emplace_back(csr::shm::construct_csr_matrix_single_input<PsiT_Matrix>(
                                        OrbMat[1](OrbMat.extension(1),{0,NAEB}),
                                        1e-8,'H',gTG.Node()));

    hdf_archive dummy;
    auto TG = TaskGroup_(gTG,std::string("DummyTG"),1,gTG.getTotalCores());
    auto HOps(ham.getHamiltonianOperations(false,true,WTYPE,PsiT,1e-6,1e-6,TG,TG,dummy));

    // Calculates Overlap, G
    //SlaterDetOperations SDet( SlaterDetOperations_shared<ComplexType>(NMO,NAEA) );
    SlaterDetOperations_shared<ComplexType> SDet(NMO,NAEA);

    shmCMatrix G({NEL,NMO},shared_allocator<ComplexType>{TG.TG_local()});
    ComplexType Ovlp;
    SDet.MixedDensityMatrix(PsiT[0],OrbMat[0],
        G.sliced(0,NAEA),std::addressof(Ovlp),true);
    if(WTYPE==COLLINEAR) {
      ComplexType Ovlp_;  
      SDet.MixedDensityMatrix(PsiT[1],OrbMat[1](OrbMat.extension(1),{0,NAEB}),
        G.sliced(NAEA,NAEA+NAEB),std::addressof(Ovlp_),true);
      Ovlp *= Ovlp_;
    }  
    REQUIRE( real(Ovlp) == Approx(1.0) );
    REQUIRE( imag(Ovlp) == Approx(0.0) );

    boost::multi::array<ComplexType,2> Eloc({1,3});
    boost::multi::array_ref<ComplexType,2> Gw(std::addressof(*G.origin()),{NEL*NMO,1});
    HOps.energy(Eloc,Gw,0,TG.getCoreID()==0);
    Eloc[0][0] = ( TG.Node() += Eloc[0][0] );
    Eloc[0][1] = ( TG.Node() += Eloc[0][1] );
    Eloc[0][2] = ( TG.Node() += Eloc[0][2] );
    if(std::abs(file_data.E0+file_data.E1)>1e-8) {
      REQUIRE( real(Eloc[0][0]) == Approx(real(file_data.E0+file_data.E1)) );
      REQUIRE( imag(Eloc[0][0]) == Approx(imag(file_data.E0+file_data.E1)) );
    } else {
      app_log()<<" E1: " <<setprecision(12) <<Eloc[0][0] <<std::endl;
    }
    if(std::abs(file_data.E2)>1e-8) {
      REQUIRE( real(Eloc[0][1]) == Approx(real(file_data.E2)));
      REQUIRE( imag(Eloc[0][1]) == Approx(imag(file_data.E2)));
    } else {
      app_log()<<" EJ: " <<setprecision(12) <<Eloc[0][2] <<std::endl;
      app_log()<<" EXX: " <<setprecision(12) <<Eloc[0][1] <<std::endl;
    }

    double sqrtdt = std::sqrt(0.01);
    auto nCV = HOps.local_number_of_cholesky_vectors();

    shmCMatrix X({nCV,1},shared_allocator<ComplexType>{TG.TG_local()});
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

    shmCMatrix vHS({NMO*NMO,1},shared_allocator<ComplexType>{TG.TG_local()});
    TG.local_barrier();
    Timer.reset("Generic");
    Timer.start("Generic");
    HOps.vHS(X,vHS,sqrtdt);
    TG.local_barrier();
    Timer.stop("Generic");
    app_log()<<" Time in vHS: " <<Timer.total("Generic") <<std::endl;
    ComplexType Vsum=0;
    for(int i=0; i<vHS.size(); i++)
        Vsum += vHS[i][0];
    if(std::abs(file_data.Vsum)>1e-8) {
      REQUIRE( real(Vsum) == Approx(real(file_data.Vsum)) );
      REQUIRE( imag(Vsum) == Approx(imag(file_data.Vsum)) );
    } else {
      app_log()<<" Vsum: " <<setprecision(12) <<Vsum <<std::endl;
    }
  }
}

TEST_CASE("test_thc_simple_serial", "[hamiltonian_operations]")
{
  OHMMS::Controller->initialize(0, NULL);
  auto world = boost::mpi3::environment::get_world_instance();

  if(not file_exists("./thc.h5") ||
     not file_exists("./wfn_thc.dat") ) {
    app_log()<<" Skipping test_thc_simple test. thc.h5 and ./wfn_thc.dat files not found. \n";
  } else {

    // Global Task Group
    afqmc::GlobalTaskGroup gTG(world);

    auto file_data = read_test_results_from_hdf<ValueType>("./thc.h5");
    int NMO=file_data.NMO;
    int NAEA=file_data.NAEA;
    int NAEB=file_data.NAEB;

    std::map<std::string,AFQMCInfo> InfoMap;
    InfoMap.insert ( std::pair<std::string,AFQMCInfo>("info0",AFQMCInfo{"info0",NMO,NAEA,NAEB}) );
    HamiltonianFactory HamFac(InfoMap);
    const char *xml_block =
"<Hamiltonian name=\"ham0\" type=\"THC\" info=\"info0\"> \
    <parameter name=\"filetype\">hdf5</parameter> \
    <parameter name=\"filename\">./thc.h5</parameter> \
    <parameter name=\"cutoff_decomposition\">1e-5</parameter> \
  </Hamiltonian> \
";
    Libxml2Document doc;
    bool okay = doc.parseFromString(xml_block);
    REQUIRE(okay);
    std::string ham_name("ham0");
    HamFac.push(ham_name,doc.getRoot());

    Hamiltonian& ham = HamFac.getHamiltonian(gTG,ham_name);

    using shm_Alloc = boost::mpi3::intranode::allocator<ComplexType>;
    using shmCMatrix = ComplexMatrix<shared_allocator<ComplexType>>;
    boost::multi::array<ComplexType,3> OrbMat;
    int walker_type = readWfn(std::string("./wfn_thc.dat"),OrbMat,NMO,NAEA,NAEB);
    int NEL = (walker_type==0)?NAEA:(NAEA+NAEB);
    WALKER_TYPES WTYPE = CLOSED;
    if(walker_type==1) WTYPE = COLLINEAR;
    if(walker_type==2) WTYPE = NONCOLLINEAR;

    std::vector<PsiT_Matrix> PsiT;
    PsiT.reserve(2);
    PsiT.emplace_back(csr::shm::construct_csr_matrix_single_input<PsiT_Matrix>(
                                        OrbMat[0],1e-8,'H',gTG.Node()));
    if(WTYPE==COLLINEAR)
      PsiT.emplace_back(csr::shm::construct_csr_matrix_single_input<PsiT_Matrix>(
                                        OrbMat[1](OrbMat.extension(1),{0,NAEB}),
                                        1e-8,'H',gTG.Node()));

    hdf_archive dummy;
    auto TG = TaskGroup_(gTG,std::string("DummyTG"),1,1);
    auto HOps(ham.getHamiltonianOperations(false,false,WTYPE,PsiT,1e-6,1e-6,TG,TG,dummy));

    // Calculates Overlap, G
    //SlaterDetOperations SDet( SlaterDetOperations_shared<ComplexType>(NMO,NAEA) );
    SlaterDetOperations_shared<ComplexType> SDet(NMO,NAEA);

    int nw=1;

    shmCMatrix Gbuff({nw,NEL*NMO},shared_allocator<ComplexType>{TG.TG_local()});
    boost::multi::array_ref<ComplexType,2> G(std::addressof(*Gbuff.origin()),{NEL,NMO});
    ComplexType Ovlp;
    SDet.MixedDensityMatrix(PsiT[0],OrbMat[0],G,std::addressof(Ovlp),true);
    REQUIRE( real(Ovlp) == Approx(1.0) );
    REQUIRE( imag(Ovlp) == Approx(0.0) );

    boost::multi::array_ref<ComplexType,2> Gw(std::addressof(*Gbuff.origin()),{nw,NEL*NMO});
    boost::multi::array<ComplexType,2> Eloc({nw,3});

    if(TG.Node().root())
      for(int i=1; i<nw; i++)
        Gw[i] = Gw[0];

    Timer.reset("Generic");
    Timer.start("Generic");
    HOps.energy(Eloc,Gw,0,true);
    TG.local_barrier();
    Timer.stop("Generic");
    app_log()<<" Time in energy: " <<Timer.total("Generic") <<std::endl;
    if(std::abs(file_data.E0+file_data.E1)>1e-8) {
      REQUIRE( real(Eloc[0][0]) == Approx(real(file_data.E0+file_data.E1)) );
      REQUIRE( imag(Eloc[0][0]) == Approx(imag(file_data.E0+file_data.E1)) );
    } else {
      app_log()<<" E1: " <<setprecision(12) <<Eloc[0][0] <<std::endl;
    }
    if(std::abs(file_data.E2)>1e-8) {
      auto E_ = Eloc[0][1]+Eloc[0][2];
      REQUIRE( real(E_) == Approx(real(file_data.E2)) );
      REQUIRE( imag(E_) == Approx(imag(file_data.E2)) );
    } else {
      app_log()<<" EXX, EJ: " <<setprecision(12) <<Eloc[0][1] <<" " <<Eloc[0][2] <<std::endl;
    }

    double sqrtdt = std::sqrt(0.01);
    auto nCV = HOps.local_number_of_cholesky_vectors();

    shmCMatrix X({nCV,nw},shared_allocator<ComplexType>{TG.TG_local()});
    Timer.reset("Generic");
    Timer.start("Generic");
    HOps.vbias(Gw,X,sqrtdt);
    TG.local_barrier();
    Timer.stop("Generic");
    app_log()<<" Time in vbias: " <<Timer.total("Generic") <<std::endl;
    ComplexType Xsum=0;
    for(int i=0; i<X.size(); i++)
        Xsum += X[i][0];
    if(std::abs(file_data.Xsum)>1e-8) {
      REQUIRE( real(Xsum) == Approx(real(file_data.Xsum)) );
      REQUIRE( imag(Xsum) == Approx(imag(file_data.Xsum)) );
    } else {
      app_log()<<" Xsum: " <<setprecision(12) <<Xsum <<std::endl;
    }

    shmCMatrix vHS({nw,NMO*NMO},shared_allocator<ComplexType>{TG.TG_local()});
    // doing twice to get reasonable timing estimate
    HOps.vHS(X,vHS,sqrtdt);
    TG.local_barrier();
    Timer.reset("Generic");
    Timer.start("Generic");
    HOps.vHS(X,vHS,sqrtdt);
    TG.local_barrier();
    Timer.stop("Generic");
    app_log()<<" Time in vHS: " <<Timer.total("Generic") <<std::endl;
    ComplexType Vsum=0;
    for(int i=0; i<vHS.shape()[1]; i++)
        Vsum += vHS[0][i];
    if(std::abs(file_data.Vsum)>1e-8) {
      REQUIRE( real(Vsum) == Approx(real(file_data.Vsum)) );
      REQUIRE( imag(Vsum) == Approx(imag(file_data.Vsum)) );
    } else {
      app_log()<<" Vsum: " <<setprecision(12) <<Vsum <<std::endl;
    }


  }

}

TEST_CASE("test_thc_simple_shared", "[hamiltonian_operations]")
{
  OHMMS::Controller->initialize(0, NULL);
  auto world = boost::mpi3::environment::get_world_instance();

  if(not file_exists("./thc.h5") ||
     not file_exists("./wfn_thc.dat") ) {
    app_log()<<" Skipping test_thc_simple test. thc.h5 and ./wfn_thc.dat files not found. \n";
  } else {

    // Global Task Group
    afqmc::GlobalTaskGroup gTG(world);

    auto file_data = read_test_results_from_hdf<ValueType>("./thc.h5");
    int NMO=file_data.NMO;
    int NAEA=file_data.NAEA;
    int NAEB=file_data.NAEB;

    std::map<std::string,AFQMCInfo> InfoMap;
    InfoMap.insert ( std::pair<std::string,AFQMCInfo>("info0",AFQMCInfo{"info0",NMO,NAEA,NAEB}) );
    HamiltonianFactory HamFac(InfoMap);
    const char *xml_block =
"<Hamiltonian name=\"ham0\" type=\"THC\" info=\"info0\"> \
    <parameter name=\"filetype\">hdf5</parameter> \
    <parameter name=\"filename\">./thc.h5</parameter> \
  </Hamiltonian> \
";
    Libxml2Document doc;
    bool okay = doc.parseFromString(xml_block);
    REQUIRE(okay);
    std::string ham_name("ham0");
    HamFac.push(ham_name,doc.getRoot());

    Hamiltonian& ham = HamFac.getHamiltonian(gTG,ham_name);

    using shm_Alloc = boost::mpi3::intranode::allocator<ComplexType>;
    using shmCMatrix = ComplexMatrix<shared_allocator<ComplexType>>;
    boost::multi::array<ComplexType,3> OrbMat;
    int walker_type = readWfn(std::string("./wfn_thc.dat"),OrbMat,NMO,NAEA,NAEB);
    int NEL = (walker_type==0)?NAEA:(NAEA+NAEB);
    WALKER_TYPES WTYPE = CLOSED;
    if(walker_type==1) WTYPE = COLLINEAR;
    if(walker_type==2) WTYPE = NONCOLLINEAR;

    std::vector<PsiT_Matrix> PsiT;
    PsiT.reserve(2);
    PsiT.emplace_back(csr::shm::construct_csr_matrix_single_input<PsiT_Matrix>(
                                        OrbMat[0],1e-8,'H',gTG.Node()));
    if(WTYPE==COLLINEAR)
      PsiT.emplace_back(csr::shm::construct_csr_matrix_single_input<PsiT_Matrix>(
                                        OrbMat[1](OrbMat.extension(1),{0,NAEB}),
                                        1e-8,'H',gTG.Node()));

    int ncores = gTG.getTotalCores();
    if(file_exists("ncores.txt")) {
      ifstream in("ncores.txt");
      in>>ncores;
      in.close();
    }

    hdf_archive dummy;
    auto TG = TaskGroup_(gTG,std::string("DummyTG"),1,ncores);
    auto HOps(ham.getHamiltonianOperations(false,false,WTYPE,PsiT,1e-6,1e-6,TG,TG,dummy));

    // Calculates Overlap, G
    //SlaterDetOperations SDet( SlaterDetOperations_shared<ComplexType>(NMO,NAEA) );
    SlaterDetOperations_shared<ComplexType> SDet(NMO,NAEA);

    shmCMatrix G({NEL,NMO},shared_allocator<ComplexType>{TG.TG_local()});
    ComplexType Ovlp;
    SDet.MixedDensityMatrix(PsiT[0],OrbMat[0],G,std::addressof(Ovlp),true);
    REQUIRE( real(Ovlp) == Approx(1.0) );
    REQUIRE( imag(Ovlp) == Approx(0.0) );

    boost::multi::array_ref<ComplexType,2> Gw(std::addressof(*G.origin()),{1,NEL*NMO});
    boost::multi::array<ComplexType,2> Eloc({1,3});
    Timer.reset("Generic");
    Timer.start("Generic");
    HOps.energy(Eloc,Gw,0,TG.getCoreID()==0);
    TG.local_barrier();
    Timer.stop("Generic");
    app_log()<<" Time in energy: " <<Timer.total("Generic") <<std::endl;
    Eloc[0][0] = ( TG.Node() += Eloc[0][0] );
    Eloc[0][1] = ( TG.Node() += Eloc[0][1] );
    Eloc[0][2] = ( TG.Node() += Eloc[0][2] );
    if(std::abs(file_data.E0+file_data.E1)>1e-8) {
      REQUIRE( real(Eloc[0][0]) == Approx(real(file_data.E0+file_data.E1)) );
      REQUIRE( imag(Eloc[0][0]) == Approx(imag(file_data.E0+file_data.E1)) );
    } else {
      app_log()<<" E1: " <<setprecision(12) <<Eloc[0][0] <<std::endl;
    }
    if(std::abs(file_data.E2)>1e-8) {
      auto E_ = Eloc[0][1]+Eloc[0][2];
      REQUIRE( real(E_) == Approx(real(file_data.E2)) );
      REQUIRE( imag(E_) == Approx(imag(file_data.E2)) );
    } else {
      app_log()<<" EXX, EJ: " <<setprecision(12) <<Eloc[0][1] <<" " <<Eloc[0][2] <<std::endl;
    }

    double sqrtdt = std::sqrt(0.01);
    auto nCV = HOps.local_number_of_cholesky_vectors();

    shmCMatrix X({nCV,1},shared_allocator<ComplexType>{TG.TG_local()});
    Timer.reset("Generic");
    Timer.start("Generic");
    HOps.vbias(Gw,X,sqrtdt);
    TG.local_barrier();
    Timer.stop("Generic");
    app_log()<<" Time in vbias: " <<Timer.total("Generic") <<std::endl;
    ComplexType Xsum=0,X2sum=0;
    for(int i=0; i<X.size(); i++)
        Xsum += X[i][0];
    for(int i=0; i<X.size(); i++)
        X2sum += X[i][0]*X[i][0];
    if(std::abs(file_data.Xsum)>1e-8) {
      REQUIRE( real(Xsum) == Approx(real(file_data.Xsum)) );
      REQUIRE( imag(Xsum) == Approx(imag(file_data.Xsum)) );
    } else {
      app_log()<<" Xsum, EJ: " <<setprecision(12) <<Xsum <<" " <<X2sum/0.01/2.0 <<std::endl;
    }

    shmCMatrix vHS({1,NMO*NMO},shared_allocator<ComplexType>{TG.TG_local()});
    // doing twice to get reasonable timing estimate
    HOps.vHS(X,vHS,sqrtdt);
    TG.local_barrier();
    Timer.reset("Generic");
    Timer.start("Generic");
    HOps.vHS(X,vHS,sqrtdt);
    TG.local_barrier();
    Timer.stop("Generic");
    app_log()<<" Time in vHS: " <<Timer.total("Generic") <<std::endl;
    ComplexType Vsum=0;
    for(int i=0; i<vHS.shape()[1]; i++)
        Vsum += vHS[0][i];
    if(std::abs(file_data.Vsum)>1e-8) {
      REQUIRE( real(Vsum) == Approx(real(file_data.Vsum)) );
      REQUIRE( imag(Vsum) == Approx(imag(file_data.Vsum)) );
    } else {
      app_log()<<" Vsum: " <<setprecision(12) <<Vsum <<std::endl;
    }


  }

}

TEST_CASE("test_thc_shared_testLuv", "[hamiltonian_operations]")
{
  OHMMS::Controller->initialize(0, NULL);
  auto world = boost::mpi3::environment::get_world_instance();

  if(not file_exists("./thc.h5") ||
     not file_exists("./wfn_thc.dat") ) {
    app_log()<<" Skipping test_thc_simple test. thc.h5 and ./wfn.dat files not found. \n";
  } else {

    // Global Task Group
    afqmc::GlobalTaskGroup gTG(world);

    auto file_data = read_test_results_from_hdf<ValueType>("./thc.h5");
    int NMO=file_data.NMO;
    int NAEA=file_data.NAEA;
    int NAEB=file_data.NAEB;

    std::map<std::string,AFQMCInfo> InfoMap;
    InfoMap.insert ( std::pair<std::string,AFQMCInfo>("info0",AFQMCInfo{"info0",NMO,NAEA,NAEB}) );
    HamiltonianFactory HamFac(InfoMap);
    const char *xml_block =
"<Hamiltonian name=\"ham0\" type=\"THC\" info=\"info0\"> \
    <parameter name=\"filetype\">hdf5</parameter> \
    <parameter name=\"filename\">./thc.h5</parameter> \
    <parameter name=\"useHalfRotatedMuv\">no</parameter> \
  </Hamiltonian> \
";
    Libxml2Document doc;
    bool okay = doc.parseFromString(xml_block);
    REQUIRE(okay);
    std::string ham_name("ham0");
    HamFac.push(ham_name,doc.getRoot());

    Hamiltonian& ham = HamFac.getHamiltonian(gTG,ham_name);

    using shm_Alloc = boost::mpi3::intranode::allocator<ComplexType>;
    using shmCMatrix = ComplexMatrix<shared_allocator<ComplexType>>;
    boost::multi::array<ComplexType,3> OrbMat;
    int walker_type = readWfn(std::string("./wfn_thc.dat"),OrbMat,NMO,NAEA,NAEB);
    int NEL = (walker_type==0)?NAEA:(NAEA+NAEB);
    WALKER_TYPES WTYPE = CLOSED;
    if(walker_type==1) WTYPE = COLLINEAR;
    if(walker_type==2) WTYPE = NONCOLLINEAR;

    std::vector<PsiT_Matrix> PsiT;
    PsiT.reserve(2);
    PsiT.emplace_back(csr::shm::construct_csr_matrix_single_input<PsiT_Matrix>(
                                        OrbMat[0],1e-8,'H',gTG.Node()));
    if(WTYPE==COLLINEAR)
      PsiT.emplace_back(csr::shm::construct_csr_matrix_single_input<PsiT_Matrix>(
                                        OrbMat[1](OrbMat.extension(1),{0,NAEB}),
                                        1e-8,'H',gTG.Node()));

    int ncores = gTG.getTotalCores();
    if(file_exists("ncores.txt")) {
      ifstream in("ncores.txt");
      in>>ncores;
      in.close();
    }

    hdf_archive dummy;
    auto TG = TaskGroup_(gTG,std::string("DummyTG"),1,ncores);
    // NOTE: This will force the replacement of HalfRotatedLuv by Luv to test the energy of the
    //       non-rotated factorization
    THCHamiltonian& thcHam = boost::get<THCHamiltonian>(ham);
    auto HOps(thcHam.getHamiltonianOperations(false,false,WTYPE,PsiT,1e-6,1e-6,TG,TG,dummy));

    // Calculates Overlap, G
    //SlaterDetOperations SDet( SlaterDetOperations_shared<ComplexType>(NMO,NAEA) );
    SlaterDetOperations_shared<ComplexType> SDet(NMO,NAEA);

    shmCMatrix G({NEL,NMO},shared_allocator<ComplexType>{TG.TG_local()});
    ComplexType Ovlp;
    SDet.MixedDensityMatrix(PsiT[0],OrbMat[0],G,std::addressof(Ovlp),true);
    REQUIRE( real(Ovlp) == Approx(1.0) );
    REQUIRE( imag(Ovlp) == Approx(0.0) );

    boost::multi::array_ref<ComplexType,2> Gw(std::addressof(*G.origin()),{1,NEL*NMO});
    boost::multi::array<ComplexType,2> Eloc({1,3});
    Timer.reset("Generic");
    Timer.start("Generic");
    HOps.energy(Eloc,Gw,0,TG.getCoreID()==0);
    TG.local_barrier();
    Timer.stop("Generic");
    app_log()<<" Time in energy: " <<Timer.total("Generic") <<std::endl;
    Eloc[0][0] = ( TG.Node() += Eloc[0][0] );
    Eloc[0][1] = ( TG.Node() += Eloc[0][1] );
    Eloc[0][2] = ( TG.Node() += Eloc[0][2] );
    if(std::abs(file_data.E0+file_data.E1)>1e-8) {
      REQUIRE( real(Eloc[0][0]) == Approx(real(file_data.E0+file_data.E1)) );
      REQUIRE( imag(Eloc[0][0]) == Approx(imag(file_data.E0+file_data.E1)) );
    } else {
      app_log()<<" E1: " <<setprecision(12) <<Eloc[0][0] <<std::endl;
    }
    if(std::abs(file_data.E2)>1e-8) {
      auto E_ = Eloc[0][1]+Eloc[0][2];
      REQUIRE( real(E_) == Approx(real(file_data.E2)) );
      REQUIRE( imag(E_) == Approx(imag(file_data.E2)) );
    } else {
      app_log()<<" EXX, EJ: " <<setprecision(12) <<Eloc[0][1] <<" " <<Eloc[0][2] <<std::endl;
    }

    double sqrtdt = std::sqrt(0.01);
    auto nCV = HOps.local_number_of_cholesky_vectors();

    shmCMatrix X({nCV,1},shared_allocator<ComplexType>{TG.TG_local()});
    Timer.reset("Generic");
    Timer.start("Generic");
    HOps.vbias(Gw,X,sqrtdt);
    TG.local_barrier();
    Timer.stop("Generic");
    app_log()<<" Time in vbias: " <<Timer.total("Generic") <<std::endl;
    ComplexType Xsum=0,X2sum=0;
    for(int i=0; i<X.size(); i++)
        Xsum += X[i][0];
    for(int i=0; i<X.size(); i++)
        X2sum += X[i][0]*X[i][0];
    if(std::abs(file_data.Xsum)>1e-8) {
      REQUIRE( real(Xsum) == Approx(real(file_data.Xsum)) );
      REQUIRE( imag(Xsum) == Approx(imag(file_data.Xsum)) );
    } else {
      app_log()<<" Xsum, EJ: " <<setprecision(12) <<Xsum <<" " <<X2sum/0.01/2.0 <<std::endl;
    }

    shmCMatrix vHS({1,NMO*NMO},shared_allocator<ComplexType>{TG.TG_local()});
    // doing twice to get reasonable timing estimate
    HOps.vHS(X,vHS,sqrtdt);
    TG.local_barrier();
    Timer.reset("Generic");
    Timer.start("Generic");
    HOps.vHS(X,vHS,sqrtdt);
    TG.local_barrier();
    Timer.stop("Generic");
    app_log()<<" Time in vHS: " <<Timer.total("Generic") <<std::endl;
    ComplexType Vsum=0;
    for(int i=0; i<vHS.shape()[1]; i++)
        Vsum += vHS[0][i];
    if(std::abs(file_data.Vsum)>1e-8) {
      REQUIRE( real(Vsum) == Approx(real(file_data.Vsum)) );
      REQUIRE( imag(Vsum) == Approx(imag(file_data.Vsum)) );
    } else {
      app_log()<<" Vsum: " <<setprecision(12) <<Vsum <<std::endl;
    }


  }

}

}
