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

#include "AFQMC/Hamiltonians/HamiltonianFactory.h"
#include "AFQMC/Hamiltonians/Hamiltonian.hpp"
#include "AFQMC/Utilities/myTimer.h"
#include "AFQMC/Utilities/readWfn.h"
#include "AFQMC/SlaterDeterminantOperations/SlaterDetOperations.hpp"
#include "AFQMC/Utilities/test_utils.hpp"

#include "AFQMC/Matrix/csr_matrix_construct.hpp"
#include "AFQMC/Numerics/ma_blas.hpp"
#include "AFQMC/Matrix/mpi3_SHMBuffer.hpp"

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

TEST_CASE("ham_factory_factorized_closed_pure", "[hamiltonian_factory]")
{
  OHMMS::Controller->initialize(0, NULL);
  auto world = boost::mpi3::environment::get_world_instance();

  if(not file_exists("./afqmc.h5") ||
     not file_exists("./wfn.dat") ) {
    app_log()<<" Skipping ham_factory_factorized_closed_pure text. afqmc.h5 or wfn.dat files not found. \n";
  } else {

    // Global Task Group
    afqmc::GlobalTaskGroup gTG(world);

    auto file_data = read_test_results_from_hdf<ValueType>("./afqmc.h5");
    int NMO=file_data.NMO;
    int NAEA=file_data.NAEA;
    int NAEB=file_data.NAEB;
    REQUIRE(NAEA==NAEB);

    std::map<std::string,AFQMCInfo> InfoMap;
    InfoMap.insert ( std::pair<std::string,AFQMCInfo>("info0",AFQMCInfo{"info0",NMO,NAEA,NAEB}) );
    HamiltonianFactory HamFac(InfoMap);
    const char *xml_block =
"<Hamiltonian name=\"ham0\" type=\"SparseGeneral\" info=\"info0\"> \
    <parameter name=\"filetype\">hdf5</parameter> \
    <parameter name=\"filename\">./afqmc.h5</parameter> \
    <parameter name=\"cutoff_decomposition\">1e-5</parameter> \
  </Hamiltonian> \
";
    Libxml2Document doc;
    bool okay = doc.parseFromString(xml_block);
    REQUIRE(okay);
    std::string ham_name("ham0");
    HamFac.push(ham_name,doc.getRoot());

    Hamiltonian& ham_ = HamFac.getHamiltonian(gTG,ham_name);
    FactorizedSparseHamiltonian& ham = boost::get<FactorizedSparseHamiltonian>(ham_);

    // build HamiltonianOperations
    if(file_exists("./wfn.dat")) {

        using shm_Alloc = boost::mpi3::intranode::allocator<ComplexType>;
        boost::multi::array<ComplexType,3> OrbMat;
        int wfn_type = readWfn(std::string("./wfn.dat"),OrbMat,NMO,NAEA,NAEB);
        if(wfn_type != 0)
            APP_ABORT(" Error: ham_factory_factorized_closed_pure expects a RHF wave function.\n");

        PsiT_Matrix TrialWfn = csr::shm::construct_csr_matrix_single_input<PsiT_Matrix>(
                                        OrbMat[0],1e-8,'H',gTG.Node());

        // Calculates Overlap, G
        SlaterDetOperations<ComplexType> SDet(NMO,NAEA);

        boost::multi::array<ComplexType,2> G({NAEA,NMO});
        boost::multi::array_ref<ComplexType,1> G0(G.origin(),extensions<1u>{NMO*NAEA});
        auto Ovlp = SDet.MixedDensityMatrix(TrialWfn,OrbMat[0],G,true);
        REQUIRE( real(Ovlp*Ovlp) == Approx(1.0) );
        REQUIRE( imag(Ovlp*Ovlp) == Approx(0.0) );

        //assert(gTG.getTotalNodes()>1);
        boost::multi::array<ComplexType,1> hij = ham.halfRotatedHij(CLOSED,
                                                                   std::addressof(TrialWfn),
                                                                   std::addressof(TrialWfn));
        ComplexType E1 = ma::dot(hij,G0);

        // this will be slightly dependent on cutoff, so be carefull
        if(std::abs(file_data.E1)>1e-8) {
          REQUIRE( real(E1) == Approx(real(file_data.E1)) );
          REQUIRE( imag(E1) == Approx(imag(file_data.E1)) );
        } else {
          app_log()<<" E1: " <<setprecision(12) <<E1 <<std::endl;
        }

        using array_ = std::array<std::size_t,4>;
        auto TG = TaskGroup_(gTG,std::string("DummyTG"),1,gTG.getTotalCores());
        std::size_t zero(0);

        std::map<IndexType,std::pair<bool,IndexType>> occ_a;
        for(int i=0; i<NMO; i++) occ_a[i] = {false,0};
        auto nel = TrialWfn.shape()[0];
        REQUIRE( nel == NAEA );
        for(std::size_t i=0; i<nel; i++) {
          if(TrialWfn.num_non_zero_elements(i) != 1)
            APP_ABORT(" Error: TrialWfn is not of pureSD form: Too many nnz. \n");
          auto val = TrialWfn.non_zero_values_data(i)[0];
          if(std::abs(val-ValueType(1.0)) > 1e-8)
            APP_ABORT(" Error: TrialWfn is not of pureSD form: Coeff != 1.0. \n");
          auto orb = TrialWfn.non_zero_indices2_data(i)[0];
          occ_a[orb] = {true,IndexType(i)};
        }

        // V2 uses std::size_t to store pointers_begin/end.
        // TODO (FDM): addCoulomb=false not yet implemented.
        bool addCoulomb = true;
        auto V2(ham.generateHijkl(CLOSED,addCoulomb,TG,occ_a,occ_a,1e-5));
        REQUIRE(V2.shape()[0] == NAEA*NMO);
        REQUIRE(V2.shape()[0] == V2.shape()[1]);

        // to get rid of narrowing warnings
        auto V2view(V2[array_{zero,V2.shape()[0],zero,V2.shape()[1]}]);
        REQUIRE(V2view.shape()[0] == NAEA*NMO);
        REQUIRE(V2view.shape()[0] == V2view.shape()[1]);
        boost::multi::array<ComplexType,1> V0(extensions<1u>{V2.shape()[0]});
        ma::product(V2view,G0,V0);
        ComplexType E2 = 0.5*ma::dot(G0,V0);
        if(std::abs(file_data.E2)>1e-8) {
          //REQUIRE( real(E2) == Approx(real(file_data.E2)));
          //REQUIRE( imag(E2) == Approx(imag(file_data.E2)));
        } else {
          app_log()<<" E2: " <<setprecision(12) <<E2 <<std::endl;
        }

        // Spvn
        double sqrtdt = std::sqrt(0.01);
        boost::multi::array<ComplexType,2> vn0({NMO,NMO});
        auto Spvn(ham.calculateHSPotentials(1e-6,TG,vn0));
        auto SpvnT(csr::shm::transpose(Spvn));
        auto rotSpvnT(sparse_rotate::halfRotateCholeskyMatrixForBias(CLOSED,TG,
                                std::addressof(TrialWfn),
                                std::addressof(TrialWfn),
                                Spvn,1e-6));

        auto Spvnview(Spvn[array_{zero,Spvn.shape()[0],zero,Spvn.shape()[1]}]);
        auto SpvnTview(SpvnT[array_{zero,SpvnT.shape()[0],zero,SpvnT.shape()[1]}]);
        auto rotSpvnTview(rotSpvnT[array_{zero,rotSpvnT.shape()[0],zero,rotSpvnT.shape()[1]}]);

        boost::multi::array<ComplexType,3> GM({1,NMO,NMO});
        boost::multi::array_ref<ComplexType,1> G0M(GM.origin(),extensions<1u>{NMO*NMO});
        Ovlp = SDet.MixedDensityMatrix(TrialWfn,OrbMat[0],
            GM[0],false);
        REQUIRE( real(Ovlp*Ovlp) == Approx(1.0) );
        REQUIRE( imag(Ovlp*Ovlp) == Approx(0.0) );

        boost::multi::array<ComplexType,1> X(extensions<1u>{Spvn.shape()[1]});
        using ma::T;
        ma::product(2.0*sqrtdt,T(Spvnview),
                    G0M.sliced(0,NMO*NMO),0.,X);
        ComplexType Xsum=0;
        for(int i=0; i<X.size(); i++)
            Xsum += X[i];
        if(std::abs(file_data.Xsum)>1e-8) {
          REQUIRE( real(Xsum) == Approx(real(file_data.Xsum)) );
          REQUIRE( imag(Xsum) == Approx(imag(file_data.Xsum)) );
        } else {
          app_log()<<" Xsum: " <<setprecision(12) <<Xsum <<std::endl;
        }

        ma::product(2.0*sqrtdt,rotSpvnTview,G0,0.,X);
        Xsum=0;
        for(int i=0; i<X.size(); i++)
            Xsum += X[i];
        if(std::abs(file_data.Xsum)>1e-8) {
          REQUIRE( real(Xsum) == Approx(real(file_data.Xsum)) );
          REQUIRE( imag(Xsum) == Approx(imag(file_data.Xsum)) );
        } else {
          app_log()<<" Xsum: " <<setprecision(12) <<Xsum <<std::endl;
        }

        ma::product(2.0*sqrtdt,SpvnTview,
                    G0M.sliced(0,NMO*NMO),0.,X);
        Xsum=0;
        for(int i=0; i<X.size(); i++)
            Xsum += X[i];
        if(std::abs(file_data.Xsum)>1e-8) {
          REQUIRE( real(Xsum) == Approx(real(file_data.Xsum)) );
          REQUIRE( imag(Xsum) == Approx(imag(file_data.Xsum)) );
        } else {
          app_log()<<" Xsum: " <<setprecision(12) <<Xsum <<std::endl;
        }

        boost::multi::array<ComplexType,1> vHS(extensions<1u>{NMO*NMO});
        ma::product(sqrtdt,Spvnview,X,0.,vHS);
        ComplexType Vsum=0;
        for(int i=0; i<vHS.size(); i++)
            Vsum += vHS[i];
        if(std::abs(file_data.Vsum)>1e-8) {
          REQUIRE( real(Vsum) == Approx(real(file_data.Vsum)) );
          REQUIRE( imag(Vsum) == Approx(imag(file_data.Vsum)) );
        } else {
          app_log()<<" Vsum: " <<setprecision(12) <<Vsum <<std::endl;
        }

    } else {
      app_error()<<" Skipping test. Missing wfn file: wfn.dat. \n";
    }

  }

}

TEST_CASE("ham_factory_factorized_collinear_with_rotation", "[hamiltonian_factory]")
{
  OHMMS::Controller->initialize(0, NULL);
  auto world = boost::mpi3::environment::get_world_instance();

  if(not file_exists("./afqmc_collinear.h5") ||
     not file_exists("./wfn_collinear.dat") ) {
    app_log()<<" Skipping ham_factory_factorized_collinear_with_rotation text. afqmc_collinear.h5 or wfn_collinear.dat files not found. \n";
  } else {

    // Global Task Group
    afqmc::GlobalTaskGroup gTG(world);

    auto file_data = read_test_results_from_hdf<ValueType>("./afqmc_collinear.h5");
    int NMO=file_data.NMO;
    int NAEA=file_data.NAEA;
    int NAEB=file_data.NAEB;

    std::map<std::string,AFQMCInfo> InfoMap;
    InfoMap.insert ( std::pair<std::string,AFQMCInfo>("info0",AFQMCInfo{"info0",NMO,NAEA,NAEB}) );
    HamiltonianFactory HamFac(InfoMap);
    const char *xml_block =
"<Hamiltonian name=\"ham0\" type=\"SparseGeneral\" info=\"info0\"> \
    <parameter name=\"filetype\">hdf5</parameter> \
    <parameter name=\"filename\">./afqmc_collinear.h5</parameter> \
    <parameter name=\"cutoff_decomposition\">1e-5</parameter> \
  </Hamiltonian> \
";
    Libxml2Document doc;
    bool okay = doc.parseFromString(xml_block);
    REQUIRE(okay);
    std::string ham_name("ham0");
    HamFac.push(ham_name,doc.getRoot());

    Hamiltonian& ham_ = HamFac.getHamiltonian(gTG,ham_name);
    FactorizedSparseHamiltonian& ham = boost::get<FactorizedSparseHamiltonian>(ham_);

    // build HamiltonianOperations
    if(file_exists("./wfn_collinear.dat")) {

        using shm_Alloc = boost::mpi3::intranode::allocator<ComplexType>;
        boost::multi::array<ComplexType,3> OrbMat;
        readWfn(std::string("./wfn_collinear.dat"),OrbMat,NMO,NAEA,NAEB);

        std::pair<PsiT_Matrix,PsiT_Matrix> TrialWfn =
                std::make_pair(csr::shm::construct_csr_matrix_single_input<PsiT_Matrix>(
                                        OrbMat[0],1e-8,'H',gTG.Node()),
                               csr::shm::construct_csr_matrix_single_input<PsiT_Matrix>(
                                        OrbMat[1],1e-8,'H',gTG.Node())
                              );

        // Calculates Overlap, G
        SlaterDetOperations<ComplexType> SDet(NMO,NAEA);

        boost::multi::array<ComplexType,2> G({NAEA+NAEB,NMO});
        boost::multi::array_ref<ComplexType,1> G0(G.origin(),extensions<1u>{NMO*(NAEA+NAEB)});
        auto Ovlp = SDet.MixedDensityMatrix(TrialWfn.first,OrbMat[0],
            G.sliced(0,NAEA),true);
        Ovlp *= SDet.MixedDensityMatrix(TrialWfn.second,OrbMat[1],
            G.sliced(NAEA,NAEA+NAEB),true);
        REQUIRE( real(Ovlp) == Approx(1.0) );
        REQUIRE( imag(Ovlp) == Approx(0.0) );

        //assert(gTG.getTotalNodes()>1);
        boost::multi::array<ComplexType,1> hij = ham.halfRotatedHij(COLLINEAR,
                                                                   std::addressof(TrialWfn.first),
                                                                   std::addressof(TrialWfn.second));
        ComplexType E1 = ma::dot(hij,G0);

        // this will be slightly dependent on cutoff, so be carefull
        if(std::abs(file_data.E1)>1e-8) {
          REQUIRE( real(E1) == Approx(real(file_data.E1)) );
          REQUIRE( imag(E1) == Approx(imag(file_data.E1)) );
        } else {
          app_log()<<" E1: " <<setprecision(12) <<E1 <<std::endl;
        }

        using array_ = std::array<std::size_t,4>;
        auto TG = TaskGroup_(gTG,std::string("DummyTG"),1,gTG.getTotalCores());
        std::size_t zero(0);

        // V2 uses std::size_t to store pointers_begin/end.
        bool addCoulomb = true;
        auto V2(ham.halfRotatedHijkl(COLLINEAR,addCoulomb,TG,std::addressof(TrialWfn.first),
                                          std::addressof(TrialWfn.second),1e-5));
        REQUIRE(V2.shape()[0] == (NAEA+NAEB)*NMO);
        REQUIRE(V2.shape()[0] == V2.shape()[1]);

        // to get rid of narrowing warnings
        auto V2view(V2[array_{zero,V2.shape()[0],zero,V2.shape()[1]}]);
        REQUIRE(V2view.shape()[0] == (NAEA+NAEB)*NMO);
        REQUIRE(V2view.shape()[0] == V2view.shape()[1]);
        boost::multi::array<ComplexType,1> V0(extensions<1u>{V2.shape()[0]});
        ma::product(V2view,G0,V0);
        ComplexType E2 = 0.5*ma::dot(G0,V0);
        if(std::abs(file_data.E2)>1e-8) {
          REQUIRE( real(E2) == Approx(real(file_data.E2)));
          REQUIRE( imag(E2) == Approx(imag(file_data.E2)));
        } else {
          app_log()<<" E2: " <<setprecision(12) <<E2 <<std::endl;
        }

        // Spvn
        double sqrtdt = std::sqrt(0.01);
        boost::multi::array<ComplexType,2> vn0({NMO,NMO});
        auto Spvn(ham.calculateHSPotentials(1e-6,TG,vn0));
        auto SpvnT(csr::shm::transpose(Spvn));
        auto rotSpvnT(sparse_rotate::halfRotateCholeskyMatrixForBias(COLLINEAR,TG,
                                std::addressof(TrialWfn.first),
                                std::addressof(TrialWfn.second),
                                Spvn,1e-6));

        auto Spvnview(Spvn[array_{zero,Spvn.shape()[0],zero,Spvn.shape()[1]}]);
        auto SpvnTview(SpvnT[array_{zero,SpvnT.shape()[0],zero,SpvnT.shape()[1]}]);
        auto rotSpvnTview(rotSpvnT[array_{zero,rotSpvnT.shape()[0],zero,rotSpvnT.shape()[1]}]);

        boost::multi::array<ComplexType,3> GM({2,NMO,NMO});
        boost::multi::array_ref<ComplexType,1> G0M(GM.origin(),extensions<1u>{2*NMO*NMO});
        Ovlp = SDet.MixedDensityMatrix(TrialWfn.first,OrbMat[0],
            GM[0],false);
        Ovlp *= SDet.MixedDensityMatrix(TrialWfn.second,OrbMat[1],
            GM[1],false);
        REQUIRE( real(Ovlp) == Approx(1.0) );
        REQUIRE( imag(Ovlp) == Approx(0.0) );

        boost::multi::array<ComplexType,1> X(extensions<1u>{Spvn.shape()[1]});
        using ma::T;
        ma::product(sqrtdt,T(Spvnview),
                    G0M.sliced(0,NMO*NMO),0.,X);
        ma::product(sqrtdt,T(Spvnview),
                    G0M.sliced(NMO*NMO,2*NMO*NMO),1.,X);
        ComplexType Xsum=0;
        for(int i=0; i<X.size(); i++)
            Xsum += X[i];
        if(std::abs(file_data.Xsum)>1e-8) {
          REQUIRE( real(Xsum) == Approx(real(file_data.Xsum)) );
          REQUIRE( imag(Xsum) == Approx(imag(file_data.Xsum)) );
        } else {
          app_log()<<" Xsum: " <<setprecision(12) <<Xsum <<std::endl;
        }

        ma::product(sqrtdt,rotSpvnTview,G0,0.,X);
        Xsum=0;
        for(int i=0; i<X.size(); i++)
            Xsum += X[i];
        if(std::abs(file_data.Xsum)>1e-8) {
          REQUIRE( real(Xsum) == Approx(real(file_data.Xsum)) );
          REQUIRE( imag(Xsum) == Approx(imag(file_data.Xsum)) );
        } else {
          app_log()<<" Xsum: " <<setprecision(12) <<Xsum <<std::endl;
        }

        ma::product(sqrtdt,SpvnTview,
                    G0M.sliced(0,NMO*NMO),0.,X);
        ma::product(sqrtdt,SpvnTview,
                    G0M.sliced(NMO*NMO,2*NMO*NMO),1.,X);
        Xsum=0;
        for(int i=0; i<X.size(); i++)
            Xsum += X[i];
        if(std::abs(file_data.Xsum)>1e-8) {
          REQUIRE( real(Xsum) == Approx(real(file_data.Xsum)) );
          REQUIRE( imag(Xsum) == Approx(imag(file_data.Xsum)) );
        } else {
          app_log()<<" Xsum: " <<setprecision(12) <<Xsum <<std::endl;
        }

        boost::multi::array<ComplexType,1> vHS(extensions<1u>{NMO*NMO});
        ma::product(sqrtdt,Spvnview,X,0.,vHS);
        ComplexType Vsum=0;
        for(int i=0; i<vHS.size(); i++)
            Vsum += vHS[i];
        if(std::abs(file_data.Vsum)>1e-8) {
          REQUIRE( real(Vsum) == Approx(real(file_data.Vsum)) );
          REQUIRE( imag(Vsum) == Approx(imag(file_data.Vsum)) );
        } else {
          app_log()<<" Vsum: " <<setprecision(12) <<Vsum <<std::endl;
        }

    } else {
      app_error()<<" Skipping test of HamiltonianOperations. Missing wfn file: wfn_collinear.dat. \n";
    }

  }

}

TEST_CASE("ham_factory_dist_ham_factorized_collinear_with_rotation", "[hamiltonian_factory]")
{
  OHMMS::Controller->initialize(0, NULL);
  auto world = boost::mpi3::environment::get_world_instance();

  if(not file_exists("./afqmc_collinear.h5") ||
     not file_exists("./wfn_collinear.dat") ) {
    app_log()<<" Skipping ham_factory_factorized_collinear_with_rotation text. afqmc_collinear.h5 or wfn_collinear.dat files not found. \n";
  } else {

    // Global Task Group
    afqmc::GlobalTaskGroup gTG(world);

    auto file_data = read_test_results_from_hdf<ValueType>("./afqmc_collinear.h5");
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
    <parameter name=\"filename\">./afqmc_collinear.h5</parameter> \
    <parameter name=\"cutoff_decomposition\">1e-5</parameter> \
  </Hamiltonian> \
";
    Libxml2Document doc;
    bool okay = doc.parseFromString(xml_block);
    REQUIRE(okay);
    std::string ham_name("ham0");
    HamFac.push(ham_name,doc.getRoot());

    Hamiltonian& ham_ = HamFac.getHamiltonian(gTG,ham_name);
    FactorizedSparseHamiltonian& ham = boost::get<FactorizedSparseHamiltonian>(ham_);

    // build HamiltonianOperations
    if(file_exists("./wfn_collinear.dat")) {

        // new
        using shm_Alloc = boost::mpi3::intranode::allocator<ComplexType>;
        boost::multi::array<ComplexType,3> OrbMat;
        readWfn(std::string("./wfn_collinear.dat"),OrbMat,NMO,NAEA,NAEB);

        std::pair<PsiT_Matrix,PsiT_Matrix> TrialWfn =
                std::make_pair(csr::shm::construct_csr_matrix_single_input<PsiT_Matrix>(
                                        OrbMat[0],1e-8,'H',gTG.Node()),
                               csr::shm::construct_csr_matrix_single_input<PsiT_Matrix>(
                                        OrbMat[1],1e-8,'H',gTG.Node())
                              );

        // Calculates Overlap, G
        SlaterDetOperations<ComplexType> SDet(NMO,NAEA);

        boost::multi::array<ComplexType,2> G({NAEA+NAEB,NMO});
        boost::multi::array_ref<ComplexType,1> G0(G.origin(),extensions<1u>{NMO*(NAEA+NAEB)});
        auto Ovlp = SDet.MixedDensityMatrix(TrialWfn.first,OrbMat[0],
            G.sliced(0,NAEA),true);
        Ovlp *= SDet.MixedDensityMatrix(TrialWfn.second,OrbMat[1],
            G.sliced(NAEA,NAEA+NAEB),true);
        REQUIRE( real(Ovlp) == Approx(1.0) );
        REQUIRE( imag(Ovlp) == Approx(0.0) );

        //assert(gTG.getTotalNodes()>1);
        boost::multi::array<ComplexType,1> hij = ham.halfRotatedHij(COLLINEAR,
                                                                   std::addressof(TrialWfn.first),
                                                                   std::addressof(TrialWfn.second));
        ComplexType E1 = ma::dot(hij,G0);

        // this will be slightly dependent on cutoff, so be carefull
        if(std::abs(file_data.E1)>1e-8) {
          REQUIRE( real(E1) == Approx(real(file_data.E1)) );
          REQUIRE( imag(E1) == Approx(imag(file_data.E1)) );
        } else {
          app_log()<<" E1: " <<setprecision(12) <<E1 <<std::endl;
        }

        using array_ = std::array<std::size_t,4>;
        auto TG = TaskGroup_(gTG,std::string("DummyTG"),gTG.getTotalNodes(),gTG.getTotalCores());
        std::size_t zero(0);

        // V2 uses std::size_t to store pointers_begin/end.
        bool addCoulomb = true;
        auto V2(ham.halfRotatedHijkl(COLLINEAR,true,TG,std::addressof(TrialWfn.first),
                                          std::addressof(TrialWfn.second),1e-5));
        REQUIRE(V2.shape()[0] == (NAEA+NAEB)*NMO);
        REQUIRE(V2.shape()[0] == V2.shape()[1]);

        // to get rid of narrowing warnings
        auto V2view(V2[array_{zero,V2.shape()[0],zero,V2.shape()[1]}]);
        REQUIRE(V2view.shape()[0] == (NAEA+NAEB)*NMO);
        REQUIRE(V2view.shape()[0] == V2view.shape()[1]);
        boost::multi::array<ComplexType,1> V0(extensions<1u>{V2.shape()[0]});
        ma::product(V2view,G0,V0);
        ComplexType E2 = 0.5*ma::dot(G0,V0);
        E2 = ( TG.Cores() += E2 );
        if(std::abs(file_data.E2)>1e-8) {
          REQUIRE( real(E2) == Approx(real(file_data.E2)));
          REQUIRE( imag(E2) == Approx(imag(file_data.E2)));
        } else {
          app_log()<<" E2: " <<setprecision(12) <<E2 <<std::endl;
        }

        // Spvn
        double sqrtdt = std::sqrt(0.01);
        boost::multi::array<ComplexType,2> vn0({NMO,NMO});
        auto Spvn(ham.calculateHSPotentials(1e-6,TG,vn0));
        auto SpvnT(csr::shm::transpose(Spvn));
        auto rotSpvnT(sparse_rotate::halfRotateCholeskyMatrixForBias(COLLINEAR,TG,
                                std::addressof(TrialWfn.first),
                                std::addressof(TrialWfn.second),
                                Spvn,1e-6));

        auto Spvnview(Spvn[array_{zero,Spvn.shape()[0],zero,Spvn.shape()[1]}]);
        auto SpvnTview(SpvnT[array_{zero,SpvnT.shape()[0],zero,SpvnT.shape()[1]}]);
        auto rotSpvnTview(rotSpvnT[array_{zero,rotSpvnT.shape()[0],zero,rotSpvnT.shape()[1]}]);

        boost::multi::array<ComplexType,3> GM({2,NMO,NMO});
        boost::multi::array_ref<ComplexType,1> G0M(GM.origin(),extensions<1u>{2*NMO*NMO});
        Ovlp = SDet.MixedDensityMatrix(TrialWfn.first,OrbMat[0],
            GM[0],false);
        Ovlp *= SDet.MixedDensityMatrix(TrialWfn.second,OrbMat[1],
            GM[1],false);
        REQUIRE( real(Ovlp) == Approx(1.0) );
        REQUIRE( imag(Ovlp) == Approx(0.0) );

        boost::multi::array<ComplexType,1> X(extensions<1u>{Spvn.shape()[1]});
        using ma::T;
        ma::product(sqrtdt,T(Spvnview),
                    G0M.sliced(0,NMO*NMO),0.,X);
        ma::product(sqrtdt,T(Spvnview),
                    G0M.sliced(NMO*NMO,2*NMO*NMO),1.,X);
        ComplexType Xsum=0;
        for(int i=0; i<X.size(); i++)
            Xsum += X[i];
        Xsum = ( TG.Cores() += Xsum );
        if(std::abs(file_data.Xsum)>1e-8) {
          REQUIRE( real(Xsum) == Approx(real(file_data.Xsum)) );
          REQUIRE( imag(Xsum) == Approx(imag(file_data.Xsum)) );
        } else {
          app_log()<<" Xsum: " <<setprecision(12) <<Xsum <<std::endl;
        }

        ma::product(sqrtdt,rotSpvnTview,G0,0.,X);
        Xsum=0;
        for(int i=0; i<X.size(); i++)
            Xsum += X[i];
        Xsum = ( TG.Cores() += Xsum );
        if(std::abs(file_data.Xsum)>1e-8) {
          REQUIRE( real(Xsum) == Approx(real(file_data.Xsum)) );
          REQUIRE( imag(Xsum) == Approx(imag(file_data.Xsum)) );
        } else {
          app_log()<<" Xsum: " <<Xsum <<setprecision(12) <<std::endl;
        }

        ma::product(sqrtdt,SpvnTview,
                    G0M.sliced(0,NMO*NMO),0.,X);
        ma::product(sqrtdt,SpvnTview,
                    G0M.sliced(NMO*NMO,2*NMO*NMO),1.,X);
        Xsum=0;
        for(int i=0; i<X.size(); i++)
            Xsum += X[i];
        Xsum = ( TG.Cores() += Xsum );
        if(std::abs(file_data.Xsum)>1e-8) {
          REQUIRE( real(Xsum) == Approx(real(file_data.Xsum)) );
          REQUIRE( imag(Xsum) == Approx(imag(file_data.Xsum)) );
        } else {
          app_log()<<" Xsum: " <<setprecision(12) <<Xsum <<std::endl;
        }

        boost::multi::array<ComplexType,1> vHS(extensions<1u>{NMO*NMO});
        ma::product(sqrtdt,Spvnview,X,0.,vHS);
        ComplexType Vsum=0;
        for(int i=0; i<vHS.size(); i++)
            Vsum += vHS[i];
        Vsum = ( TG.Cores() += Vsum );
        if(std::abs(file_data.Vsum)>1e-8) {
          REQUIRE( real(Vsum) == Approx(real(file_data.Vsum)) );
          REQUIRE( imag(Vsum) == Approx(imag(file_data.Vsum)) );
        } else {
          app_log()<<" Vsum: " <<setprecision(12) <<Vsum <<std::endl;
        }

    } else {
      app_error()<<" Skipping test of HamiltonianOperations. Missing wfn file: wfn_collinear.dat. \n";
    }

  }

}

TEST_CASE("ham_generation_timing_hdf", "[hamiltonian_factory]")
{
  OHMMS::Controller->initialize(0, NULL);
  auto world = boost::mpi3::environment::get_world_instance();

  if(not file_exists("./afqmc_timing.h5")) {
    app_log()<<" Skipping ham_fac_timing text. afqmc_timing.h5 file not found. \n";
  } else {

    // Global Task Group
    afqmc::GlobalTaskGroup gTG(world);

    int NMO,NAEA,NAEB;
    std::tie(NMO,NAEA,NAEB) = read_info_from_hdf("./afqmc_timing.h5");

    std::map<std::string,AFQMCInfo> InfoMap;
    InfoMap.insert ( std::pair<std::string,AFQMCInfo>("info0",AFQMCInfo{"info0",NMO,NAEA,NAEB}) );
    HamiltonianFactory HamFac(InfoMap);

    const char *xml_block =
"<Hamiltonian name=\"ham0\" type=\"SparseGeneral\" info=\"info0\"> \
    <parameter name=\"filetype\">hdf5</parameter> \
    <parameter name=\"version\">new</parameter> \
    <parameter name=\"filename\">./afqmc_timing.h5</parameter> \
    <parameter name=\"cutoff_decomposition\">1e-5</parameter> \
  </Hamiltonian> \
";
    Libxml2Document doc;
    bool okay = doc.parseFromString(xml_block);
    REQUIRE(okay);
    std::string ham_name("ham0");
    HamFac.push(ham_name,doc.getRoot());

    myTimer Timer_;
    Timer_.start("GenTest0");
    Hamiltonian& ham = HamFac.getHamiltonian(gTG,ham_name);
    Timer_.stop("GenTest0");
    app_log()<<"\n*********************************************************************\n"
             <<" Time to create hamiltonian in ham_generation_timing_hdf: "
             <<Timer_.total("GenTest0")
             <<"\n*********************************************************************\n\n";

  }

}

}
