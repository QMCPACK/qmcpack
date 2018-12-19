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
#include "Utilities/RandomGenerator.h"
#include "Utilities/SimpleRandom.h"

#undef APP_ABORT
#define APP_ABORT(x) {std::cout << x <<std::endl; exit(0);}

#include <sys/stat.h>
#include <unistd.h>
#include <stdio.h>
#include <string>
#include <vector>
#include <complex>
#include <iomanip>

#include "AFQMC/Utilities/test_utils.hpp"

#include "boost/multi_array.hpp"
#include "AFQMC/Hamiltonians/HamiltonianFactory.h"
#include "AFQMC/Hamiltonians/Hamiltonian.hpp"
#include "AFQMC/Utilities/myTimer.h"
#include "AFQMC/SlaterDeterminantOperations/SlaterDetOperations.hpp"
#include "AFQMC/Wavefunctions/WavefunctionFactory.h"
#include "AFQMC/Walkers/WalkerSet.hpp"

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

using boost::extents;
using boost::indices;
using range_t = boost::multi_array_types::index_range;

namespace qmcplusplus
{


using namespace afqmc;

TEST_CASE("wfn_fac_sdet", "[wavefunction_factory]")
{

  if(not file_exists("./afqmc.h5") ||
     not file_exists("./wfn.dat") ) {
    app_log()<<" Skipping wfn_fac_collinear_sdet text. afqmc.h5 and ./wfn.dat files not found. \n";
  } else {

    // mpi3
    communicator& world = OHMMS::Controller->comm;

    // Global Task Group
    GlobalTaskGroup gTG(world);

    auto file_data = read_test_results_from_hdf<ValueType>("./afqmc.h5");
    int NMO=file_data.NMO;
    int NAEA=file_data.NAEA;
    int NAEB=file_data.NAEB;
    WALKER_TYPES type = afqmc::getWalkerType("wfn.dat");

    std::map<std::string,AFQMCInfo> InfoMap;
    InfoMap.insert ( std::pair<std::string,AFQMCInfo>("info0",AFQMCInfo{"info0",NMO,NAEA,NAEB}) );
    HamiltonianFactory HamFac(InfoMap);
    const char *ham_xml_block =
"<Hamiltonian name=\"ham0\" type=\"SparseGeneral\" info=\"info0\"> \
    <parameter name=\"filetype\">hdf5</parameter> \
    <parameter name=\"filename\">./afqmc.h5</parameter> \
    <parameter name=\"cutoff_decomposition\">1e-5</parameter> \
  </Hamiltonian> \
";
    Libxml2Document doc;
    bool okay = doc.parseFromString(ham_xml_block);
    REQUIRE(okay);
    std::string ham_name("ham0");
    HamFac.push(ham_name,doc.getRoot());
    Hamiltonian& ham = HamFac.getHamiltonian(gTG,ham_name);


    auto TG = TaskGroup_(gTG,std::string("WfnTG"),1,gTG.getTotalCores());
    int nwalk = 11; // choose prime number to force non-trivial splits in shared routines
    RandomGenerator_t rng;

const char *wlk_xml_block_closed =
"<WalkerSet name=\"wset0\">  \
  <parameter name=\"walker_type\">closed</parameter>  \
</WalkerSet> \
";
const char *wlk_xml_block_coll =
"<WalkerSet name=\"wset0\">  \
  <parameter name=\"walker_type\">collinear</parameter>  \
</WalkerSet> \
";
const char *wlk_xml_block_noncol =
"<WalkerSet name=\"wset0\">  \
  <parameter name=\"walker_type\">noncollinear</parameter>  \
</WalkerSet> \
";

    const char *wlk_xml_block = ( (type==CLOSED)?(wlk_xml_block_closed): 
                                  (type==COLLINEAR?wlk_xml_block_coll:wlk_xml_block_noncol) );     

    Libxml2Document doc3;
    okay = doc3.parseFromString(wlk_xml_block);
    REQUIRE(okay);

    const char *wfn_xml_block =
"<Wavefunction name=\"wfn0\" info=\"info0\"> \
      <parameter name=\"filetype\">ascii</parameter> \
      <parameter name=\"filename\">./wfn.dat</parameter> \
      <parameter name=\"cutoff\">1e-6</parameter> \
      <parameter name=\"restart_file\">dummy.h5</parameter> \
  </Wavefunction> \
";
    Libxml2Document doc2;
    okay = doc2.parseFromString(wfn_xml_block);
    REQUIRE(okay);
    std::string wfn_name("wfn0");
    WavefunctionFactory WfnFac(InfoMap); 
    WfnFac.push(wfn_name,doc2.getRoot());
    Wavefunction& wfn = WfnFac.getWavefunction(TG,TG,wfn_name,type,&ham,1e-6,nwalk); 

    WalkerSet wset(TG,doc3.getRoot(),InfoMap["info0"],&rng);
    auto initial_guess = WfnFac.getInitialGuess(wfn_name);
    REQUIRE(initial_guess.shape()[0]==2);
    REQUIRE(initial_guess.shape()[1]==NMO);
    REQUIRE(initial_guess.shape()[2]==NAEA);
    
    if(type == COLLINEAR) 
        wset.resize(nwalk,initial_guess[0],
                         initial_guess[1][indices[range_t()][range_t(0,NAEB)]]);
    else
        wset.resize(nwalk,initial_guess[0],
                         initial_guess[0]);

    wfn.Overlap(wset);
    for(auto it = wset.begin(); it!=wset.end(); ++it) {
      REQUIRE(real(it->overlap()) == Approx(1.0));
      REQUIRE(imag(it->overlap()) == Approx(0.0));
    }

    using shm_Alloc = boost::mpi3::intranode::allocator<ComplexType>;
    using SHM_Buffer = mpi3_SHMBuffer<ComplexType>;

    wfn.Energy(wset);
    if(std::abs(file_data.E0+file_data.E1+file_data.E2)>1e-8) {
      for(auto it = wset.begin(); it!=wset.end(); ++it) {
        REQUIRE( real(it->E1()) == Approx(real(file_data.E0+file_data.E1)));
        REQUIRE( real(it->EXX()+it->EJ()) == Approx(real(file_data.E2)));
        REQUIRE( imag(it->energy()) == Approx(imag(file_data.E0+file_data.E1+file_data.E2)));
      }
    } else {
      app_log()<<" E: " <<wset[0].energy() <<std::endl; 
      app_log()<<" E0+E1: " <<wset[0].E1() <<std::endl; 
      app_log()<<" EJ: " <<wset[0].EJ() <<std::endl; 
      app_log()<<" EXX: " <<wset[0].EXX() <<std::endl; 
    }

    auto size_of_G = wfn.size_of_G_for_vbias();
    SHM_Buffer Gbuff(TG.TG_local(),nwalk*size_of_G);
    int Gdim1 = (wfn.transposed_G_for_vbias()?nwalk:size_of_G);
    int Gdim2 = (wfn.transposed_G_for_vbias()?size_of_G:nwalk);
    boost::multi_array_ref<ComplexType,2> G(Gbuff.data(),extents[Gdim1][Gdim2]);
    wfn.MixedDensityMatrix_for_vbias(wset,G);

    double sqrtdt = std::sqrt(0.01);
    auto nCV = wfn.local_number_of_cholesky_vectors();
    SHM_Buffer Xbuff(TG.TG_local(),nCV*nwalk);
    boost::multi_array_ref<ComplexType,2> X(Xbuff.data(),extents[nCV][nwalk]);
    wfn.vbias(G,X,sqrtdt);
    ComplexType Xsum=0;
    if(std::abs(file_data.Xsum)>1e-8) { 
      for(int n=0; n<nwalk; n++) {
        Xsum=0;
        for(int i=0; i<X.shape()[0]; i++)
          Xsum += X[i][n];
        REQUIRE( real(Xsum) == Approx(real(file_data.Xsum)) );
        REQUIRE( imag(Xsum) == Approx(imag(file_data.Xsum)) );
      }
    } else {
      Xsum=0;
      for(int i=0; i<X.shape()[0]; i++)
        Xsum += X[i][0];
      app_log()<<" Xsum: " <<setprecision(12) <<Xsum <<std::endl;
    }    

    SHM_Buffer vHSbuff(TG.TG_local(),NMO*NMO*nwalk);
    int vdim1 = (wfn.transposed_vHS()?nwalk:NMO*NMO);
    int vdim2 = (wfn.transposed_vHS()?NMO*NMO:nwalk);
    boost::multi_array_ref<ComplexType,2> vHS(vHSbuff.data(),extents[vdim1][vdim2]);
    wfn.vHS(X,vHS,sqrtdt);
    TG.local_barrier();
    ComplexType Vsum=0;
    if(std::abs(file_data.Vsum)>1e-8) { 
      for(int n=0; n<nwalk; n++) {
        Vsum=0;
        for(int i=0; i<vHS.shape()[0]; i++)
          Vsum += vHS[i][n];
        REQUIRE( real(Vsum) == Approx(real(file_data.Vsum)) );
        REQUIRE( imag(Vsum) == Approx(imag(file_data.Vsum)) );
      }
    } else {
      Vsum=0;
      for(int i=0; i<vHS.shape()[0]; i++)
        Vsum += vHS[i][0];
      app_log()<<" Vsum: " <<setprecision(12) <<Vsum <<std::endl;
    }

  // Restarting Wavefunction from file
    const char *wfn_xml_block_restart =
"<Wavefunction name=\"wfn1\" info=\"info0\"> \
      <parameter name=\"filetype\">hdf5</parameter> \
      <parameter name=\"filename\">./dummy.h5</parameter> \
      <parameter name=\"cutoff\">1e-6</parameter> \
  </Wavefunction> \
";
    Libxml2Document doc4;
    okay = doc4.parseFromString(wfn_xml_block_restart);
    REQUIRE(okay);
    wfn_name = "wfn1";
    WfnFac.push(wfn_name,doc4.getRoot());
    Wavefunction& wfn2 = WfnFac.getWavefunction(TG,TG,wfn_name,type,nullptr,1e-6,nwalk);

    WalkerSet wset2(TG,doc3.getRoot(),InfoMap["info0"],&rng);
    //auto initial_guess = WfnFac.getInitialGuess(wfn_name);
    REQUIRE(initial_guess.shape()[0]==2);
    REQUIRE(initial_guess.shape()[1]==NMO);
    REQUIRE(initial_guess.shape()[2]==NAEA);

    if(type == COLLINEAR)
        wset2.resize(nwalk,initial_guess[0],
                         initial_guess[1][indices[range_t()][range_t(0,NAEB)]]);
    else
        wset2.resize(nwalk,initial_guess[0],
                         initial_guess[0]);

    wfn2.Overlap(wset2);
    for(auto it = wset2.begin(); it!=wset2.end(); ++it) {
      REQUIRE(real(it->overlap()) == Approx(1.0));
      REQUIRE(imag(it->overlap()) == Approx(0.0));
    }

    wfn2.Energy(wset2);
    if(std::abs(file_data.E0+file_data.E1+file_data.E2)>1e-8) {
      for(auto it = wset2.begin(); it!=wset2.end(); ++it) {
        REQUIRE( real(it->E1()) == Approx(real(file_data.E0+file_data.E1)));
        REQUIRE( real(it->EXX()+it->EJ()) == Approx(real(file_data.E2)));
        REQUIRE( imag(it->energy()) == Approx(imag(file_data.E0+file_data.E1+file_data.E2)));
      }
    } else {
      app_log()<<" E0+E1: " <<setprecision(12) <<wset[0].E1() <<std::endl; 
      app_log()<<" EJ: " <<setprecision(12) <<wset[0].EJ() <<std::endl; 
      app_log()<<" EXX: " <<setprecision(12) <<wset[0].EXX() <<std::endl; 
    }

    REQUIRE(size_of_G == wfn2.size_of_G_for_vbias());
    wfn2.MixedDensityMatrix_for_vbias(wset2,G);

    REQUIRE( nCV == wfn2.local_number_of_cholesky_vectors());
    wfn2.vbias(G,X,sqrtdt);
    Xsum=0;
    if(std::abs(file_data.Xsum)>1e-8) {
      for(int n=0; n<nwalk; n++) {
        Xsum=0;
        for(int i=0; i<X.shape()[0]; i++)
          Xsum += X[i][n];
        REQUIRE( real(Xsum) == Approx(real(file_data.Xsum)) );
        REQUIRE( imag(Xsum) == Approx(imag(file_data.Xsum)) );
      }
    } else {
      Xsum=0;
      for(int i=0; i<X.shape()[0]; i++)
        Xsum += X[i][0];
      app_log()<<" Xsum: " <<setprecision(12) <<Xsum <<std::endl;
    }

    wfn2.vHS(X,vHS,sqrtdt);
    TG.local_barrier();
    Vsum=0;
    if(std::abs(file_data.Vsum)>1e-8) {
      for(int n=0; n<nwalk; n++) {
        Vsum=0;
        for(int i=0; i<vHS.shape()[0]; i++)
          Vsum += vHS[i][n];
        REQUIRE( real(Vsum) == Approx(real(file_data.Vsum)) );
        REQUIRE( imag(Vsum) == Approx(imag(file_data.Vsum)) );
      }
    } else {
      Vsum=0;
      for(int i=0; i<vHS.shape()[0]; i++)
        Vsum += vHS[i][0];
      app_log()<<" Vsum: " <<setprecision(12) <<Vsum <<std::endl;
    }

    TG.Global().barrier();
    // remove temporary file
    if(TG.Node().root())
      remove("dummy.h5");
  }

}

TEST_CASE("wfn_fac_sdet_distributed", "[wavefunction_factory]")
{

  if(not file_exists("./afqmc.h5") ||
     not file_exists("./wfn.dat") ) {
    app_log()<<" Skipping wfn_fac_sdet_distributed text. afqmc.h5 and ./wfn.dat files not found. \n";
  } else {

    // mpi3
    communicator& world = OHMMS::Controller->comm;

    // Global Task Group
    GlobalTaskGroup gTG(world);

    auto file_data = read_test_results_from_hdf<ValueType>("./afqmc.h5");
    int NMO=file_data.NMO;
    int NAEA=file_data.NAEA;
    int NAEB=file_data.NAEB;
    WALKER_TYPES type = afqmc::getWalkerType("wfn.dat");

    std::map<std::string,AFQMCInfo> InfoMap;
    InfoMap.insert ( std::pair<std::string,AFQMCInfo>("info0",AFQMCInfo{"info0",NMO,NAEA,NAEB}) );
    HamiltonianFactory HamFac(InfoMap);
    const char *ham_xml_block =
"<Hamiltonian name=\"ham0\" type=\"SparseGeneral\" info=\"info0\"> \
    <parameter name=\"filetype\">hdf5</parameter> \
    <parameter name=\"filename\">./afqmc.h5</parameter> \
    <parameter name=\"cutoff_decomposition\">1e-5</parameter> \
  </Hamiltonian> \
";
    Libxml2Document doc;
    bool okay = doc.parseFromString(ham_xml_block);
    REQUIRE(okay);
    std::string ham_name("ham0");
    HamFac.push(ham_name,doc.getRoot());
    Hamiltonian& ham = HamFac.getHamiltonian(gTG,ham_name);


    auto TG = TaskGroup_(gTG,std::string("WfnTG"),1,gTG.getTotalCores());
    auto TGwfn = TaskGroup_(gTG,std::string("WfnTG"),gTG.getTotalNodes(),gTG.getTotalCores());
    int nwalk = 11; // choose prime number to force non-trivial splits in shared routines
    RandomGenerator_t rng;

const char *wlk_xml_block_closed =
"<WalkerSet name=\"wset0\">  \
  <parameter name=\"walker_type\">closed</parameter>  \
</WalkerSet> \
";
const char *wlk_xml_block_coll =
"<WalkerSet name=\"wset0\">  \
  <parameter name=\"walker_type\">collinear</parameter>  \
</WalkerSet> \
";
const char *wlk_xml_block_noncol =
"<WalkerSet name=\"wset0\">  \
  <parameter name=\"walker_type\">noncollinear</parameter>  \
</WalkerSet> \
";

    const char *wlk_xml_block = ( (type==CLOSED)?(wlk_xml_block_closed):
                                  (type==COLLINEAR?wlk_xml_block_coll:wlk_xml_block_noncol) );

    Libxml2Document doc3;
    okay = doc3.parseFromString(wlk_xml_block);
    REQUIRE(okay);

    const char *wfn_xml_block =
"<Wavefunction name=\"wfn0\" info=\"info0\"> \
      <parameter name=\"filetype\">ascii</parameter> \
      <parameter name=\"filename\">./wfn.dat</parameter> \
      <parameter name=\"cutoff\">1e-6</parameter> \
      <parameter name=\"restart_file\">dummy.h5</parameter> \
  </Wavefunction> \
";
    Libxml2Document doc2;
    okay = doc2.parseFromString(wfn_xml_block);
    REQUIRE(okay);
    std::string wfn_name("wfn0");
    WavefunctionFactory WfnFac(InfoMap); 
    WfnFac.push(wfn_name,doc2.getRoot());
    Wavefunction& wfn = WfnFac.getWavefunction(TGwfn,TGwfn,wfn_name,type,&ham,1e-6,nwalk); 

    WalkerSet wset(TG,doc3.getRoot(),InfoMap["info0"],&rng);
    auto initial_guess = WfnFac.getInitialGuess(wfn_name);
    REQUIRE(initial_guess.shape()[0]==2);
    REQUIRE(initial_guess.shape()[1]==NMO);
    REQUIRE(initial_guess.shape()[2]==NAEA);

    if(type == COLLINEAR)
        wset.resize(nwalk,initial_guess[0],
                         initial_guess[1][indices[range_t()][range_t(0,NAEB)]]);
    else
        wset.resize(nwalk,initial_guess[0],
                         initial_guess[0]);

    wfn.Overlap(wset);
    for(auto it = wset.begin(); it!=wset.end(); ++it) {
      REQUIRE(real(it->overlap()) == Approx(1.0));
      REQUIRE(imag(it->overlap()) == Approx(0.0));
    }

    using shm_Alloc = boost::mpi3::intranode::allocator<ComplexType>;
    using SHM_Buffer = mpi3_SHMBuffer<ComplexType>;
 
    wfn.Energy(wset);
    if(std::abs(file_data.E0+file_data.E1+file_data.E2)>1e-8) {
      for(auto it = wset.begin(); it!=wset.end(); ++it) {
        REQUIRE( real(it->E1()) == Approx(real(file_data.E0+file_data.E1)));
        REQUIRE( real(it->EXX()+it->EJ()) == Approx(real(file_data.E2)));
        REQUIRE( imag(it->energy()) == Approx(imag(file_data.E0+file_data.E1+file_data.E2)));
      }
    } else {
      app_log()<<" E0+E1: " <<setprecision(12) <<wset[0].E1() <<std::endl;
      app_log()<<" EJ: " <<setprecision(12) <<wset[0].EJ() <<std::endl;
      app_log()<<" EXX: " <<setprecision(12) <<wset[0].EXX() <<std::endl;
    }

    auto size_of_G = wfn.size_of_G_for_vbias();
    SHM_Buffer Gbuff(TG.TG_local(),nwalk*size_of_G);
    int Gdim1 = (wfn.transposed_G_for_vbias()?nwalk:size_of_G);
    int Gdim2 = (wfn.transposed_G_for_vbias()?size_of_G:nwalk);
    boost::multi_array_ref<ComplexType,2> G(Gbuff.data(),extents[Gdim1][Gdim2]);
    wfn.MixedDensityMatrix_for_vbias(wset,G);

    double sqrtdt = std::sqrt(0.01);
    auto nCV = wfn.local_number_of_cholesky_vectors();
    SHM_Buffer Xbuff(TG.TG_local(),nCV*nwalk);
    boost::multi_array_ref<ComplexType,2> X(Xbuff.data(),extents[nCV][nwalk]);
    wfn.vbias(G,X,sqrtdt);

    ComplexType Xsum=0;
    if(std::abs(file_data.Xsum)>1e-8) {
      for(int n=0; n<nwalk; n++) {
        Xsum=0;
        if(TGwfn.TG_local().root()) 
          for(int i=0; i<X.shape()[0]; i++) 
            Xsum += X[i][n];
        Xsum = TGwfn.TG().all_reduce_value(Xsum,std::plus<>());
        REQUIRE( real(Xsum) == Approx(real(file_data.Xsum)) );
        REQUIRE( imag(Xsum) == Approx(imag(file_data.Xsum)) );
      }
    } else {
      Xsum=0;
      if(TGwfn.TG_local().root()) 
        for(int i=0; i<X.shape()[0]; i++)
          Xsum += X[i][0];
      Xsum = TGwfn.TG().all_reduce_value(Xsum,std::plus<>());
      app_log()<<" Xsum: " <<setprecision(12) <<Xsum <<std::endl;
    }

    // vbias must be reduced if false
    if(not wfn.distribution_over_cholesky_vectors()) {
      boost::multi_array<ComplexType,2> T(extents[nCV][nwalk]);
      if(TGwfn.TG_local().root())
        std::copy_n(X.origin(),X.num_elements(),T.origin());
      else
        std::fill_n(T.origin(),T.num_elements(),ComplexType(0.0,0.0));
      TGwfn.TG().all_reduce_in_place_n(T.origin(),T.num_elements(),std::plus<>());
      if(TGwfn.TG_local().root())
        std::copy_n(T.origin(),T.num_elements(),X.origin());
      TGwfn.TG_local().barrier();
    }

    SHM_Buffer vHSbuff(TG.TG_local(),NMO*NMO*nwalk);
    int vdim1 = (wfn.transposed_vHS()?nwalk:NMO*NMO);
    int vdim2 = (wfn.transposed_vHS()?NMO*NMO:nwalk);
    boost::multi_array_ref<ComplexType,2> vHS(vHSbuff.data(),extents[vdim1][vdim2]);
    wfn.vHS(X,vHS,sqrtdt);
    TG.local_barrier();
    ComplexType Vsum=0;
    if(std::abs(file_data.Vsum)>1e-8) {
      for(int n=0; n<nwalk; n++) {
        Vsum=0;
        if(TGwfn.TG_local().root()) 
          for(int i=0; i<vHS.shape()[0]; i++)
            Vsum += vHS[i][n];
        Vsum = TGwfn.TG().all_reduce_value(Vsum,std::plus<>());
        REQUIRE( real(Vsum) == Approx(real(file_data.Vsum)) );
        REQUIRE( imag(Vsum) == Approx(imag(file_data.Vsum)) );
      }
    } else {
      Vsum=0;
      if(TGwfn.TG_local().root()) 
        for(int i=0; i<vHS.shape()[0]; i++)
          Vsum += vHS[i][0];
      Vsum = TGwfn.TG().all_reduce_value(Vsum,std::plus<>());
      app_log()<<" Vsum: " <<setprecision(12) <<Vsum <<std::endl;
    }

  // Restarting Wavefunction from file
    const char *wfn_xml_block_restart =
"<Wavefunction name=\"wfn1\" info=\"info0\"> \
      <parameter name=\"filetype\">hdf5</parameter> \
      <parameter name=\"filename\">./dummy.h5</parameter> \
      <parameter name=\"cutoff\">1e-6</parameter> \
  </Wavefunction> \
";
    Libxml2Document doc4;
    okay = doc4.parseFromString(wfn_xml_block_restart);
    REQUIRE(okay);
    wfn_name = "wfn1";
    WfnFac.push(wfn_name,doc4.getRoot());
    Wavefunction& wfn2 = WfnFac.getWavefunction(TGwfn,TGwfn,wfn_name,type,nullptr,1e-8,nwalk);

    WalkerSet wset2(TG,doc3.getRoot(),InfoMap["info0"],&rng);
    //auto initial_guess = WfnFac.getInitialGuess(wfn_name);
    REQUIRE(initial_guess.shape()[0]==2);
    REQUIRE(initial_guess.shape()[1]==NMO);
    REQUIRE(initial_guess.shape()[2]==NAEA);

    if(type == COLLINEAR)
        wset2.resize(nwalk,initial_guess[0],
                         initial_guess[1][indices[range_t()][range_t(0,NAEB)]]);
    else
        wset2.resize(nwalk,initial_guess[0],
                         initial_guess[0]);

    wfn2.Overlap(wset2);
    for(auto it = wset2.begin(); it!=wset2.end(); ++it) {
      REQUIRE(real(it->overlap()) == Approx(1.0));
      REQUIRE(imag(it->overlap()) == Approx(0.0));
    }

    wfn2.Energy(wset2);
    if(std::abs(file_data.E0+file_data.E1+file_data.E2)>1e-8) {
      for(auto it = wset2.begin(); it!=wset2.end(); ++it) {
        REQUIRE( real(it->E1()) == Approx(real(file_data.E0+file_data.E1)));
        REQUIRE( real(it->EXX()+it->EJ()) == Approx(real(file_data.E2)));
        REQUIRE( imag(it->energy()) == Approx(imag(file_data.E0+file_data.E1+file_data.E2)));
      }
    } else {
      app_log()<<" E: " <<wset[0].energy() <<std::endl;
      app_log()<<" E0+E1: " <<wset[0].E1() <<std::endl;
      app_log()<<" EJ: " <<wset[0].EJ() <<std::endl;
      app_log()<<" EXX: " <<wset[0].EXX() <<std::endl;
    }

    REQUIRE(size_of_G == wfn2.size_of_G_for_vbias());
    wfn2.MixedDensityMatrix_for_vbias(wset2,G);

    nCV = wfn2.local_number_of_cholesky_vectors();
    Xbuff.resize(nCV*nwalk); 
    boost::multi_array_ref<ComplexType,2> X2(Xbuff.data(),extents[nCV][nwalk]);
    wfn2.vbias(G,X2,sqrtdt);
    Xsum=0;
    if(std::abs(file_data.Xsum)>1e-8) {
      for(int n=0; n<nwalk; n++) {
        Xsum=0;
        if(TGwfn.TG_local().root())
          for(int i=0; i<X2.shape()[0]; i++)
            Xsum += X2[i][n];
        Xsum = TGwfn.TG().all_reduce_value(Xsum,std::plus<>());
        REQUIRE( real(Xsum) == Approx(real(file_data.Xsum)) );
        REQUIRE( imag(Xsum) == Approx(imag(file_data.Xsum)) );
      }
    } else {
      Xsum=0;
      if(TGwfn.TG_local().root())
        for(int i=0; i<X2.shape()[0]; i++)
          Xsum += X2[i][0];
      Xsum = TGwfn.TG().all_reduce_value(Xsum,std::plus<>());
      app_log()<<" Xsum: " <<setprecision(12) <<Xsum <<std::endl;
    }

    // vbias must be reduced if false
    if(not wfn.distribution_over_cholesky_vectors()) {
      boost::multi_array<ComplexType,2> T(extents[nCV][nwalk]);
      if(TGwfn.TG_local().root())
        std::copy_n(X2.origin(),X2.num_elements(),T.origin());
      else
        std::fill_n(T.origin(),T.num_elements(),ComplexType(0.0,0.0));
      TGwfn.TG().all_reduce_in_place_n(T.origin(),T.num_elements(),std::plus<>());
      if(TGwfn.TG_local().root())
        std::copy_n(T.origin(),T.num_elements(),X.origin());
      TGwfn.TG_local().barrier();
    }

    wfn2.vHS(X2,vHS,sqrtdt);
    TG.local_barrier();
    Vsum=0;
    if(std::abs(file_data.Vsum)>1e-8) {
      for(int n=0; n<nwalk; n++) {
        Vsum=0;
        if(TGwfn.TG_local().root())
          for(int i=0; i<vHS.shape()[0]; i++)
            Vsum += vHS[i][n];
        Vsum = TGwfn.TG().all_reduce_value(Vsum,std::plus<>());
        REQUIRE( real(Vsum) == Approx(real(file_data.Vsum)) );
        REQUIRE( imag(Vsum) == Approx(imag(file_data.Vsum)) );
      }
    } else {
      Vsum=0;
      if(TGwfn.TG_local().root())
        for(int i=0; i<vHS.shape()[0]; i++)
          Vsum += vHS[i][0];
      Vsum = TGwfn.TG().all_reduce_value(Vsum,std::plus<>());
      app_log()<<" Vsum: " <<setprecision(12) <<Vsum <<std::endl;
    }

    TG.Global().barrier();
    // remove temporary file
    if(TG.Node().root())
      remove("dummy.h5");


  }
}

TEST_CASE("wfn_fac_collinear_multidet", "[wavefunction_factory]")
{

  if(not file_exists("./afqmc_msd.h5") ||
     not file_exists("./wfn_msd.dat") ) {
    app_log()<<" Skipping wfn_fac_collinear_multidet text. afqmc_msd.h5 and ./wfn_msd.dat files not found. \n";
  } else {

    // mpi3
    communicator& world = OHMMS::Controller->comm;

    // Global Task Group
    GlobalTaskGroup gTG(world);

    auto file_data = read_test_results_from_hdf<ValueType>("./afqmc_msd.h5");
    int NMO=file_data.NMO;
    int NAEA=file_data.NAEA;
    int NAEB=file_data.NAEB;

    std::map<std::string,AFQMCInfo> InfoMap;
    InfoMap.insert ( std::pair<std::string,AFQMCInfo>("info0",AFQMCInfo{"info0",NMO,NAEA,NAEB}) );
    HamiltonianFactory HamFac(InfoMap);
    const char *ham_xml_block =
"<Hamiltonian name=\"ham0\" type=\"SparseGeneral\" info=\"info0\"> \
    <parameter name=\"filetype\">hdf5</parameter> \
    <parameter name=\"filename\">./afqmc_msd.h5</parameter> \
    <parameter name=\"cutoff_decomposition\">1e-5</parameter> \
  </Hamiltonian> \
";
    Libxml2Document doc;
    bool okay = doc.parseFromString(ham_xml_block);
    REQUIRE(okay);
    std::string ham_name("ham0");
    HamFac.push(ham_name,doc.getRoot());
    Hamiltonian& ham = HamFac.getHamiltonian(gTG,ham_name);


    auto TG = TaskGroup_(gTG,std::string("WfnTG"),1,gTG.getTotalCores());
    int nwalk = 11; // choose prime number to force non-trivial splits in shared routines
    RandomGenerator_t rng;

const char *wlk_xml_block =
"<WalkerSet name=\"wset0\">  \
  <parameter name=\"walker_type\">collinear</parameter>  \
</WalkerSet> \
";
    Libxml2Document doc3;
    okay = doc3.parseFromString(wlk_xml_block);
    REQUIRE(okay);

    const char *wfn_xml_block =
"<Wavefunction name=\"wfn0\" info=\"info0\"> \
      <parameter name=\"filetype\">ascii</parameter> \
      <parameter name=\"filename\">./wfn_msd.dat</parameter> \
      <parameter name=\"cutoff\">1e-6</parameter> \
  </Wavefunction> \
";
    Libxml2Document doc2;
    okay = doc2.parseFromString(wfn_xml_block);
    REQUIRE(okay);
    std::string wfn_name("wfn0");
    WavefunctionFactory WfnFac(InfoMap); 
    WfnFac.push(wfn_name,doc2.getRoot());
    Wavefunction& wfn = WfnFac.getWavefunction(TG,TG,wfn_name,COLLINEAR,&ham,1e-6,nwalk); 

    WalkerSet wset(TG,doc3.getRoot(),InfoMap["info0"],&rng);
    auto initial_guess = WfnFac.getInitialGuess(wfn_name);
    REQUIRE(initial_guess.shape()[0]==2);
    REQUIRE(initial_guess.shape()[1]==NMO);
    REQUIRE(initial_guess.shape()[2]==NAEA);
    wset.resize(nwalk,initial_guess[0],
                         initial_guess[1][indices[range_t()][range_t(0,NAEB)]]);

    // no guarantee that overlap is 1.0
    wfn.Overlap(wset);

    using shm_Alloc = boost::mpi3::intranode::allocator<ComplexType>;
    using SHM_Buffer = mpi3_SHMBuffer<ComplexType>;

    wfn.Energy(wset);
    if(std::abs(file_data.E0+file_data.E1+file_data.E2)>1e-8) {
      for(auto it = wset.begin(); it!=wset.end(); ++it) {
        REQUIRE( real(it->E1()) == Approx(real(file_data.E0+file_data.E1)));
        REQUIRE( real(it->EXX()+it->EJ()) == Approx(real(file_data.E2)));
        REQUIRE( imag(it->energy()) == Approx(imag(file_data.E0+file_data.E1+file_data.E2)));
      }
    } else {
      app_log()<<" E: " <<wset[0].E1() <<" " <<wset[0].EXX()+wset[0].EJ() <<std::endl;
    }

    auto size_of_G = wfn.size_of_G_for_vbias();
    SHM_Buffer Gbuff(TG.TG_local(),nwalk*size_of_G);
    int Gdim1 = (wfn.transposed_G_for_vbias()?nwalk:size_of_G);
    int Gdim2 = (wfn.transposed_G_for_vbias()?size_of_G:nwalk);
    boost::multi_array_ref<ComplexType,2> G(Gbuff.data(),extents[Gdim1][Gdim2]);
    wfn.MixedDensityMatrix_for_vbias(wset,G);

    double sqrtdt = std::sqrt(0.01);
    auto nCV = wfn.local_number_of_cholesky_vectors();
    SHM_Buffer Xbuff(TG.TG_local(),nCV*nwalk);
    boost::multi_array_ref<ComplexType,2> X(Xbuff.data(),extents[nCV][nwalk]);
    wfn.vbias(G,X,sqrtdt);
    ComplexType Xsum=0;
    if(std::abs(file_data.Xsum)>1e-8) {
      for(int n=0; n<nwalk; n++) {
        Xsum=0;
        for(int i=0; i<X.shape()[0]; i++)
          Xsum += X[i][n];
        REQUIRE( real(Xsum) == Approx(real(file_data.Xsum)) );
        REQUIRE( imag(Xsum) == Approx(imag(file_data.Xsum)) );
      }
    } else {
      Xsum=0;
      for(int i=0; i<X.shape()[0]; i++)
        Xsum += X[i][0];
      app_log()<<" Xsum: " <<setprecision(12) <<Xsum <<std::endl;
    }


    SHM_Buffer vHSbuff(TG.TG_local(),NMO*NMO*nwalk);
    int vdim1 = (wfn.transposed_vHS()?nwalk:NMO*NMO);
    int vdim2 = (wfn.transposed_vHS()?NMO*NMO:nwalk);
    boost::multi_array_ref<ComplexType,2> vHS(vHSbuff.data(),extents[vdim1][vdim2]);
    wfn.vHS(X,vHS,sqrtdt);
    TG.local_barrier();
    ComplexType Vsum=0;
    if(std::abs(file_data.Vsum)>1e-8) {
      for(int n=0; n<nwalk; n++) {
        Vsum=0;
        for(int i=0; i<vHS.shape()[0]; i++)
          Vsum += vHS[i][n];
        REQUIRE( real(Vsum) == Approx(real(file_data.Vsum)) );
        REQUIRE( imag(Vsum) == Approx(imag(file_data.Vsum)) );
      }
    } else {
      Vsum=0;
      for(int i=0; i<vHS.shape()[0]; i++)
        Vsum += vHS[i][0];
      app_log()<<" Vsum: " <<setprecision(12) <<Vsum <<std::endl;
    }
  }
}

TEST_CASE("wfn_fac_collinear_multidet_distributed", "[wavefunction_factory]")
{

  if(not file_exists("./afqmc_msd.h5") ||
     not file_exists("./wfn_msd.dat") ) {
    app_log()<<" Skipping wfn_fac_collinear_multidet text. afqmc_msd.h5 and ./wfn_msd.dat files not found. \n";
  } else {

    // mpi3
    communicator& world = OHMMS::Controller->comm;

    // Global Task Group
    GlobalTaskGroup gTG(world);

    auto file_data = read_test_results_from_hdf<ValueType>("./afqmc_msd.h5");
    int NMO=file_data.NMO;
    int NAEA=file_data.NAEA;
    int NAEB=file_data.NAEB;

    std::map<std::string,AFQMCInfo> InfoMap;
    InfoMap.insert ( std::pair<std::string,AFQMCInfo>("info0",AFQMCInfo{"info0",NMO,NAEA,NAEB}) );
    HamiltonianFactory HamFac(InfoMap);
    const char *ham_xml_block =
"<Hamiltonian name=\"ham0\" type=\"SparseGeneral\" info=\"info0\"> \
    <parameter name=\"filetype\">hdf5</parameter> \
    <parameter name=\"filename\">./afqmc_msd.h5</parameter> \
    <parameter name=\"cutoff_decomposition\">1e-5</parameter> \
  </Hamiltonian> \
";
    Libxml2Document doc;
    bool okay = doc.parseFromString(ham_xml_block);
    REQUIRE(okay);
    std::string ham_name("ham0");
    HamFac.push(ham_name,doc.getRoot());
    Hamiltonian& ham = HamFac.getHamiltonian(gTG,ham_name);


    auto TG = TaskGroup_(gTG,std::string("TG"),1,gTG.getTotalCores());
    auto TGwfn = TaskGroup_(gTG,std::string("WfnTG"),gTG.getTotalNodes(),gTG.getTotalCores());
    int nwalk = 11; // choose prime number to force non-trivial splits in shared routines
    RandomGenerator_t rng;

const char *wlk_xml_block =
"<WalkerSet name=\"wset0\">  \
  <parameter name=\"walker_type\">collinear</parameter>  \
</WalkerSet> \
";
    Libxml2Document doc3;
    okay = doc3.parseFromString(wlk_xml_block);
    REQUIRE(okay);

    const char *wfn_xml_block =
"<Wavefunction name=\"wfn0\" info=\"info0\"> \
      <parameter name=\"filetype\">ascii</parameter> \
      <parameter name=\"filename\">./wfn_msd.dat</parameter> \
      <parameter name=\"cutoff\">1e-6</parameter> \
  </Wavefunction> \
";
    Libxml2Document doc2;
    okay = doc2.parseFromString(wfn_xml_block);
    REQUIRE(okay);
    std::string wfn_name("wfn0");
    WavefunctionFactory WfnFac(InfoMap); 
    WfnFac.push(wfn_name,doc2.getRoot());
    Wavefunction& wfn = WfnFac.getWavefunction(TG,TGwfn,wfn_name,COLLINEAR,&ham,1e-6,nwalk); 

    WalkerSet wset(TG,doc3.getRoot(),InfoMap["info0"],&rng);
    auto initial_guess = WfnFac.getInitialGuess(wfn_name);
    REQUIRE(initial_guess.shape()[0]==2);
    REQUIRE(initial_guess.shape()[1]==NMO);
    REQUIRE(initial_guess.shape()[2]==NAEA);
    wset.resize(nwalk,initial_guess[0],
                         initial_guess[1][indices[range_t()][range_t(0,NAEB)]]);

    // no guarantee that overlap is 1.0
    wfn.Overlap(wset);

    using shm_Alloc = boost::mpi3::intranode::allocator<ComplexType>;
    using SHM_Buffer = mpi3_SHMBuffer<ComplexType>;

    wfn.Energy(wset);
    if(std::abs(file_data.E0+file_data.E1+file_data.E2)>1e-8) {
      for(auto it = wset.begin(); it!=wset.end(); ++it) {
        REQUIRE( real(it->E1()) == Approx(real(file_data.E0+file_data.E1)));
        REQUIRE( real(it->EXX()+it->EJ()) == Approx(real(file_data.E2)));
        REQUIRE( imag(it->energy()) == Approx(imag(file_data.E0+file_data.E1+file_data.E2)));
      }
    } else {
      app_log()<<" E: " <<wset[0].energy() <<std::endl;
    }

    auto size_of_G = wfn.size_of_G_for_vbias();
    SHM_Buffer Gbuff(TG.TG_local(),nwalk*size_of_G);
    int Gdim1 = (wfn.transposed_G_for_vbias()?nwalk:size_of_G);
    int Gdim2 = (wfn.transposed_G_for_vbias()?size_of_G:nwalk);
    boost::multi_array_ref<ComplexType,2> G(Gbuff.data(),extents[Gdim1][Gdim2]);
    wfn.MixedDensityMatrix_for_vbias(wset,G);

    double sqrtdt = std::sqrt(0.01);
    auto nCV = wfn.local_number_of_cholesky_vectors();
    SHM_Buffer Xbuff(TG.TG_local(),nCV*nwalk);
    boost::multi_array_ref<ComplexType,2> X(Xbuff.data(),extents[nCV][nwalk]);
    wfn.vbias(G,X,sqrtdt);
    ComplexType Xsum=0;
    if(std::abs(file_data.Xsum)>1e-8) {
      for(int n=0; n<nwalk; n++) {
        Xsum=0;
        for(int i=0; i<X.shape()[0]; i++)
          Xsum += X[i][n];
        REQUIRE( real(Xsum) == Approx(real(file_data.Xsum)) );
        REQUIRE( imag(Xsum) == Approx(imag(file_data.Xsum)) );
      }
    } else {
      Xsum=0;
      for(int i=0; i<X.shape()[0]; i++)
        Xsum += X[i][0];
      app_log()<<" Xsum: " <<setprecision(12) <<Xsum <<std::endl;
    }


    SHM_Buffer vHSbuff(TG.TG_local(),NMO*NMO*nwalk);
    int vdim1 = (wfn.transposed_vHS()?nwalk:NMO*NMO);
    int vdim2 = (wfn.transposed_vHS()?NMO*NMO:nwalk);
    boost::multi_array_ref<ComplexType,2> vHS(vHSbuff.data(),extents[vdim1][vdim2]);
    wfn.vHS(X,vHS,sqrtdt);
    TG.local_barrier();
    ComplexType Vsum=0;
    if(std::abs(file_data.Vsum)>1e-8) {
      for(int n=0; n<nwalk; n++) {
        Vsum=0;
        for(int i=0; i<vHS.shape()[0]; i++)
          Vsum += vHS[i][n];
        REQUIRE( real(Vsum) == Approx(real(file_data.Vsum)) );
        REQUIRE( imag(Vsum) == Approx(imag(file_data.Vsum)) );
      }
    } else {
      Vsum=0;
      for(int i=0; i<vHS.shape()[0]; i++)
        Vsum += vHS[i][0];
      app_log()<<" Vsum: " <<setprecision(12) <<Vsum <<std::endl;
    }
  }
}

}
