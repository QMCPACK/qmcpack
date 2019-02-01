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
#include "Utilities/Timer.h"

#undef APP_ABORT
#define APP_ABORT(x) {std::cout << x <<std::endl; exit(0);}

#include <sys/stat.h>
#include <unistd.h>
#include <stdio.h>
#include <string>
#include <vector>
#include <complex>
#include <iomanip>
#include <random>

#include "AFQMC/Utilities/test_utils.hpp"

#include "AFQMC/Hamiltonians/HamiltonianFactory.h"
#include "AFQMC/Hamiltonians/Hamiltonian.hpp"
#include "AFQMC/Utilities/myTimer.h"
#include "AFQMC/SlaterDeterminantOperations/SlaterDetOperations.hpp"
#include "AFQMC/Wavefunctions/WavefunctionFactory.h"
#include "AFQMC/Walkers/WalkerSet.hpp"

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

TEST_CASE("wfn_fac_sdet", "[wavefunction_factory]")
{
  OHMMS::Controller->initialize(0, NULL);
  auto world = boost::mpi3::environment::get_world_instance();

  if(not file_exists("./afqmc.h5") ||
     not file_exists("./wfn.dat") ) {
    app_log()<<" Skipping wfn_fac_collinear_sdet text. afqmc.h5 and ./wfn.dat files not found. \n";
  } else {

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
"<Hamiltonian name=\"ham0\" info=\"info0\"> \
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


    //auto TG = TaskGroup_(gTG,std::string("WfnTG"),1,1);
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

    //for(int nw=1; nw<2; nw*=2)
    {
      //nwalk=nw;
      WalkerSet wset(TG,doc3.getRoot(),InfoMap["info0"],&rng);
      auto initial_guess = WfnFac.getInitialGuess(wfn_name);
      REQUIRE(initial_guess.shape()[0]==2);
      REQUIRE(initial_guess.shape()[1]==NMO);
      REQUIRE(initial_guess.shape()[2]==NAEA);

      if(type == COLLINEAR)
        wset.resize(nwalk,initial_guess[0],
                         initial_guess[1](initial_guess.extension(1),{0,NAEB}));
      else
        wset.resize(nwalk,initial_guess[0],
                         initial_guess[0]);

      wfn.Overlap(wset);
      for(auto it = wset.begin(); it!=wset.end(); ++it) {
        REQUIRE(real(it->overlap()) == Approx(1.0));
        REQUIRE(imag(it->overlap()) == Approx(0.0));
      }

      using shmCMatrix = boost::multi::array<ComplexType,2,shared_allocator<ComplexType>>;

      qmcplusplus::Timer Time;
      double t1;
      Time.restart();
      wfn.Energy(wset);
      TG.local_barrier();
      t1=Time.elapsed();
      if(std::abs(file_data.E0+file_data.E1+file_data.E2)>1e-8) {
        for(auto it = wset.begin(); it!=wset.end(); ++it) {
          REQUIRE( real(it->E1()) == Approx(real(file_data.E0+file_data.E1)));
          REQUIRE( real(it->EXX()+it->EJ()) == Approx(real(file_data.E2)));
          REQUIRE( imag(it->energy()) == Approx(imag(file_data.E0+file_data.E1+file_data.E2)));
        }
      } else {
        app_log()<<" E: " <<setprecision(12) <<wset[0].energy() <<" Time: " <<t1 <<std::endl;
        app_log()<<" E0+E1: " <<setprecision(12) <<wset[0].E1() <<std::endl;
        app_log()<<" EJ: " <<setprecision(12) <<wset[0].EJ() <<std::endl;
        app_log()<<" EXX: " <<setprecision(12) <<wset[0].EXX() <<std::endl;
      }

      auto size_of_G = wfn.size_of_G_for_vbias();
      int Gdim1 = (wfn.transposed_G_for_vbias()?nwalk:size_of_G);
      int Gdim2 = (wfn.transposed_G_for_vbias()?size_of_G:nwalk);
      shmCMatrix G({Gdim1,Gdim2},shared_allocator<ComplexType>{TG.TG_local()});
      wfn.MixedDensityMatrix_for_vbias(wset,G);

      double sqrtdt = std::sqrt(0.01);
      auto nCV = wfn.local_number_of_cholesky_vectors();
      shmCMatrix X({nCV,nwalk},shared_allocator<ComplexType>{TG.TG_local()});
      Time.restart();
      wfn.vbias(G,X,sqrtdt);
      TG.local_barrier();
      t1=Time.elapsed();
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
        ComplexType Xsum2=0;
        for(int i=0; i<X.shape()[0]; i++) {
          Xsum += X[i][0];
          Xsum2 += 0.5*X[i][0]*X[i][0];
        }
        app_log()<<" Xsum: " <<setprecision(12) <<Xsum <<" Time: " <<t1 <<std::endl;
        app_log()<<" Xsum2 (EJ): " <<setprecision(12) <<Xsum2/sqrtdt/sqrtdt <<std::endl;
      }

      int vdim1 = (wfn.transposed_vHS()?nwalk:NMO*NMO);
      int vdim2 = (wfn.transposed_vHS()?NMO*NMO:nwalk);
      shmCMatrix vHS({vdim1,vdim2},shared_allocator<ComplexType>{TG.TG_local()});
      Time.restart();
      wfn.vHS(X,vHS,sqrtdt);
      TG.local_barrier();
      t1=Time.elapsed();
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
        if(wfn.transposed_vHS()) {
          for(int i=0; i<vHS.shape()[1]; i++)
            Vsum += vHS[0][i];
        } else {
          for(int i=0; i<vHS.shape()[0]; i++)
            Vsum += vHS[i][0];
        }
        app_log()<<" Vsum: " <<setprecision(12) <<Vsum <<" Time: " <<t1 <<std::endl;
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
                         initial_guess[1](initial_guess.extension(1),{0,NAEB}));
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
        ComplexType Xsum2(0.0);
        for(int i=0; i<X.shape()[0]; i++) {
          Xsum += X[i][0];
          Xsum2 += 0.5*X[i][0]*X[i][0];
        }
        app_log()<<" Xsum: " <<setprecision(12) <<Xsum <<std::endl;
        app_log()<<" Xsum2 (EJ): " <<setprecision(12) <<Xsum2/sqrtdt/sqrtdt <<std::endl;
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
}

TEST_CASE("wfn_fac_sdet_distributed", "[wavefunction_factory]")
{
  OHMMS::Controller->initialize(0, NULL);
  auto world = boost::mpi3::environment::get_world_instance();

  if(not file_exists("./afqmc.h5") ||
     not file_exists("./wfn.dat") ) {
    app_log()<<" Skipping wfn_fac_sdet_distributed text. afqmc.h5 and ./wfn.dat files not found. \n";
  } else {

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
                         initial_guess[1](initial_guess.extension(1),{0,NAEB}));
    else
        wset.resize(nwalk,initial_guess[0],
                         initial_guess[0]);

    wfn.Overlap(wset);
    for(auto it = wset.begin(); it!=wset.end(); ++it) {
      REQUIRE(real(it->overlap()) == Approx(1.0));
      REQUIRE(imag(it->overlap()) == Approx(0.0));
    }

    using shmCMatrix = boost::multi::array<ComplexType,2,shared_allocator<ComplexType>>;

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
    int Gdim1 = (wfn.transposed_G_for_vbias()?nwalk:size_of_G);
    int Gdim2 = (wfn.transposed_G_for_vbias()?size_of_G:nwalk);
    shmCMatrix G({Gdim1,Gdim2},shared_allocator<ComplexType>{TG.TG_local()});
    wfn.MixedDensityMatrix_for_vbias(wset,G);

    double sqrtdt = std::sqrt(0.01);
    auto nCV = wfn.local_number_of_cholesky_vectors();
    shmCMatrix X({nCV,nwalk},shared_allocator<ComplexType>{TG.TG_local()});
    wfn.vbias(G,X,sqrtdt);

    ComplexType Xsum=0;
    if(std::abs(file_data.Xsum)>1e-8) {
      for(int n=0; n<nwalk; n++) {
        Xsum=0;
        if(TGwfn.TG_local().root())
          for(int i=0; i<X.shape()[0]; i++)
            Xsum += X[i][n];
        Xsum = ( TGwfn.TG() += Xsum );
        REQUIRE( real(Xsum) == Approx(real(file_data.Xsum)) );
        REQUIRE( imag(Xsum) == Approx(imag(file_data.Xsum)) );
      }
    } else {
      Xsum=0;
      if(TGwfn.TG_local().root())
        for(int i=0; i<X.shape()[0]; i++)
          Xsum += X[i][0];
      Xsum = ( TGwfn.TG() += Xsum );
      app_log()<<" Xsum: " <<setprecision(12) <<Xsum <<std::endl;
    }

    // vbias must be reduced if false
    if(not wfn.distribution_over_cholesky_vectors()) {
      boost::multi::array<ComplexType,2> T({nCV,nwalk});
      if(TGwfn.TG_local().root())
        std::copy_n(X.origin(),X.num_elements(),T.origin());
      else
        std::fill_n(T.origin(),T.num_elements(),ComplexType(0.0,0.0));
      TGwfn.TG().all_reduce_in_place_n(std::addressof(*T.origin()),T.num_elements(),std::plus<>());
      if(TGwfn.TG_local().root())
        std::copy_n(T.origin(),T.num_elements(),X.origin());
      TGwfn.TG_local().barrier();
    }

    int vdim1 = (wfn.transposed_vHS()?nwalk:NMO*NMO);
    int vdim2 = (wfn.transposed_vHS()?NMO*NMO:nwalk);
    shmCMatrix vHS({vdim1,vdim2},shared_allocator<ComplexType>{TG.TG_local()});
    wfn.vHS(X,vHS,sqrtdt);
    TG.local_barrier();
    ComplexType Vsum=0;
    if(std::abs(file_data.Vsum)>1e-8) {
      for(int n=0; n<nwalk; n++) {
        Vsum=0;
        if(TGwfn.TG_local().root())
          for(int i=0; i<vHS.shape()[0]; i++)
            Vsum += vHS[i][n];
        Vsum = ( TGwfn.TG() += Vsum );
        REQUIRE( real(Vsum) == Approx(real(file_data.Vsum)) );
        REQUIRE( imag(Vsum) == Approx(imag(file_data.Vsum)) );
      }
    } else {
      Vsum=0;
      if(TGwfn.TG_local().root())
        for(int i=0; i<vHS.shape()[0]; i++)
          Vsum += vHS[i][0];
      Vsum = ( TGwfn.TG() += Vsum );
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
                         initial_guess[1](initial_guess.extension(1),{0,NAEB}));
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
    boost::multi::array_ref<ComplexType,2> X2(std::addressof(*X.origin()),{nCV,nwalk});
    wfn2.vbias(G,X2,sqrtdt);
    Xsum=0;
    if(std::abs(file_data.Xsum)>1e-8) {
      for(int n=0; n<nwalk; n++) {
        Xsum=0;
        if(TGwfn.TG_local().root())
          for(int i=0; i<X2.shape()[0]; i++)
            Xsum += X2[i][n];
        Xsum = ( TGwfn.TG() += Xsum );
        REQUIRE( real(Xsum) == Approx(real(file_data.Xsum)) );
        REQUIRE( imag(Xsum) == Approx(imag(file_data.Xsum)) );
      }
    } else {
      Xsum=0;
      if(TGwfn.TG_local().root())
        for(int i=0; i<X2.shape()[0]; i++)
          Xsum += X2[i][0];
      Xsum = ( TGwfn.TG() += Xsum );
      app_log()<<" Xsum: " <<setprecision(12) <<Xsum <<std::endl;
    }

    // vbias must be reduced if false
    if(not wfn.distribution_over_cholesky_vectors()) {
      boost::multi::array<ComplexType,2> T({nCV,nwalk});
      if(TGwfn.TG_local().root())
        std::copy_n(X2.origin(),X2.num_elements(),T.origin());
      else
        std::fill_n(T.origin(),T.num_elements(),ComplexType(0.0,0.0));
      TGwfn.TG().all_reduce_in_place_n(std::addressof(*T.origin()),T.num_elements(),std::plus<>());
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
        Vsum = ( TGwfn.TG() += Vsum );
        REQUIRE( real(Vsum) == Approx(real(file_data.Vsum)) );
        REQUIRE( imag(Vsum) == Approx(imag(file_data.Vsum)) );
      }
    } else {
      Vsum=0;
      if(TGwfn.TG_local().root())
        for(int i=0; i<vHS.shape()[0]; i++)
          Vsum += vHS[i][0];
      Vsum = ( TGwfn.TG() += Vsum );
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
  OHMMS::Controller->initialize(0, NULL);
  auto world = boost::mpi3::environment::get_world_instance();

  if(not file_exists("./afqmc_msd.h5") ||
     not file_exists("./wfn_msd.dat") ) {
    app_log()<<" Skipping wfn_fac_collinear_multidet text. afqmc_msd.h5 and ./wfn_msd.dat files not found. \n";
  } else {

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
                         initial_guess[1](initial_guess.extension(1),{0,NAEB}));

    // no guarantee that overlap is 1.0
    wfn.Overlap(wset);

    using shmCMatrix = boost::multi::array<ComplexType,2,shared_allocator<ComplexType>>;

    wfn.Energy(wset);
    if(std::abs(file_data.E0+file_data.E1+file_data.E2)>1e-8) {
      for(auto it = wset.begin(); it!=wset.end(); ++it) {
        REQUIRE( real(it->E1()) == Approx(real(file_data.E0+file_data.E1)));
        REQUIRE( real(it->EXX()+it->EJ()) == Approx(real(file_data.E2)));
        REQUIRE( imag(it->energy()) == Approx(imag(file_data.E0+file_data.E1+file_data.E2)));
      }
    } else {
      app_log()<<" E: " <<wset[0].E1() <<" " <<wset[0].EXX() <<" " <<wset[0].EJ() <<std::endl;
    }

    auto size_of_G = wfn.size_of_G_for_vbias();
    int Gdim1 = (wfn.transposed_G_for_vbias()?nwalk:size_of_G);
    int Gdim2 = (wfn.transposed_G_for_vbias()?size_of_G:nwalk);
    shmCMatrix G({Gdim1,Gdim2},shared_allocator<ComplexType>{TG.TG_local()});
    wfn.MixedDensityMatrix_for_vbias(wset,G);

    double sqrtdt = std::sqrt(0.01);
    auto nCV = wfn.local_number_of_cholesky_vectors();
    shmCMatrix X({nCV,nwalk},shared_allocator<ComplexType>{TG.TG_local()});
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


    int vdim1 = (wfn.transposed_vHS()?nwalk:NMO*NMO);
    int vdim2 = (wfn.transposed_vHS()?NMO*NMO:nwalk);
    shmCMatrix vHS({vdim1,vdim2},shared_allocator<ComplexType>{TG.TG_local()});
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
  OHMMS::Controller->initialize(0, NULL);
  auto world = boost::mpi3::environment::get_world_instance();

  if(not file_exists("./afqmc_msd.h5") ||
     not file_exists("./wfn_msd.dat") ) {
    app_log()<<" Skipping wfn_fac_collinear_multidet text. afqmc_msd.h5 and ./wfn_msd.dat files not found. \n";
  } else {

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
                         initial_guess[1](initial_guess.extension(1),{0,NAEB}));

    // no guarantee that overlap is 1.0
    wfn.Overlap(wset);

    using shmCMatrix = boost::multi::array<ComplexType,2,shared_allocator<ComplexType>>;

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
    int Gdim1 = (wfn.transposed_G_for_vbias()?nwalk:size_of_G);
    int Gdim2 = (wfn.transposed_G_for_vbias()?size_of_G:nwalk);
    shmCMatrix G({Gdim1,Gdim2},shared_allocator<ComplexType>{TG.TG_local()});
    wfn.MixedDensityMatrix_for_vbias(wset,G);

    double sqrtdt = std::sqrt(0.01);
    auto nCV = wfn.local_number_of_cholesky_vectors();
    shmCMatrix X({nCV,nwalk},shared_allocator<ComplexType>{TG.TG_local()});
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


    int vdim1 = (wfn.transposed_vHS()?nwalk:NMO*NMO);
    int vdim2 = (wfn.transposed_vHS()?NMO*NMO:nwalk);
    shmCMatrix vHS({vdim1,vdim2},shared_allocator<ComplexType>{TG.TG_local()});
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

TEST_CASE("wfn_fac_collinear_phmsd", "[wavefunction_factory]")
{
  OHMMS::Controller->initialize(0, NULL);
  auto world = boost::mpi3::environment::get_world_instance();

  if(not file_exists("./afqmc_phmsd.h5") ||
     not file_exists("./wfn_phmsd.dat") ) {
    app_log()<<" Skipping wfn_fac_collinear_phmsd text. afqmc_msd.h5 and ./wfn_msd.dat files not found. \n";
  } else {

    // Global Task Group
    GlobalTaskGroup gTG(world);

    auto file_data = read_test_results_from_hdf<ValueType>("./afqmc_phmsd.h5");
    int NMO=file_data.NMO;
    int NAEA=file_data.NAEA;
    int NAEB=file_data.NAEB;

    std::map<std::string,AFQMCInfo> InfoMap;
    InfoMap.insert ( std::pair<std::string,AFQMCInfo>("info0",AFQMCInfo{"info0",NMO,NAEA,NAEB}) );
    HamiltonianFactory HamFac(InfoMap);
    const char *ham_xml_block =
"<Hamiltonian name=\"ham0\" info=\"info0\"> \
    <parameter name=\"filetype\">hdf5</parameter> \
    <parameter name=\"filename\">./afqmc_phmsd.h5</parameter> \
    <parameter name=\"cutoff_decomposition\">1e-5</parameter> \
    <parameter name=\"useHalfRotatedMuv\">no</parameter> \
  </Hamiltonian> \
";
    Libxml2Document doc;
    bool okay = doc.parseFromString(ham_xml_block);
    REQUIRE(okay);
    std::string ham_name("ham0");
    HamFac.push(ham_name,doc.getRoot());
    Hamiltonian& ham = HamFac.getHamiltonian(gTG,ham_name);


    auto TG = TaskGroup_(gTG,std::string("WfnTG"),1,gTG.getTotalCores());
    //auto TG = TaskGroup_(gTG,std::string("WfnTG"),1,1);
    //int nwalk = 1; // choose prime number to force non-trivial splits in shared routines
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
"<Wavefunction name=\"wfn0\" type=\"phmsd\" info=\"info0\"> \
      <parameter name=\"filetype\">ascii</parameter> \
      <parameter name=\"filename\">./wfn_phmsd.dat</parameter> \
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

    const char *wfn_xml_block2 =
"<Wavefunction name=\"wfn1\" type=\"nomsd\" info=\"info0\"> \
      <parameter name=\"filetype\">ascii</parameter> \
      <parameter name=\"filename\">./wfn_phmsd.dat</parameter> \
      <parameter name=\"cutoff\">1e-6</parameter> \
  </Wavefunction> \
";

//#define __compare__
#ifdef __compare__
    Libxml2Document doc2_;
    okay = doc2_.parseFromString(wfn_xml_block2);
    REQUIRE(okay);
    std::string wfn2_name("wfn1");
    WfnFac.push(wfn2_name,doc2_.getRoot());
    Wavefunction& nomsd = WfnFac.getWavefunction(TG,TG,wfn2_name,COLLINEAR,&ham,1e-6,nwalk);

#endif

    WalkerSet wset(TG,doc3.getRoot(),InfoMap["info0"],&rng);
    auto initial_guess = WfnFac.getInitialGuess(wfn_name);
    REQUIRE(initial_guess.shape()[0]==2);
    REQUIRE(initial_guess.shape()[1]==NMO);
    REQUIRE(initial_guess.shape()[2]==NAEA);


    std::default_random_engine generator;
    std::uniform_real_distribution<double> distribution(-0.005,0.005);
    for(int i=0; i<NMO; i++) {
      for(int j=0; j<NAEA; j++)
        initial_guess[0][i][j] += distribution(generator);
      for(int j=0; j<NAEB; j++)
        initial_guess[1][i][j] = initial_guess[0][i][j];
        //initial_guess[1][i][j] += distribution(generator);
    }
    wset.resize(nwalk,initial_guess[0],
                         initial_guess[1](initial_guess.extension(1),{0,NAEB}));
    qmcplusplus::Timer Time;
    // no guarantee that overlap is 1.0
    double t1;
#ifdef __compare__
    Time.restart();
    nomsd.Overlap(wset);
    t1=Time.elapsed();
    app_log()<<" NOMSD Overlap: " <<setprecision(12) <<wset[0].overlap() <<" " <<t1 <<std::endl;
#endif
    Time.restart();
    wfn.Overlap(wset);
    t1=Time.elapsed();
    app_log()<<" PHMSD Overlap: " <<setprecision(12) <<wset[0].overlap() <<" " <<t1 <<std::endl;
    for(int i=1; i<nwalk; i++) {
      REQUIRE( real(wset[0].overlap()) == Approx(real(wset[i].overlap())));
      REQUIRE( imag(wset[0].overlap()) == Approx(imag(wset[i].overlap())));
    }

    using shmCMatrix = boost::multi::array<ComplexType,2,shared_allocator<ComplexType>>;

    Time.restart();
    wfn.Energy(wset);
    t1=Time.elapsed();
    app_log()<<" PHMSD E: " <<setprecision(12) <<wset[0].energy() <<" "
             <<wset[0].E1() <<" " <<wset[0].EXX() <<" " <<wset[0].EJ() <<" " <<t1 <<std::endl;
    for(int i=1; i<nwalk; i++) {
      REQUIRE( real(wset[0].E1()) == Approx(real(wset[i].E1())));
      REQUIRE( imag(wset[0].E1()) == Approx(imag(wset[i].E1())));
      REQUIRE( real(wset[0].EJ()) == Approx(real(wset[i].EJ())));
      REQUIRE( imag(wset[0].EJ()) == Approx(imag(wset[i].EJ())));
      REQUIRE( real(wset[0].EXX()) == Approx(real(wset[i].EXX())));
      REQUIRE( imag(wset[0].EXX()) == Approx(imag(wset[i].EXX())));
    }
#ifdef __compare__
    Time.restart();
    nomsd.Energy(wset);
    t1=Time.elapsed();
    app_log()<<" NOMSD E: " <<setprecision(12) <<wset[0].energy() <<" "
             <<wset[0].E1() <<" " <<wset[0].EXX() <<" " <<wset[0].EJ() <<" " <<t1  <<std::endl;
     
      shmCMatrix Gph({2*NMO*NMO,nwalk},shared_allocator<ComplexType>{TG.TG_local()});
      shmCMatrix Gno({2*NMO*NMO,nwalk},shared_allocator<ComplexType>{TG.TG_local()});
      wfn.MixedDensityMatrix(wset,Gph,false,false);
      nomsd.MixedDensityMatrix(wset,Gno,false,false);
      std::cout<<" Comparing G \n";
      for(int i=0; i<NMO; i++)
       for(int j=0; j<NMO; j++)
        if(std::abs(Gph[i*NMO+j][0]-Gno[i*NMO+j][0]) > 1e-8)
          std::cout<<i <<" " <<j <<" " <<Gph[i*NMO+j][0] <<" " <<Gno[i*NMO+j][0] <<" "
                   <<std::abs(Gph[i*NMO+j][0]-Gno[i*NMO+j][0]) <<std::endl;
#endif

    auto size_of_G = wfn.size_of_G_for_vbias();
    int Gdim1 = (wfn.transposed_G_for_vbias()?nwalk:size_of_G);
    int Gdim2 = (wfn.transposed_G_for_vbias()?size_of_G:nwalk);
    shmCMatrix G({Gdim1,Gdim2},shared_allocator<ComplexType>{TG.TG_local()});
    wfn.MixedDensityMatrix_for_vbias(wset,G);
/*
std::cout<<" G: \n";
if(wfn.transposed_G_for_vbias())
 for(int i=0; i<size_of_G; i++)
  std::cout<<i <<" " <<G[0][i] <<"\n";
else
 for(int i=0; i<size_of_G; i++)
  std::cout<<i <<" " <<G[i][0] <<"\n";
*/
    double sqrtdt = std::sqrt(0.01);
    auto nCV = wfn.local_number_of_cholesky_vectors();
    shmCMatrix X({nCV,nwalk},shared_allocator<ComplexType>{TG.TG_local()});
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
      ComplexType Xsum2=0;
      for(int i=0; i<X.shape()[0]; i++) {
        Xsum += X[i][0];
        Xsum2 += 0.5*X[i][0]*X[i][0];
      }
      app_log()<<" Xsum: " <<setprecision(12) <<Xsum <<std::endl;
      app_log()<<" Xsum2 (EJ): " <<setprecision(12) <<Xsum2/sqrtdt/sqrtdt <<std::endl;
      for(int n=1; n<nwalk; n++) {
        ComplexType Xsum_=0;
        for(int i=0; i<X.shape()[0]; i++)
          Xsum_ += X[i][n];
        REQUIRE( real(Xsum) == Approx(real(Xsum_)) );
        REQUIRE( imag(Xsum) == Approx(imag(Xsum_)) );
      }
    }

    int vdim1 = (wfn.transposed_vHS()?nwalk:NMO*NMO);
    int vdim2 = (wfn.transposed_vHS()?NMO*NMO:nwalk);
    shmCMatrix vHS({vdim1,vdim2},shared_allocator<ComplexType>{TG.TG_local()});
    wfn.vHS(X,vHS,sqrtdt);
    TG.local_barrier();
    ComplexType Vsum=0;
    if(std::abs(file_data.Vsum)>1e-8) {
      for(int n=0; n<nwalk; n++) {
        Vsum=0;
        if(wfn.transposed_vHS())
          for(int i=0; i<vHS.shape()[1]; i++)
            Vsum += vHS[n][i];
        else
          for(int i=0; i<vHS.shape()[0]; i++)
            Vsum += vHS[i][n];
        REQUIRE( real(Vsum) == Approx(real(file_data.Vsum)) );
        REQUIRE( imag(Vsum) == Approx(imag(file_data.Vsum)) );
      }
    } else {
      Vsum=0;
      if(wfn.transposed_vHS())
        for(int i=0; i<vHS.shape()[1]; i++)
          Vsum += vHS[0][i];
      else
        for(int i=0; i<vHS.shape()[0]; i++)
          Vsum += vHS[i][0];
      app_log()<<" Vsum: " <<setprecision(12) <<Vsum <<std::endl;
      for(int n=1; n<nwalk; n++) {
        ComplexType Vsum_=0;
        if(wfn.transposed_vHS())
          for(int i=0; i<vHS.shape()[1]; i++)
            Vsum_ += vHS[n][i];
        else
          for(int i=0; i<vHS.shape()[0]; i++)
            Vsum_ += vHS[i][n];
        REQUIRE( real(Vsum) == Approx(real(Vsum_)) );
        REQUIRE( imag(Vsum) == Approx(imag(Vsum_)) );
      }
    }

    boost::multi::array<ComplexType,1> vMF(extensions<1u>{nCV});
    wfn.vMF(vMF);
    ComplexType vMFsum=0;
    {
      vMFsum=0;
      for(int i=0; i<vMF.shape()[0]; i++)
        vMFsum += vMF[i];
      app_log()<<" vMFsum: " <<setprecision(12) <<vMFsum <<std::endl;
    }

#ifdef __compare__
    assert(nCV == nomsd.local_number_of_cholesky_vectors());

      auto size_of_G2 = nomsd.size_of_G_for_vbias();
      int Gdim1_ = (nomsd.transposed_G_for_vbias()?nwalk:size_of_G2);
      int Gdim2_ = (nomsd.transposed_G_for_vbias()?size_of_G2:nwalk);
      shmCMatrix G_({Gdim1_,Gdim2_},shared_allocator<ComplexType>{TG.TG_local()});
      nomsd.MixedDensityMatrix_for_vbias(wset,G_);
/*
      std::cout<<" Comparing G \n";
      for(int i=0; i<NMO; i++)
       for(int j=0; j<NMO; j++) {
        if(std::abs(Gph[i*NMO+j][0]-G_[i*NMO+j][0]) > 1e-8)
          std::cout<<i <<" - " <<j <<" " <<Gph[i*NMO+j][0] <<" " <<G_[i*NMO+j][0] <<" "
                   <<std::abs(Gph[i*NMO+j][0]-G_[i*NMO+j][0]) <<std::endl;
        if(std::abs(G_[i*NMO+j][0]-Gno[i*NMO+j][0]) > 1e-8)
          std::cout<<i <<" -- " <<j <<" " <<G_[i*NMO+j][0] <<" " <<Gno[i*NMO+j][0] <<" "
                   <<std::abs(G_[i*NMO+j][0]-Gno[i*NMO+j][0]) <<std::endl;
       }
*/
      boost::multi::array_ref<ComplexType,2> X2(std::addressof(*X.origin())+nCV*nwalk,{nCV,nwalk});
      nomsd.vbias(G_,X2,sqrtdt);
      Xsum=0;
      ComplexType Xsum2(0.0);
      for(int i=0; i<X2.shape()[0]; i++) {
        Xsum += X2[i][0];
        Xsum2 += 0.5*X2[i][0]*X2[i][0];
      }
      app_log()<<" Xsum (NOMSD): " <<setprecision(12) <<Xsum <<std::endl;
      app_log()<<" Xsum2 (EJ): " <<setprecision(12) <<Xsum2/sqrtdt/sqrtdt <<std::endl;

    int vdim1_ = (nomsd.transposed_vHS()?nwalk:NMO*NMO);
    int vdim2_ = (nomsd.transposed_vHS()?NMO*NMO:nwalk);
    shmCMatrix vHS_({vdim1_,vdim2_},shared_allocator<ComplexType>{TG.TG_local()});
    nomsd.vHS(X2,vHS_,sqrtdt);
    TG.local_barrier();
    Vsum=0;
    if(std::abs(file_data.Vsum)>1e-8) {
      for(int n=0; n<nwalk; n++) {
        Vsum=0;
        if(nomsd.transposed_vHS())
          for(int i=0; i<vHS_.shape()[1]; i++)
            Vsum += vHS_[0][i];
        else
          for(int i=0; i<vHS_.shape()[0]; i++)
            Vsum += vHS_[i][0];
        REQUIRE( real(Vsum) == Approx(real(file_data.Vsum)) );
        REQUIRE( imag(Vsum) == Approx(imag(file_data.Vsum)) );
      }
    } else {
      Vsum=0;
      if(nomsd.transposed_vHS())
        for(int i=0; i<vHS_.shape()[1]; i++)
          Vsum += vHS_[0][i];
      else
        for(int i=0; i<vHS_.shape()[0]; i++)
          Vsum += vHS_[i][0];
      app_log()<<" Vsum: " <<setprecision(12) <<Vsum <<std::endl;
    }

    boost::multi::array<ComplexType,1> vMF2(extensions<1u>{nCV});
    nomsd.vMF(vMF2);
    vMFsum=0;
    {
      vMFsum=0;
      app_log()<<" vMF: " <<std::endl;
      for(int i=0; i<vMF2.shape()[0]; i++) {
        vMFsum += vMF2[i];
//        if(std::abs(vMF[i]-vMF2[i]) > 1e-8)
//          app_log()<<i <<": " <<setprecision(12) <<vMF[i] <<" " <<vMF2[i] <<" " <<std::abs(vMF[i]-vMF2[i]) <<std::endl;
      }
      app_log()<<" vMFsum: " <<setprecision(12) <<vMFsum <<std::endl;
    }

#endif
  }
}

}
