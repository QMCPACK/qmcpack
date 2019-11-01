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

#undef NDEBUG

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
//#include "AFQMC/Utilities/myTimer.h"
#include "Utilities/Timer.h"
#include "AFQMC/Hamiltonians/HamiltonianFactory.h"
#include "AFQMC/Hamiltonians/Hamiltonian.hpp"
#include "AFQMC/Wavefunctions/WavefunctionFactory.h"
#include "AFQMC/Wavefunctions/Wavefunction.hpp"
#include "AFQMC/Propagators/PropagatorFactory.h"
#include "AFQMC/Propagators/Propagator.hpp"
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

void propg_fac_shared(boost::mpi3::communicator & world)
{
  if(not file_exists(UTEST_HAMIL) ||
     not file_exists(UTEST_WFN) ) {
    app_log()<<" Skipping ham_ops_basic_serial. Hamiltonian or wavefunction file not found. \n";
    app_log()<<" Run unit test with --hamil /path/to/hamil.h5 and --wfn /path/to/wfn.dat.\n";
  } else {

    TimerManager.set_timer_threshold(timer_level_coarse);
    setup_timers(AFQMCTimers, AFQMCTimerNames,timer_level_coarse);

    // Global Task Group
    afqmc::GlobalTaskGroup gTG(world);

    int NMO,NAEA,NAEB;
    std::tie(NMO,NAEA,NAEB) = read_info_from_hdf(UTEST_HAMIL);

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
    const char *ham_xml_block = hamil_xml.c_str();
    Libxml2Document doc;
    bool okay = doc.parseFromString(ham_xml_block);
    REQUIRE(okay);
    std::string ham_name("ham0");
    HamFac.push(ham_name,doc.getRoot());
    Hamiltonian& ham = HamFac.getHamiltonian(gTG,ham_name);

    auto TG = TaskGroup_(gTG,std::string("WfnTG"),1,gTG.getTotalCores());
    int nwalk = 11; // choose prime number to force non-trivial splits in shared routines
    RandomGenerator_t rng;

    WALKER_TYPES type = afqmc::getWalkerType(UTEST_WFN);
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

    // Create unique restart filename to avoid issues with running tests in parallel
    // through ctest.
    std::string restart_file = create_test_hdf(UTEST_WFN, UTEST_HAMIL);
    std::string wfn_xml =
"<Wavefunction name=\"wfn0\" info=\"info0\"> \
      <parameter name=\"filetype\">ascii</parameter> \
      <parameter name=\"filename\">"+UTEST_WFN+"</parameter> \
      <parameter name=\"cutoff\">1e-6</parameter> \
      <parameter name=\"restart_file\">"+restart_file+"</parameter> \
  </Wavefunction> \
";
    const char *wfn_xml_block = wfn_xml.c_str();
    Libxml2Document doc2;
    okay = doc2.parseFromString(wfn_xml_block);
    REQUIRE(okay);
    std::string wfn_name("wfn0");
    WavefunctionFactory WfnFac(InfoMap);
    WfnFac.push(wfn_name,doc2.getRoot());
    Wavefunction& wfn = WfnFac.getWavefunction(TG,TG,wfn_name,type,&ham,1e-6,nwalk);

    WalkerSet wset(TG,doc3.getRoot(),InfoMap["info0"],&rng);
    auto initial_guess = WfnFac.getInitialGuess(wfn_name);
    REQUIRE(initial_guess.size(0)==2);
    REQUIRE(initial_guess.size(1)==NMO);
    REQUIRE(initial_guess.size(2)==NAEA);
    wset.resize(nwalk,initial_guess[0],initial_guess[0]);
//                         initial_guess[1](XXX.extension(0),{0,NAEB}));

const char *propg_xml_block =
"<Propagator name=\"prop0\">  \
</Propagator> \
";
    Libxml2Document doc4;
    okay = doc4.parseFromString(propg_xml_block);
    REQUIRE(okay);
    std::string prop_name("prop0");
    PropagatorFactory PropgFac(InfoMap);
    PropgFac.push(prop_name,doc4.getRoot());
    Propagator& prop = PropgFac.getPropagator(TG,prop_name,wfn,&rng);

std::cout<<setprecision(12);
    wfn.Energy(wset);
    {
      ComplexType eav=0,ov=0;
      for(auto it=wset.begin(); it!=wset.end(); ++it) {
        eav += *it->weight()*(it->energy());
        ov += *it->weight();
      }
      app_log()<<" Initial Energy: " <<(eav/ov).real() <<std::endl;
    }
    double tot_time=0;
    RealType dt=0.01;
    RealType Eshift=std::abs(ComplexType(*wset[0].overlap()));
    for(int i=0; i<4; i++) {
      prop.Propagate(2,wset,Eshift,dt,1);
      wfn.Energy(wset);
      ComplexType eav=0,ov=0;
      for(auto it=wset.begin(); it!=wset.end(); ++it) {
        eav += *it->weight()*(it->energy());
        ov += *it->weight();
      }
      tot_time+=2*dt;
      app_log()<<" -- " <<i <<" " <<tot_time <<" " <<(eav/ov).real() <<std::endl;
      wfn.Orthogonalize(wset,true);
    }
    for(int i=0; i<4; i++) {
      prop.Propagate(4,wset,Eshift,dt,1);
      wfn.Energy(wset);
      ComplexType eav=0,ov=0;
      for(auto it=wset.begin(); it!=wset.end(); ++it) {
        eav += *it->weight()*(it->energy());
        ov += *it->weight();
      }
      tot_time+=4*dt;
      app_log()<<" -- " <<i <<" " <<tot_time <<" " <<(eav/ov).real() <<std::endl;
      wfn.Orthogonalize(wset,true);
    }

    for(int i=0; i<4; i++) {
      prop.Propagate(4,wset,Eshift,dt,2);
      wfn.Energy(wset);
      ComplexType eav=0,ov=0;
      for(auto it=wset.begin(); it!=wset.end(); ++it) {
        eav += *it->weight()*(it->energy());
        ov += *it->weight();
      }
      tot_time+=4*dt;
      app_log()<<" -- " <<i <<" " <<tot_time <<" " <<(eav/ov).real() <<std::endl;
      wfn.Orthogonalize(wset,true);
    }
    for(int i=0; i<4; i++) {
      prop.Propagate(5,wset,Eshift,2*dt,2);
      wfn.Energy(wset);
      ComplexType eav=0,ov=0;
      for(auto it=wset.begin(); it!=wset.end(); ++it) {
        eav += *it->weight()*(it->energy());
        ov += *it->weight();
      }
      tot_time+=5*2*dt;
      app_log()<<" -- " <<i <<" " <<tot_time <<" " <<(eav/ov).real() <<std::endl;
      wfn.Orthogonalize(wset,true);
    }

    TimerManager.print(nullptr);

  }
}

void propg_fac_distributed(boost::mpi3::communicator & world, int ngrp)
{

  if(not file_exists(UTEST_HAMIL) ||
     not file_exists(UTEST_WFN) ) {
    app_log()<<" Skipping ham_ops_basic_serial. Hamiltonian or wavefunction file not found. \n";
    app_log()<<" Run unit test with --hamil /path/to/hamil.h5 and --wfn /path/to/wfn.dat.\n";
  } else {

    TimerManager.set_timer_threshold(timer_level_coarse);
    setup_timers(AFQMCTimers, AFQMCTimerNames,timer_level_coarse);

    // Global Task Group
    afqmc::GlobalTaskGroup gTG(world);

    int NMO,NAEA,NAEB;
    std::tie(NMO,NAEA,NAEB) = read_info_from_hdf(UTEST_HAMIL);

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
    const char *ham_xml_block = hamil_xml.c_str();
    Libxml2Document doc;
    bool okay = doc.parseFromString(ham_xml_block);
    REQUIRE(okay);
    std::string ham_name("ham0");
    HamFac.push(ham_name,doc.getRoot());
    Hamiltonian& ham = HamFac.getHamiltonian(gTG,ham_name);

    auto TG = TaskGroup_(gTG,std::string("WfnTG"),1,gTG.getTotalCores());
    auto TGprop = TaskGroup_(gTG,std::string("WfnTG"),ngrp,gTG.getTotalCores());
    //int nwalk = 4; // choose prime number to force non-trivial splits in shared routines
    int nwalk = 11; // choose prime number to force non-trivial splits in shared routines
    RandomGenerator_t rng;
    qmcplusplus::Timer Time;

    WALKER_TYPES type = afqmc::getWalkerType(UTEST_WFN);
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

    std::string restart_file = create_test_hdf(UTEST_WFN, UTEST_HAMIL);
    std::string wfn_xml =
"<Wavefunction name=\"wfn0\" info=\"info0\"> \
      <parameter name=\"filetype\">ascii</parameter> \
      <parameter name=\"filename\">"+UTEST_WFN+"</parameter> \
      <parameter name=\"cutoff\">1e-6</parameter> \
      <parameter name=\"restart_file\">"+restart_file+"</parameter> \
  </Wavefunction> \
";
    const char *wfn_xml_block = wfn_xml.c_str();
    Libxml2Document doc2;
    okay = doc2.parseFromString(wfn_xml_block);
    REQUIRE(okay);
    std::string wfn_name("wfn0");
    WavefunctionFactory WfnFac(InfoMap);
    WfnFac.push(wfn_name,doc2.getRoot());
    Wavefunction& wfn = WfnFac.getWavefunction(TGprop,TGprop,wfn_name,type,&ham,1e-6,nwalk);

    WalkerSet wset(TG,doc3.getRoot(),InfoMap["info0"],&rng);
    auto initial_guess = WfnFac.getInitialGuess(wfn_name);
    REQUIRE(initial_guess.size(0)==2);
    REQUIRE(initial_guess.size(1)==NMO);
    REQUIRE(initial_guess.size(2)==NAEA);
    wset.resize(nwalk,initial_guess[0],initial_guess[0]);

const char *propg_xml_block0 =
"<Propagator name=\"prop0\">  \
      <parameter name=\"nnodes\">";
const char *propg_xml_block1 =
"</parameter> \
</Propagator> \
";
    std::string str_ = std::string("<Propagator name=\"prop0\"> <parameter name=\"nnodes\">") +
                       std::to_string(gTG.getTotalNodes()) +
                       std::string("</parameter> </Propagator>");
    Libxml2Document doc4;
    okay = doc4.parseFromString(str_.c_str());
    REQUIRE(okay);
    std::string prop_name("prop0");
    PropagatorFactory PropgFac(InfoMap);
    PropgFac.push(prop_name,doc4.getRoot());
    Propagator& prop = PropgFac.getPropagator(TGprop,prop_name,wfn,&rng);
    wfn.Energy(wset);
    {
      ComplexType eav=0,ov=0;
      for(auto it=wset.begin(); it!=wset.end(); ++it) {
        eav += *it->weight()*(it->energy());
        ov += *it->weight();
      }
      app_log()<<" Initial Energy: " <<(eav/ov).real() <<std::endl;
    }
    double tot_time=0,t1;
    RealType dt=0.01;
    RealType Eshift=std::abs(ComplexType(*wset[0].overlap()));
    Time.restart();
    for(int i=0; i<4; i++) {
      prop.Propagate(2,wset,Eshift,dt,1);
      wfn.Energy(wset);
      ComplexType eav=0,ov=0;
      for(auto it=wset.begin(); it!=wset.end(); ++it) {
        eav += *it->weight()*(it->energy());
        ov += *it->weight();
      }
      tot_time+=2*dt;
      wfn.Orthogonalize(wset,true);
      t1=Time.elapsed();
      app_log()<<" -- " <<i <<" " <<tot_time <<" " <<(eav/ov).real() <<" Time: " <<t1 <<std::endl;
    }
    for(int i=0; i<4; i++) {
      prop.Propagate(4,wset,Eshift,dt,1);
      wfn.Energy(wset);
      ComplexType eav=0,ov=0;
      for(auto it=wset.begin(); it!=wset.end(); ++it) {
        eav += *it->weight()*(it->energy());
        ov += *it->weight();
      }
      tot_time+=4*dt;
      wfn.Orthogonalize(wset,true);
      t1=Time.elapsed();
      app_log()<<" -- " <<i <<" " <<tot_time <<" " <<(eav/ov).real() <<" Time: " <<t1 <<std::endl;
    }

    for(int i=0; i<4; i++) {
      prop.Propagate(4,wset,Eshift,dt,2);
      wfn.Energy(wset);
      ComplexType eav=0,ov=0;
      for(auto it=wset.begin(); it!=wset.end(); ++it) {
        eav += *it->weight()*(it->energy());
        ov += *it->weight();
      }
      tot_time+=4*dt;
      wfn.Orthogonalize(wset,true);
      t1=Time.elapsed();
      app_log()<<" -- " <<i <<" " <<tot_time <<" " <<(eav/ov).real() <<" Time: " <<t1 <<std::endl;
    }
    for(int i=0; i<4; i++) {
      prop.Propagate(5,wset,Eshift,2*dt,2);
      wfn.Energy(wset);
      ComplexType eav=0,ov=0;
      for(auto it=wset.begin(); it!=wset.end(); ++it) {
        eav += *it->weight()*(it->energy());
        ov += *it->weight();
      }
      tot_time+=5*2*dt;
      wfn.Orthogonalize(wset,true);
      t1=Time.elapsed();
      app_log()<<" -- " <<i <<" " <<tot_time <<" " <<(eav/ov).real() <<" Time: " <<t1 <<std::endl;
    }

    TimerManager.print(nullptr);
  }
}

TEST_CASE("propg_fac_shared", "[propagator_factory]")
{
  OHMMS::Controller->initialize(0, NULL);
  auto world = boost::mpi3::environment::get_world_instance();
  if(not world.root()) infoLog.pause();

#ifdef ENABLE_CUDA
  auto node = world.split_shared(world.rank());
  qmc_cuda::CUDA_INIT(node);
#endif

  propg_fac_shared(world);

}

TEST_CASE("propg_fac_distributed", "[propagator_factory]")
{
  OHMMS::Controller->initialize(0, NULL);
  auto world = boost::mpi3::environment::get_world_instance();
  if(not world.root()) infoLog.pause();

#ifdef ENABLE_CUDA
  auto node = world.split_shared(world.rank());
  int ngrp(world.size());
  qmc_cuda::CUDA_INIT(node);
#else
  auto node = world.split_shared(world.rank());
  int ngrp(world.size()/node.size());
#endif

  propg_fac_distributed(world,ngrp);

}

}
