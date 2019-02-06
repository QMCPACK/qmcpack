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
#include "AFQMC/Utilities/myTimer.h"
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

template<class Allocator>
void propg_fac_shared(boost::mpi3::communicator & world)
{
  using pointer = device_ptr<ComplexType>;

  if(not file_exists("./afqmc.h5") ||
     not file_exists("./wfn.dat") ) {
    app_log()<<" Skipping propg_fac_shared. afqmc.h5 and ./wfn.dat files not found. \n";
  } else {

    TimerManager.set_timer_threshold(timer_level_coarse);
    setup_timers(AFQMCTimers, AFQMCTimerNames,timer_level_coarse);

    // Global Task Group
    afqmc::GlobalTaskGroup gTG(world);

    int NMO,NAEA,NAEB;
    std::tie(NMO,NAEA,NAEB) = read_info_from_hdf("./afqmc.h5");

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

    auto TG = TaskGroup_(gTG,std::string("WfnTG"),1,gTG.getTotalCores());
    int nwalk = 11; // choose prime number to force non-trivial splits in shared routines
    RandomGenerator_t rng;

    Allocator alloc_(make_localTG_allocator<ComplexType>(TG));

const char *wlk_xml_block =
"<WalkerSet name=\"wset0\">  \
  <parameter name=\"walker_type\">closed</parameter>  \
</WalkerSet> \
";
    Libxml2Document doc3;
    okay = doc3.parseFromString(wlk_xml_block);
    REQUIRE(okay);

    const char *wfn_xml_block =
"<Wavefunction name=\"wfn0\" info=\"info0\"> \
      <parameter name=\"filetype\">ascii</parameter> \
      <parameter name=\"filename\">./wfn.dat</parameter> \
      <parameter name=\"cutoff\">1e-6</parameter> \
  </Wavefunction> \
";
    Libxml2Document doc2;
    okay = doc2.parseFromString(wfn_xml_block);
    REQUIRE(okay);
    std::string wfn_name("wfn0");
    WavefunctionFactory WfnFac(InfoMap);
    WfnFac.push(wfn_name,doc2.getRoot());
    Wavefunction& wfn = WfnFac.getWavefunction(TG,TG,wfn_name,CLOSED,&ham,1e-6,nwalk);

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

TEST_CASE("propg_fac_distributed", "[propagator_factory]")
{
  OHMMS::Controller->initialize(0, NULL);
  auto world = boost::mpi3::environment::get_world_instance();

  if(not file_exists("./afqmc.h5") ||
     not file_exists("./wfn.dat") ) {
    app_log()<<" Skipping propg_fac_shared. afqmc.h5 and ./wfn.dat files not found. \n";
  } else {

    TimerManager.set_timer_threshold(timer_level_coarse);
    setup_timers(AFQMCTimers, AFQMCTimerNames,timer_level_coarse);

    // Global Task Group
    afqmc::GlobalTaskGroup gTG(world);

    int NMO,NAEA,NAEB;
    std::tie(NMO,NAEA,NAEB) = read_info_from_hdf("./afqmc.h5");

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
    auto TGprop = TaskGroup_(gTG,std::string("WfnTG"),gTG.getTotalNodes(),gTG.getTotalCores());
    int nwalk = 11; // choose prime number to force non-trivial splits in shared routines
    RandomGenerator_t rng;

const char *wlk_xml_block =
"<WalkerSet name=\"wset0\">  \
  <parameter name=\"walker_type\">closed</parameter>  \
</WalkerSet> \
";
    Libxml2Document doc3;
    okay = doc3.parseFromString(wlk_xml_block);
    REQUIRE(okay);

    const char *wfn_xml_block =
"<Wavefunction name=\"wfn0\" info=\"info0\"> \
      <parameter name=\"filetype\">ascii</parameter> \
      <parameter name=\"filename\">./wfn.dat</parameter> \
      <parameter name=\"cutoff\">1e-6</parameter> \
  </Wavefunction> \
";
    Libxml2Document doc2;
    okay = doc2.parseFromString(wfn_xml_block);
    REQUIRE(okay);
    std::string wfn_name("wfn0");
    WavefunctionFactory WfnFac(InfoMap);
    WfnFac.push(wfn_name,doc2.getRoot());
    Wavefunction& wfn = WfnFac.getWavefunction(TGprop,TG,wfn_name,CLOSED,&ham,1e-6,nwalk);

    WalkerSet wset(TG,doc3.getRoot(),InfoMap["info0"],&rng);
    auto initial_guess = WfnFac.getInitialGuess(wfn_name);
    REQUIRE(initial_guess.size(0)==2);
    REQUIRE(initial_guess.size(1)==NMO);
    REQUIRE(initial_guess.size(2)==NAEA);
    wset.resize(nwalk,initial_guess[0],initial_guess[0]);
//                         initial_guess[1](XXX.extension(0),{0,NAEB}));

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

  }
}

TEST_CASE("propg_fac_shared", "[propagator_factory]")
{
  OHMMS::Controller->initialize(0, NULL);
  auto world = boost::mpi3::environment::get_world_instance();

#ifdef QMC_CUDA
  auto node = world.split_shared(world.rank());

  qmc_cuda::CUDA_INIT(node);
  using Alloc = qmc_cuda::cuda_gpu_allocator<ComplexType>;
#else
  using Alloc = shared_allocator<ComplexType>;
#endif

  propg_fac_shared<Alloc>(world);

}

}
