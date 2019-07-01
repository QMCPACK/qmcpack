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
#include <Utilities/NewTimer.h>
#include "Utilities/Timer.h"
#include "Utilities/OutputManager.h"

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
void test_read_phmsd(boost::mpi3::communicator & world)
{

  using pointer = device_ptr<ComplexType>;

  if( not file_exists(UTEST_WFN) ) {
    app_log()<<" Skipping read_phmsd. Wavefunction file not found. \n";
    app_log()<<" Run unit test with --wfn /path/to/wfn.dat.\n";
  } else {

    // Global Task Group
    GlobalTaskGroup gTG(world);
    auto TG = TaskGroup_(gTG,std::string("WfnTG"),1,gTG.getTotalCores());
    auto TGwfn = TaskGroup_(gTG,std::string("WfnTG"),1,gTG.getTotalCores());
    Allocator alloc_(make_localTG_allocator<ComplexType>(TG));

    int NMO;
    int NAEA;
    int NAEB;
    std::tie(NMO,NAEA,NAEB) = read_info_from_wfn(UTEST_WFN, "PHMSD");
    hdf_archive dump;
    if(!dump.open(UTEST_WFN, H5F_ACC_RDONLY)) {
      app_error()<<"Error reading wavefunction file.\n";
    }
    if(!dump.push("Wavefunction",false)) {
      app_error()<<" Error in WavefunctionFactory: Group Wavefunction not found. \n";
      APP_ABORT("");
    }
    if(!dump.push("PHMSD",false)) {
      app_error()<<" Error in WavefunctionFactory: Group PHMSD not found. \n";
      APP_ABORT("");
    }
    int ndets_to_read = -1;
    std::string wfn_type;
    WALKER_TYPES walker_type = COLLINEAR;
    std::vector<PsiT_Matrix> PsiT_MO; // read_ph_wavefunction returns the full MO matrix
    ph_excitations<int,ComplexType> abij = read_ph_wavefunction_hdf(dump,ndets_to_read,walker_type,
                                                                    gTG.Node(),NMO,NAEA,NAEB,PsiT_MO,
                                                                    wfn_type);
    std::vector<int> buff(ndets_to_read*(NAEA+NAEB));
    boost::multi::array_ref<int,2> occs(to_address(buff.data()), {ndets_to_read, NAEA+NAEB});
    if(!dump.readEntry(buff, "occs"))
      APP_ABORT("Error reading occs array.\n");
    std::vector<ComplexType> ci_coeffs(ndets_to_read);
    if(!dump.readEntry(ci_coeffs, "ci_coeffs"))
      APP_ABORT("Error reading occs array.\n");
    using std::get;
    auto cit = abij.configurations_begin();
    std::vector<int> configa(NAEA), configb(NAEB);
    // Is it fortuitous that the order of determinants is the same?
    for(int nd = 0; nd < ndets_to_read; nd++, ++cit) {
      int alpha_ix = get<0>(*cit);
      int beta_ix = get<1>(*cit);
      auto ci = get<2>(*cit);
      abij.get_configuration(0, alpha_ix, configa);
      abij.get_configuration(1, beta_ix, configb);
      std::sort(configa.begin(),configa.end());
      std::sort(configb.begin(),configb.end());
      for(int i = 0; i < NAEA; i++) {
        REQUIRE(configa[i]==occs[nd][i]);
      }
      for(int i = 0; i < NAEB; i++) {
        REQUIRE(configb[i]==occs[nd][i+NAEA]+NMO);
      }
      REQUIRE(std::abs(ci_coeffs[nd])==std::abs(ci));
    }
    // Check sign of permutation.
    REQUIRE(abij.number_of_configurations() == ndets_to_read);
  }

}

void getBasicWavefunction(std::vector<int>& occs, std::vector<ComplexType>& coeffs, int NEL)
{
    hdf_archive dump;
    if(!dump.open(UTEST_WFN, H5F_ACC_RDONLY)) {
      app_error()<<"Error reading wavefunction file.\n";
    }
    if(!dump.push("Wavefunction",false)) {
      app_error()<<" Error in WavefunctionFactory: Group Wavefunction not found. \n";
      APP_ABORT("");
    }
    if(!dump.push("PHMSD",false)) {
      app_error()<<" Error in WavefunctionFactory: Group PHMSD not found. \n";
      APP_ABORT("");
    }
    std::vector<int> Idata(5);
    if(!dump.readEntry(Idata, "dims"))
      APP_ABORT("Errro reading dims array\n");
    int ndets = Idata[4];
    occs.resize(ndets*NEL);
    if(!dump.readEntry(occs, "occs"))
      APP_ABORT("Error reading occs array.\n");
    std::vector<ComplexType> ci_coeffs(ndets);
    if(!dump.readEntry(coeffs, "ci_coeffs"))
      APP_ABORT("Error reading occs array.\n");
}

// Construct PsiT^{dagger}
void getSlaterMatrix(boost::multi::array<ComplexType,2>& SM, boost::multi::array_ref<int,1>& occs, int NEL)
{
  using std::fill_n;
  fill_n(SM.origin(),SM.num_elements(),ComplexType(0.0));
  for(int i=0; i<NEL; i++)
    SM[i][occs[i]] = ComplexType(1.0);
}

// Very dumb.
int getExcitation(std::vector<int>& occi, std::vector<int> occj, std::vector<int>& excitation, int& perm, bool print=false)
{
  std::vector<int> from_orb, to_orb;
  // Work out which orbitals are excited from / to.
  std::set_difference(occj.begin(), occj.end(),
                      occi.begin(), occi.end(),
                      std::inserter(from_orb, from_orb.begin()));
  std::set_difference(occi.begin(), occi.end(),
                      occj.begin(), occj.end(),
                      std::inserter(to_orb, to_orb.begin()));
  int nexcit = from_orb.size();
  if(nexcit <= 2) {
    for(int i = 0; i < from_orb.size(); i++)
      excitation.push_back(from_orb[i]);
    for(int i = 0; i < to_orb.size(); i++)
      excitation.push_back(to_orb[i]);
    int nperm = 0;
    int nmove = 0;
    for(auto o : from_orb) {
      auto it = std::find(occj.begin(), occj.end(), o);
      int loc = std::distance(occj.begin(), it);
      nperm += occj.size() - loc - 1 + nmove;
      nmove += 1;
    }
    nmove = 0;
    for(auto o : to_orb) {
      auto it = std::find(occi.begin(), occi.end(), o);
      int loc = std::distance(occi.begin(), it);
      nperm += occi.size() - loc - 1 + nmove;
      nmove += 1;
    }
    if(nperm % 2 == 1)
      perm = -1;
  }
  return nexcit;
}

int decodeSpinOrbital(int spinOrb, int& spin)
{
  spin = spinOrb%2==0 ? 0 : 1;
  return spin ? (spinOrb-1) / 2 : spinOrb / 2;
}

ComplexType slaterCondon0(Hamiltonian& ham, std::vector<int>& det)
{
  ComplexType oneBody = ComplexType(0.0), twoBody = ComplexType(0.0);
  auto H1 = ham.getH1();
  int spini, spinj;
  for(auto i : det) {
    int oi = decodeSpinOrbital(i, spini);
    oneBody += H1[oi][oi];
    for(auto j : det) {
      int oj = decodeSpinOrbital(j, spinj);
      twoBody += ham.H(oi,oj,oi,oj);
      if(spini == spinj) twoBody -= ham.H(oi,oj,oj,oi);
    }
  }
  return oneBody + 0.5 * twoBody;
}

ComplexType slaterCondon1(Hamiltonian& ham, std::vector<int>& excit, std::vector<int>& det)
{
  int spini, spina;
  int oi = decodeSpinOrbital(excit[0], spini);
  int oa = decodeSpinOrbital(excit[1], spina);
  auto H1 = ham.getH1();
  ComplexType oneBody = H1[oi][oa];
  ComplexType twoBody = ComplexType(0.0);
  for(auto j : det) {
    int spinj;
    int oj = decodeSpinOrbital(j, spinj);
    if(j != excit[0]) {
      twoBody += ham.H(oi,oj,oa,oj);
      if(spini == spinj)
        twoBody -= ham.H(oi,oj,oj,oa);
    }
  }
  return oneBody + twoBody;
}

ComplexType slaterCondon2(Hamiltonian& ham, std::vector<int>& excit)
{
  ComplexType twoBody = ComplexType(0.0);
  int spini, spinj, spina, spinb;
  int oi = decodeSpinOrbital(excit[0], spini);
  int oj = decodeSpinOrbital(excit[1], spinj);
  int oa = decodeSpinOrbital(excit[2], spina);
  int ob = decodeSpinOrbital(excit[3], spinb);
  if(spini == spina)
    twoBody = ham.H(oi,oj,oa,ob);
  if(spini == spinb)
    twoBody -= ham.H(oi,oj,ob,oa);
  return twoBody;
}

void createDeterminant(boost::multi::array_ref<int,1>& occa, boost::multi::array_ref<int,1>& occb, std::vector<int>& det)
{
  for(auto i : occa)
    det.push_back(2*i);
  for(auto i : occb)
    det.push_back(2*i+1);
  std::sort(det.begin(), det.end());
}

void computeVariationalEnergy(Wavefunction& wfn, boost::multi::array_ref<int,2>& occs, Hamiltonian& ham, int NAEA, int NAEB)
{
  int ndets = occs.size(0);
  boost::multi::array<ComplexType,2> H({ndets,ndets});
  for(int idet = 0; idet < ndets; idet++) {
    boost::multi::array_ref<int,1> occi_a(occs[idet].origin(), {NAEA});
    boost::multi::array_ref<int,1> occi_b(occs[idet].origin()+NAEA, {NAEB});
    std::vector<int> deti;
    createDeterminant(occi_a, occi_b, deti);
    for(int jdet = 0; jdet < ndets; jdet++) {
      boost::multi::array_ref<int,1> occj_a(occs[jdet].origin(), {NAEA});
      boost::multi::array_ref<int,1> occj_b(occs[jdet].origin()+NAEA, {NAEB});
      std::vector<int> detj;
      createDeterminant(occj_a, occj_b, detj);
      int perm = 1;
      std::vector<int> excit;
      int nexcit = getExcitation(deti, detj, excit, perm);
      // Compute <Di|H|Dj>
      if(nexcit == 0) {
        H[idet][jdet] = ComplexType(perm)*slaterCondon0(ham, detj);
      } else if(nexcit == 1) {
        H[idet][jdet] = ComplexType(perm)*slaterCondon1(ham, excit, detj);
      } else if(nexcit == 2) {
        H[idet][jdet] = ComplexType(perm)*slaterCondon2(ham, excit);
      } else {
        H[idet][jdet] = ComplexType(0.0);
      }
    }
  }
  using TVec = boost::multi::array<RealType,1>;
  using TMat = boost::multi::array<ComplexType,2>;
  using eigSys = std::pair<TVec,TMat>;
  eigSys Sol = ma::symEig<TVec,TMat>(H);
  std::cout << "Variational Energy: " << Sol.first[0] << std::endl;
}

ComplexType contractOneBody(std::vector<int>& deti, std::vector<int>& detj, std::vector<int>& excit, int perm, boost::multi::array_ref<ComplexType,2>)
{
}

void computeMeanFieldShift(Wavefunction& wfn, boost::multi::array_ref<int,2>& occs, Hamiltonian& ham, int NAEA, int NAEB)
{
  int ndets = occs.size(0);
  boost::multi::array<ComplexType,2> H({ndets,ndets});
  ComplexType numer = ComplexType(0.0);
  ComplexType denom = ComplexType(0.0);
  for(int idet = 0; idet < ndets; idet++) {
    boost::multi::array_ref<int,1> occi_a(occs[idet].origin(), {NAEA});
    boost::multi::array_ref<int,1> occi_b(occs[idet].origin()+NAEA, {NAEB});
    std::vector<int> deti;
    createDeterminant(occi_a, occi_b, deti);
    for(int jdet = 0; jdet < ndets; jdet++) {
      boost::multi::array_ref<int,1> occj_a(occs[jdet].origin(), {NAEA});
      boost::multi::array_ref<int,1> occj_b(occs[jdet].origin()+NAEA, {NAEB});
      std::vector<int> detj;
      createDeterminant(occj_a, occj_b, detj);
      int perm = 1;
      std::vector<int> excit;
      int nexcit = getExcitation(deti, detj, excit, perm);
    }
  }
}

//template<class Allocator, class SDet_Type>
template<class Allocator>
//void test_phmsd(boost::mpi3::communicator & world, Allocator alloc)
void test_phmsd(boost::mpi3::communicator& world)
{

  using pointer = device_ptr<ComplexType>;

  if( not file_exists(UTEST_WFN) || not file_exists(UTEST_HAMIL)) {
    app_log()<<" Skipping read_phmsd. Wavefunction and/or Hamiltonian file not found. \n";
    app_log()<<" Run unit test with --wfn /path/to/wfn.h5 and --hamil /path/to/hamil.h5.\n";
  } else {

    // Global Task Group
    GlobalTaskGroup gTG(world);
    auto TG = TaskGroup_(gTG,std::string("WfnTG"),1,gTG.getTotalCores());
    auto TGwfn = TaskGroup_(gTG,std::string("WfnTG"),1,gTG.getTotalCores());
    Allocator alloc_(make_localTG_allocator<ComplexType>(TG));

    int NMO;
    int NAEA;
    int NAEB;
    std::tie(NMO,NAEA,NAEB) = read_info_from_wfn(UTEST_WFN, "PHMSD");
    // Test overlap.
    //wfn.Overlap(wset);
    WALKER_TYPES type = afqmc::getWalkerTypeHDF5(UTEST_WFN, "PHMSD");
    std::map<std::string,AFQMCInfo> InfoMap;
    InfoMap.insert ( std::pair<std::string,AFQMCInfo>("info0",AFQMCInfo{"info0",NMO,NAEA,NAEB}) );
    HamiltonianFactory HamFac(InfoMap);
    std::string hamil_xml =
  "<Hamiltonian name=\"ham0\" info=\"info0\"> \
  <parameter name=\"filetype\">hdf5</parameter> \
  <parameter name=\"filename\">"+UTEST_HAMIL+"</parameter> \
  <parameter name=\"cutoff_decomposition\">1e-12</parameter> \
  <parameter name=\"cutoff_1bar\">1e-12</parameter> \
  </Hamiltonian> \
  ";
    const char *ham_xml_block = hamil_xml.c_str();
    Libxml2Document doc;
    bool okay = doc.parseFromString(ham_xml_block);
    REQUIRE(okay);
    std::string ham_name("ham0");
    HamFac.push(ham_name,doc.getRoot());
    Hamiltonian& ham = HamFac.getHamiltonian(gTG,ham_name);

    std::string wfn_xml =
  "<Wavefunction name=\"wfn0\" info=\"info0\" type=\"phmsd\"> \
      <parameter name=\"filetype\">hdf5</parameter> \
      <parameter name=\"filename\">"+UTEST_WFN+"</parameter> \
      <parameter name=\"cutoff\">1e-6</parameter> \
  </Wavefunction> \
  ";
    const char *wfn_xml_block = wfn_xml.c_str();
    Libxml2Document doc2;
    okay = doc2.parseFromString(wfn_xml_block);
    REQUIRE(okay);
    std::string wfn_name("wfn0");
    WavefunctionFactory WfnFac(InfoMap);
    WfnFac.push(wfn_name,doc2.getRoot());
    int nwalk = 1;
    Wavefunction& wfn = WfnFac.getWavefunction(TGwfn,TGwfn,wfn_name,type,&ham,1e-6,nwalk);
    const char *wlk_xml_block =
    "<WalkerSet name=\"wset0\">  \
      <parameter name=\"walker_type\">collinear</parameter>  \
    </WalkerSet> \
    ";
    Libxml2Document doc3;
    okay = doc3.parseFromString(wlk_xml_block);
    REQUIRE(okay);
    RandomGenerator_t rng;
    WalkerSet wset(TG,doc3.getRoot(),InfoMap["info0"],&rng);
    auto initial_guess = WfnFac.getInitialGuess(wfn_name);
    REQUIRE(initial_guess.size(0)==2);
    REQUIRE(initial_guess.size(1)==NMO);
    REQUIRE(initial_guess.size(2)==NAEA);

    wset.resize(nwalk,initial_guess[0],
                initial_guess[1](initial_guess.extension(1),{0,NAEB}));
    // 1. Test Overlap Explicitly
    // 1.a Get raw occupancies and coefficients from file.
    std::vector<ComplexType> coeffs;
    std::vector<int> buff;
    getBasicWavefunction(buff, coeffs, NAEA+NAEB);
    int ndets = coeffs.size();
    boost::multi::array_ref<int,2> occs(buff.data(), {ndets,NAEA+NAEB});
    // 1.b Compute overlap of trial wavefunction compotents.
    boost::multi::array<ComplexType,2> Orbs({NMO,NMO});
    for(int i=0; i<NMO; i++)
      Orbs[i][i] = ComplexType(1.0);
    boost::multi::array<ComplexType,2> TrialA({NAEB,NMO}), TrialB({NAEB,NMO});
    auto sdet = wfn.getSlaterDetOperations();
    ComplexType ovlp;
    ComplexType ovlp_sum = ComplexType(0.0);
    ComplexType logovlp;
    for(int idet=0; idet<coeffs.size(); idet++) {
      // Construct slater matrix from given set of occupied orbitals.
      boost::multi::array_ref<int,1> oa(occs[idet].origin(), {NAEA});
      getSlaterMatrix(TrialA, oa, NAEA);
      boost::multi::array_ref<int,1> ob(occs[idet].origin()+NAEA, {NAEB});
      getSlaterMatrix(TrialB, ob, NAEB);
      ComplexType ovlpa = sdet->Overlap(TrialA, wset[0].SlaterMatrix(Alpha), logovlp);
      ComplexType ovlpb = sdet->Overlap(TrialB, wset[0].SlaterMatrix(Beta), logovlp);
      ovlp_sum += ma::conj(coeffs[idet])*ovlpa*ovlpb;
    }
    wfn.Overlap(wset);
    std::cout.precision(16);
    for(auto it = wset.begin(); it!=wset.end(); ++it) {
      REQUIRE(real(*it->overlap()) == Approx(real(ovlp_sum)));
      REQUIRE(imag(*it->overlap()) == Approx(imag(ovlp_sum)));
    }
    auto nCV = wfn.local_number_of_cholesky_vectors();
    boost::multi::array<ComplexType,1> vMF(iextensions<1u>{nCV});
    wfn.vMF(vMF);
    wfn.Energy(wset);
    computeVariationalEnergy(wfn, occs, ham, NAEA, NAEB);
    for(int i=0; i < vMF.size(); i++) {
      std::cout << vMF[i] << std::endl;
    }
    std::cout << wset[0].energy() << std::endl;
  }
}

TEST_CASE("test_read_phmsd", "[test_read_phmsd]")
{
  OHMMS::Controller->initialize(0, NULL);
  auto world = boost::mpi3::environment::get_world_instance();
  if(not world.root()) infoLog.pause();

#ifdef ENABLE_CUDA
  auto node = world.split_shared(world.rank());
  qmc_cuda::CUDA_INIT(node);
  using Alloc = qmc_cuda::cuda_gpu_allocator<ComplexType>;
#else
  using Alloc = shared_allocator<ComplexType>;
#endif

  test_read_phmsd<Alloc>(world);
}

TEST_CASE("test_phmsd", "[read_phmsd]")
{
  OHMMS::Controller->initialize(0, NULL);
  auto world = boost::mpi3::environment::get_world_instance();
  if(not world.root()) infoLog.pause();

#ifdef ENABLE_CUDA
  auto node = world.split_shared(world.rank());
  qmc_cuda::CUDA_INIT(node);
  using Alloc = qmc_cuda::cuda_gpu_allocator<ComplexType>;
#else
  using Alloc = shared_allocator<ComplexType>;
#endif

  //test_phmsd<Alloc,SlaterDetOperations_serial<Alloc>>(world);
  test_phmsd<Alloc>(world);
}

}
