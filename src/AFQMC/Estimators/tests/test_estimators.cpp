//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Miguel Morales, moralessilva2@llnl.gov, Lawrence Livermore National Laboratory
//
// File created by: Fionn Malone, malone14@llnl.gov, Lawrence Livermore National Laboratory
//////////////////////////////////////////////////////////////////////////////////////

#include "Message/catch_mpi_main.hpp"

#include "Configuration.h"

#include "OhmmsData/Libxml2Doc.h"
#include "OhmmsApp/ProjectData.h"
#include "io/hdf_archive.h"

#undef APP_ABORT
#define APP_ABORT(x) {std::cout << x <<std::endl; exit(0);}

#include <stdio.h>
#include <string>

#include "AFQMC/config.h"
#include "AFQMC/Matrix/tests/matrix_helpers.h"
#include "AFQMC/Hamiltonians/HamiltonianFactory.h"
#include "AFQMC/Hamiltonians/Hamiltonian.hpp"
#include "AFQMC/Wavefunctions/WavefunctionFactory.h"
#include "AFQMC/Wavefunctions/Wavefunction.hpp"
#include "AFQMC/Walkers/WalkerSet.hpp"
#include "AFQMC/Estimators/EstimatorBase.h"
#include "AFQMC/Estimators/BackPropagatedEstimator.h"
#include "AFQMC/Utilities/test_utils.hpp"

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

TEST_CASE("back_propagated_rdm", "[estimators]")
{
  OHMMS::Controller->initialize(0, NULL);
  auto world = boost::mpi3::environment::get_world_instance();

  if(not file_exists("./wfn.dat") ) {
    app_log()<<" Skipping back_propagated_rdm. ./wfn.dat files not found. \n";
  } else {

    // Global Task Group
    afqmc::GlobalTaskGroup gTG(world);

    auto file_data = read_test_results_from_hdf<ValueType>("./afqmc.h5");
    int NMO = file_data.NMO;
    int NAEA = file_data.NAEA;
    int NAEB = file_data.NAEB;

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

const char *wfn_xml_block =
"<Wavefunction name=\"wfn0\" type=\"NOMSD\" info=\"info0\"> \
  <parameter name=\"filetype\">ascii</parameter> \
  <parameter name=\"filename\">./wfn.dat</parameter> \
  <parameter name=\"cutoff\">1e-6</parameter> \
</Wavefunction> \
";
    WALKER_TYPES walker_type = COLLINEAR;
    Libxml2Document doc2;
    okay = doc2.parseFromString(wfn_xml_block);
    REQUIRE(okay);
    std::string wfn_name("wfn0");
    WavefunctionFactory WfnFac(InfoMap);
    WfnFac.push(wfn_name,doc2.getRoot());
    Wavefunction& wfn = WfnFac.getWavefunction(TG,TG,wfn_name,walker_type,&ham,1e-6,nwalk);

const char *wlk_xml_block =
"<WalkerSet name=\"wset0\">  \
  <parameter name=\"walker_type\">collinear</parameter>  \
  <parameter name=\"back_propagation_steps\">10</parameter>  \
</WalkerSet> \
";
    Libxml2Document doc3;
    okay = doc3.parseFromString(wlk_xml_block);
    REQUIRE(okay);


    WalkerSet wset(TG,doc3.getRoot(),InfoMap["info0"],&rng);
    auto initial_guess = WfnFac.getInitialGuess(wfn_name);
    REQUIRE(initial_guess.shape()[0]==2);
    REQUIRE(initial_guess.shape()[1]==NMO);
    REQUIRE(initial_guess.shape()[2]==NAEA);
    wset.resize(nwalk,initial_guess[0],initial_guess[0]);
    using EstimPtr = std::shared_ptr<EstimatorBase>;
    std::vector<EstimPtr> estimators;
    const char *est_xml_block =
"<Estimator name=\"back_propagation\"> \
      <parameter name=\"block_size\">2</parameter> \
  </Estimator> \
";
    Libxml2Document doc4;
    okay = doc4.parseFromString(est_xml_block);
    REQUIRE(okay);
    bool impsamp = true;
    Wavefunction* wfn0 = &wfn;
    estimators.emplace_back(static_cast<EstimPtr>(std::make_shared<BackPropagatedEstimator>(TG,InfoMap["info0"],"none",doc4.getRoot(),
                                                                                            walker_type,*wfn0,impsamp)));
    hdf_archive dump;
    std::ofstream out;
    dump.create("estimates.h5");
    dump.open("estimates.h5");
    for (int iblock = 0; iblock < 10; iblock++) {
      estimators[0]->accumulate_block(wset);
      estimators[0]->print(out,dump,wset);
    }
    dump.close();
    std::vector<ComplexType> read_data(2*NMO*NMO);
    hdf_archive reader;
    // Read from a particular block.
    if(!reader.open("estimates.h5",H5F_ACC_RDONLY)) {
      app_error()<<" Error opening estimates.h5. \n";
      APP_ABORT("");
    }
    reader.read(read_data, "BackPropagated/one_rdm_4");
    reader.close();
    boost::multi::array_ref<ComplexType,3> BPRDM(read_data.data(), {2,NMO,NMO});
    ComplexType trace = ComplexType(0.0);
    for(int i = 0; i < NMO; i++) trace += BPRDM[0][i][i] + BPRDM[1][i][i];
    REQUIRE(trace.real()==(NAEA+NAEB));
    // Test the RDM. Since no back propagation has been performed the RDM should be
    // identical to the mixed estimate.
    boost::multi::array<ComplexType,3> OrbMat;
    readWfn(std::string("./wfn.dat"),OrbMat,NMO,NAEA,NAEB);
    SlaterDetOperations<ComplexType> SDet(NMO,NAEA);
    boost::multi::array<ComplexType,2> G;
    G.reextent({NMO,NMO});
    auto Ovlp = SDet.MixedDensityMatrix_noHerm(OrbMat[0],OrbMat[0],G);
    verify_approx(G, BPRDM[0]);
    Ovlp = SDet.MixedDensityMatrix_noHerm(OrbMat[1],OrbMat[1],G);
    verify_approx(G, BPRDM[1]);

  }
}

}
