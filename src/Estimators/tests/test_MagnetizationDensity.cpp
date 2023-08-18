//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2022 QMCPACK developers
//
// File developed by: Raymond Clay, rclay@sandia.gov, Sandia National Laboratories
//                    Cody A. Melton, cmelton@sandia.gov, Sandia National Laboratories
//
// File created by: Raymond Clay, rclay@sandia.gov, Sandia National Laboratories
//////////////////////////////////////////////////////////////////////////////////////


#include "catch.hpp"
#include "Estimators/MagnetizationDensityInput.h"
#include "ValidMagnetizationDensityInput.h"
#include "Configuration.h"
//#include "QMCHamiltonians/MagDensityEstimator.h"
//for wavefunction
#include "OhmmsData/Libxml2Doc.h"
#include "QMCWaveFunctions/WaveFunctionComponent.h"
#include "QMCWaveFunctions/TrialWaveFunction.h"
#include "QMCWaveFunctions/Fermion/DiracDeterminant.h"
#include "QMCWaveFunctions/Fermion/SlaterDet.h"
#include "QMCWaveFunctions/SpinorSet.h"
#include "Particle/Walker.h"
#include "QMCHamiltonians/OperatorBase.h"
//for nonlocal moves
#include "QMCWaveFunctions/tests/ConstantSPOSet.h"
//for Hamiltonian manipulations.
#include "Particle/ParticleSet.h"
#include "LongRange/EwaldHandler3D.h"
#include "QMCWaveFunctions/WaveFunctionTypes.hpp"

#include "Estimators/MagnetizationDensity.h"
#include "EstimatorTesting.h"

#include "Utilities/RuntimeOptions.h"

using std::string;


namespace qmcplusplus
{
namespace testing
{
class MagnetizationDensityTests
{
  using WF       = WaveFunctionTypes<QMCTraits::ValueType, QMCTraits::FullPrecValueType>;
  using Position = QMCTraits::PosType;
  using Data     = MagnetizationDensity::Data;
  using Real     = WF::Real;
  using Value    = WF::Value;

public:
  void testCopyConstructor(const MagnetizationDensity& magdens)
  {
    MagnetizationDensity magdens2(magdens);
    CHECK(magdens.nsamples_ == magdens2.nsamples_);
    CHECK(magdens.integrator_ == magdens2.integrator_);
    CHECK(magdens.data_locality_ == magdens2.data_locality_);
    for (int idim = 0; idim < OHMMS_DIM; idim++)
    {
      CHECK(magdens.gdims_[idim] == magdens2.gdims_[idim]);
      CHECK(magdens.grid_[idim] == magdens2.grid_[idim]);
      CHECK(magdens.rcorner_[idim] == Approx(magdens2.rcorner_[idim]));
      CHECK(magdens.center_[idim] == Approx(magdens2.center_[idim]));
    }
  }

  void testData(const MagnetizationDensity& magdens, const Data& data)
  {
    for (size_t i = 0; i < data.size(); i++)
      CHECK(magdens.data_[i] == Approx(data[i]));
  }

  void testIntegrationFunctions(const MagnetizationDensity& magdens)
  {
    //For xgrid, we use dx=1.0 from 0 to 8.
    //For ygrid, we use f(x)=x^3.
    std::vector<Value> ygrid = {0, 1, 8, 27, 64, 125, 216, 343, 512};
    Value result_simpsons(0.0);
    result_simpsons = magdens.integrateBySimpsonsRule(ygrid, Real(1.0));

    CHECK(std::real(result_simpsons) == Approx(Real(1024)));
  }

  void testGrids(const MagnetizationDensity& magdens)
  {
    int npoints = 5;
    std::vector<Real> xgrid(npoints);
    Real start = 0.0;
    Real stop  = 1.5;

    Real delta = (stop - start) / Real(npoints - 1);

    magdens.generateUniformGrid(xgrid, start, stop);
    for (int i = 0; i < npoints; i++)
      CHECK(xgrid[i] == Approx(start + i * delta));

    FakeRandom rng;
    magdens.generateRandomGrid(xgrid, rng, start, stop);

    for (int i = 0; i < npoints; i++)
    {
      bool ok = (xgrid[i] >= start) && (xgrid[i] <= stop);
      CHECK(ok);
    }
  }
  int computeBinAccessor(const MagnetizationDensity& magdens, const Position& r, const int spin_index)
  {
    return magdens.computeBin(r, spin_index);
  }
};
} //namespace testing


#ifdef QMC_COMPLEX
TEST_CASE("MagnetizationDensity::MagnetizationDensity(SPInput, Lattice, SpeciesSet)", "[estimators]")
{
  using namespace testing;
  Libxml2Document doc;
  auto input_xml = magdensity::valid_mag_density_input_sections[magdensity::Inputs::valid_magdensity_input];
  bool okay      = doc.parseFromString(input_xml);
  REQUIRE(okay);
  xmlNodePtr node = doc.getRoot();
  MagnetizationDensityInput mdi(node);

  auto lattice = testing::makeTestLattice();

  MagnetizationDensity magdens(std::move(mdi), lattice);
  MagnetizationDensityTests magdenstest;
  magdenstest.testCopyConstructor(magdens);
}

TEST_CASE("MagnetizationDensity::spawnCrowdClone()", "[estimators]")
{
  using namespace testing;
  Libxml2Document doc;
  auto input_xml = magdensity::valid_mag_density_input_sections[magdensity::Inputs::valid_magdensity_input];
  bool okay      = doc.parseFromString(input_xml);
  REQUIRE(okay);
  xmlNodePtr node = doc.getRoot();
  MagnetizationDensityInput mdi(node);

  auto lattice = testing::makeTestLattice();

  MagnetizationDensity original(std::move(mdi), lattice);
  auto clone = original.spawnCrowdClone();
  REQUIRE(clone != nullptr);
  REQUIRE(clone.get() != &original);
  REQUIRE(dynamic_cast<decltype(&original)>(clone.get()) != nullptr);
}
TEST_CASE("MagnetizationDensity::integrals", "[estimators]")
{
  using namespace testing;
  Libxml2Document doc;
  auto input_xml = magdensity::valid_mag_density_input_sections[magdensity::Inputs::valid_magdensity_input];
  bool okay      = doc.parseFromString(input_xml);
  REQUIRE(okay);
  xmlNodePtr node = doc.getRoot();
  MagnetizationDensityInput mdi(node);

  auto lattice = testing::makeTestLattice();

  MagnetizationDensity magdens(std::move(mdi), lattice);

  MagnetizationDensityTests magdenstest;
  magdenstest.testIntegrationFunctions(magdens);
  magdenstest.testGrids(magdens);
}

TEST_CASE("MagnetizationDensity::gridAssignment", "[estimators]")
{
  using namespace testing;
  int nbintest                   = 8;
  ParticleSet::ParticlePos Rtest = {{0.92832101, 0, 1.77067004},           //bin 0
                                    {-2.32013211, 0, 5.31201011},          //bin 1
                                    {3.48086859, 1.61996772, 1.77067004},  //bin 2
                                    {0.23241546, 1.61996772, 5.31201011},  //bin 3
                                    {3.48086859, -1.61996772, 1.77067004}, //bin 4
                                    {0.23241546, -1.61996772, 5.31201011}, //bin 5
                                    {6.03341616, 0, 1.77067004},           //bin 6
                                    {2.78496304, 0, 5.31201011}};          //bin 7

  ParticleSet::ParticleLayout lattice;
  lattice.R(0, 0) = 5.10509515;
  lattice.R(0, 1) = -3.23993545;
  lattice.R(0, 2) = 0.00000000;
  lattice.R(1, 0) = 5.10509515;
  lattice.R(1, 1) = 3.23993545;
  lattice.R(1, 2) = 0.00000000;
  lattice.R(2, 0) = -6.49690625;
  lattice.R(2, 1) = 0.00000000;
  lattice.R(2, 2) = 7.08268015;

  lattice.BoxBConds     = true; //periodic
  lattice.LR_dim_cutoff = 15;
  lattice.reset();
  //Shamelessly stealing this from test_einset.cpp.

  //Now to construct the input.  See ValidMagnetizationDensity.h, item 4 for the actual input.
  auto mag_input_xml = testing::magdensity::valid_mag_density_input_sections
      [testing::magdensity::Inputs::valid_magdensity_input_unittest];
  Libxml2Document doc;
  bool okay = doc.parseFromString(mag_input_xml);
  REQUIRE(okay);
  xmlNodePtr node = doc.getRoot();

  MagnetizationDensityInput maginput(node);

  maginput.calculateDerivedParameters(lattice);
  MagnetizationDensity magdensity(std::move(maginput), lattice);

  MagnetizationDensityTests magdenstest;
  for (int ibin = 0; ibin < nbintest; ibin++)
  {
    int mag_bin_x = magdenstest.computeBinAccessor(magdensity, Rtest[ibin], 0);
    int mag_bin_y = magdenstest.computeBinAccessor(magdensity, Rtest[ibin], 1);
    int mag_bin_z = magdenstest.computeBinAccessor(magdensity, Rtest[ibin], 2);

    //Here we go from the flattened spatial grid layout to the grid+spin layout
    //that we'll dump in data_.
    int test_bin_x = OHMMS_DIM * ibin + 0;
    int test_bin_y = OHMMS_DIM * ibin + 1;
    int test_bin_z = OHMMS_DIM * ibin + 2;

    CHECK(mag_bin_x == test_bin_x);
    CHECK(mag_bin_y == test_bin_y);
    CHECK(mag_bin_z == test_bin_z);
  }
}

TEST_CASE("MagnetizationDensity::integralAPI", "[estimators]")
{
  using namespace testing;
  using WF    = WaveFunctionTypes<QMCTraits::ValueType, QMCTraits::FullPrecValueType>;
  using Real  = WF::Real;
  using Value = WF::Value;

  Libxml2Document doc_simpsons;
  Libxml2Document doc_mc;
  auto input_xml_simpsons = magdensity::valid_mag_density_input_sections[magdensity::Inputs::valid_magdensity_input];
  auto input_xml_mc       = magdensity::valid_mag_density_input_sections[magdensity::Inputs::valid_magdensity_input_dr];
  bool okay_simpsons      = doc_simpsons.parseFromString(input_xml_simpsons);
  bool okay_mc            = doc_mc.parseFromString(input_xml_mc);
  REQUIRE(okay_simpsons);
  REQUIRE(okay_mc);
  xmlNodePtr node_simpsons = doc_simpsons.getRoot();
  xmlNodePtr node_mc       = doc_mc.getRoot();
  MagnetizationDensityInput mdi_simpsons(node_simpsons);
  MagnetizationDensityInput mdi_mc(node_mc);

  auto lattice = testing::makeTestLattice();

  MagnetizationDensity magdens_simpsons(std::move(mdi_simpsons), lattice);
  MagnetizationDensity magdens_mc(std::move(mdi_mc), lattice);


  std::vector<Value> ygrid = {0, 1, 8, 27, 64, 125, 216, 343, 512};
  Value result_simpsons(0.0);
  Value result_mc(0.0);


  result_simpsons = magdens_simpsons.integrateMagnetizationDensity(ygrid);
  result_mc       = magdens_mc.integrateMagnetizationDensity(ygrid);

  CHECK(std::real(result_simpsons) ==
        Approx(Real(16.2539682539)));                  //From scipy.integrate.  Note, this input file has 64 samples.
                                                       //Since I only use 9 entries, this integral is internally treated
                                                       //as from [0 to 2pi*8/63]
  CHECK(std::real(result_mc) == Approx(Real(10.125))); //Divide sum(ygrid) by nsamples=128
}

TEST_CASE("MagnetizationDensity::IntegrationTest", "[estimators]")
{
  app_log() << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n";
  app_log() << "!!!!   Evaluate MagDensity   !!!!\n";
  app_log() << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n";

  //For now, do a small square case.
  const int nelec       = 2;
  const int norb        = 2;
  using WF              = WaveFunctionTypes<QMCTraits::ValueType, QMCTraits::FullPrecValueType>;
  using Real            = WF::Real;
  using Value           = WF::Value;
  using Grad            = WF::Grad;
  using ValueVector     = Vector<Value>;
  using GradVector      = Vector<Grad>;
  using ValueMatrix     = Matrix<Value>;
  using PropertySetType = OperatorBase::PropertySetType;
  using MCPWalker       = Walker<QMCTraits, PtclOnLatticeTraits>;
  using Data            = MagnetizationDensity::Data;
  using GradMatrix      = Matrix<Grad>;
  using namespace testing;

  // O2 test example from pwscf non-collinear calculation.
  ParticleSet::ParticleLayout lattice;
  lattice.R(0, 0) = 5.10509515;
  lattice.R(0, 1) = -3.23993545;
  lattice.R(0, 2) = 0.00000000;
  lattice.R(1, 0) = 5.10509515;
  lattice.R(1, 1) = 3.23993545;
  lattice.R(1, 2) = 0.00000000;
  lattice.R(2, 0) = -6.49690625;
  lattice.R(2, 1) = 0.00000000;
  lattice.R(2, 2) = 7.08268015;

  lattice.BoxBConds     = true; //periodic
  lattice.LR_dim_cutoff = 15;
  lattice.reset();
  //Shamelessly stealing this from test_einset.cpp.

  //Now to construct the input.  See ValidMagnetizationDensity.h, item 4 for the actual input.
  auto mag_input_xml = testing::magdensity::valid_mag_density_input_sections
      [testing::magdensity::Inputs::valid_magdensity_input_unittest];
  Libxml2Document doc;
  bool okay = doc.parseFromString(mag_input_xml);
  REQUIRE(okay);
  xmlNodePtr node = doc.getRoot();

  MagnetizationDensityInput maginput(node);

  maginput.calculateDerivedParameters(lattice);

  const SimulationCell simulation_cell(lattice);
  ParticleSet elec_(simulation_cell);

  elec_.setName("elec");
  elec_.create({2});

  elec_.R[0]     = {5, 0, 0};
  elec_.R[1]     = {2.22798, 0, 4.249609};
  elec_.spins[0] = 1.9;
  elec_.spins[1] = 2.5410;
  elec_.setSpinor(true);

  SpeciesSet& tspecies       = elec_.getSpeciesSet();
  int upIdx                  = tspecies.addSpecies("u");
  int chargeIdx              = tspecies.addAttribute("charge");
  tspecies(chargeIdx, upIdx) = -1;

  elec_.resetGroups();
  elec_.update();
  // </steal>


  //The values for the spinor SPO (up and down channels) were taken from a python script.
  //Since we are considering only a frozen electron configuration, we lock the spatial parts
  //of the spinor down to test the spin integration.
  //
  //The spin state of these two orbitals corresponds to theta=45, phi=45, so
  //[cos(theta/2), e^(i*phi)*sin(theta/2)].  This yields an expectation value of
  //[0.5,0.5,1/sqrt(2)] per electron when done analytically.
  ValueVector uprow0{Value(0.92387953, 0), Value(0.92387953, 0.)};
  ValueVector dnrow0{Value(0.27059805, 0.27059805), Value(0.27059805, 0.27059805)};
  ValueVector uprow1{Value(0.29131988, 0.87674747), Value(0.81078057, 0.44293144)};
  ValueVector dnrow1{Value(-0.17146777, 0.342119), Value(0.10774051, 0.36720375)};

  ValueMatrix mup, mdn;
  mup.resize(nelec, norb);
  mdn.resize(nelec, norb);

  for (int iorb = 0; iorb < norb; iorb++)
  {
    mup(0, iorb) = uprow0[iorb];
    mdn(0, iorb) = dnrow0[iorb];
    mup(1, iorb) = uprow1[iorb];
    mdn(1, iorb) = dnrow1[iorb];
  }
  auto spo_up = std::make_unique<ConstantSPOSet>("ConstantUpSet", nelec, norb);
  auto spo_dn = std::make_unique<ConstantSPOSet>("ConstantDnSet", nelec, norb);

  spo_up->setRefVals(mup);
  spo_dn->setRefVals(mdn);
  auto spinor_set = std::make_unique<SpinorSet>("ConstSpinorSet");
  spinor_set->set_spos(std::move(spo_up), std::move(spo_dn));

  auto dd = std::make_unique<DiracDeterminant<>>(std::move(spinor_set), 0, nelec);

  std::vector<std::unique_ptr<DiracDeterminantBase>> dirac_dets;
  dirac_dets.push_back(std::move(dd));
  auto sd = std::make_unique<SlaterDet>(elec_, std::move(dirac_dets));

  RuntimeOptions runtime_options;
  TrialWaveFunction psi(runtime_options);
  psi.addComponent(std::move(sd));

  //psi.addComponent(std::move(dd));
  psi.evaluateLog(elec_); //make sure all intermediates are calculated and initialized.

  MagnetizationDensity magdensity(std::move(maginput), lattice);

  //Now to create a fake walker with property set.
  PropertySetType Observables;
  MCPWalker my_walker(nelec);
  my_walker.Weight = 1.0;

  elec_.saveWalker(my_walker);

  //Now to create the crowd related quantities:
  int nwalkers = 2;
  std::vector<MCPWalker> walkers(nwalkers);
  std::vector<ParticleSet> psets{elec_, elec_};

  psets[1].R[0][0] = 0;
  psets[1].R[0][1] = 0.5;
  psets[1].R[0][2] = 0;
  psets[1].R[1][0] = 1;
  psets[1].R[1][1] = 3;
  psets[1].R[1][2] = 1;

  std::vector<UPtr<TrialWaveFunction>> twfcs(2);
  for (int iw = 0; iw < nwalkers; iw++)
    twfcs[iw] = psi.makeClone(psets[iw]);

  auto updateWalker = [](auto& walker, auto& pset_target, auto& trial_wavefunction) {
    pset_target.update(true);
    pset_target.donePbyP();
    trial_wavefunction.evaluateLog(pset_target);
    pset_target.saveWalker(walker);
  };
  for (int iw = 0; iw < nwalkers; iw++)
    updateWalker(walkers[iw], psets[iw], *(twfcs[iw]));

  auto ref_walkers(makeRefVector<MCPWalker>(walkers));
  auto ref_psets(makeRefVector<ParticleSet>(psets));
  auto ref_twfcs(convertUPtrToRefVector(twfcs));

  FakeRandom rng;
  magdensity.accumulate(ref_walkers, ref_psets, ref_twfcs, rng);

  //Now the reference data
  //
  //Note.  These reference values are derived from an independently written python code.
  //At this spin configuration, the observable should return these values.  When the python code
  //performs the integration over all spin variables, only then will one agree with analytic results.
  //More details found in spinor_45deg_magdensity_test.py
  int datsize = 24;
  Data ref_data(datsize, 0);

  //Spinors yield the same value when integrated, but the choice of bin is different between walker 1 and 2.
  ref_data[0] = -0.97448154;
  ref_data[1] = -0.37462387;
  ref_data[2] = 2.36817514;

  ref_data[12] = -0.97448154;
  ref_data[13] = -0.37462387;
  ref_data[14] = 2.36817514;

  //Spinors yield the same value when integrated, but the choice of bin is different.
  ref_data[18] = 1.20557377;
  ref_data[19] = -0.60536469;
  ref_data[20] = 0.98980165;

  ref_data[21] = 1.20557377;
  ref_data[22] = -0.60536469;
  ref_data[23] = 0.98980165;


  MagnetizationDensityTests magdenstest;
  magdenstest.testData(magdensity, ref_data);
}
#endif
} // namespace qmcplusplus
