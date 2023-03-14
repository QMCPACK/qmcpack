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

using std::string;

namespace qmcplusplus
{
#ifdef QMC_COMPLEX
TEST_CASE("MagDensity", "[hamiltonian]")
{
  app_log() << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n";
  app_log() << "!!!!   Evaluate MagDensity   !!!!\n";
  app_log() << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n";

  //For now, do a small square case.
  const int nelec   = 2;
  const int norb    = 2;
  using WF          = WaveFunctionTypes<QMCTraits::ValueType, QMCTraits::FullPrecValueType>;
  using Real        = WF::Real;
  using Value       = WF::Value;
  using Grad        = WF::Grad;
  using ValueVector = Vector<Value>;
  using GradVector  = Vector<Grad>;
  using ValueMatrix = Matrix<Value>;
  using GradMatrix  = Matrix<Grad>;

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
  auto mag_input_xml = testing::magdensity::valid_mag_density_input_sections[testing::magdensity::Inputs::valid_magdensity_input_unittest];
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

  elec_.R[0][0] = 5;
  elec_.R[0][1] = 0;
  elec_.R[0][2] = 0;
  elec_.R[1][0] = 2.22798;
  elec_.R[1][1] = 0;
  elec_.R[1][2] = 4.249609;

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
  auto spo_up = std::make_unique<ConstantSPOSet>("ConstantUpSet",nelec,norb);
  auto spo_dn = std::make_unique<ConstantSPOSet>("ConstantDnSet",nelec,norb);

  spo_up->setRefVals(mup);
  spo_dn->setRefVals(mdn);
  auto spinor_set = std::make_unique<SpinorSet>("ConstSpinorSet");
  spinor_set->set_spos(std::move(spo_up), std::move(spo_dn));

  auto dd = std::make_unique<DiracDeterminant<>>(std::move(spinor_set), 0, nelec);

  TrialWaveFunction psi;
  psi.addComponent(std::move(dd));
  psi.evaluateLog(elec_); //make sure all intermediates are calculated and initialized.
   
  MagnetizationDensity magdensity(std::move(maginput),lattice);
  /*
  MagDensityEstimator magdensity(elec_, psi);

  const char* magtxt = "<tmp> \
  <estimator name=\"MagDensity\" type=\"magdensity\" delta=\"0.5 0.5 0.5\" spin_integral=\"simpson\" nsamples=\"9\"/> \
  </tmp> \
  ";

  Libxml2Document doc;
  bool okay = doc.parseFromString(magtxt);
  REQUIRE(okay);

  xmlNodePtr root = doc.getRoot();

  xmlNodePtr magxml = xmlFirstElementChild(root);
  magdensity.put(magxml);
  //Now to create a fake walker with property set.
  PropertySetType Observables;
  MCPWalker my_walker(nelec);
  my_walker.Weight = 1.0;

  elec_.saveWalker(my_walker);
  //register the observable and evaluate.
  magdensity.setHistories(my_walker);
  magdensity.addObservables(Observables, elec_.Collectables);
  magdensity.evaluate(elec_);
  app_log() << "Collectables = \n";
  for (int i = 0; i < elec_.Collectables.size(); i++)
    app_log() << i << " " << elec_.Collectables[i] << std::endl;
  */
  //Note.  These reference values are derived from an independently written python code.
  //At this spin configuration, the observable should return these values.  When the python code
  //performs the integration over all spin variables, only then will one agree with analytic results.
  /*
  //Spin of first grid point.
  REQUIRE(elec_.Collectables[0] == Approx(-0.97448154));
  REQUIRE(elec_.Collectables[1] == Approx(-0.37462387));
  REQUIRE(elec_.Collectables[2] == Approx(2.36817514));
  //Spin of last grid point.
  REQUIRE(elec_.Collectables[21] == Approx(1.20557377));
  REQUIRE(elec_.Collectables[22] == Approx(-0.60536469));
  REQUIRE(elec_.Collectables[23] == Approx(0.98980165));

  //All other grid points should have zero spin, because there are no electrons there.
  for (int i = 3; i < 21; i++)
  {
    REQUIRE(elec_.Collectables[i] == Approx(0.0));
  }

  //This concludes the correctness test of the estimator.  Now to test cloning:
  elec_.Collectables[0] = 0;
  elec_.Collectables[1] = 0;
  elec_.Collectables[2] = 0;

  elec_.Collectables[21] = 0;
  elec_.Collectables[22] = 0;
  elec_.Collectables[23] = 0;

  std::unique_ptr<OperatorBase> magclone = magdensity.makeClone(elec_, psi);
  magdensity.evaluate(elec_);

  //Spin of first grid point.
  REQUIRE(elec_.Collectables[0] == Approx(-0.97448154));
  REQUIRE(elec_.Collectables[1] == Approx(-0.37462387));
  REQUIRE(elec_.Collectables[2] == Approx(2.36817514));
  //Spin of last grid point.
  REQUIRE(elec_.Collectables[21] == Approx(1.20557377));
  REQUIRE(elec_.Collectables[22] == Approx(-0.60536469));
  REQUIRE(elec_.Collectables[23] == Approx(0.98980165));

  //All other grid points should have zero spin, because there are no electrons there.
  for (int i = 3; i < 21; i++)
  {
    REQUIRE(elec_.Collectables[i] == Approx(0.0));
  }
  */
}
#endif
} // namespace qmcplusplus
