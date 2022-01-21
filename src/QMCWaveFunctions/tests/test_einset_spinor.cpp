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

#include "OhmmsData/Libxml2Doc.h"
#include "OhmmsPETE/OhmmsMatrix.h"
#include "Particle/ParticleSet.h"
#include "Particle/ParticleSetPool.h"
#include "QMCWaveFunctions/WaveFunctionComponent.h"
#include "QMCWaveFunctions/SPOSetBuilderFactory.h"
#include "Utilities/ResourceCollection.h"

#include <stdio.h>
#include <string>
#include <limits>

using std::string;

namespace qmcplusplus
{
//Now we test the spinor set with Einspline orbitals from HDF.
#ifdef QMC_COMPLEX
TEST_CASE("Einspline SpinorSet from HDF", "[wavefunction]")
{
  app_log() << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n";
  app_log() << "!!!!!  Einspline SpinorSet from HDF   !!!!!\n";
  app_log() << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n";

  using ValueType = SPOSet::ValueType;
  using RealType  = SPOSet::RealType;
  Communicate* c = OHMMS::Controller;

  ParticleSet::ParticleLayout_t lattice;
  // O2 test example from pwscf non-collinear calculation.
  lattice.R(0, 0) = 5.10509515;
  lattice.R(0, 1) = -3.23993545;
  lattice.R(0, 2) = 0.00000000;
  lattice.R(1, 0) = 5.10509515;
  lattice.R(1, 1) = 3.23993545;
  lattice.R(1, 2) = 0.00000000;
  lattice.R(2, 0) = -6.49690625;
  lattice.R(2, 1) = 0.00000000;
  lattice.R(2, 2) = 7.08268015;

  ParticleSetPool ptcl = ParticleSetPool(c);
  ptcl.setSimulationCell(lattice);
  auto ions_uptr = std::make_unique<ParticleSet>(ptcl.getSimulationCell());
  auto elec_uptr = std::make_unique<ParticleSet>(ptcl.getSimulationCell());
  ParticleSet& ions_(*ions_uptr);
  ParticleSet& elec_(*elec_uptr);

  ions_.setName("ion");
  ptcl.addParticleSet(std::move(ions_uptr));
  ions_.create(2);

  ions_.R[0][0] = 0.00000000;
  ions_.R[0][1] = 0.00000000;
  ions_.R[0][2] = 1.08659253;
  ions_.R[1][0] = 0.00000000;
  ions_.R[1][1] = 0.00000000;
  ions_.R[1][2] = -1.08659253;

  elec_.setName("elec");
  ptcl.addParticleSet(std::move(elec_uptr));
  elec_.create(3);
  elec_.R[0][0] = 0.1;
  elec_.R[0][1] = -0.3;
  elec_.R[0][2] = 1.0;
  elec_.R[1][0] = -0.1;
  elec_.R[1][1] = 0.3;
  elec_.R[1][2] = 1.0;
  elec_.R[2][0] = 0.1;
  elec_.R[2][1] = 0.2;
  elec_.R[2][2] = 0.3;

  elec_.spins[0] = 0.0;
  elec_.spins[1] = 0.2;
  elec_.spins[2] = 0.4;
  elec_.setSpinor(true);

  SpeciesSet& tspecies       = elec_.getSpeciesSet();
  int upIdx                  = tspecies.addSpecies("u");
  int chargeIdx              = tspecies.addAttribute("charge");
  tspecies(chargeIdx, upIdx) = -1;

  elec_.update();
  ions_.update();


  const char* particles = "<tmp> \
   <sposet_builder name=\"A\" type=\"einspline\" href=\"o2_45deg_spins.pwscf.h5\" tilematrix=\"1 0 0 0 1 0 0 0 1\" twistnum=\"0\" source=\"ion\" size=\"3\" precision=\"float\" meshfactor=\"4.0\"> \
     <sposet name=\"myspo\" size=\"3\"> \
       <occupation mode=\"ground\"/> \
     </sposet> \
   </sposet_builder> \
   </tmp> \
";

  Libxml2Document doc;
  bool okay = doc.parseFromString(particles);
  CHECK(okay);

  xmlNodePtr root = doc.getRoot();

  xmlNodePtr ein1 = xmlFirstElementChild(root);

  SPOSetBuilderFactory fac(c, elec_, ptcl.getPool());
  auto& builder = fac.createSPOSetBuilder(ein1);

  SPOSet* spo = builder.createSPOSet(ein1);
  CHECK(spo);

  SPOSet::ValueMatrix_t psiM(elec_.R.size(), spo->getOrbitalSetSize());
  SPOSet::GradMatrix_t dpsiM(elec_.R.size(), spo->getOrbitalSetSize());
  SPOSet::ValueMatrix_t dspsiM(elec_.R.size(), spo->getOrbitalSetSize()); //spin gradient
  SPOSet::ValueMatrix_t d2psiM(elec_.R.size(), spo->getOrbitalSetSize());

  //These are the reference values computed from a spin-polarized calculation,
  //with the assumption that the coefficients for phi^\uparrow
  SPOSet::ValueMatrix_t psiM_up(elec_.R.size(), spo->getOrbitalSetSize());
  SPOSet::ValueMatrix_t psiM_down(elec_.R.size(), spo->getOrbitalSetSize());
  SPOSet::ValueMatrix_t psiM_ref(elec_.R.size(), spo->getOrbitalSetSize());
  SPOSet::ValueMatrix_t dspsiM_ref(elec_.R.size(), spo->getOrbitalSetSize());

  SPOSet::GradMatrix_t dpsiM_up(elec_.R.size(), spo->getOrbitalSetSize());
  SPOSet::GradMatrix_t dpsiM_down(elec_.R.size(), spo->getOrbitalSetSize());
  SPOSet::GradMatrix_t dpsiM_ref(elec_.R.size(), spo->getOrbitalSetSize());

  SPOSet::ValueMatrix_t d2psiM_up(elec_.R.size(), spo->getOrbitalSetSize());
  SPOSet::ValueMatrix_t d2psiM_down(elec_.R.size(), spo->getOrbitalSetSize());
  SPOSet::ValueMatrix_t d2psiM_ref(elec_.R.size(), spo->getOrbitalSetSize());


  //These reference values were generated as follows:
  // 1.) Non-Collinear O2 calculation in PBC's performed using Quantum Espresso.
  // 2.) Spinor wavefunction converted to HDF5 using convertpw4qmc tool.  Mainly, this places the up channel in the spin_0 slot, and spin down in spin_1.
  // 3.) Up and down values for spinor components are directly evaluated using plane waves to make reference values.

  //reference values, elec 0, orbital 0
  psiM_up[0][0]       = ValueType(3.0546848269585873, 2.6698880339914073);
  dpsiM_up[0][0][0]   = ValueType(-0.1851673255644419, -0.1618419361786101);
  dpsiM_up[0][0][1]   = ValueType(0.60078848166567, 0.5251078699998748);
  dpsiM_up[0][0][2]   = ValueType(-4.882727805715862, -4.267653505749376);
  d2psiM_up[0][0]     = ValueType(-5.949965529898391, -5.200452764159391);
  psiM_down[0][0]     = ValueType(1.26529305723931, 1.105896575232331);
  dpsiM_down[0][0][0] = ValueType(-0.07669912992155427, -0.0670361472429191);
  dpsiM_down[0][0][1] = ValueType(0.24885435007410817, 0.21750534167815624);
  dpsiM_down[0][0][2] = ValueType(-2.0224938024371926, -1.7677088591422807);
  d2psiM_down[0][0]   = ValueType(-2.4645607715652593, -2.1540811306541467);

  //reference values, elec 0, orbital 1
  psiM_up[0][1]       = ValueType(-0.05877385066521865, -1.652862402896465);
  dpsiM_up[0][1][0]   = ValueType(0.003806703582748908, 0.1070679097294273);
  dpsiM_up[0][1][1]   = ValueType(-0.01231680454381807, -0.34637011425999203);
  dpsiM_up[0][1][2]   = ValueType(0.09976715349787516, 2.8057021736885446);
  d2psiM_up[0][1]     = ValueType(0.11243204605091028, 3.1618388472842645);
  psiM_down[0][1]     = ValueType(0.14190658959478397, 3.9903718588991324);
  dpsiM_down[0][1][0] = ValueType(-0.009192765914265412, -0.258482966548934);
  dpsiM_down[0][1][1] = ValueType(0.0297380860827861, 0.8362132640454849);
  dpsiM_down[0][1][2] = ValueType(-0.24088367433889216, -6.773578493991127);
  d2psiM_down[0][1]   = ValueType(-0.27145824651495754, -7.633377309590959);

  //reference values, elec 0, orbital 2
  psiM_up[0][2]       = ValueType(1.7341610599591415, 3.26452640029962);
  dpsiM_up[0][2][0]   = ValueType(0.02410454300409052, 0.04537650582197769);
  dpsiM_up[0][2][1]   = ValueType(-0.0812208525894339, -0.15289674365675554);
  dpsiM_up[0][2][2]   = ValueType(2.056759046129918, 3.871811527443685);
  d2psiM_up[0][2]     = ValueType(-0.3715079628152589, -0.6993565364031098);
  psiM_down[0][2]     = ValueType(0.7183159092489255, 1.3522053467852482);
  dpsiM_down[0][2][0] = ValueType(0.009984891703495969, 0.01879614056538452);
  dpsiM_down[0][2][1] = ValueType(-0.033643334635896874, -0.06333185576262795);
  dpsiM_down[0][2][2] = ValueType(0.8519414150210098, 1.603750121892103);
  d2psiM_down[0][2]   = ValueType(-0.15388532808689237, -0.2896807896155573);


  //reference values, elec 1, orbital 0
  psiM_up[1][0]       = ValueType(3.0526148989188244, 2.6680787492636187);
  dpsiM_up[1][0][0]   = ValueType(0.20449301174627027, 0.17873282170446286);
  dpsiM_up[1][0][1]   = ValueType(-0.6096780888298439, -0.5328779559603193);
  dpsiM_up[1][0][2]   = ValueType(-4.885040183718155, -4.269674660852541);
  d2psiM_up[1][0]     = ValueType(-5.875072106235885, -5.134991421765417);
  psiM_down[1][0]     = ValueType(1.2644358252460446, 1.1051472909800415);
  dpsiM_down[1][0][0] = ValueType(0.0847027693565092, 0.0740326785381825);
  dpsiM_down[1][0][1] = ValueType(-0.25253678670740615, -0.22072360765802454);
  dpsiM_down[1][0][2] = ValueType(-2.023451052801436, -1.768545254958754);
  d2psiM_down[1][0]   = ValueType(-2.433542859102463, -2.126969850545346);

  //reference values, elec 1, orbital 1
  psiM_up[1][1]       = ValueType(-0.05872929134760467, -1.6516107719315123);
  dpsiM_up[1][1][0]   = ValueType(-0.0042225364734192325, -0.11876835593196035);
  dpsiM_up[1][1][1]   = ValueType(0.012491965861615007, 0.35129150754532346);
  dpsiM_up[1][1][2]   = ValueType(0.09980846579193113, 2.806855260627992);
  d2psiM_up[1][1]     = ValueType(0.11086616211845124, 3.1178291585160025);
  psiM_down[1][1]     = ValueType(0.14179908178203693, 3.9873502499791);
  dpsiM_down[1][1][0] = ValueType(0.010197668920898767, 0.2867312658960351);
  dpsiM_down[1][1][1] = ValueType(-0.030160592987572725, -0.8480940968707702);
  dpsiM_down[1][1][2] = ValueType(-0.24098310461934494, -6.776362513721667);
  d2psiM_down[1][1]   = ValueType(-0.26767894399782044, -7.527130337670782);

  //reference values, elec 1, orbital 2
  psiM_up[1][2]       = ValueType(1.7338733833257558, 3.263984881354726);
  dpsiM_up[1][2][0]   = ValueType(-0.02165584086872901, -0.0407670903699481);
  dpsiM_up[1][2][1]   = ValueType(0.08288083949305346, 0.15602164581174188);
  dpsiM_up[1][2][2]   = ValueType(2.0621151061966456, 3.881894235760205);
  d2psiM_up[1][2]     = ValueType(-0.3566890854259599, -0.6714605501817572);
  psiM_down[1][2]     = ValueType(0.7181968101586865, 1.3519810682722548);
  dpsiM_down[1][2][0] = ValueType(-0.00897085696147509, -0.01688677968381685);
  dpsiM_down[1][2][1] = ValueType(0.03433096009876233, 0.06462627360298884);
  dpsiM_down[1][2][2] = ValueType(0.8541600134552085, 1.6079267140278);
  d2psiM_down[1][2]   = ValueType(-0.1477482607270697, -0.2781267037329471);


  //reference values, elec 2, orbital 0
  psiM_up[2][0]       = ValueType(4.774972481925916, 4.173472045741891);
  dpsiM_up[2][0][0]   = ValueType(-0.9457258313862555, -0.8265932416325468);
  dpsiM_up[2][0][1]   = ValueType(-1.899691502097708, -1.6603886891267976);
  dpsiM_up[2][0][2]   = ValueType(0.7771301258291673, 0.6792356362625406);
  d2psiM_up[2][0]     = ValueType(-20.40945585069239, -17.83848971293892);
  psiM_down[2][0]     = ValueType(1.9778599493693798, 1.7286974948556129);
  dpsiM_down[2][0][0] = ValueType(-0.3917329590645693, -0.3423838452062018);
  dpsiM_down[2][0][1] = ValueType(-0.786878588427934, -0.6877513351266712);
  dpsiM_down[2][0][2] = ValueType(0.3218978428488249, 0.281346363058232);
  d2psiM_down[2][0]   = ValueType(-8.45387947982597, -7.388894117402044);

  //reference values, elec 2, orbital 1
  psiM_up[2][1]       = ValueType(-0.095146182382511, -2.6757440636563343);
  dpsiM_up[2][1][0]   = ValueType(0.01912387482485274, 0.53780199541144);
  dpsiM_up[2][1][1]   = ValueType(0.03838799057297392, 1.0795586887258484);
  dpsiM_up[2][1][2]   = ValueType(-0.013683016882420245, -0.38479709783829663);
  d2psiM_up[2][1]     = ValueType(0.41702609987278866, 11.727776988772089);
  psiM_down[2][1]     = ValueType(0.22972652344682393, 6.459831671158625);
  dpsiM_down[2][1][0] = ValueType(-0.046172654628017486, -1.2983723819140731);
  dpsiM_down[2][1][1] = ValueType(-0.09268554947869961, -2.606290867185097);
  dpsiM_down[2][1][2] = ValueType(0.03303644631176311, 0.9289838072933512);
  d2psiM_down[2][1]   = ValueType(-1.006891760427076, -28.313415815931304);

  //reference values, elec 2, orbital 2
  psiM_up[2][2]       = ValueType(0.5573944761518197, 1.0492847452220198);
  dpsiM_up[2][2][0]   = ValueType(0.06369314000215545, 0.11990123738728313);
  dpsiM_up[2][2][1]   = ValueType(0.1265010825081423, 0.23813600324436654);
  dpsiM_up[2][2][2]   = ValueType(1.6373025933118952, 3.082192181081695);
  d2psiM_up[2][2]     = ValueType(-0.8650133588132842, -1.6283710474610622);
  psiM_down[2][2]     = ValueType(0.23088102734804172, 0.43462602449526555);
  dpsiM_down[2][2][0] = ValueType(0.02638262680056288, 0.0496645026299118);
  dpsiM_down[2][2][1] = ValueType(0.052398844163716374, 0.0986386889965079);
  dpsiM_down[2][2][2] = ValueType(0.6781955750070273, 1.276680192776833);
  d2psiM_down[2][2]   = ValueType(-0.3583001666840082, -0.6744889374649511);


  RealType h  = 0.001;
  RealType h2 = 0.1;
  for (unsigned int iat = 0; iat < 3; iat++)
  {
    RealType s = elec_.spins[iat];
    RealType coss(0.0), sins(0.0);

    coss = std::cos(s);
    sins = std::sin(s);

    ValueType eis(coss, sins);
    ValueType emis(coss, -sins);
    ValueType eye(0, 1.0);
    //Using the reference values for the up and down channels invdividually, we build the total reference spinor value
    //consistent with the current spin value of particle iat.
    psiM_ref[iat][0] = eis * psiM_up[iat][0] + emis * psiM_down[iat][0];
    psiM_ref[iat][1] = eis * psiM_up[iat][1] + emis * psiM_down[iat][1];
    psiM_ref[iat][2] = eis * psiM_up[iat][2] + emis * psiM_down[iat][2];

    dspsiM_ref[iat][0] = eye * (eis * psiM_up[iat][0] - emis * psiM_down[iat][0]);
    dspsiM_ref[iat][1] = eye * (eis * psiM_up[iat][1] - emis * psiM_down[iat][1]);
    dspsiM_ref[iat][2] = eye * (eis * psiM_up[iat][2] - emis * psiM_down[iat][2]);

    dpsiM_ref[iat][0] = eis * dpsiM_up[iat][0] + emis * dpsiM_down[iat][0];
    dpsiM_ref[iat][1] = eis * dpsiM_up[iat][1] + emis * dpsiM_down[iat][1];
    dpsiM_ref[iat][2] = eis * dpsiM_up[iat][2] + emis * dpsiM_down[iat][2];

    d2psiM_ref[iat][0] = eis * d2psiM_up[iat][0] + emis * d2psiM_down[iat][0];
    d2psiM_ref[iat][1] = eis * d2psiM_up[iat][1] + emis * d2psiM_down[iat][1];
    d2psiM_ref[iat][2] = eis * d2psiM_up[iat][2] + emis * d2psiM_down[iat][2];
  }

  //OK.  Going to test evaluate_notranspose with spin.
  spo->evaluate_notranspose(elec_, 0, elec_.R.size(), psiM, dpsiM, d2psiM);

  for (unsigned int iat = 0; iat < 3; iat++)
  {
    CHECK(psiM[iat][0] == ComplexApprox(psiM_ref[iat][0]).epsilon(h));
    CHECK(psiM[iat][1] == ComplexApprox(psiM_ref[iat][1]).epsilon(h));
    CHECK(psiM[iat][2] == ComplexApprox(psiM_ref[iat][2]).epsilon(h));

    CHECK(dpsiM[iat][0][0] == ComplexApprox(dpsiM_ref[iat][0][0]).epsilon(h));
    CHECK(dpsiM[iat][0][1] == ComplexApprox(dpsiM_ref[iat][0][1]).epsilon(h));
    CHECK(dpsiM[iat][0][2] == ComplexApprox(dpsiM_ref[iat][0][2]).epsilon(h));

    CHECK(dpsiM[iat][1][0] == ComplexApprox(dpsiM_ref[iat][1][0]).epsilon(h));
    CHECK(dpsiM[iat][1][1] == ComplexApprox(dpsiM_ref[iat][1][1]).epsilon(h));
    CHECK(dpsiM[iat][1][2] == ComplexApprox(dpsiM_ref[iat][1][2]).epsilon(h));

    CHECK(dpsiM[iat][2][0] == ComplexApprox(dpsiM_ref[iat][2][0]).epsilon(h));
    CHECK(dpsiM[iat][2][1] == ComplexApprox(dpsiM_ref[iat][2][1]).epsilon(h));
    CHECK(dpsiM[iat][2][2] == ComplexApprox(dpsiM_ref[iat][2][2]).epsilon(h));


    CHECK(d2psiM[iat][0] == ComplexApprox(d2psiM_ref[iat][0]).epsilon(h2));
    CHECK(d2psiM[iat][1] == ComplexApprox(d2psiM_ref[iat][1]).epsilon(h2));
    CHECK(d2psiM[iat][2] == ComplexApprox(d2psiM_ref[iat][2]).epsilon(h2));
  }

  //Now we're going to test evaluateValue and evaluateVGL:

  int OrbitalSetSize = spo->getOrbitalSetSize();
  //temporary arrays for holding the values of the up and down channels respectively.
  SPOSet::ValueVector_t psi_work;

  //temporary arrays for holding the gradients of the up and down channels respectively.
  SPOSet::GradVector_t dpsi_work;

  //temporary arrays for holding the laplacians of the up and down channels respectively.
  SPOSet::ValueVector_t d2psi_work;

  //temporary arrays for holding the spin gradient
  SPOSet::ValueVector_t dspsi_work;


  psi_work.resize(OrbitalSetSize);

  dpsi_work.resize(OrbitalSetSize);

  d2psi_work.resize(OrbitalSetSize);

  dspsi_work.resize(OrbitalSetSize);

  //We worked hard to generate nice reference data above.  Let's generate a test for evaluateV
  //and evaluateVGL by perturbing the electronic configuration by dR, and then make
  //single particle moves that bring it back to our original R reference values.

  //Our perturbation vector.
  ParticleSet::ParticlePos_t dR;
  dR.resize(3);
  //Values chosen based on divine inspiration.  Not important.
  dR[0][0] = 0.1;
  dR[0][1] = 0.2;
  dR[0][2] = 0.1;
  dR[1][0] = -0.1;
  dR[1][1] = -0.05;
  dR[1][2] = 0.05;
  dR[2][0] = -0.1;
  dR[2][1] = 0.0;
  dR[2][2] = 0.0;

  //The new R of our perturbed particle set.
  ParticleSet::ParticlePos_t Rnew;
  Rnew.resize(3);
  Rnew    = elec_.R + dR;
  elec_.R = Rnew;
  elec_.update();

  //Now we test evaluateValue()
  for (unsigned int iat = 0; iat < 3; iat++)
  {
    psi_work = 0.0;
    elec_.makeMove(iat, -dR[iat], false);
    spo->evaluateValue(elec_, iat, psi_work);

    CHECK(psi_work[0] == ComplexApprox(psiM_ref[iat][0]).epsilon(h));
    CHECK(psi_work[1] == ComplexApprox(psiM_ref[iat][1]).epsilon(h));
    CHECK(psi_work[2] == ComplexApprox(psiM_ref[iat][2]).epsilon(h));
    elec_.rejectMove(iat);
  }

  //Now we test evaluateVGL()
  for (unsigned int iat = 0; iat < 3; iat++)
  {
    psi_work   = 0.0;
    dpsi_work  = 0.0;
    d2psi_work = 0.0;
    dspsi_work = 0.0;

    elec_.makeMove(iat, -dR[iat], false);
    spo->evaluateVGL_spin(elec_, iat, psi_work, dpsi_work, d2psi_work, dspsi_work);

    CHECK(psi_work[0] == ComplexApprox(psiM_ref[iat][0]).epsilon(h));
    CHECK(psi_work[1] == ComplexApprox(psiM_ref[iat][1]).epsilon(h));
    CHECK(psi_work[2] == ComplexApprox(psiM_ref[iat][2]).epsilon(h));

    CHECK(dpsi_work[0][0] == ComplexApprox(dpsiM_ref[iat][0][0]).epsilon(h));
    CHECK(dpsi_work[0][1] == ComplexApprox(dpsiM_ref[iat][0][1]).epsilon(h));
    CHECK(dpsi_work[0][2] == ComplexApprox(dpsiM_ref[iat][0][2]).epsilon(h));

    CHECK(dpsi_work[1][0] == ComplexApprox(dpsiM_ref[iat][1][0]).epsilon(h));
    CHECK(dpsi_work[1][1] == ComplexApprox(dpsiM_ref[iat][1][1]).epsilon(h));
    CHECK(dpsi_work[1][2] == ComplexApprox(dpsiM_ref[iat][1][2]).epsilon(h));

    CHECK(dpsi_work[2][0] == ComplexApprox(dpsiM_ref[iat][2][0]).epsilon(h));
    CHECK(dpsi_work[2][1] == ComplexApprox(dpsiM_ref[iat][2][1]).epsilon(h));
    CHECK(dpsi_work[2][2] == ComplexApprox(dpsiM_ref[iat][2][2]).epsilon(h));

    CHECK(d2psi_work[0] == ComplexApprox(d2psiM_ref[iat][0]).epsilon(h2));
    CHECK(d2psi_work[1] == ComplexApprox(d2psiM_ref[iat][1]).epsilon(h2));
    CHECK(d2psi_work[2] == ComplexApprox(d2psiM_ref[iat][2]).epsilon(h2));

    CHECK(dspsi_work[0] == ComplexApprox(dspsiM_ref[iat][0]).epsilon(h));
    CHECK(dspsi_work[1] == ComplexApprox(dspsiM_ref[iat][1]).epsilon(h));
    CHECK(dspsi_work[2] == ComplexApprox(dspsiM_ref[iat][2]).epsilon(h));

    elec_.rejectMove(iat);
  }

  //Now we test evaluateSpin:

  for (unsigned int iat = 0; iat < 3; iat++)
  {
    psi_work   = 0.0;
    dspsi_work = 0.0;

    elec_.makeMove(iat, -dR[iat], false);
    spo->evaluate_spin(elec_, iat, psi_work, dspsi_work);

    CHECK(psi_work[0] == ComplexApprox(psiM_ref[iat][0]).epsilon(h));
    CHECK(psi_work[1] == ComplexApprox(psiM_ref[iat][1]).epsilon(h));
    CHECK(psi_work[2] == ComplexApprox(psiM_ref[iat][2]).epsilon(h));

    CHECK(dspsi_work[0] == ComplexApprox(dspsiM_ref[iat][0]).epsilon(h));
    CHECK(dspsi_work[1] == ComplexApprox(dspsiM_ref[iat][1]).epsilon(h));
    CHECK(dspsi_work[2] == ComplexApprox(dspsiM_ref[iat][2]).epsilon(h));

    elec_.rejectMove(iat);
  }


  // test batched interface
  // first move elec_ back to original positions for reference
  Rnew    = elec_.R - dR;
  elec_.R = Rnew;
  elec_.update();

  //now create second walker, with permuted particle positions
  ParticleSet elec_2(elec_);
  // permute electrons
  elec_2.R[0]     = elec_.R[1];
  elec_2.R[1]     = elec_.R[2];
  elec_2.R[2]     = elec_.R[0];
  elec_2.spins[0] = elec_.spins[1];
  elec_2.spins[1] = elec_.spins[2];
  elec_2.spins[2] = elec_.spins[0];

  ResourceCollection pset_res("test_pset_res");
  elec_.createResource(pset_res);

  RefVectorWithLeader<ParticleSet> p_list(elec_);
  p_list.push_back(elec_);
  p_list.push_back(elec_2);

  ResourceCollectionTeamLock<ParticleSet> mw_pset_lock(pset_res, p_list);

  //update all walkers
  elec_.mw_update(p_list);

  std::unique_ptr<SPOSet> spo_2(spo->makeClone());
  RefVectorWithLeader<SPOSet> spo_list(*spo);
  spo_list.push_back(*spo);
  spo_list.push_back(*spo_2);

  SPOSet::ValueMatrix_t psiM_2(elec_.R.size(), spo->getOrbitalSetSize());
  SPOSet::GradMatrix_t dpsiM_2(elec_.R.size(), spo->getOrbitalSetSize());
  SPOSet::ValueMatrix_t d2psiM_2(elec_.R.size(), spo->getOrbitalSetSize());

  RefVector<SPOSet::ValueMatrix_t> logdet_list;
  RefVector<SPOSet::GradMatrix_t> dlogdet_list;
  RefVector<SPOSet::ValueMatrix_t> d2logdet_list;

  logdet_list.push_back(psiM);
  logdet_list.push_back(psiM_2);
  dlogdet_list.push_back(dpsiM);
  dlogdet_list.push_back(dpsiM_2);
  d2logdet_list.push_back(d2psiM);
  d2logdet_list.push_back(d2psiM_2);

  spo->mw_evaluate_notranspose(spo_list, p_list, 0, 3, logdet_list, dlogdet_list, d2logdet_list);
  for (unsigned int iat = 0; iat < 3; iat++)
  {
    //walker 0
    CHECK(logdet_list[0].get()[iat][0] == ComplexApprox(psiM_ref[iat][0]).epsilon(h));
    CHECK(logdet_list[0].get()[iat][1] == ComplexApprox(psiM_ref[iat][1]).epsilon(h));
    CHECK(logdet_list[0].get()[iat][2] == ComplexApprox(psiM_ref[iat][2]).epsilon(h));
    CHECK(dlogdet_list[0].get()[iat][0][0] == ComplexApprox(dpsiM_ref[iat][0][0]).epsilon(h));
    CHECK(dlogdet_list[0].get()[iat][0][1] == ComplexApprox(dpsiM_ref[iat][0][1]).epsilon(h));
    CHECK(dlogdet_list[0].get()[iat][0][2] == ComplexApprox(dpsiM_ref[iat][0][2]).epsilon(h));
    CHECK(dlogdet_list[0].get()[iat][1][0] == ComplexApprox(dpsiM_ref[iat][1][0]).epsilon(h));
    CHECK(dlogdet_list[0].get()[iat][1][1] == ComplexApprox(dpsiM_ref[iat][1][1]).epsilon(h));
    CHECK(dlogdet_list[0].get()[iat][1][2] == ComplexApprox(dpsiM_ref[iat][1][2]).epsilon(h));
    CHECK(dlogdet_list[0].get()[iat][2][0] == ComplexApprox(dpsiM_ref[iat][2][0]).epsilon(h));
    CHECK(dlogdet_list[0].get()[iat][2][1] == ComplexApprox(dpsiM_ref[iat][2][1]).epsilon(h));
    CHECK(dlogdet_list[0].get()[iat][2][2] == ComplexApprox(dpsiM_ref[iat][2][2]).epsilon(h));
    CHECK(d2logdet_list[0].get()[iat][0] == ComplexApprox(d2psiM_ref[iat][0]).epsilon(h2));
    CHECK(d2logdet_list[0].get()[iat][1] == ComplexApprox(d2psiM_ref[iat][1]).epsilon(h2));
    CHECK(d2logdet_list[0].get()[iat][2] == ComplexApprox(d2psiM_ref[iat][2]).epsilon(h2));

    //walker 1, permuted from reference
    CHECK(logdet_list[1].get()[iat][0] == ComplexApprox(psiM_ref[(iat + 1) % 3][0]).epsilon(h));
    CHECK(logdet_list[1].get()[iat][1] == ComplexApprox(psiM_ref[(iat + 1) % 3][1]).epsilon(h));
    CHECK(logdet_list[1].get()[iat][2] == ComplexApprox(psiM_ref[(iat + 1) % 3][2]).epsilon(h));
    CHECK(dlogdet_list[1].get()[iat][0][0] == ComplexApprox(dpsiM_ref[(iat + 1) % 3][0][0]).epsilon(h));
    CHECK(dlogdet_list[1].get()[iat][0][1] == ComplexApprox(dpsiM_ref[(iat + 1) % 3][0][1]).epsilon(h));
    CHECK(dlogdet_list[1].get()[iat][0][2] == ComplexApprox(dpsiM_ref[(iat + 1) % 3][0][2]).epsilon(h));
    CHECK(dlogdet_list[1].get()[iat][1][0] == ComplexApprox(dpsiM_ref[(iat + 1) % 3][1][0]).epsilon(h));
    CHECK(dlogdet_list[1].get()[iat][1][1] == ComplexApprox(dpsiM_ref[(iat + 1) % 3][1][1]).epsilon(h));
    CHECK(dlogdet_list[1].get()[iat][1][2] == ComplexApprox(dpsiM_ref[(iat + 1) % 3][1][2]).epsilon(h));
    CHECK(dlogdet_list[1].get()[iat][2][0] == ComplexApprox(dpsiM_ref[(iat + 1) % 3][2][0]).epsilon(h));
    CHECK(dlogdet_list[1].get()[iat][2][1] == ComplexApprox(dpsiM_ref[(iat + 1) % 3][2][1]).epsilon(h));
    CHECK(dlogdet_list[1].get()[iat][2][2] == ComplexApprox(dpsiM_ref[(iat + 1) % 3][2][2]).epsilon(h));
    CHECK(d2logdet_list[1].get()[iat][0] == ComplexApprox(d2psiM_ref[(iat + 1) % 3][0]).epsilon(h2));
    CHECK(d2logdet_list[1].get()[iat][1] == ComplexApprox(d2psiM_ref[(iat + 1) % 3][1]).epsilon(h2));
    CHECK(d2logdet_list[1].get()[iat][2] == ComplexApprox(d2psiM_ref[(iat + 1) % 3][2]).epsilon(h2));
  }

  //first, lets displace all the electrons in each walker.
  for (int iat = 0; iat < 3; iat++)
  {
    std::vector<ParticleSet::SingleParticlePos_t> displs = {dR[iat], dR[iat]};
    elec_.mw_makeMove(p_list, iat, displs);
    std::vector<bool> accept = {true, true};
    elec_.mw_accept_rejectMove(p_list, iat, accept);
  }
  elec_.mw_update(p_list);

  SPOSet::ValueVector_t psi_work_2(OrbitalSetSize);
  SPOSet::GradVector_t dpsi_work_2(OrbitalSetSize);
  SPOSet::ValueVector_t d2psi_work_2(OrbitalSetSize);
  SPOSet::ValueVector_t dspsi_work_2(OrbitalSetSize);

  RefVector<SPOSet::ValueVector_t> psi_v_list   = {psi_work, psi_work_2};
  RefVector<SPOSet::GradVector_t> dpsi_v_list   = {dpsi_work, dpsi_work_2};
  RefVector<SPOSet::ValueVector_t> d2psi_v_list = {d2psi_work, d2psi_work_2};
  RefVector<SPOSet::ValueVector_t> dspsi_v_list = {dspsi_work, dspsi_work_2};
  //check mw_evaluateVGLWithSpin
  for (int iat = 0; iat < 3; iat++)
  {
    //reset values to zero, updates the ref vectors to zero as well
    psi_work     = 0.0;
    dpsi_work    = 0.0;
    d2psi_work   = 0.0;
    dspsi_work   = 0.0;
    psi_work_2   = 0.0;
    dpsi_work_2  = 0.0;
    d2psi_work_2 = 0.0;
    dspsi_work_2 = 0.0;

    std::vector<ParticleSet::SingleParticlePos_t> displs = {-dR[iat], -dR[iat]};
    elec_.mw_makeMove(p_list, iat, displs);
    spo->mw_evaluateVGLWithSpin(spo_list, p_list, iat, psi_v_list, dpsi_v_list, d2psi_v_list, dspsi_v_list);
    //walker 0
    CHECK(psi_v_list[0].get()[0] == ComplexApprox(psiM_ref[iat][0]).epsilon(h));
    CHECK(psi_v_list[0].get()[1] == ComplexApprox(psiM_ref[iat][1]).epsilon(h));
    CHECK(psi_v_list[0].get()[2] == ComplexApprox(psiM_ref[iat][2]).epsilon(h));

    CHECK(dpsi_v_list[0].get()[0][0] == ComplexApprox(dpsiM_ref[iat][0][0]).epsilon(h));
    CHECK(dpsi_v_list[0].get()[0][1] == ComplexApprox(dpsiM_ref[iat][0][1]).epsilon(h));
    CHECK(dpsi_v_list[0].get()[0][2] == ComplexApprox(dpsiM_ref[iat][0][2]).epsilon(h));

    CHECK(dpsi_v_list[0].get()[1][0] == ComplexApprox(dpsiM_ref[iat][1][0]).epsilon(h));
    CHECK(dpsi_v_list[0].get()[1][1] == ComplexApprox(dpsiM_ref[iat][1][1]).epsilon(h));
    CHECK(dpsi_v_list[0].get()[1][2] == ComplexApprox(dpsiM_ref[iat][1][2]).epsilon(h));

    CHECK(dpsi_v_list[0].get()[2][0] == ComplexApprox(dpsiM_ref[iat][2][0]).epsilon(h));
    CHECK(dpsi_v_list[0].get()[2][1] == ComplexApprox(dpsiM_ref[iat][2][1]).epsilon(h));
    CHECK(dpsi_v_list[0].get()[2][2] == ComplexApprox(dpsiM_ref[iat][2][2]).epsilon(h));

    CHECK(d2psi_v_list[0].get()[0] == ComplexApprox(d2psiM_ref[iat][0]).epsilon(h2));
    CHECK(d2psi_v_list[0].get()[1] == ComplexApprox(d2psiM_ref[iat][1]).epsilon(h2));
    CHECK(d2psi_v_list[0].get()[2] == ComplexApprox(d2psiM_ref[iat][2]).epsilon(h2));

    CHECK(dspsi_v_list[0].get()[0] == ComplexApprox(dspsiM_ref[iat][0]).epsilon(h));
    CHECK(dspsi_v_list[0].get()[1] == ComplexApprox(dspsiM_ref[iat][1]).epsilon(h));
    CHECK(dspsi_v_list[0].get()[2] == ComplexApprox(dspsiM_ref[iat][2]).epsilon(h));

    //walker 1, permuted from reference
    CHECK(psi_v_list[1].get()[0] == ComplexApprox(psiM_ref[(iat + 1) % 3][0]).epsilon(h));
    CHECK(psi_v_list[1].get()[1] == ComplexApprox(psiM_ref[(iat + 1) % 3][1]).epsilon(h));
    CHECK(psi_v_list[1].get()[2] == ComplexApprox(psiM_ref[(iat + 1) % 3][2]).epsilon(h));

    CHECK(dpsi_v_list[1].get()[0][0] == ComplexApprox(dpsiM_ref[(iat + 1) % 3][0][0]).epsilon(h));
    CHECK(dpsi_v_list[1].get()[0][1] == ComplexApprox(dpsiM_ref[(iat + 1) % 3][0][1]).epsilon(h));
    CHECK(dpsi_v_list[1].get()[0][2] == ComplexApprox(dpsiM_ref[(iat + 1) % 3][0][2]).epsilon(h));

    CHECK(dpsi_v_list[1].get()[1][0] == ComplexApprox(dpsiM_ref[(iat + 1) % 3][1][0]).epsilon(h));
    CHECK(dpsi_v_list[1].get()[1][1] == ComplexApprox(dpsiM_ref[(iat + 1) % 3][1][1]).epsilon(h));
    CHECK(dpsi_v_list[1].get()[1][2] == ComplexApprox(dpsiM_ref[(iat + 1) % 3][1][2]).epsilon(h));

    CHECK(dpsi_v_list[1].get()[2][0] == ComplexApprox(dpsiM_ref[(iat + 1) % 3][2][0]).epsilon(h));
    CHECK(dpsi_v_list[1].get()[2][1] == ComplexApprox(dpsiM_ref[(iat + 1) % 3][2][1]).epsilon(h));
    CHECK(dpsi_v_list[1].get()[2][2] == ComplexApprox(dpsiM_ref[(iat + 1) % 3][2][2]).epsilon(h));

    CHECK(d2psi_v_list[1].get()[0] == ComplexApprox(d2psiM_ref[(iat + 1) % 3][0]).epsilon(h2));
    CHECK(d2psi_v_list[1].get()[1] == ComplexApprox(d2psiM_ref[(iat + 1) % 3][1]).epsilon(h2));
    CHECK(d2psi_v_list[1].get()[2] == ComplexApprox(d2psiM_ref[(iat + 1) % 3][2]).epsilon(h2));

    CHECK(dspsi_v_list[1].get()[0] == ComplexApprox(dspsiM_ref[(iat + 1) % 3][0]).epsilon(h));
    CHECK(dspsi_v_list[1].get()[1] == ComplexApprox(dspsiM_ref[(iat + 1) % 3][1]).epsilon(h));
    CHECK(dspsi_v_list[1].get()[2] == ComplexApprox(dspsiM_ref[(iat + 1) % 3][2]).epsilon(h));

    std::vector<bool> accept = {false, false};
    elec_.mw_accept_rejectMove(p_list, iat, accept);
  }
}

#endif //QMC_COMPLEX


} // namespace qmcplusplus
