//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Mark Dewing, markdewing@gmail.com, University of Illinois at Urbana-Champaign
//
// File created by: Mark Dewing, markdewing@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


#include "catch.hpp"

#include <cstdio>
#include <string>
#include <limits>

#include "OhmmsData/Libxml2Doc.h"
#include "OhmmsPETE/OhmmsMatrix.h"
#include "DistanceTable.h"
#include "Particle/ParticleSet.h"
#include "Particle/ParticleSetPool.h"
#include "QMCWaveFunctions/WaveFunctionComponent.h"
#include "BsplineFactory/EinsplineSetBuilder.h"
#include "BsplineFactory/EinsplineSpinorSetBuilder.h"
#include <ResourceCollection.h>

using std::string;

namespace qmcplusplus
{
TEST_CASE("Hybridrep SPO from HDF diamond_1x1x1", "[wavefunction]")
{
  Communicate* c = OHMMS::Controller;

  ParticleSet::ParticleLayout lattice;
  // diamondC_1x1x1
  lattice.R = {3.37316115, 3.37316115, 0.0, 0.0, 3.37316115, 3.37316115, 3.37316115, 0.0, 3.37316115};

  ParticleSetPool ptcl = ParticleSetPool(c);
  ptcl.setSimulationCell(lattice);
  auto ions_uptr = std::make_unique<ParticleSet>(ptcl.getSimulationCell());
  auto elec_uptr = std::make_unique<ParticleSet>(ptcl.getSimulationCell());
  ParticleSet& ions_(*ions_uptr);
  ParticleSet& elec_(*elec_uptr);

  ions_.setName("ion");
  ptcl.addParticleSet(std::move(ions_uptr));
  ions_.create({2});
  ions_.R[0]              = {0.0, 0.0, 0.0};
  ions_.R[1]              = {1.68658058, 1.68658058, 1.68658058};
  SpeciesSet& ion_species = ions_.getSpeciesSet();
  int C_Idx               = ion_species.addSpecies("C");
  int C_chargeIdx         = ion_species.addAttribute("charge");
  int cutoffIdx           = ion_species.addAttribute("cutoff_radius");
  int lmaxIdx             = ion_species.addAttribute("lmax");

  ion_species(C_chargeIdx, C_Idx) = 4;
  ion_species(cutoffIdx, C_Idx)   = 0.9;
  ion_species(lmaxIdx, C_Idx)     = 3;

  elec_.setName("elec");
  ptcl.addParticleSet(std::move(elec_uptr));
  elec_.create({3});
  elec_.R[0]                 = {0.4, 0.0, 0.0};
  elec_.R[1]                 = {0.0, 1.0, 0.0};
  elec_.R[2]                 = {0.0, 0.5, 0.5};
  SpeciesSet& tspecies       = elec_.getSpeciesSet();
  int upIdx                  = tspecies.addSpecies("u");
  int chargeIdx              = tspecies.addAttribute("charge");
  tspecies(chargeIdx, upIdx) = -1;

  //diamondC_1x1x1
  const char* particles = R"(<tmp>
<determinantset type="einspline" href="diamondC_1x1x1.pwscf.h5" tilematrix="1 0 0 0 1 0 0 0 1" twistnum="0" source="ion" meshfactor="1.0" precision="float" size="4" hybridrep="yes"/>
</tmp>
)";

  Libxml2Document doc;
  bool okay = doc.parseFromString(particles);
  REQUIRE(okay);

  xmlNodePtr root = doc.getRoot();

  xmlNodePtr ein1 = xmlFirstElementChild(root);

  EinsplineSetBuilder einSet(elec_, ptcl.getPool(), c, ein1);
  auto spo = einSet.createSPOSetFromXML(ein1);
  REQUIRE(spo);

  ions_.update();
  elec_.update();

  // for vgl
  SPOSet::ValueMatrix psiM(elec_.R.size(), spo->getOrbitalSetSize());
  SPOSet::GradMatrix dpsiM(elec_.R.size(), spo->getOrbitalSetSize());
  SPOSet::ValueMatrix d2psiM(elec_.R.size(), spo->getOrbitalSetSize());
  spo->evaluate_notranspose(elec_, 0, elec_.R.size(), psiM, dpsiM, d2psiM);

  // due to the different ordering of bands skip the tests on CUDA+Real builds
  // checking evaluations, reference values are not independently generated.
  // electron 0
  // value
  CHECK(std::real(psiM[0][0]) == Approx(-0.6594096422));
  CHECK(std::real(psiM[0][1]) == Approx(-1.3352056742));
  // grad
  CHECK(std::real(dpsiM[0][0][0]) == Approx(-0.8762991428));
  CHECK(std::real(dpsiM[0][0][1]) == Approx(0.0000000044));
  CHECK(std::real(dpsiM[0][0][2]) == Approx(0.0000000044));
  CHECK(std::real(dpsiM[0][1][0]) == Approx(-0.8603816628));
  CHECK(std::real(dpsiM[0][1][1]) == Approx(4.3501935005));
  CHECK(std::real(dpsiM[0][1][2]) == Approx(-0.6386129856));
  // lapl
  CHECK(std::real(d2psiM[0][0]) == Approx(-4.1090884209));
  CHECK(std::real(d2psiM[0][1]) == Approx(22.3851032257).epsilon(3e-5));

  // electron 1
  // value
  CHECK(std::real(psiM[1][0]) == Approx(-0.8886948824));
  CHECK(std::real(psiM[1][1]) == Approx(1.4194120169));
  // grad
  CHECK(std::real(dpsiM[1][0][0]) == Approx(-0.0000183403));
  CHECK(std::real(dpsiM[1][0][1]) == Approx(0.1655139178));
  CHECK(std::real(dpsiM[1][0][2]) == Approx(-0.0000193077));
  CHECK(std::real(dpsiM[1][1][0]) == Approx(-1.3131694794));
  CHECK(std::real(dpsiM[1][1][1]) == Approx(-1.1174004078));
  CHECK(std::real(dpsiM[1][1][2]) == Approx(-0.8462534547));
  // lapl
  CHECK(std::real(d2psiM[1][0]) == Approx(1.3313053846));
  CHECK(std::real(d2psiM[1][1]) == Approx(-4.712583065));

  // electron 2
  // value
  CHECK(std::real(psiM[2][0]) == Approx(-0.8694150782));
  CHECK(std::real(psiM[2][1]) == Approx(0.780789871));
  // grad
  CHECK(std::real(dpsiM[2][0][0]) == Approx(-0.0602699547));
  CHECK(std::real(dpsiM[2][0][1]) == Approx(-0.2770632381));
  CHECK(std::real(dpsiM[2][0][2]) == Approx(-0.2770632488));
  CHECK(std::real(dpsiM[2][1][0]) == Approx(-1.5982062406));
  CHECK(std::real(dpsiM[2][1][1]) == Approx(1.0020672904));
  CHECK(std::real(dpsiM[2][1][2]) == Approx(-1.9794520201));
  // lapl
  CHECK(std::real(d2psiM[2][0]) == Approx(1.1232769428).epsilon(5e-5));
  CHECK(std::real(d2psiM[2][1]) == Approx(-4.9779265738).epsilon(3e-5));
}

TEST_CASE("Hybridrep SPO from HDF diamond_2x1x1", "[wavefunction]")
{
  Communicate* c = OHMMS::Controller;

  ParticleSet::ParticleLayout lattice;
  // diamondC_2x1x1
  lattice.R = {6.7463223, 6.7463223, 0.0, 0.0, 3.37316115, 3.37316115, 3.37316115, 0.0, 3.37316115};

  ParticleSetPool ptcl = ParticleSetPool(c);
  ptcl.setSimulationCell(lattice);
  auto ions_uptr = std::make_unique<ParticleSet>(ptcl.getSimulationCell());
  auto elec_uptr = std::make_unique<ParticleSet>(ptcl.getSimulationCell());
  ParticleSet& ions_(*ions_uptr);
  ParticleSet& elec_(*elec_uptr);

  ions_.setName("ion");
  ptcl.addParticleSet(std::move(ions_uptr));
  ions_.create({4});
  ions_.R[0]              = {0.0, 0.0, 0.0};
  ions_.R[1]              = {1.68658058, 1.68658058, 1.68658058};
  ions_.R[2]              = {3.37316115, 3.37316115, 0.0};
  ions_.R[3]              = {5.05974173, 5.05974173, 1.68658058};
  SpeciesSet& ion_species = ions_.getSpeciesSet();
  int C_Idx               = ion_species.addSpecies("C");
  int C_chargeIdx         = ion_species.addAttribute("charge");
  int cutoffIdx           = ion_species.addAttribute("cutoff_radius");
  int lmaxIdx             = ion_species.addAttribute("lmax");

  ion_species(C_chargeIdx, C_Idx) = 4;
  ion_species(cutoffIdx, C_Idx)   = 0.9;
  ion_species(lmaxIdx, C_Idx)     = 3;

  elec_.setName("elec");
  ptcl.addParticleSet(std::move(elec_uptr));
  elec_.create({3});
  elec_.R[0]                 = {0.4, 0.0, 0.0};
  elec_.R[1]                 = {0.0, 1.0, 0.0};
  elec_.R[2]                 = {0.0, 0.5, 0.5};
  SpeciesSet& tspecies       = elec_.getSpeciesSet();
  int upIdx                  = tspecies.addSpecies("u");
  int chargeIdx              = tspecies.addAttribute("charge");
  tspecies(chargeIdx, upIdx) = -1;

  //diamondC_2x1x1
  const char* particles = R"(<tmp>
<determinantset type="einspline" href="diamondC_2x1x1.pwscf.h5" tilematrix="2 0 0 0 1 0 0 0 1" twistnum="0" source="ion" meshfactor="1.0" precision="float" size="4" hybridrep="yes"/>
</tmp>
)";

  Libxml2Document doc;
  bool okay = doc.parseFromString(particles);
  REQUIRE(okay);

  xmlNodePtr root = doc.getRoot();

  xmlNodePtr ein1 = xmlFirstElementChild(root);

  EinsplineSetBuilder einSet(elec_, ptcl.getPool(), c, ein1);
  auto spo = einSet.createSPOSetFromXML(ein1);
  REQUIRE(spo);

  ions_.update();
  elec_.update();

  ParticleSet::RealType r;
  ParticleSet::PosType dr;
  elec_.getDistTable(0).get_first_neighbor(0, r, dr, false);
  app_log() << std::setprecision(14) << "check r^2 against dr^2. "
            << "r = " << r << " dr = " << dr << std::endl;
  app_log() << "abs(r^2 - dr^2) = " << std::abs(r * r - dot(dr, dr))
            << " epsilon = " << std::numeric_limits<double>::epsilon() << std::endl;
#if defined(MIXED_PRECISION)
  REQUIRE(std::abs(r * r - dot(dr, dr)) < std::numeric_limits<double>::epsilon() * 1e8);
#else
  REQUIRE(std::abs(r * r - dot(dr, dr)) < std::numeric_limits<double>::epsilon());
#endif

  // for vgl
  SPOSet::ValueMatrix psiM(elec_.R.size(), spo->getOrbitalSetSize());
  SPOSet::GradMatrix dpsiM(elec_.R.size(), spo->getOrbitalSetSize());
  SPOSet::ValueMatrix d2psiM(elec_.R.size(), spo->getOrbitalSetSize());
  spo->evaluate_notranspose(elec_, 0, elec_.R.size(), psiM, dpsiM, d2psiM);

#if defined(QMC_COMPLEX)
  using ValueType = SPOSet::ValueType;
  // electron 0
  // value
  CHECK(psiM[0][0] == ComplexApprox(ValueType{0.6776432991, 0.6776432991}));
  CHECK(psiM[0][1] == ComplexApprox(ValueType{1.0759553909, 1.0762499571}));
  // grad
  CHECK(dpsiM[0][0][0] == ComplexApprox(ValueType{0.8782411218, 0.878241539}));
  CHECK(dpsiM[0][0][1] == ComplexApprox(ValueType{0.004904394, 0.0049043936}));
  CHECK(dpsiM[0][0][2] == ComplexApprox(ValueType{-0.0049044029, -0.0049044029}));
  CHECK(dpsiM[0][1][0] == ComplexApprox(ValueType{1.1041458845, 1.1067043543}));
  CHECK(dpsiM[0][1][1] == ComplexApprox(ValueType{0.6333346963, 0.6384321451}));
  CHECK(dpsiM[0][1][2] == ComplexApprox(ValueType{-0.6333346963, -0.6384321451}));
  // lapl
  CHECK(d2psiM[0][0] == ComplexApprox(ValueType{4.0779185295, 4.0779790878}).epsilon(1e-4));
  CHECK(d2psiM[0][1] == ComplexApprox(ValueType{-0.7860302329, -0.7897151113}).epsilon(1e-4));

  // electron 1
  // value
  CHECK(psiM[1][0] == ComplexApprox(ValueType{0.9008999467, 0.9008999467}));
  CHECK(psiM[1][1] == ComplexApprox(ValueType{1.2383049726, 1.2383049726}));
  // grad
  CHECK(dpsiM[1][0][0] == ComplexApprox(ValueType{0.0025820041, 0.0025820041}));
  CHECK(dpsiM[1][0][1] == ComplexApprox(ValueType{-0.1880052537, -0.1880052537}));
  CHECK(dpsiM[1][0][2] == ComplexApprox(ValueType{-0.0025404284, -0.0025404284}));
  CHECK(dpsiM[1][1][0] == ComplexApprox(ValueType{0.1069662273, 0.1069453433}));
  CHECK(dpsiM[1][1][1] == ComplexApprox(ValueType{-0.4364597797, -0.43649593}));
  CHECK(dpsiM[1][1][2] == ComplexApprox(ValueType{-0.106951952, -0.1069145575}));
  // lapl
  CHECK(d2psiM[1][0] == ComplexApprox(ValueType{-1.3757134676, -1.3757134676}));
  CHECK(d2psiM[1][1] == ComplexApprox(ValueType{-2.4803137779, -2.4919104576}));

  // electron 2
  // value
  CHECK(psiM[2][0] == ComplexApprox(ValueType{0.8841851282, 0.8841851282}));
  CHECK(psiM[2][1] == ComplexApprox(ValueType{1.0331713017, 1.0291547321}));
  // grad
  CHECK(dpsiM[2][0][0] == ComplexApprox(ValueType{0.0688574613, 0.0688574613}));
  CHECK(dpsiM[2][0][1] == ComplexApprox(ValueType{0.2735091889, 0.2735091889}));
  CHECK(dpsiM[2][0][2] == ComplexApprox(ValueType{0.2666900514, 0.2666900514}));
  CHECK(dpsiM[2][1][0] == ComplexApprox(ValueType{0.5398793935, 0.5376840012}));
  CHECK(dpsiM[2][1][1] == ComplexApprox(ValueType{0.7525391523, 0.7550587239}));
  CHECK(dpsiM[2][1][2] == ComplexApprox(ValueType{-0.1224437827, -0.1228616516}));
  // lapl
  CHECK(d2psiM[2][0] == ComplexApprox(ValueType{-1.1799273657, -1.2044699918}).epsilon(1e-4));
  CHECK(d2psiM[2][1] == ComplexApprox(ValueType{-2.0339757673, -1.885562226}).epsilon(1e-4));
#endif

  // test batched interfaces
  ParticleSet elec_2(elec_);
  // interchange positions
  elec_2.R[0] = elec_.R[1];
  elec_2.R[1] = elec_.R[0];
  elec_2.update();
  RefVectorWithLeader<ParticleSet> p_list(elec_);
  p_list.push_back(elec_);
  p_list.push_back(elec_2);

  std::unique_ptr<SPOSet> spo_2(spo->makeClone());
  RefVectorWithLeader<SPOSet> spo_list(*spo);
  spo_list.push_back(*spo);
  spo_list.push_back(*spo_2);

  ResourceCollection pset_res("test_pset_res");
  ResourceCollection spo_res("test_spo_res");

  elec_.createResource(pset_res);
  spo->createResource(spo_res);

  ResourceCollectionTeamLock<ParticleSet> mw_pset_lock(pset_res, p_list);
  ResourceCollectionTeamLock<SPOSet> mw_sposet_lock(spo_res, spo_list);

  SPOSet::ValueVector psi(spo->getOrbitalSetSize());
  SPOSet::GradVector dpsi(spo->getOrbitalSetSize());
  SPOSet::ValueVector d2psi(spo->getOrbitalSetSize());
  SPOSet::ValueVector psi_2(spo->getOrbitalSetSize());
  SPOSet::GradVector dpsi_2(spo->getOrbitalSetSize());
  SPOSet::ValueVector d2psi_2(spo->getOrbitalSetSize());

  RefVector<SPOSet::ValueVector> psi_v_list;
  RefVector<SPOSet::GradVector> dpsi_v_list;
  RefVector<SPOSet::ValueVector> d2psi_v_list;

  psi_v_list.push_back(psi);
  psi_v_list.push_back(psi_2);
  dpsi_v_list.push_back(dpsi);
  dpsi_v_list.push_back(dpsi_2);
  d2psi_v_list.push_back(d2psi);
  d2psi_v_list.push_back(d2psi_2);

  spo->mw_evaluateVGL(spo_list, p_list, 0, psi_v_list, dpsi_v_list, d2psi_v_list);
#if defined(QMC_COMPLEX)
  using ValueType = SPOSet::ValueType;
  // value
  CHECK(psi_v_list[1].get()[0] == ComplexApprox(ValueType{0.9008999467, 0.9008999467}));
  CHECK(psi_v_list[1].get()[1] == ComplexApprox(ValueType{1.2383049726, 1.2383049726}));
  // grad
  CHECK(dpsi_v_list[1].get()[0][0] == ComplexApprox(ValueType{0.0025820041, 0.0025820041}));
  CHECK(dpsi_v_list[1].get()[0][1] == ComplexApprox(ValueType{-0.1880052537, -0.1880052537}));
  CHECK(dpsi_v_list[1].get()[0][2] == ComplexApprox(ValueType{-0.0025404284, -0.0025404284}));
  CHECK(dpsi_v_list[1].get()[1][0] == ComplexApprox(ValueType{0.1069662273, 0.1069453433}));
  CHECK(dpsi_v_list[1].get()[1][1] == ComplexApprox(ValueType{-0.4364597797, -0.43649593}));
  CHECK(dpsi_v_list[1].get()[1][2] == ComplexApprox(ValueType{-0.106951952, -0.1069145575}));
  // lapl
  CHECK(d2psi_v_list[1].get()[0] == ComplexApprox(ValueType{-1.3757134676, -1.3757134676}));
  CHECK(d2psi_v_list[1].get()[1] == ComplexApprox(ValueType{-2.4803137779, -2.4919104576}));
#endif
}

} // namespace qmcplusplus
