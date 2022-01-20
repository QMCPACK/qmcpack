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

#include "OhmmsData/Libxml2Doc.h"
#include "OhmmsPETE/OhmmsMatrix.h"
#include "DistanceTable.h"
#include "Particle/ParticleSet.h"
#include "Particle/ParticleSetPool.h"
#include "QMCWaveFunctions/WaveFunctionComponent.h"
#include "QMCWaveFunctions/EinsplineSetBuilder.h"
#include "QMCWaveFunctions/EinsplineSpinorSetBuilder.h"

#include <stdio.h>
#include <string>
#include <limits>

using std::string;

namespace qmcplusplus
{
#if !defined(QMC_CUDA) && !defined(ENABLE_OFFLOAD)
TEST_CASE("Hybridrep SPO from HDF diamond_1x1x1", "[wavefunction]")
{
  Communicate* c = OHMMS::Controller;

  ParticleSet::ParticleLayout_t lattice;
  // diamondC_1x1x1
  lattice.R(0, 0) = 3.37316115;
  lattice.R(0, 1) = 3.37316115;
  lattice.R(0, 2) = 0.0;
  lattice.R(1, 0) = 0.0;
  lattice.R(1, 1) = 3.37316115;
  lattice.R(1, 2) = 3.37316115;
  lattice.R(2, 0) = 3.37316115;
  lattice.R(2, 1) = 0.0;
  lattice.R(2, 2) = 3.37316115;

  ParticleSetPool ptcl = ParticleSetPool(c);
  ptcl.setSimulationCell(lattice);
  auto ions_uptr = std::make_unique<ParticleSet>(ptcl.getSimulationCell());
  auto elec_uptr = std::make_unique<ParticleSet>(ptcl.getSimulationCell());
  ParticleSet& ions_(*ions_uptr);
  ParticleSet& elec_(*elec_uptr);

  ions_.setName("ion");
  ptcl.addParticleSet(std::move(ions_uptr));
  ions_.create(2);
  ions_.R[0][0] = 0.0;
  ions_.R[0][1] = 0.0;
  ions_.R[0][2] = 0.0;
  ions_.R[1][0] = 1.68658058;
  ions_.R[1][1] = 1.68658058;
  ions_.R[1][2] = 1.68658058;

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
  elec_.create(2);
  elec_.R[0][0] = 0.4;
  elec_.R[0][1] = 0.0;
  elec_.R[0][2] = 0.0;
  elec_.R[1][0] = 0.0;
  elec_.R[1][1] = 1.0;
  elec_.R[1][2] = 0.0;

  SpeciesSet& tspecies       = elec_.getSpeciesSet();
  int upIdx                  = tspecies.addSpecies("u");
  int chargeIdx              = tspecies.addAttribute("charge");
  tspecies(chargeIdx, upIdx) = -1;

  //diamondC_1x1x1
  const char* particles = "<tmp> \
<determinantset type=\"einspline\" href=\"diamondC_1x1x1.pwscf.h5\" tilematrix=\"1 0 0 0 1 0 0 0 1\" twistnum=\"0\" source=\"ion\" meshfactor=\"1.0\" precision=\"float\" size=\"4\" hybridrep=\"yes\"/> \
</tmp> \
";

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
  SPOSet::ValueMatrix_t psiM(elec_.R.size(), spo->getOrbitalSetSize());
  SPOSet::GradMatrix_t dpsiM(elec_.R.size(), spo->getOrbitalSetSize());
  SPOSet::ValueMatrix_t d2psiM(elec_.R.size(), spo->getOrbitalSetSize());
  spo->evaluate_notranspose(elec_, 0, elec_.R.size(), psiM, dpsiM, d2psiM);

#if !defined(QMC_CUDA) || defined(QMC_COMPLEX)
  // due to the different ordering of bands skip the tests on CUDA+Real builds
  // checking evaluations, reference values are not independently generated.
  // electron 0
  // value
  REQUIRE(std::real(psiM[0][0]) == Approx(-0.6594096422));
  REQUIRE(std::real(psiM[0][1]) == Approx(-1.3352056742));
  // grad
  REQUIRE(std::real(dpsiM[0][0][0]) == Approx(-0.8762991428));
  REQUIRE(std::real(dpsiM[0][0][1]) == Approx(0.0000000044));
  REQUIRE(std::real(dpsiM[0][0][2]) == Approx(0.0000000044));
  REQUIRE(std::real(dpsiM[0][1][0]) == Approx(-0.8603816628));
  REQUIRE(std::real(dpsiM[0][1][1]) == Approx(4.3501935005));
  REQUIRE(std::real(dpsiM[0][1][2]) == Approx(-0.6386129856));
  // lapl
  REQUIRE(std::real(d2psiM[0][0]) == Approx(-4.1090884209));
#if defined(MIXED_PRECISION)
  REQUIRE(std::real(d2psiM[0][1]) == Approx(22.3851032257).epsilon(3e-5));
#else
  REQUIRE(std::real(d2psiM[0][1]) == Approx(22.3851032257));
#endif

  // electron 1
  // value
  REQUIRE(std::real(psiM[1][0]) == Approx(-0.8886948824));
  REQUIRE(std::real(psiM[1][1]) == Approx(1.4194120169));
  // grad
  REQUIRE(std::real(dpsiM[1][0][0]) == Approx(-0.0000183403));
  REQUIRE(std::real(dpsiM[1][0][1]) == Approx(0.1655139178));
  REQUIRE(std::real(dpsiM[1][0][2]) == Approx(-0.0000193077));
  REQUIRE(std::real(dpsiM[1][1][0]) == Approx(-1.3131694794));
  REQUIRE(std::real(dpsiM[1][1][1]) == Approx(-1.1174004078));
  REQUIRE(std::real(dpsiM[1][1][2]) == Approx(-0.8462534547));
  // lapl
  REQUIRE(std::real(d2psiM[1][0]) == Approx(1.3313053846));
  REQUIRE(std::real(d2psiM[1][1]) == Approx(-4.712583065));
#endif
}

TEST_CASE("Hybridrep SPO from HDF diamond_2x1x1", "[wavefunction]")
{
  Communicate* c = OHMMS::Controller;

  ParticleSet::ParticleLayout_t lattice;
  // diamondC_2x1x1
  lattice.R(0, 0) = 6.7463223;
  lattice.R(0, 1) = 6.7463223;
  lattice.R(0, 2) = 0.0;
  lattice.R(1, 0) = 0.0;
  lattice.R(1, 1) = 3.37316115;
  lattice.R(1, 2) = 3.37316115;
  lattice.R(2, 0) = 3.37316115;
  lattice.R(2, 1) = 0.0;
  lattice.R(2, 2) = 3.37316115;

  ParticleSetPool ptcl = ParticleSetPool(c);
  ptcl.setSimulationCell(lattice);
  auto ions_uptr = std::make_unique<ParticleSet>(ptcl.getSimulationCell());
  auto elec_uptr = std::make_unique<ParticleSet>(ptcl.getSimulationCell());
  ParticleSet& ions_(*ions_uptr);
  ParticleSet& elec_(*elec_uptr);

  ions_.setName("ion");
  ptcl.addParticleSet(std::move(ions_uptr));
  ions_.create(4);
  ions_.R[0][0] = 0.0;
  ions_.R[0][1] = 0.0;
  ions_.R[0][2] = 0.0;
  ions_.R[1][0] = 1.68658058;
  ions_.R[1][1] = 1.68658058;
  ions_.R[1][2] = 1.68658058;
  ions_.R[2][0] = 3.37316115;
  ions_.R[2][1] = 3.37316115;
  ions_.R[2][2] = 0.0;
  ions_.R[3][0] = 5.05974173;
  ions_.R[3][1] = 5.05974173;
  ions_.R[3][2] = 1.68658058;

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
  elec_.create(2);
  elec_.R[0][0] = 0.4;
  elec_.R[0][1] = 0.0;
  elec_.R[0][2] = 0.0;
  elec_.R[1][0] = 0.0;
  elec_.R[1][1] = 1.0;
  elec_.R[1][2] = 0.0;

  SpeciesSet& tspecies       = elec_.getSpeciesSet();
  int upIdx                  = tspecies.addSpecies("u");
  int chargeIdx              = tspecies.addAttribute("charge");
  tspecies(chargeIdx, upIdx) = -1;

  //diamondC_2x1x1
  const char* particles = "<tmp> \
<determinantset type=\"einspline\" href=\"diamondC_2x1x1.pwscf.h5\" tilematrix=\"2 0 0 0 1 0 0 0 1\" twistnum=\"0\" source=\"ion\" meshfactor=\"1.0\" precision=\"float\" size=\"4\" hybridrep=\"yes\"/> \
</tmp> \
";

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
  std::cout << std::setprecision(14) << "check r^2 against dr^2. "
            << "r = " << r << " dr = " << dr << std::endl;
  std::cout << "abs(r^2 - dr^2) = " << std::abs(r * r - dot(dr, dr))
            << " epsilon = " << std::numeric_limits<double>::epsilon() << std::endl;
#if defined(MIXED_PRECISION)
  REQUIRE(std::abs(r * r - dot(dr, dr)) < std::numeric_limits<double>::epsilon() * 1e8);
#else
  REQUIRE(std::abs(r * r - dot(dr, dr)) < std::numeric_limits<double>::epsilon());
#endif

  // for vgl
  SPOSet::ValueMatrix_t psiM(elec_.R.size(), spo->getOrbitalSetSize());
  SPOSet::GradMatrix_t dpsiM(elec_.R.size(), spo->getOrbitalSetSize());
  SPOSet::ValueMatrix_t d2psiM(elec_.R.size(), spo->getOrbitalSetSize());
  spo->evaluate_notranspose(elec_, 0, elec_.R.size(), psiM, dpsiM, d2psiM);

#if !defined(QMC_CUDA) || defined(QMC_COMPLEX)
  // real part
  // due to the different ordering of bands skip the tests on CUDA+Real builds
  // checking evaluations, reference values are not independently generated.

  // electron 0
  // value
  REQUIRE(std::real(psiM[0][0]) == Approx(0.6776432991));
  REQUIRE(std::real(psiM[0][1]) == Approx(1.0759553909));
  // grad
  REQUIRE(std::real(dpsiM[0][0][0]) == Approx(0.8782411218));
  REQUIRE(std::real(dpsiM[0][0][1]) == Approx(0.004904394));
  REQUIRE(std::real(dpsiM[0][0][2]) == Approx(-0.0049044029));
  REQUIRE(std::real(dpsiM[0][1][0]) == Approx(1.1041458845));
  REQUIRE(std::real(dpsiM[0][1][1]) == Approx(0.6333346963));
  REQUIRE(std::real(dpsiM[0][1][2]) == Approx(-0.6333346963));
  // lapl
  REQUIRE(std::real(d2psiM[0][0]) == Approx(4.0779185295).epsilon(1e-4));
  REQUIRE(std::real(d2psiM[0][1]) == Approx(-0.7860302329).epsilon(1e-4));

  // electron 1
  // value
  REQUIRE(std::real(psiM[1][0]) == Approx(0.9008999467));
  REQUIRE(std::real(psiM[1][1]) == Approx(1.2383049726));
  // grad
  REQUIRE(std::real(dpsiM[1][0][0]) == Approx(0.0025820041));
  REQUIRE(std::real(dpsiM[1][0][1]) == Approx(-0.1880052537));
  REQUIRE(std::real(dpsiM[1][0][2]) == Approx(-0.0025404284));
  REQUIRE(std::real(dpsiM[1][1][0]) == Approx(0.1069662273));
  REQUIRE(std::real(dpsiM[1][1][1]) == Approx(-0.4364597797));
  REQUIRE(std::real(dpsiM[1][1][2]) == Approx(-0.106951952));
  // lapl
  REQUIRE(std::real(d2psiM[1][0]) == Approx(-1.3757134676));
  REQUIRE(std::real(d2psiM[1][1]) == Approx(-2.4803137779));
#endif

#if defined(QMC_COMPLEX)
  // imaginary part
  // electron 0
  // value
  REQUIRE(std::imag(psiM[0][0]) == Approx(0.6776432991));
  REQUIRE(std::imag(psiM[0][1]) == Approx(1.0762499571));
  // grad
  REQUIRE(std::imag(dpsiM[0][0][0]) == Approx(0.878241539));
  REQUIRE(std::imag(dpsiM[0][0][1]) == Approx(0.0049043936));
  REQUIRE(std::imag(dpsiM[0][0][2]) == Approx(-0.0049044029));
  REQUIRE(std::imag(dpsiM[0][1][0]) == Approx(1.1067043543));
  REQUIRE(std::imag(dpsiM[0][1][1]) == Approx(0.6384321451));
  REQUIRE(std::imag(dpsiM[0][1][2]) == Approx(-0.6384321451));
  // lapl
  REQUIRE(std::imag(d2psiM[0][0]) == Approx(4.0779790878).epsilon(1e-4));
  REQUIRE(std::imag(d2psiM[0][1]) == Approx(-0.7897151113).epsilon(1e-4));
  // electron 1
  // value
  REQUIRE(std::imag(psiM[1][0]) == Approx(0.9008999467));
  REQUIRE(std::imag(psiM[1][1]) == Approx(1.2383049726));
  // grad
  REQUIRE(std::imag(dpsiM[1][0][0]) == Approx(0.0025820041));
  REQUIRE(std::imag(dpsiM[1][0][1]) == Approx(-0.1880052537));
  REQUIRE(std::imag(dpsiM[1][0][2]) == Approx(-0.0025404284));
  REQUIRE(std::imag(dpsiM[1][1][0]) == Approx(0.1069453433));
  REQUIRE(std::imag(dpsiM[1][1][1]) == Approx(-0.43649593));
  REQUIRE(std::imag(dpsiM[1][1][2]) == Approx(-0.1069145575));
  // lapl
  REQUIRE(std::imag(d2psiM[1][0]) == Approx(-1.3757134676));
  REQUIRE(std::imag(d2psiM[1][1]) == Approx(-2.4919104576));
#endif

  // test batched interfaces

  ParticleSet elec_2(elec_);
  // interchange positions
  elec_2.R[0] = elec_.R[1];
  elec_2.R[1] = elec_.R[0];
  elec_2.update();
  RefVectorWithLeader<ParticleSet> P_list(elec_);
  P_list.push_back(elec_);
  P_list.push_back(elec_2);

  std::unique_ptr<SPOSet> spo_2(spo->makeClone());
  RefVectorWithLeader<SPOSet> spo_list(*spo);
  spo_list.push_back(*spo);
  spo_list.push_back(*spo_2);

  SPOSet::ValueVector_t psi(spo->getOrbitalSetSize());
  SPOSet::GradVector_t dpsi(spo->getOrbitalSetSize());
  SPOSet::ValueVector_t d2psi(spo->getOrbitalSetSize());
  SPOSet::ValueVector_t psi_2(spo->getOrbitalSetSize());
  SPOSet::GradVector_t dpsi_2(spo->getOrbitalSetSize());
  SPOSet::ValueVector_t d2psi_2(spo->getOrbitalSetSize());

  RefVector<SPOSet::ValueVector_t> psi_v_list;
  RefVector<SPOSet::GradVector_t> dpsi_v_list;
  RefVector<SPOSet::ValueVector_t> d2psi_v_list;

  psi_v_list.push_back(psi);
  psi_v_list.push_back(psi_2);
  dpsi_v_list.push_back(dpsi);
  dpsi_v_list.push_back(dpsi_2);
  d2psi_v_list.push_back(d2psi);
  d2psi_v_list.push_back(d2psi_2);

  spo->mw_evaluateVGL(spo_list, P_list, 0, psi_v_list, dpsi_v_list, d2psi_v_list);
#if !defined(QMC_CUDA) || defined(QMC_COMPLEX)
  // real part
  // due to the different ordering of bands skip the tests on CUDA+Real builds
  // checking evaluations, reference values are not independently generated.
  // value
  REQUIRE(std::real(psi_v_list[1].get()[0]) == Approx(0.9008999467));
  REQUIRE(std::real(psi_v_list[1].get()[1]) == Approx(1.2383049726));
  // grad
  REQUIRE(std::real(dpsi_v_list[1].get()[0][0]) == Approx(0.0025820041));
  REQUIRE(std::real(dpsi_v_list[1].get()[0][1]) == Approx(-0.1880052537));
  REQUIRE(std::real(dpsi_v_list[1].get()[0][2]) == Approx(-0.0025404284));
  REQUIRE(std::real(dpsi_v_list[1].get()[1][0]) == Approx(0.1069662273));
  REQUIRE(std::real(dpsi_v_list[1].get()[1][1]) == Approx(-0.4364597797));
  REQUIRE(std::real(dpsi_v_list[1].get()[1][2]) == Approx(-0.106951952));
  // lapl
  REQUIRE(std::real(d2psi_v_list[1].get()[0]) == Approx(-1.3757134676));
  REQUIRE(std::real(d2psi_v_list[1].get()[1]) == Approx(-2.4803137779));
#endif

#if defined(QMC_COMPLEX)
  // imaginary part
  // value
  REQUIRE(std::imag(psi_v_list[1].get()[0]) == Approx(0.9008999467));
  REQUIRE(std::imag(psi_v_list[1].get()[1]) == Approx(1.2383049726));
  // grad
  REQUIRE(std::imag(dpsi_v_list[1].get()[0][0]) == Approx(0.0025820041));
  REQUIRE(std::imag(dpsi_v_list[1].get()[0][1]) == Approx(-0.1880052537));
  REQUIRE(std::imag(dpsi_v_list[1].get()[0][2]) == Approx(-0.0025404284));
  REQUIRE(std::imag(dpsi_v_list[1].get()[1][0]) == Approx(0.1069453433));
  REQUIRE(std::imag(dpsi_v_list[1].get()[1][1]) == Approx(-0.43649593));
  REQUIRE(std::imag(dpsi_v_list[1].get()[1][2]) == Approx(-0.1069145575));
  // lapl
  REQUIRE(std::imag(d2psi_v_list[1].get()[0]) == Approx(-1.3757134676));
  REQUIRE(std::imag(d2psi_v_list[1].get()[1]) == Approx(-2.4919104576));
#endif
}
#endif

} // namespace qmcplusplus
