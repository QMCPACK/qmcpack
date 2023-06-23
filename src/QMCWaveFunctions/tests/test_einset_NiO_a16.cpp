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
#include "Particle/ParticleSet.h"
#include "Particle/ParticleSetPool.h"
#include "QMCWaveFunctions/WaveFunctionComponent.h"
#include "QMCWaveFunctions/EinsplineSetBuilder.h"
#include "QMCWaveFunctions/EinsplineSpinorSetBuilder.h"
#include <ResourceCollection.h>

#include <stdio.h>
#include <string>
#include <limits>

using std::string;

namespace qmcplusplus
{

TEST_CASE("Einspline SPO from HDF NiO a16 97 electrons", "[wavefunction]")
{
  Communicate* c = OHMMS::Controller;

  ParticleSet::ParticleLayout lattice;
  lattice.R = {3.94055, 3.94055, 7.8811, 3.94055, 3.94055, -7.8811, -7.8811, 7.8811, 0};

  ParticleSetPool ptcl = ParticleSetPool(c);
  ptcl.setSimulationCell(lattice);
  auto ions_uptr = std::make_unique<ParticleSet>(ptcl.getSimulationCell());
  auto elec_uptr = std::make_unique<ParticleSet>(ptcl.getSimulationCell());
  ParticleSet& ions_(*ions_uptr);
  ParticleSet& elec_(*elec_uptr);

  ions_.setName("ion");
  ptcl.addParticleSet(std::move(ions_uptr));
  ions_.create({8, 8});
  ions_.R[0]  = {0.75, 0.25, 0};
  ions_.R[1]  = {0.75, 0.25, 0.5};
  ions_.R[2]  = {0.75, 0.75, 0.25};
  ions_.R[3]  = {0.75, 0.75, 0.75};
  ions_.R[4]  = {0.25, 0.75, 0};
  ions_.R[5]  = {0.25, 0.25, 0.25};
  ions_.R[6]  = {0.25, 0.75, 0.5};
  ions_.R[7]  = {0.25, 0.25, 0.75};
  ions_.R[8]  = {0, 0, 0};
  ions_.R[9]  = {0, 0, 0.5};
  ions_.R[10] = {0, 0.5, 0.25};
  ions_.R[11] = {0, 0.5, 0.75};
  ions_.R[12] = {0.5, 0.5, 0};
  ions_.R[13] = {0.5, 0, 0.25};
  ions_.R[14] = {0.5, 0.5, 0.5};
  ions_.R[15] = {0.5, 0, 0.75};

  // convert to cartesian
  ions_.R.setUnit(PosUnit::Lattice);
  ions_.convert2Cart(ions_.R);

  // printout coordinates
  SpeciesSet& ispecies = ions_.getSpeciesSet();
  int Oidx             = ispecies.addSpecies("O");
  int Niidx            = ispecies.addSpecies("Ni");
  ions_.print(app_log());

  elec_.setName("elec");
  ptcl.addParticleSet(std::move(elec_uptr));
  elec_.create({97});
  elec_.R[0] = {0.0, 0.0, 0.0};
  elec_.R[1] = {0.0, 1.0, 0.0};
  elec_.R[2] = {0.0, 1.1, 0.0};
  elec_.R[3] = {0.0, 1.2, 0.0};
  elec_.R[4] = {0.0, 1.3, 0.0};

  SpeciesSet& tspecies       = elec_.getSpeciesSet();
  int upIdx                  = tspecies.addSpecies("u");
  int chargeIdx              = tspecies.addAttribute("charge");
  tspecies(chargeIdx, upIdx) = -1;

  //diamondC_2x1x1
  const char* particles = R"(<tmp>
<determinantset type="einspline" href="NiO-fcc-supertwist111-supershift000-S4.h5" tilematrix="1 0 0 0 1 1 0 2 -2" twistnum="0" source="ion" meshfactor="1.0" precision="float" size="97" gpu="omptarget"/>
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

  // for vgl
  SPOSet::ValueMatrix psiM(elec_.R.size(), spo->getOrbitalSetSize());
  SPOSet::GradMatrix dpsiM(elec_.R.size(), spo->getOrbitalSetSize());
  SPOSet::ValueMatrix d2psiM(elec_.R.size(), spo->getOrbitalSetSize());
  spo->evaluate_notranspose(elec_, 0, elec_.R.size(), psiM, dpsiM, d2psiM);

#if !defined(QMC_COMPLEX)
  // real part
  // due to the different ordering of bands skip the tests on CUDA+Real builds
  // checking evaluations, reference values are not independently generated.
  // value
  CHECK(std::real(psiM[1][0]) == Approx(-2.6693785191));
  CHECK(std::real(psiM[1][1]) == Approx(-2.669519186));
  // grad
  CHECK(std::real(dpsiM[1][0][0]) == Approx(-0.0002679008));
  CHECK(std::real(dpsiM[1][0][1]) == Approx(7.622625351));
  CHECK(std::real(dpsiM[1][0][2]) == Approx(-0.0003001692));
  CHECK(std::real(dpsiM[1][1][0]) == Approx(-0.0002825252));
  CHECK(std::real(dpsiM[1][1][1]) == Approx(7.6229834557));
  CHECK(std::real(dpsiM[1][1][2]) == Approx(-0.0002339484));
  // lapl
  CHECK(std::real(d2psiM[1][0]) == Approx(-2.6130394936).epsilon(0.0001));
  CHECK(std::real(d2psiM[1][1]) == Approx(-2.6067698002).epsilon(0.0001));
#endif

#if defined(QMC_COMPLEX)
  // imaginary part
  // value
  CHECK(std::imag(psiM[1][0]) == Approx(2.6693942547));
  CHECK(std::imag(psiM[1][1]) == Approx(-2.6693785191));
  // grad
  CHECK(std::imag(dpsiM[1][0][0]) == Approx(-0.000179037));
  CHECK(std::imag(dpsiM[1][0][1]) == Approx(-7.6221408844));
  CHECK(std::imag(dpsiM[1][0][2]) == Approx(-0.0000232533));
  CHECK(std::imag(dpsiM[1][1][0]) == Approx(-0.0002679008));
  CHECK(std::imag(dpsiM[1][1][1]) == Approx(7.622625351));
  CHECK(std::imag(dpsiM[1][1][2]) == Approx(-0.0003001692));
  // lapl
  CHECK(std::imag(d2psiM[1][0]) == Approx(2.6158618927).epsilon(0.0001));
  CHECK(std::imag(d2psiM[1][1]) == Approx(-2.6130394936).epsilon(0.0001));
#endif

  // test batched interfaces

  ParticleSet elec_2(elec_);
  // interchange positions
  elec_2.R[0] = elec_.R[1];
  elec_2.R[1] = elec_.R[0];
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
#if !defined(QMC_COMPLEX)
  // real part
  // due to the different ordering of bands skip the tests on CUDA+Real builds
  // checking evaluations, reference values are not independently generated.
  // value
  CHECK(std::real(psi_v_list[1].get()[0]) == Approx(-2.6693785191));
  CHECK(std::real(psi_v_list[1].get()[1]) == Approx(-2.669519186));
  // grad
  CHECK(std::real(dpsi_v_list[1].get()[0][0]) == Approx(-0.0002679008));
  CHECK(std::real(dpsi_v_list[1].get()[0][1]) == Approx(7.622625351));
  CHECK(std::real(dpsi_v_list[1].get()[0][2]) == Approx(-0.0003001692));
  CHECK(std::real(dpsi_v_list[1].get()[1][0]) == Approx(-0.0002825252));
  CHECK(std::real(dpsi_v_list[1].get()[1][1]) == Approx(7.6229834557));
  CHECK(std::real(dpsi_v_list[1].get()[1][2]) == Approx(-0.0002339484));
  // lapl
  CHECK(std::real(d2psi_v_list[1].get()[0]) == Approx(-2.6130394936).epsilon(0.0001));
  CHECK(std::real(d2psi_v_list[1].get()[1]) == Approx(-2.6067698002).epsilon(0.0001));
#endif

#if defined(QMC_COMPLEX)
  // imaginary part
  // value
  CHECK(std::imag(psi_v_list[1].get()[0]) == Approx(2.6693942547));
  CHECK(std::imag(psi_v_list[1].get()[1]) == Approx(-2.6693785191));
  // grad
  CHECK(std::imag(dpsi_v_list[1].get()[0][0]) == Approx(-0.000179037));
  CHECK(std::imag(dpsi_v_list[1].get()[0][1]) == Approx(-7.6221408844));
  CHECK(std::imag(dpsi_v_list[1].get()[0][2]) == Approx(-0.0000232533));
  CHECK(std::imag(dpsi_v_list[1].get()[1][0]) == Approx(-0.0002679008));
  CHECK(std::imag(dpsi_v_list[1].get()[1][1]) == Approx(7.622625351));
  CHECK(std::imag(dpsi_v_list[1].get()[1][2]) == Approx(-0.0003001692));
  // lapl
  CHECK(std::imag(d2psi_v_list[1].get()[0]) == Approx(2.6158618927).epsilon(0.0001));
  CHECK(std::imag(d2psi_v_list[1].get()[1]) == Approx(-2.6130394936).epsilon(0.0001));
#endif

  const size_t nw = 2;
  std::vector<SPOSet::ValueType> ratio_v(nw);
  std::vector<SPOSet::GradType> grads_v(nw);

  Vector<SPOSet::ValueType, OffloadPinnedAllocator<SPOSet::ValueType>> inv_row(5);
  inv_row[0] = 0.1;
  inv_row[1] = 0.2;
  inv_row[2] = 0.3;
  inv_row[3] = 0.4;
  inv_row[4] = 0.5;
  inv_row.updateTo();

  std::vector<const SPOSet::ValueType*> inv_row_ptr(nw, inv_row.device_data());

  SPOSet::OffloadMWVGLArray phi_vgl_v;
  phi_vgl_v.resize(QMCTraits::DIM_VGL, nw, 5);
  spo->mw_evaluateVGLandDetRatioGrads(spo_list, p_list, 0, inv_row_ptr, phi_vgl_v, ratio_v, grads_v);
  phi_vgl_v.updateFrom();
#if !defined(QMC_COMPLEX)
  CHECK(std::real(ratio_v[0]) == Approx(-0.4838374162));
  CHECK(std::real(grads_v[0][0]) == Approx(-24.6573209338));
  CHECK(std::real(grads_v[0][1]) == Approx(24.6573187288));
  CHECK(std::real(grads_v[0][2]) == Approx(-0.0000002709));
  CHECK(std::real(phi_vgl_v(0, 0, 0)) == Approx(-1.6125482321));
  CHECK(std::real(ratio_v[1]) == Approx(-2.4885484445));
  CHECK(std::real(grads_v[1][0]) == Approx(-0.6732621026));
  CHECK(std::real(grads_v[1][1]) == Approx(-2.6037152796));
  CHECK(std::real(grads_v[1][2]) == Approx(-0.0001673779));
  CHECK(std::real(phi_vgl_v(0, 1, 0)) == Approx(-2.6693785191));
#else
  CHECK(std::imag(ratio_v[0]) == Approx(0.0005530113));
  CHECK(std::imag(grads_v[0][0]) == Approx(-0.000357729));
  CHECK(std::imag(grads_v[0][1]) == Approx(-0.0005462062));
  CHECK(std::imag(grads_v[0][2]) == Approx(0.0004516766));
  CHECK(std::imag(phi_vgl_v(0, 0, 0)) == Approx(1.6124026775));
  CHECK(std::imag(ratio_v[1]) == Approx(0.0009635784));
  CHECK(std::imag(grads_v[1][0]) == Approx(0.0001453869));
  CHECK(std::imag(grads_v[1][1]) == Approx(-0.0000331457));
  CHECK(std::imag(grads_v[1][2]) == Approx(0.0000295465));
  CHECK(std::imag(phi_vgl_v(0, 1, 0)) == Approx(2.6693942547));
#endif
}
} // namespace qmcplusplus
