//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2024 QMCPACK developers.
//
// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////


#include "catch.hpp"

#include "OhmmsData/Libxml2Doc.h"
#include "OhmmsPETE/OhmmsMatrix.h"
#include "Particle/ParticleSet.h"
#include "Particle/ParticleSetPool.h"
#include "QMCWaveFunctions/WaveFunctionComponent.h"
#include "BsplineFactory/EinsplineSetBuilder.h"
#include "BsplineFactory/EinsplineSpinorSetBuilder.h"
#include <ResourceCollection.h>

#include <stdio.h>
#include <string>
#include <limits>
#include <regex>

using std::string;

namespace qmcplusplus
{
void test_einset_LiH_x(bool use_offload)
{
  Communicate* c = OHMMS::Controller;

  ParticleSet::ParticleLayout lattice;
  lattice.R = {-3.55, 0.0, 3.55, 0.0, 3.55, 3.55, -3.55, 3.55, 0.0};

  ParticleSetPool ptcl = ParticleSetPool(c);
  ptcl.setSimulationCell(lattice);
  auto ions_uptr = std::make_unique<ParticleSet>(ptcl.getSimulationCell());
  auto elec_uptr = std::make_unique<ParticleSet>(ptcl.getSimulationCell());
  ParticleSet& ions_(*ions_uptr);
  ParticleSet& elec_(*elec_uptr);

  ions_.setName("ion");
  ptcl.addParticleSet(std::move(ions_uptr));
  ions_.create({1, 1});
  ions_.R[0] = {0.0, 0.0, 0.0};
  ions_.R[1] = {3.55, 3.55, 3.55};

  elec_.create({1, 1});
  elec_.setName("elec");
  ptcl.addParticleSet(std::move(elec_uptr));
  elec_.create({2});
  elec_.R[0] = {0.0, 0.0, 0.0};
  elec_.R[1] = {0.0, 1.0, 0.0};

  SpeciesSet& tspecies       = elec_.getSpeciesSet();
  int upIdx                  = tspecies.addSpecies("u");
  int chargeIdx              = tspecies.addAttribute("charge");
  tspecies(chargeIdx, upIdx) = -1;

  // LiH
  std::string spo_xml = R"XML(
<sposet_collection type="einspline" href="LiH-x.pwscf.h5" tilematrix="1 0 0 0 1 0 0 0 1" source="ion" precision="float" gpu="omptarget" twist="0.5 0.0 0.0">
  <sposet name="updet" size="6">
    <occupation mode="ground"/>
  </sposet>
</sposet_collection>
)XML";

  if (!use_offload)
    spo_xml = std::regex_replace(spo_xml, std::regex("omptarget"), "no");

  Libxml2Document doc;
  bool okay = doc.parseFromString(spo_xml);
  REQUIRE(okay);

  xmlNodePtr root = doc.getRoot();

  xmlNodePtr ein1 = xmlFirstElementChild(root);

  EinsplineSetBuilder einSet(elec_, ptcl.getPool(), c, root);
  auto spo = einSet.createSPOSetFromXML(ein1);
  REQUIRE(spo);

  // Test the case where the number of psi values is not the orbital set size
  // or the number of electrons/2
  const int psi_size = 3;
  SPOSet::ValueMatrix psiM(elec_.R.size(), psi_size);
  SPOSet::GradMatrix dpsiM(elec_.R.size(), psi_size);
  SPOSet::ValueMatrix d2psiM(elec_.R.size(), psi_size);
  spo->evaluate_notranspose(elec_, 0, elec_.R.size(), psiM, dpsiM, d2psiM);

  // value
  CHECK(std::real(psiM[0][0]) == Approx(12.3543100357));
  CHECK(std::real(psiM[0][1]) == Approx(0.0));
#if defined(QMC_COMPLEX)
  CHECK(std::real(psiM[1][0]) == Approx(1.1857329607));
  CHECK(std::real(psiM[1][1]) == Approx(-0.4717386365));
  // grad
  CHECK(std::real(dpsiM[1][0][0]) == Approx(-0.00828881));
  CHECK(std::real(dpsiM[1][0][1]) == Approx(-2.8782308102));
  CHECK(std::real(dpsiM[1][0][2]) == Approx(0.0082882112));
  CHECK(std::real(dpsiM[1][1][0]) == Approx(-0.4088457525));
  CHECK(std::real(dpsiM[1][1][1]) == Approx(-0.2475463897));
  CHECK(std::real(dpsiM[1][1][2]) == Approx(0.4088463187));
  // lapl
  CHECK(std::real(d2psiM[1][0]) == Approx(1.7295608521));
  CHECK(std::real(d2psiM[1][1]) == Approx(0.7432643771));
#else
  CHECK(std::real(psiM[1][0]) == Approx(-1.1857329607));
  CHECK(std::real(psiM[1][1]) == Approx(0.4717386365));
  // grad
  CHECK(std::real(dpsiM[1][0][0]) == Approx(0.0083514443));
  CHECK(std::real(dpsiM[1][0][1]) == Approx(2.8783009052));
  CHECK(std::real(dpsiM[1][0][2]) == Approx(-0.0083516147));
  CHECK(std::real(dpsiM[1][1][0]) == Approx(0.4088440537));
  CHECK(std::real(dpsiM[1][1][1]) == Approx(0.2475452274));
  CHECK(std::real(dpsiM[1][1][2]) == Approx(-0.408844173));
  // lapl
  CHECK(std::real(d2psiM[1][0]) == Approx(-1.7209997177));
  CHECK(std::real(d2psiM[1][1]) == Approx(-0.7427335978));
#endif
  // Test the batched interface
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

  const int ne = elec_.R.size();
  SPOSet::ValueMatrix psi(ne, psi_size);
  SPOSet::GradMatrix dpsi(ne, psi_size);
  SPOSet::ValueMatrix d2psi(ne, psi_size);
  SPOSet::ValueMatrix psi_2(ne, psi_size);
  SPOSet::GradMatrix dpsi_2(ne, psi_size);
  SPOSet::ValueMatrix d2psi_2(ne, psi_size);

  RefVector<SPOSet::ValueMatrix> psi_v_list;
  RefVector<SPOSet::GradMatrix> dpsi_v_list;
  RefVector<SPOSet::ValueMatrix> d2psi_v_list;

  psi_v_list.push_back(psi);
  psi_v_list.push_back(psi_2);
  dpsi_v_list.push_back(dpsi);
  dpsi_v_list.push_back(dpsi_2);
  d2psi_v_list.push_back(d2psi);
  d2psi_v_list.push_back(d2psi_2);
  spo->mw_evaluate_notranspose(spo_list, p_list, 0, elec_.R.size(), psi_v_list, dpsi_v_list, d2psi_v_list);

  CHECK(std::real(psi_v_list[0].get()[0][0]) == Approx(12.3543100357));
#if defined(QMC_COMPLEX)
  CHECK(std::real(psi_v_list[0].get()[1][0]) == Approx(1.1857329607));
  // Second particle set had particle positions flipped
  CHECK(std::real(psi_v_list[1].get()[0][0]) == Approx(1.1857329607));
#else
  CHECK(std::real(psi_v_list[0].get()[1][0]) == Approx(-1.1857329607));
  // Second particle set had particle positions flipped
  CHECK(std::real(psi_v_list[1].get()[0][0]) == Approx(-1.1857329607));
#endif
  CHECK(std::real(psi_v_list[1].get()[1][0]) == Approx(12.3543100357));

  const size_t nw = 2;
  std::vector<SPOSet::ValueType> ratio_v(nw);
  std::vector<SPOSet::GradType> grads_v(nw);

  Vector<SPOSet::ValueType, OffloadPinnedAllocator<SPOSet::ValueType>> inv_row(5);
  inv_row = {0.1, 0.2, 0.3, 0.4, 0.5};
  inv_row.updateTo();

  std::vector<const SPOSet::ValueType*> inv_row_ptr(nw, use_offload ? inv_row.device_data() : inv_row.data());

  SPOSet::OffloadMWVGLArray phi_vgl_v;
  phi_vgl_v.resize(QMCTraits::DIM_VGL, nw, 5);
  spo->mw_evaluateVGLandDetRatioGrads(spo_list, p_list, 0, inv_row_ptr, phi_vgl_v, ratio_v, grads_v);
#if defined(QMC_COMPLEX)
  CHECK(ratio_v[0] == ComplexApprox(std::complex{0.111643,0.111644}));
  CHECK(grads_v[0][0] == ComplexApprox(std::complex{-3.73874,-22.8802}));
  CHECK(grads_v[0][1] == ComplexApprox(std::complex{-7.20916,10.4173}));
  CHECK(grads_v[0][2] == ComplexApprox(std::complex{-10.9493,0.954407}));
  CHECK(ratio_v[1] == ComplexApprox(std::complex{-0.494879,0.363164}));
  CHECK(grads_v[1][0] == ComplexApprox(std::complex{-1.85816,0.899594}));
  CHECK(grads_v[1][1] == ComplexApprox(std::complex{0.24621,-0.28615}));
  CHECK(grads_v[1][2] == ComplexApprox(std::complex{0.147884,0.67006}));
#else
  CHECK(std::real(ratio_v[0]) == Approx(0.1116431436));
  CHECK(std::real(grads_v[0][0]) == Approx(-8.4771142931));
  CHECK(std::real(grads_v[0][1]) == Approx(-45.2452890181));
  CHECK(std::real(grads_v[0][2]) == Approx(15.7149783314));
  CHECK(std::real(ratio_v[1]) == Approx(0.4948832706));
  CHECK(std::real(grads_v[1][0]) == Approx(-1.1981267596));
  CHECK(std::real(grads_v[1][1]) == Approx(0.0360827888));
  CHECK(std::real(grads_v[1][2]) == Approx(0.639729227));
#endif
}

TEST_CASE("Einspline SPO from HDF LiH-x", "[wavefunction]")
{
  test_einset_LiH_x(true);
  test_einset_LiH_x(false);
}
} // namespace qmcplusplus
