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
#include "QMCHamiltonians/NLPPJob.h"
#include "DistanceTable.h"
#include <ResourceCollection.h>

#include <stdio.h>
#include <string>
#include <limits>
#include <regex>

using std::string;

namespace qmcplusplus
{
template<typename T>
using OffloadVector = Vector<T, OffloadPinnedAllocator<T>>;

void test_einset_LiH_x(bool use_offload)
{
  Communicate* c = OHMMS::Controller;

  Lattice lattice;
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
  ions_.update();

  elec_.create({1, 1});
  elec_.setName("elec");
  ptcl.addParticleSet(std::move(elec_uptr));
  elec_.create({2});
  elec_.R[0]                = {0.0, 0.0, 0.0};
  elec_.R[1]                = {0.0, 1.0, 0.0};
  const auto ei_table_index = elec_.addTable(ions_);
  elec_.update();

  SpeciesSet& tspecies       = elec_.getSpeciesSet();
  int upIdx                  = tspecies.addSpecies("u");
  int chargeIdx              = tspecies.addAttribute("charge");
  tspecies(chargeIdx, upIdx) = -1;

  // LiH
  std::string spo_xml = R"XML(
<sposet_collection type="einspline" href="LiH-x.pwscf.h5" tilematrix="1 0 0 0 1 0 0 0 1" source="ion" precision="float" gpu="omptarget" twist="0.5 0.0 0.0">
  <sposet name="updet" size="6"/>
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
  CHECK(std::real(d2psiM[1][0]) == Approx(1.7295608521).epsilon(4e-5));
  CHECK(std::real(d2psiM[1][1]) == Approx(0.7432643771).epsilon(2e-5));
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
  const size_t nw = 2;
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
  ParticleSet::mw_update(p_list);

  const size_t norb = 5;
  //Test SplineR2R/C2C::mw_evaluate_notranspose
  {
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
  }
  //Test SplineR2R/C2C::mw_evaluateVGLandDetRatioGrads
  {
    std::vector<SPOSet::ValueType> ratio_v(nw);
    std::vector<SPOSet::GradType> grads_v(nw);

    OffloadVector<SPOSet::ValueType> inv_row(norb);
    inv_row = {0.1, 0.2, 0.3, 0.4, 0.5};
    inv_row.updateTo();

    std::vector<const SPOSet::ValueType*> inv_row_ptr(nw, use_offload ? inv_row.device_data() : inv_row.data());

    SPOSet::OffloadMWVGLArray phi_vgl_v;
    phi_vgl_v.resize(QMCTraits::DIM_VGL, nw, norb);
    spo->mw_evaluateVGLandDetRatioGrads(spo_list, p_list, 0, inv_row_ptr, phi_vgl_v, ratio_v, grads_v);
#if defined(QMC_COMPLEX)
    CHECK(ratio_v[0] == ComplexApprox(std::complex{0.111643, 0.111644}));
    CHECK(grads_v[0][0] == ComplexApprox(std::complex{-3.73874, -22.8802}));
    CHECK(grads_v[0][1] == ComplexApprox(std::complex{-7.20916, 10.4173}));
    CHECK(grads_v[0][2] == ComplexApprox(std::complex{-10.9493, 0.954407}));
    CHECK(ratio_v[1] == ComplexApprox(std::complex{-0.494879, 0.363164}));
    CHECK(grads_v[1][0] == ComplexApprox(std::complex{-1.85816, 0.899594}));
    CHECK(grads_v[1][1] == ComplexApprox(std::complex{0.24621, -0.28615}));
    CHECK(grads_v[1][2] == ComplexApprox(std::complex{0.147884, 0.67006}));
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

  //Test SplineR2R/C2C::mw_evaluateDetRatios")
  {
    // make VPs
    const size_t nvp_                  = 4;
    const size_t nvp_2                 = 3;
    const std::vector<size_t> nvp_list = {nvp_, nvp_2};
    VirtualParticleSet VP_(elec_);
    VirtualParticleSet VP_2(elec_2);

    // move VPs
    std::vector<ParticleSet::SingleParticlePos> newpos_vp_(nvp_);
    std::vector<ParticleSet::SingleParticlePos> newpos_vp_2(nvp_2);
    for (int i = 0; i < nvp_; i++)
    {
      newpos_vp_[i][0] = 0.1 * i;
      newpos_vp_[i][1] = 0.2 * i;
      newpos_vp_[i][2] = 0.3 * i;
    }
    for (int i = 0; i < nvp_2; i++)
    {
      newpos_vp_2[i][0] = 0.2 * i;
      newpos_vp_2[i][1] = 0.3 * i;
      newpos_vp_2[i][2] = 0.4 * i;
    }

    const auto& ei_table_  = elec_.getDistTableAB(ei_table_index);
    const auto& ei_table_2 = elec_2.getDistTableAB(ei_table_index);
    // make virtual move of elec 0, reference ion 1
    NLPPJob<SPOSet::RealType> job_(1, 0, ei_table_.getDistances()[0][1], -ei_table_.getDisplacements()[0][1]);
    // make virtual move of elec 1, reference ion 0
    NLPPJob<SPOSet::RealType> job_2(0, 1, ei_table_2.getDistances()[1][0], -ei_table_2.getDisplacements()[1][0]);

    // make VP refvec
    RefVectorWithLeader<VirtualParticleSet> vp_list(VP_, {VP_, VP_2});

    ResourceCollection vp_res("test_vp_res");
    VP_.createResource(vp_res);
    ResourceCollectionTeamLock<VirtualParticleSet> mw_vpset_lock(vp_res, vp_list);
    VirtualParticleSet::mw_makeMoves(vp_list, p_list, {newpos_vp_, newpos_vp_2}, {job_, job_2}, false);

    // fill invrow with dummy data for each walker
    OffloadVector<SPOSet::ValueType> psiMinv_data_0(norb), psiMinv_data_1(norb);
    SPOSet::ValueVector psiMinv_ref_0(psiMinv_data_0.data(), norb), psiMinv_ref_1(psiMinv_data_1.data(), norb);
    for (int i = 0; i < norb; i++)
    {
      psiMinv_data_0[i] = 0.1 * i;
      psiMinv_data_1[i] = 0.2 * i;
    }
    psiMinv_data_0.updateTo();
    psiMinv_data_1.updateTo();

    std::vector<const SPOSet::ValueType*> invRow_ptr_list;
    if (spo->isOMPoffload())
      invRow_ptr_list = {psiMinv_data_0.device_data(), psiMinv_data_1.device_data()};
    else
      invRow_ptr_list = {psiMinv_data_0.data(), psiMinv_data_1.data()};

    // ratios_list
    std::vector<std::vector<SPOSet::ValueType>> ratios_list(nw);
    for (size_t iw = 0; iw < nw; iw++)
      ratios_list[iw].resize(nvp_list[iw]);

    // just need dummy refvec with correct size
    SPOSet::ValueVector tmp_psi0(norb), tmp_psi1(norb);
    RefVector<SPOSet::ValueVector> tmp_psi_list{tmp_psi0, tmp_psi1};
    spo->mw_evaluateDetRatios(spo_list, RefVectorWithLeader<const VirtualParticleSet>(VP_, {VP_, VP_2}), tmp_psi_list,
                              invRow_ptr_list, ratios_list);

    std::vector<SPOSet::ValueType> ratios_ref_0(nvp_);
    std::vector<SPOSet::ValueType> ratios_ref_1(nvp_2);
    // single-walker functions for reference
    spo->evaluateDetRatios(VP_, tmp_psi0, psiMinv_ref_0, ratios_ref_0);
    spo_2->evaluateDetRatios(VP_2, tmp_psi1, psiMinv_ref_1, ratios_ref_1);

    for (int ivp = 0; ivp < nvp_; ivp++)
      CHECK(Approx(std::real(ratios_list[0][ivp])) == std::real(ratios_ref_0[ivp]));
    for (int ivp = 0; ivp < nvp_2; ivp++)
      CHECK(Approx(std::real(ratios_list[1][ivp])) == std::real(ratios_ref_1[ivp]));
    for (int ivp = 0; ivp < nvp_; ivp++)
      CHECK(Approx(std::imag(ratios_list[0][ivp])) == std::imag(ratios_ref_0[ivp]));
    for (int ivp = 0; ivp < nvp_2; ivp++)
      CHECK(Approx(std::imag(ratios_list[1][ivp])) == std::imag(ratios_ref_1[ivp]));

#if defined(QMC_COMPLEX)
    CHECK(ComplexApprox(ratios_ref_0[0]) == std::complex{-0.749192, -0.749192});
    CHECK(ComplexApprox(ratios_ref_0[1]) == std::complex{-0.613026, -0.608225});
    CHECK(ComplexApprox(ratios_ref_0[2]) == std::complex{-0.387559, -0.40694});
    CHECK(ComplexApprox(ratios_ref_0[3]) == std::complex{-0.220177, -0.299219});
    CHECK(ComplexApprox(ratios_ref_1[0]) == std::complex{-1.49838, -1.49838});
    CHECK(ComplexApprox(ratios_ref_1[1]) == std::complex{-0.87387, -1.10799});
    CHECK(ComplexApprox(ratios_ref_1[2]) == std::complex{-0.356906, -0.748977});

    CHECK(ComplexApprox(ratios_list[0][0]) == std::complex{-0.749192, -0.749192});
    CHECK(ComplexApprox(ratios_list[0][1]) == std::complex{-0.613026, -0.608225});
    CHECK(ComplexApprox(ratios_list[0][2]) == std::complex{-0.387559, -0.40694});
    CHECK(ComplexApprox(ratios_list[0][3]) == std::complex{-0.220177, -0.299219});
    CHECK(ComplexApprox(ratios_list[1][0]) == std::complex{-1.49838, -1.49838});
    CHECK(ComplexApprox(ratios_list[1][1]) == std::complex{-0.87387, -1.10799});
    CHECK(ComplexApprox(ratios_list[1][2]) == std::complex{-0.356906, -0.748977});
#else
    CHECK(ratios_ref_0[0] == Approx(-0.7491922272));
    CHECK(ratios_ref_0[1] == Approx(-0.6130255755));
    CHECK(ratios_ref_0[2] == Approx(-0.3875585041));
    CHECK(ratios_ref_0[3] == Approx(-0.2201766206));
    CHECK(ratios_ref_1[0] == Approx(-1.4983844545));
    CHECK(ratios_ref_1[1] == Approx(0.8193794169));
    CHECK(ratios_ref_1[2] == Approx(0.3709115237));

    CHECK(ratios_list[0][0] == Approx(-0.7491922272));
    CHECK(ratios_list[0][1] == Approx(-0.6130255755));
    CHECK(ratios_list[0][2] == Approx(-0.3875585041));
    CHECK(ratios_list[0][3] == Approx(-0.2201766206));
    CHECK(ratios_list[1][0] == Approx(-1.4983844545));
    CHECK(ratios_list[1][1] == Approx(0.8193794169));
    CHECK(ratios_list[1][2] == Approx(0.3709115237));
#endif
    // // print SW ref values
    // for (int ivp = 0; ivp < nvp_; ivp++)
    //   app_log() << "CHECK(Approx(std::real(ratios_ref_0[" << ivp << "])) == " << std::real(ratios_ref_0[ivp]) << ");\n";
    // for (int ivp = 0; ivp < nvp_2; ivp++)
    //   app_log() << "CHECK(Approx(std::real(ratios_ref_1[" << ivp << "])) == " << std::real(ratios_ref_1[ivp]) << ");\n";

    // // print MW ref values
    // for (int iw = 0; iw < nw; iw++)
    //   for (int ivp = 0; ivp < nvp_list[iw]; ivp++)
    //     app_log() << "CHECK(Approx(std::real(ratios_list[" << iw << "][" << ivp
    //               << "])) == " << std::real(ratios_list[iw][ivp]) << ");\n";
    // for (int iw = 0; iw < nw; iw++)
    //   for (int ivp = 0; ivp < nvp_list[iw]; ivp++)
    //     app_log() << "CHECK(Approx(std::imag(ratios_list[" << iw << "][" << ivp
    //               << "])) == " << std::imag(ratios_list[iw][ivp]) << ");\n";
  }
}

TEST_CASE("Einspline SPO from HDF LiH-x", "[wavefunction]")
{
  test_einset_LiH_x(true);
  //test_einset_LiH_x(false);
}
} // namespace qmcplusplus
