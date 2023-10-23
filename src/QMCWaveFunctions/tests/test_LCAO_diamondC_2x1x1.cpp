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
#include "LCAO/LCAOrbitalBuilder.h"
#include <ResourceCollection.h>

#include <stdio.h>
#include <string>
#include <limits>

using std::string;

namespace qmcplusplus
{

TEST_CASE("LCAO DiamondC_2x1x1", "[wavefunction]")
{
  using VT       = SPOSet::ValueType;
  Communicate* c = OHMMS::Controller;

  ParticleSet::ParticleLayout lattice;
  lattice.R = {6.7463223, 6.7463223, 0.0, 0.0, 3.37316115, 3.37316115, 3.37316115, 0.0, 3.37316115};


  ParticleSetPool ptcl = ParticleSetPool(c);
  ptcl.setSimulationCell(lattice);
  auto ions_uptr = std::make_unique<ParticleSet>(ptcl.getSimulationCell());
  auto elec_uptr = std::make_unique<ParticleSet>(ptcl.getSimulationCell());
  ParticleSet& ions_(*ions_uptr);
  ParticleSet& elec_(*elec_uptr);

  ions_.setName("ion0");
  ptcl.addParticleSet(std::move(ions_uptr));
  ions_.create({4});
  ions_.R[0]           = {0.0, 0.0, 0.0};
  ions_.R[1]           = {1.686580575, 1.686580575, 1.686580575};
  ions_.R[2]           = {3.37316115, 3.37316115, 0.0};
  ions_.R[3]           = {5.059741726, 5.059741726, 1.686580575};
  SpeciesSet& ispecies = ions_.getSpeciesSet();
  int Cidx             = ispecies.addSpecies("C");

  ions_.print(app_log());

  elec_.setName("elec");
  ptcl.addParticleSet(std::move(elec_uptr));
  elec_.create({8, 8});
  elec_.R[0] = {0.0, 1.0, 0.0};
  elec_.R[1] = {0.0, 1.1, 0.0};
  elec_.R[2] = {0.0, 1.2, 0.0};
  elec_.R[3] = {0.0, 1.3, 0.0};

  SpeciesSet& tspecies         = elec_.getSpeciesSet();
  int upIdx                    = tspecies.addSpecies("u");
  int downIdx                  = tspecies.addSpecies("d");
  int chargeIdx                = tspecies.addAttribute("charge");
  tspecies(chargeIdx, upIdx)   = -1;
  tspecies(chargeIdx, downIdx) = -1;
  // const int ei_table_index = elec_.addTable(ions_);

  // diamondC_2x1x1
  // from tests/solids/diamondC_2x1x1-Gaussian_pp/C_Diamond-Twist0.wfj.xml
  const std::string wf_xml_str = R"(
  <wavefunction name="psi0" target="e">
      <sposet_collection type="MolecularOrbital" name="LCAOBSet" source="ion0" twist="0  0  0" href="C_Diamond.h5" PBCimages="5  5  5">
      <basisset name="LCAOBSet" key="GTO" transform="no">
        <grid type="log" ri="1.e-6" rf="1.e2" npts="1001"/>
      </basisset>
      <sposet basisset="LCAOBSet" name="spo-up" size="8">
        <occupation mode="ground"/>
        <coefficient size="116" spindataset="0"/>
      </sposet>
      <sposet basisset="LCAOBSet" name="spo-dn" size="8">
        <occupation mode="ground"/>
        <coefficient size="116" spindataset="0"/>
      </sposet>
    </sposet_collection>
  </wavefunction>
  )";
  Libxml2Document doc;
  bool okay = doc.parseFromString(wf_xml_str);
  REQUIRE(okay);

  xmlNodePtr root    = doc.getRoot();
  xmlNodePtr spo_xml = xmlFirstElementChild(root);

  LCAOrbitalBuilder lcaoSet(elec_, ions_, c, spo_xml);
  auto spo = lcaoSet.createSPOSetFromXML(spo_xml);
  REQUIRE(spo);
  auto& lcao_spos = dynamic_cast<const LCAOrbitalSet&>(*spo);
  // CHECK(!lcao_spos.isIdentity());

  const int norb = spo->getOrbitalSetSize();


  // test batched interfaces
  const size_t nw = 2;

  ParticleSet elec_2(elec_);
  // interchange positions
  elec_2.R[0] = elec_.R[1];
  elec_2.R[1] = elec_.R[0];
  elec_.update();
  elec_2.update();
  RefVectorWithLeader<ParticleSet> p_list(elec_, {elec_, elec_2});


  std::unique_ptr<SPOSet> spo_2(spo->makeClone());
  RefVectorWithLeader<SPOSet> spo_list(*spo, {*spo, *spo_2});

  ResourceCollection pset_res("test_pset_res");
  ResourceCollection spo_res("test_spo_res");

  elec_.createResource(pset_res);
  spo->createResource(spo_res);

  ResourceCollectionTeamLock<ParticleSet> mw_pset_lock(pset_res, p_list);
  ResourceCollectionTeamLock<SPOSet> mw_sposet_lock(spo_res, spo_list);

  SPOSet::ValueVector psiref_0(norb);
  SPOSet::GradVector dpsiref_0(norb);
  SPOSet::ValueVector d2psiref_0(norb);
  SPOSet::ValueVector psiref_1(norb);
  SPOSet::GradVector dpsiref_1(norb);
  SPOSet::ValueVector d2psiref_1(norb);
  spo->evaluateVGL(elec_, 0, psiref_0, dpsiref_0, d2psiref_0);
  spo_2->evaluateVGL(elec_2, 0, psiref_1, dpsiref_1, d2psiref_1);

  // app_log() << "vgl_refvalues: \n" << std::setprecision(14);
  // for (int iorb = 0; iorb < 2; iorb++)
  // {
  //   app_log() << "CHECK(std::real(psiref_0[" << iorb << "]) == Approx(" << psiref_0[iorb] << "));\n";
  //   app_log() << "CHECK(std::real(d2psiref_0[" << iorb << "]) == Approx(" << d2psiref_0[iorb] << "));\n";
  //   for (int idim = 0; idim < 3; idim++)
  //     app_log() << "CHECK(std::real(dpsiref_0[" << iorb << "][" << idim << "]) == Approx(" << dpsiref_0[iorb][idim] << "));\n";
  // }
  // for (int iorb = 0; iorb < 2; iorb++)
  // {
  //   app_log() << "CHECK(std::real(psiref_1[" << iorb << "]) == Approx(" << psiref_1[iorb] << "));\n";
  //   app_log() << "CHECK(std::real(d2psiref_1[" << iorb << "]) == Approx(" << d2psiref_1[iorb] << "));\n";
  //   for (int idim = 0; idim < 3; idim++)
  //     app_log() << "CHECK(std::real(dpsiref_1[" << iorb << "][" << idim << "]) == Approx(" << dpsiref_1[iorb][idim] << "));\n";
  // }
  CHECK(std::real(psiref_0[0]) == Approx(0.23966351080363));
  CHECK(std::real(d2psiref_0[0]) == Approx(-0.52701907904562));
  CHECK(std::real(dpsiref_0[0][0]) == Approx(-0.0057821587526225));
  CHECK(std::real(dpsiref_0[0][1]) == Approx(-0.10996780053441));
  CHECK(std::real(dpsiref_0[0][2]) == Approx(0.0057818514970319));
  CHECK(std::real(psiref_0[1]) == Approx(0.26673125102484));
  CHECK(std::real(d2psiref_0[1]) == Approx(-0.56873091511411));
  CHECK(std::real(dpsiref_0[1][0]) == Approx(-1.5495607896545e-06));
  CHECK(std::real(dpsiref_0[1][1]) == Approx(-0.49161303927196));
  CHECK(std::real(dpsiref_0[1][2]) == Approx(-2.5223735120643e-07));

  CHECK(std::real(psiref_1[0]) == Approx(0.22735347471476));
  CHECK(std::real(d2psiref_1[0]) == Approx(-0.41247534585435));
  CHECK(std::real(dpsiref_1[0][0]) == Approx(-0.0064508940042537));
  CHECK(std::real(dpsiref_1[0][1]) == Approx(-0.13393939202647));
  CHECK(std::real(dpsiref_1[0][2]) == Approx(0.0064505548318947));
  CHECK(std::real(psiref_1[1]) == Approx(0.21979955506832));
  CHECK(std::real(d2psiref_1[1]) == Approx(-0.31185449762657));
  CHECK(std::real(dpsiref_1[1][0]) == Approx(-1.7377316644163e-06));
  CHECK(std::real(dpsiref_1[1][1]) == Approx(-0.44562353474568));
  CHECK(std::real(dpsiref_1[1][2]) == Approx(2.5297265181413e-07));

  SPOSet::ValueVector psi_v_1(norb);
  SPOSet::ValueVector psi_v_2(norb);
  RefVector<SPOSet::ValueVector> psi_v_list{psi_v_1, psi_v_2};
  spo->mw_evaluateValue(spo_list, p_list, 0, psi_v_list);

  SPOSet::ValueVector psi_1(norb);
  SPOSet::GradVector dpsi_1(norb);
  SPOSet::ValueVector d2psi_1(norb);
  SPOSet::ValueVector psi_2(norb);
  SPOSet::GradVector dpsi_2(norb);
  SPOSet::ValueVector d2psi_2(norb);
  RefVector<SPOSet::ValueVector> psi_list   = {psi_1, psi_2};
  RefVector<SPOSet::GradVector> dpsi_list   = {dpsi_1, dpsi_2};
  RefVector<SPOSet::ValueVector> d2psi_list = {d2psi_1, d2psi_2};
  spo->mw_evaluateVGL(spo_list, p_list, 0, psi_list, dpsi_list, d2psi_list);

  for (size_t iorb = 0; iorb < norb; iorb++)
  {
    CHECK(psi_list[0].get()[iorb] == Approx(psiref_0[iorb]));
    CHECK(psi_list[1].get()[iorb] == Approx(psiref_1[iorb]));
    CHECK(d2psi_list[0].get()[iorb] == Approx(d2psiref_0[iorb]));
    CHECK(d2psi_list[1].get()[iorb] == Approx(d2psiref_1[iorb]));
    for (size_t idim = 0; idim < SPOSet::DIM; idim++)
    {
      CHECK(dpsi_list[0].get()[iorb][idim] == Approx(dpsiref_0[iorb][idim]));
      CHECK(dpsi_list[1].get()[iorb][idim] == Approx(dpsiref_1[iorb][idim]));
    }
  }

  for (size_t iorb = 0; iorb < norb; iorb++)
    for (size_t iw = 0; iw < nw; iw++)
      CHECK(psi_v_list[iw].get()[iorb] == Approx(psi_list[iw].get()[iorb]));

  // make VPs
  const size_t nvp_                  = 4;
  const size_t nvp_2                 = 3;
  const std::vector<size_t> nvp_list = {nvp_, nvp_2};
  VirtualParticleSet VP_(elec_, nvp_);
  VirtualParticleSet VP_2(elec_2, nvp_2);

  // move VPs
  std::vector<ParticleSet::SingleParticlePos> newpos_vp_(nvp_);
  std::vector<ParticleSet::SingleParticlePos> newpos_vp_2(nvp_2);
  for (int i = 0; i < nvp_; i++)
    newpos_vp_[i] = {double(1.0 * i / nvp_), double(2.0 * i / nvp_), double(i / (2.0 * nvp_))};
  for (int i = 0; i < nvp_2; i++)
    newpos_vp_2[i] = {double(2.0 * i / nvp_2), double(i / (2.0 * nvp_2)), double(1.0 * i / nvp_2)};
  VP_.makeMoves(elec_, 0, newpos_vp_);
  VP_2.makeMoves(elec_2, 0, newpos_vp_2);

  // make VP refvec
  RefVectorWithLeader<VirtualParticleSet> vp_list(VP_, {VP_, VP_2});
  ResourceCollection vp_res("test_vp_res");
  VP_.createResource(vp_res);
  ResourceCollectionTeamLock<VirtualParticleSet> mw_vpset_lock(vp_res, vp_list);

  // fill invrow with dummy data for each walker
  std::vector<SPOSet::ValueType> psiMinv_data_(norb);
  std::vector<SPOSet::ValueType> psiMinv_data_2(norb);
  SPOSet::ValueVector psiMinv_ref_0(norb);
  SPOSet::ValueVector psiMinv_ref_1(norb);
  for (int i = 0; i < norb; i++)
  {
    psiMinv_data_[i]  = (i % 10) / 4.5 - 1.0;
    psiMinv_data_2[i] = ((i + 5) % 10) / 4.5 - 1.0;
    psiMinv_ref_0[i]  = psiMinv_data_[i];
    psiMinv_ref_1[i]  = psiMinv_data_2[i];
  }
  std::vector<const SPOSet::ValueType*> invRow_ptr_list{psiMinv_data_.data(), psiMinv_data_2.data()};

  // ratios_list
  std::vector<std::vector<SPOSet::ValueType>> ratios_list(nw);
  for (size_t iw = 0; iw < nw; iw++)
    ratios_list[iw].resize(nvp_list[iw]);

  // just need dummy refvec with correct size
  SPOSet::ValueVector tmp_psi_list(norb);
  spo->mw_evaluateDetRatios(spo_list, RefVectorWithLeader<const VirtualParticleSet>(VP_, {VP_, VP_2}),
                            RefVector<SPOSet::ValueVector>{tmp_psi_list}, invRow_ptr_list, ratios_list);

  std::vector<SPOSet::ValueType> ratios_ref_0(nvp_);
  std::vector<SPOSet::ValueType> ratios_ref_1(nvp_2);
  spo->evaluateDetRatios(VP_, tmp_psi_list, psiMinv_ref_0, ratios_ref_0);
  spo_2->evaluateDetRatios(VP_2, tmp_psi_list, psiMinv_ref_1, ratios_ref_1);
  for (int ivp = 0; ivp < nvp_; ivp++)
    CHECK(ratios_list[0][ivp] == Approx(ratios_ref_0[ivp]));
  for (int ivp = 0; ivp < nvp_2; ivp++)
    CHECK(ratios_list[1][ivp] == Approx(ratios_ref_1[ivp]));

  // app_log() << "ratios_list refvalues: \n" << std::setprecision(14);
  // for (int iw = 0; iw < nw; iw++)
  //   for (int ivp = 0; ivp < nvp_list[iw]; ivp++)
  //     app_log() << "CHECK(std::real(ratios_list[" << iw << "][" << ivp << "]) == Approx(" << ratios_list[iw][ivp]
  //               << "));\n";

  CHECK(std::real(ratios_list[0][0]) == Approx(-0.11554491049855));
  CHECK(std::real(ratios_list[0][1]) == Approx(0.19155774810121));
  CHECK(std::real(ratios_list[0][2]) == Approx(-0.16063724839636));
  CHECK(std::real(ratios_list[0][3]) == Approx(-0.37105113615831));
  CHECK(std::real(ratios_list[1][0]) == Approx(-0.12705469585615));
  CHECK(std::real(ratios_list[1][1]) == Approx(0.67930890998428));
  CHECK(std::real(ratios_list[1][2]) == Approx(0.83583922544552));
}
} // namespace qmcplusplus
