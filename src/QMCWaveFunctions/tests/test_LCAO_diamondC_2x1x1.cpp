//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2023 QMCPACK developers.
//
// File developed by: Kevin Gasperich, kgasperich@anl.gov, Argonne National Laboratory
//
// File created by: Kevin Gasperich, kgasperich@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////


#include "catch.hpp"

#include "ParticleIO/LatticeIO.h"
#include "OhmmsData/Libxml2Doc.h"
#include "OhmmsPETE/OhmmsMatrix.h"
#include "Particle/ParticleSet.h"
#include "QMCWaveFunctions/WaveFunctionComponent.h"
#include "LCAO/LCAOrbitalBuilder.h"
#include <ResourceCollection.h>
#include "QMCHamiltonians/NLPPJob.h"
#include "DistanceTable.h"

#include <stdio.h>
#include <string>
#include <limits>

using std::string;

namespace qmcplusplus
{

void test_LCAO_DiamondC_2x1x1_real()
{
  using VT       = SPOSet::ValueType;
  Communicate* c = OHMMS::Controller;

  const char* particles = R"(<simulationcell>
     <parameter name="lattice" units="bohr">
        6.7463223   6.7463223   0.0
        0.0         3.37316115  3.37316115
        3.37316115  0.0         3.37316115
     </parameter>
     <parameter name="bconds">
        p p p
     </parameter>
  </simulationcell>
)";

  Libxml2Document doc_lattice;
  REQUIRE(doc_lattice.parseFromString(particles));

  // read lattice
  ParticleSet::ParticleLayout lattice;
  LatticeParser lp(lattice);
  lp.put(doc_lattice.getRoot());
  lattice.print(app_log(), 0);

  SimulationCell simcell(lattice);
  ParticleSet ions_(simcell);

  ions_.setName("ion0");
  ions_.create({4});
  ions_.R[0]           = {0.0, 0.0, 0.0};
  ions_.R[1]           = {1.686580575, 1.686580575, 1.686580575};
  ions_.R[2]           = {3.37316115, 3.37316115, 0.0};
  ions_.R[3]           = {5.059741726, 5.059741726, 1.686580575};
  SpeciesSet& ispecies = ions_.getSpeciesSet();
  const int Cidx       = ispecies.addSpecies("C");

  ions_.print(app_log());
  ions_.update(); // propagate SoA.

  ParticleSet elec_(simcell);
  elec_.setName("elec");
  elec_.create({8, 8});
  elec_.R[0] = {0.0, 1.0, 0.0};
  elec_.R[1] = {0.0, 1.1, 0.0};
  elec_.R[2] = {0.0, 1.2, 0.0};
  elec_.R[3] = {0.0, 1.3, 0.0};

  SpeciesSet& tspecies         = elec_.getSpeciesSet();
  const int upIdx              = tspecies.addSpecies("u");
  const int downIdx            = tspecies.addSpecies("d");
  const int chargeIdx          = tspecies.addAttribute("charge");
  tspecies(chargeIdx, upIdx)   = -1;
  tspecies(chargeIdx, downIdx) = -1;
  const int ei_table_index     = elec_.addTable(ions_);

  // diamondC_2x1x1
  // from tests/solids/diamondC_2x1x1-Gaussian_pp/C_Diamond-Twist0.wfj.xml
  const std::string wf_xml_str = R"(
    <sposet_collection type="molecularorbital" name="LCAOBSet" source="ion0" transform="yes" twist="0  0  0" href="C_Diamond_2x1x1-Gaussian.h5" PBCimages="5  5  5">
      <basisset name="LCAOBSet" key="GTO" transform="yes">
        <grid type="log" ri="1.e-6" rf="1.e2" npts="1001"/>
      </basisset>
      <sposet name="spoud" size="8">
        <occupation mode="ground"/>
        <coefficient size="116" spindataset="0"/>
      </sposet>
    </sposet_collection>
  )";
  Libxml2Document doc;
  bool okay = doc.parseFromString(wf_xml_str);
  REQUIRE(okay);

  xmlNodePtr root       = doc.getRoot();
  xmlNodePtr bset_xml   = xmlFirstElementChild(root);
  xmlNodePtr sposet_xml = xmlNextElementSibling(bset_xml);

  LCAOrbitalBuilder lcaoSet(elec_, ions_, c, root);
  auto spo = lcaoSet.createSPOSetFromXML(sposet_xml);
  REQUIRE(spo);
  auto& lcao_spos = dynamic_cast<const LCAOrbitalSet&>(*spo);
  CHECK(!lcao_spos.isIdentity());

  const int norb = spo->getOrbitalSetSize();
  app_log() << "norb: " << norb << std::endl;

  // test batched interfaces with 2 walkers
  const size_t nw = 2;

  ParticleSet elec_2(elec_);
  // interchange positions
  elec_2.R[0] = elec_.R[1];
  elec_2.R[1] = elec_.R[0];

  RefVectorWithLeader<ParticleSet> p_list(elec_, {elec_, elec_2});

  std::unique_ptr<SPOSet> spo_2(spo->makeClone());
  RefVectorWithLeader<SPOSet> spo_list(*spo, {*spo, *spo_2});

  ResourceCollection pset_res("test_pset_res");
  ResourceCollection spo_res("test_spo_res");

  elec_.createResource(pset_res);
  spo->createResource(spo_res);

  ResourceCollectionTeamLock<ParticleSet> mw_pset_lock(pset_res, p_list);
  ResourceCollectionTeamLock<SPOSet> mw_sposet_lock(spo_res, spo_list);
  ParticleSet::mw_update(p_list);

  SECTION("LCAOrbitalSet::mw_evaluateVGL")
  {
    SPOSet::ValueVector psiref_0(norb);
    SPOSet::GradVector dpsiref_0(norb);
    SPOSet::ValueVector d2psiref_0(norb);
    SPOSet::ValueVector psiref_1(norb);
    SPOSet::GradVector dpsiref_1(norb);
    SPOSet::ValueVector d2psiref_1(norb);

    spo_2->evaluateVGL(elec_2, 0, psiref_1, dpsiref_1, d2psiref_1);
    spo->evaluateVGL(elec_, 0, psiref_0, dpsiref_0, d2psiref_0);
    SECTION("single-walker VGL")
    {
      CHECK(Approx(std::real(psiref_0[0])) == 0.10203360788082);
      CHECK(Approx(std::real(d2psiref_0[0])) == -0.14269632514649);
      CHECK(Approx(std::real(dpsiref_0[0][0])) == 0.00031796000022);
      CHECK(Approx(std::real(dpsiref_0[0][1])) == -0.01951178212495);
      CHECK(Approx(std::real(dpsiref_0[0][2])) == -0.00031738324507);
      CHECK(Approx(std::real(psiref_0[1])) == 0.14111282090404);
      CHECK(Approx(std::real(d2psiref_0[1])) == -0.28230884342631);
      CHECK(Approx(std::real(dpsiref_0[1][0])) == 0.01205368790709);
      CHECK(Approx(std::real(dpsiref_0[1][1])) == -0.04876640907938);
      CHECK(Approx(std::real(dpsiref_0[1][2])) == -0.01205402562448);

      CHECK(Approx(std::real(psiref_1[0])) == 0.09954792826469);
      CHECK(Approx(std::real(d2psiref_1[0])) == -0.10799468581785);
      CHECK(Approx(std::real(dpsiref_1[0][0])) == 0.00024577684629);
      CHECK(Approx(std::real(dpsiref_1[0][1])) == -0.02943596305516);
      CHECK(Approx(std::real(dpsiref_1[0][2])) == -0.00024512723822);
      CHECK(Approx(std::real(psiref_1[1])) == 0.13558495733438);
      CHECK(Approx(std::real(d2psiref_1[1])) == -0.21885179829267);
      CHECK(Approx(std::real(dpsiref_1[1][0])) == 0.00738771300147);
      CHECK(Approx(std::real(dpsiref_1[1][1])) == -0.06080141726019);
      CHECK(Approx(std::real(dpsiref_1[1][2])) == -0.00738811343748);
    }

    SPOSet::ValueVector psi_v_1(norb);
    SPOSet::ValueVector psi_v_2(norb);
    RefVector<SPOSet::ValueVector> psi_v_list{psi_v_1, psi_v_2};
    spo->mw_evaluateValue(spo_list, p_list, 0, psi_v_list);
    SECTION("multi-walker V")
    {
      CHECK(Approx(std::real(psi_v_list[0].get()[0])) == 0.10203360788082);
      CHECK(Approx(std::real(psi_v_list[0].get()[1])) == 0.14111282090404);
      CHECK(Approx(std::real(psi_v_list[1].get()[0])) == 0.09954792826469);
      CHECK(Approx(std::real(psi_v_list[1].get()[1])) == 0.13558495733438);
    }


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
    SECTION("multi-walker VGL")
    {
      CHECK(Approx(std::real(psi_list[0].get()[0])) == 0.10203360788082);
      CHECK(Approx(std::real(d2psi_list[0].get()[0])) == -0.14269632514649);
      CHECK(Approx(std::real(dpsi_list[0].get()[0][0])) == 0.00031796000022);
      CHECK(Approx(std::real(dpsi_list[0].get()[0][1])) == -0.01951178212495);
      CHECK(Approx(std::real(dpsi_list[0].get()[0][2])) == -0.00031738324507);
      CHECK(Approx(std::real(psi_list[0].get()[1])) == 0.14111282090404);
      CHECK(Approx(std::real(d2psi_list[0].get()[1])) == -0.28230884342631);
      CHECK(Approx(std::real(dpsi_list[0].get()[1][0])) == 0.01205368790709);
      CHECK(Approx(std::real(dpsi_list[0].get()[1][1])) == -0.04876640907938);
      CHECK(Approx(std::real(dpsi_list[0].get()[1][2])) == -0.01205402562448);

      CHECK(Approx(std::real(psi_list[1].get()[0])) == 0.09954792826469);
      CHECK(Approx(std::real(d2psi_list[1].get()[0])) == -0.10799468581785);
      CHECK(Approx(std::real(dpsi_list[1].get()[0][0])) == 0.00024577684629);
      CHECK(Approx(std::real(dpsi_list[1].get()[0][1])) == -0.02943596305516);
      CHECK(Approx(std::real(dpsi_list[1].get()[0][2])) == -0.00024512723822);
      CHECK(Approx(std::real(psi_list[1].get()[1])) == 0.13558495733438);
      CHECK(Approx(std::real(d2psi_list[1].get()[1])) == -0.21885179829267);
      CHECK(Approx(std::real(dpsi_list[1].get()[1][0])) == 0.00738771300147);
      CHECK(Approx(std::real(dpsi_list[1].get()[1][1])) == -0.06080141726019);
      CHECK(Approx(std::real(dpsi_list[1].get()[1][2])) == -0.00738811343748);
    }


    SECTION("compare MW/SW V/VGL")
    {
      for (size_t iorb = 0; iorb < norb; iorb++)
      {
        CHECK(Approx(std::real(psi_list[0].get()[iorb])) == std::real(psiref_0[iorb]));
        CHECK(Approx(std::real(psi_list[1].get()[iorb])) == std::real(psiref_1[iorb]));
        CHECK(Approx(std::real(d2psi_list[0].get()[iorb])) == std::real(d2psiref_0[iorb]));
        CHECK(Approx(std::real(d2psi_list[1].get()[iorb])) == std::real(d2psiref_1[iorb]));
        for (size_t idim = 0; idim < SPOSet::DIM; idim++)
        {
          CHECK(Approx(std::real(dpsi_list[0].get()[iorb][idim])) == std::real(dpsiref_0[iorb][idim]));
          CHECK(Approx(std::real(dpsi_list[1].get()[iorb][idim])) == std::real(dpsiref_1[iorb][idim]));
        }
      }
      for (size_t iorb = 0; iorb < norb; iorb++)
        for (size_t iw = 0; iw < nw; iw++)
          CHECK(Approx(std::real(psi_v_list[iw].get()[iorb])) == std::real(psi_list[iw].get()[iorb]));
#ifdef QMC_COMPLEX
      for (size_t iorb = 0; iorb < norb; iorb++)
      {
        CHECK(Approx(std::imag(psi_list[0].get()[iorb])) == std::imag(psiref_0[iorb]));
        CHECK(Approx(std::imag(psi_list[1].get()[iorb])) == std::imag(psiref_1[iorb]));
        CHECK(Approx(std::imag(d2psi_list[0].get()[iorb])) == std::imag(d2psiref_0[iorb]));
        CHECK(Approx(std::imag(d2psi_list[1].get()[iorb])) == std::imag(d2psiref_1[iorb]));
        for (size_t idim = 0; idim < SPOSet::DIM; idim++)
        {
          CHECK(Approx(std::imag(dpsi_list[0].get()[iorb][idim])) == std::imag(dpsiref_0[iorb][idim]));
          CHECK(Approx(std::imag(dpsi_list[1].get()[iorb][idim])) == std::imag(dpsiref_1[iorb][idim]));
        }
      }
      for (size_t iorb = 0; iorb < norb; iorb++)
        for (size_t iw = 0; iw < nw; iw++)
          CHECK(Approx(std::imag(psi_v_list[iw].get()[iorb])) == std::imag(psi_list[iw].get()[iorb]));
#endif
    }
    // app_log() << "vgl_refvalues: \n" << std::setprecision(14);
    // for (int iorb = 0; iorb < 2; iorb++)
    // {
    //   app_log() << "CHECK(Approx(std::real(psiref_0[" << iorb << "])) == " << psiref_0[iorb] << ");\n";
    //   app_log() << "CHECK(Approx(std::real(d2psiref_0[" << iorb << "])) == " << d2psiref_0[iorb] << ");\n";
    //   for (int idim = 0; idim < 3; idim++)
    //     app_log() << "CHECK(Approx(std::real(dpsiref_0[" << iorb << "][" << idim << "])) == " << dpsiref_0[iorb][idim] << ");\n";
    // }
    // for (int iorb = 0; iorb < 2; iorb++)
    // {
    //   app_log() << "CHECK(Approx(std::real(psiref_1[" << iorb << "])) == " << psiref_1[iorb] << ");\n";
    //   app_log() << "CHECK(Approx(std::real(d2psiref_1[" << iorb << "])) == " << d2psiref_1[iorb] << "));\n";
    //   for (int idim = 0; idim < 3; idim++)
    //     app_log() << "CHECK(Approx(std::real(dpsiref_1[" << iorb << "][" << idim << "])) == " << dpsiref_1[iorb][idim] << ");\n";
    // }
  }
  SECTION("LCAOrbitalSet::mw_evaluateDetRatios")
  {
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
    // make virtual move of elec 1, reference ion 3
    NLPPJob<SPOSet::RealType> job_2(3, 1, ei_table_2.getDistances()[1][3], -ei_table_2.getDisplacements()[1][3]);

    // make VP refvec
    RefVectorWithLeader<VirtualParticleSet> vp_list(VP_, {VP_, VP_2});

    ResourceCollection vp_res("test_vp_res");
    VP_.createResource(vp_res);
    ResourceCollectionTeamLock<VirtualParticleSet> mw_vpset_lock(vp_res, vp_list);
    VirtualParticleSet::mw_makeMoves(vp_list, p_list, {newpos_vp_, newpos_vp_2}, {job_, job_2}, false);

    // fill invrow with dummy data for each walker
    std::vector<SPOSet::ValueType> psiMinv_data_(norb);
    std::vector<SPOSet::ValueType> psiMinv_data_2(norb);
    SPOSet::ValueVector psiMinv_ref_0(norb);
    SPOSet::ValueVector psiMinv_ref_1(norb);
    for (int i = 0; i < norb; i++)
    {
      psiMinv_data_[i]  = 0.1 * i;
      psiMinv_data_2[i] = 0.2 * i;
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
    // single-walker functions for reference
    spo->evaluateDetRatios(VP_, tmp_psi_list, psiMinv_ref_0, ratios_ref_0);
    spo_2->evaluateDetRatios(VP_2, tmp_psi_list, psiMinv_ref_1, ratios_ref_1);

    for (int ivp = 0; ivp < nvp_; ivp++)
      CHECK(Approx(std::real(ratios_list[0][ivp])) == std::real(ratios_ref_0[ivp]));
    for (int ivp = 0; ivp < nvp_2; ivp++)
      CHECK(Approx(std::real(ratios_list[1][ivp])) == std::real(ratios_ref_1[ivp]));
#ifdef QMC_COMPLEX
    for (int ivp = 0; ivp < nvp_; ivp++)
    {
      CHECK(Approx(std::imag(ratios_list[0][ivp])) == std::imag(ratios_ref_0[ivp]));
      CHECK(Approx(std::imag(ratios_list[0][ivp])) == 0.0);
    }
    for (int ivp = 0; ivp < nvp_2; ivp++)
    {
      CHECK(Approx(std::imag(ratios_list[1][ivp])) == std::imag(ratios_ref_1[ivp]));
      CHECK(Approx(std::imag(ratios_list[1][ivp])) == 0.0);
    }
#endif

    CHECK(Approx(std::real(ratios_list[0][0])) == 0.19309684969511);
    CHECK(Approx(std::real(ratios_list[0][1])) == 0.19743141486366);
    CHECK(Approx(std::real(ratios_list[0][2])) == 0.17884881050205);
    CHECK(Approx(std::real(ratios_list[0][3])) == 0.15105783567230);
    CHECK(Approx(std::real(ratios_list[1][0])) == 0.38619369939021);
    CHECK(Approx(std::real(ratios_list[1][1])) == 0.38429955941922);
    CHECK(Approx(std::real(ratios_list[1][2])) == 0.32071997896196);

    CHECK(Approx(std::real(ratios_ref_0[0])) == 0.1930968497);
    CHECK(Approx(std::real(ratios_ref_0[1])) == 0.1974314149);
    CHECK(Approx(std::real(ratios_ref_0[2])) == 0.1788488105);
    CHECK(Approx(std::real(ratios_ref_0[3])) == 0.1510578357);
    CHECK(Approx(std::real(ratios_ref_1[0])) == 0.3861936994);
    CHECK(Approx(std::real(ratios_ref_1[1])) == 0.3842995594);
    CHECK(Approx(std::real(ratios_ref_1[2])) == 0.3207199790);

    //// print SW ref values
    // for (int ivp = 0; ivp < nvp_; ivp++)
    //   app_log() << "CHECK(Approx(std::real(ratios_ref_0[" << ivp << "])) == " << std::real(ratios_ref_0[ivp]) << ");\n";
    // for (int ivp = 0; ivp < nvp_2; ivp++)
    //   app_log() << "CHECK(Approx(std::real(ratios_ref_1[" << ivp << "])) == " << std::real(ratios_ref_1[ivp]) << ");\n";

    //// print MW ref values
    // for (int iw = 0; iw < nw; iw++)
    //   for (int ivp = 0; ivp < nvp_list[iw]; ivp++)
    //     app_log() << "CHECK(Approx(std::real(ratios_list[" << iw << "][" << ivp << "])) == " << std::real(ratios_list[iw][ivp]) << ");\n";
    // for (int iw = 0; iw < nw; iw++)
    //   for (int ivp = 0; ivp < nvp_list[iw]; ivp++)
    //     app_log() << "CHECK(Approx(std::imag(ratios_list[" << iw << "][" << ivp << "])) == " << std::imag(ratios_list[iw][ivp]) << ");\n";
  }
}

void test_LCAO_DiamondC_2x1x1_cplx()
{
  using VT       = SPOSet::ValueType;
  Communicate* c = OHMMS::Controller;

  const char* particles = R"(<simulationcell>
     <parameter name="lattice" units="bohr">
        6.7463223   6.7463223   0.0
        0.0         3.37316115  3.37316115
        3.37316115  0.0         3.37316115
     </parameter>
     <parameter name="bconds">
        p p p
     </parameter>
  </simulationcell>
)";

  Libxml2Document doc_lattice;
  REQUIRE(doc_lattice.parseFromString(particles));

  // read lattice
  ParticleSet::ParticleLayout lattice;
  LatticeParser lp(lattice);
  lp.put(doc_lattice.getRoot());
  lattice.print(app_log(), 0);

  SimulationCell simcell(lattice);
  ParticleSet ions_(simcell);

  ions_.setName("ion0");
  ions_.create({4});
  ions_.R[0]           = {0.0, 0.0, 0.0};
  ions_.R[1]           = {1.686580575, 1.686580575, 1.686580575};
  ions_.R[2]           = {3.37316115, 3.37316115, 0.0};
  ions_.R[3]           = {5.059741726, 5.059741726, 1.686580575};
  SpeciesSet& ispecies = ions_.getSpeciesSet();
  const int Cidx       = ispecies.addSpecies("C");

  ions_.print(app_log());
  ions_.update(); // propagate SoA.

  ParticleSet elec_(simcell);
  elec_.setName("elec");
  elec_.create({8, 8});
  elec_.R[0] = {0.0, 1.0, 0.0};
  elec_.R[1] = {0.0, 1.1, 0.0};
  elec_.R[2] = {0.0, 1.2, 0.0};
  elec_.R[3] = {0.0, 1.3, 0.0};

  SpeciesSet& tspecies         = elec_.getSpeciesSet();
  const int upIdx              = tspecies.addSpecies("u");
  const int downIdx            = tspecies.addSpecies("d");
  const int chargeIdx          = tspecies.addAttribute("charge");
  tspecies(chargeIdx, upIdx)   = -1;
  tspecies(chargeIdx, downIdx) = -1;
  const int ei_table_index     = elec_.addTable(ions_);

  // diamondC_2x1x1
  // from tests/solids/diamondC_2x1x1-Gaussian_pp_cplx/C_Diamond-tiled-cplx.wfj.xml
  const std::string wf_xml_str = R"(
    <sposet_collection type="molecularorbital" name="LCAOBSet" source="ion0" transform="yes" twist="0.07761248  0.07761248  -0.07761248" href="C_Diamond_2x1x1-Gaussian-tiled-cplx.h5" PBCimages="5  5  5">
      <basisset name="LCAOBSet" key="GTO" transform="yes">
        <grid type="log" ri="1.e-6" rf="1.e2" npts="1001"/>
      </basisset>
      <sposet name="spoud" size="8" >
        <occupation mode="ground"/>
        <coefficient size="52" spindataset="0"/>
      </sposet>
    </sposet_collection>
  )";
  Libxml2Document doc;
  bool okay = doc.parseFromString(wf_xml_str);
  REQUIRE(okay);

  xmlNodePtr root       = doc.getRoot();
  xmlNodePtr bset_xml   = xmlFirstElementChild(root);
  xmlNodePtr sposet_xml = xmlNextElementSibling(bset_xml);

  LCAOrbitalBuilder lcaoSet(elec_, ions_, c, root);
  auto spo = lcaoSet.createSPOSetFromXML(sposet_xml);
  REQUIRE(spo);
  auto& lcao_spos = dynamic_cast<const LCAOrbitalSet&>(*spo);
  CHECK(!lcao_spos.isIdentity());

  const int norb = spo->getOrbitalSetSize();
  app_log() << "norb: " << norb << std::endl;

  // test batched interfaces with 2 walkers
  const size_t nw = 2;

  ParticleSet elec_2(elec_);
  // interchange positions
  elec_2.R[0] = elec_.R[1];
  elec_2.R[1] = elec_.R[0];

  RefVectorWithLeader<ParticleSet> p_list(elec_, {elec_, elec_2});

  std::unique_ptr<SPOSet> spo_2(spo->makeClone());
  RefVectorWithLeader<SPOSet> spo_list(*spo, {*spo, *spo_2});

  ResourceCollection pset_res("test_pset_res");
  ResourceCollection spo_res("test_spo_res");

  elec_.createResource(pset_res);
  spo->createResource(spo_res);

  ResourceCollectionTeamLock<ParticleSet> mw_pset_lock(pset_res, p_list);
  ResourceCollectionTeamLock<SPOSet> mw_sposet_lock(spo_res, spo_list);
  ParticleSet::mw_update(p_list);

  SECTION("LCAOrbitalSet::mw_evaluateVGL")
  {
    SPOSet::ValueVector psiref_0(norb);
    SPOSet::GradVector dpsiref_0(norb);
    SPOSet::ValueVector d2psiref_0(norb);
    SPOSet::ValueVector psiref_1(norb);
    SPOSet::GradVector dpsiref_1(norb);
    SPOSet::ValueVector d2psiref_1(norb);

    spo_2->evaluateVGL(elec_2, 0, psiref_1, dpsiref_1, d2psiref_1);
    spo->evaluateVGL(elec_, 0, psiref_0, dpsiref_0, d2psiref_0);
    SECTION("single-walker VGL")
    {
      CHECK(Approx(std::real(psiref_0[0])) == 0.10394953298732);
      CHECK(Approx(std::real(d2psiref_0[0])) == -0.15355833615767);
      CHECK(Approx(std::real(dpsiref_0[0][0])) == -0.00066360565262);
      CHECK(Approx(std::real(dpsiref_0[0][1])) == -0.03040447635885);
      CHECK(Approx(std::real(dpsiref_0[0][2])) == 0.00066400579544);
      CHECK(Approx(std::real(psiref_0[1])) == 0.04007076321982);
      CHECK(Approx(std::real(d2psiref_0[1])) == -0.13347762878346);
      CHECK(Approx(std::real(dpsiref_0[1][0])) == 0.06427121424074);
      CHECK(Approx(std::real(dpsiref_0[1][1])) == -0.03145499248870);
      CHECK(Approx(std::real(dpsiref_0[1][2])) == 0.09703977123104);

      CHECK(Approx(std::real(psiref_1[0])) == 0.10041826765555);
      CHECK(Approx(std::real(d2psiref_1[0])) == -0.11271657756448);
      CHECK(Approx(std::real(dpsiref_1[0][0])) == -0.00056949645372);
      CHECK(Approx(std::real(dpsiref_1[0][1])) == -0.03940982127946);
      CHECK(Approx(std::real(dpsiref_1[0][2])) == 0.00056995118020);
      CHECK(Approx(std::real(psiref_1[1])) == 0.03692120161253);
      CHECK(Approx(std::real(d2psiref_1[1])) == -0.10225156633271);
      CHECK(Approx(std::real(dpsiref_1[1][0])) == 0.05862729651642);
      CHECK(Approx(std::real(dpsiref_1[1][1])) == -0.03141220770622);
      CHECK(Approx(std::real(dpsiref_1[1][2])) == 0.08454924955444);

      CHECK(Approx(std::imag(psiref_0[0])) == 0.00920567821605);
      CHECK(Approx(std::imag(d2psiref_0[0])) == -0.01897980702117);
      CHECK(Approx(std::imag(dpsiref_0[0][0])) == 0.00599757336742);
      CHECK(Approx(std::imag(dpsiref_0[0][1])) == 0.00009471083794);
      CHECK(Approx(std::imag(dpsiref_0[0][2])) == -0.00599734141981);
      CHECK(Approx(std::imag(psiref_0[1])) == -0.07297020296665);
      CHECK(Approx(std::imag(d2psiref_0[1])) == 0.24972839540365);
      CHECK(Approx(std::imag(dpsiref_0[1][0])) == -0.14785027769482);
      CHECK(Approx(std::imag(dpsiref_0[1][1])) == 0.05546078377043);
      CHECK(Approx(std::imag(dpsiref_0[1][2])) == -0.19276310593334);

      CHECK(Approx(std::imag(psiref_1[0])) == 0.00921878891696);
      CHECK(Approx(std::imag(d2psiref_1[0])) == -0.01501853597019);
      CHECK(Approx(std::imag(dpsiref_1[0][0])) == 0.00531194469490);
      CHECK(Approx(std::imag(dpsiref_1[0][1])) == 0.00017539927704);
      CHECK(Approx(std::imag(dpsiref_1[0][2])) == -0.00531168663225);
      CHECK(Approx(std::imag(psiref_1[1])) == -0.06736883166719);
      CHECK(Approx(std::imag(d2psiref_1[1])) == 0.19452200419216);
      CHECK(Approx(std::imag(dpsiref_1[1][0])) == -0.13949226666533);
      CHECK(Approx(std::imag(dpsiref_1[1][1])) == 0.05635953434634);
      CHECK(Approx(std::imag(dpsiref_1[1][2])) == -0.17227552311686);
    }

    SPOSet::ValueVector psi_v_1(norb);
    SPOSet::ValueVector psi_v_2(norb);
    RefVector<SPOSet::ValueVector> psi_v_list{psi_v_1, psi_v_2};
    spo->mw_evaluateValue(spo_list, p_list, 0, psi_v_list);
    SECTION("multi-walker V")
    {
      CHECK(Approx(std::real(psi_v_list[0].get()[0])) == 0.10394953298732);
      CHECK(Approx(std::real(psi_v_list[0].get()[1])) == 0.04007076321982);
      CHECK(Approx(std::real(psi_v_list[1].get()[0])) == 0.10041826765555);
      CHECK(Approx(std::real(psi_v_list[1].get()[1])) == 0.03692120161253);

      CHECK(Approx(std::imag(psi_v_list[0].get()[0])) == 0.00920567821605);
      CHECK(Approx(std::imag(psi_v_list[0].get()[1])) == -0.07297020296665);
      CHECK(Approx(std::imag(psi_v_list[1].get()[0])) == 0.00921878891696);
      CHECK(Approx(std::imag(psi_v_list[1].get()[1])) == -0.06736883166719);
    }


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
    SECTION("multi-walker VGL")
    {
      CHECK(Approx(std::real(psi_list[0].get()[0])) == 0.10394953298732);
      CHECK(Approx(std::real(d2psi_list[0].get()[0])) == -0.15355833615767);
      CHECK(Approx(std::real(dpsi_list[0].get()[0][0])) == -0.00066360565262);
      CHECK(Approx(std::real(dpsi_list[0].get()[0][1])) == -0.03040447635885);
      CHECK(Approx(std::real(dpsi_list[0].get()[0][2])) == 0.00066400579544);
      CHECK(Approx(std::real(psi_list[0].get()[1])) == 0.04007076321982);
      CHECK(Approx(std::real(d2psi_list[0].get()[1])) == -0.13347762878346);
      CHECK(Approx(std::real(dpsi_list[0].get()[1][0])) == 0.06427121424074);
      CHECK(Approx(std::real(dpsi_list[0].get()[1][1])) == -0.03145499248870);
      CHECK(Approx(std::real(dpsi_list[0].get()[1][2])) == 0.09703977123104);

      CHECK(Approx(std::real(psi_list[1].get()[0])) == 0.10041826765555);
      CHECK(Approx(std::real(d2psi_list[1].get()[0])) == -0.11271657756448);
      CHECK(Approx(std::real(dpsi_list[1].get()[0][0])) == -0.00056949645372);
      CHECK(Approx(std::real(dpsi_list[1].get()[0][1])) == -0.03940982127946);
      CHECK(Approx(std::real(dpsi_list[1].get()[0][2])) == 0.00056995118020);
      CHECK(Approx(std::real(psi_list[1].get()[1])) == 0.03692120161253);
      CHECK(Approx(std::real(d2psi_list[1].get()[1])) == -0.10225156633271);
      CHECK(Approx(std::real(dpsi_list[1].get()[1][0])) == 0.05862729651642);
      CHECK(Approx(std::real(dpsi_list[1].get()[1][1])) == -0.03141220770622);
      CHECK(Approx(std::real(dpsi_list[1].get()[1][2])) == 0.08454924955444);

      CHECK(Approx(std::imag(psi_list[0].get()[0])) == 0.00920567821605);
      CHECK(Approx(std::imag(d2psi_list[0].get()[0])) == -0.01897980702117);
      CHECK(Approx(std::imag(dpsi_list[0].get()[0][0])) == 0.00599757336742);
      CHECK(Approx(std::imag(dpsi_list[0].get()[0][1])) == 0.00009471083794);
      CHECK(Approx(std::imag(dpsi_list[0].get()[0][2])) == -0.00599734141981);
      CHECK(Approx(std::imag(psi_list[0].get()[1])) == -0.07297020296665);
      CHECK(Approx(std::imag(d2psi_list[0].get()[1])) == 0.24972839540365);
      CHECK(Approx(std::imag(dpsi_list[0].get()[1][0])) == -0.14785027769482);
      CHECK(Approx(std::imag(dpsi_list[0].get()[1][1])) == 0.05546078377043);
      CHECK(Approx(std::imag(dpsi_list[0].get()[1][2])) == -0.19276310593334);

      CHECK(Approx(std::imag(psi_list[1].get()[0])) == 0.00921878891696);
      CHECK(Approx(std::imag(d2psi_list[1].get()[0])) == -0.01501853597019);
      CHECK(Approx(std::imag(dpsi_list[1].get()[0][0])) == 0.00531194469490);
      CHECK(Approx(std::imag(dpsi_list[1].get()[0][1])) == 0.00017539927704);
      CHECK(Approx(std::imag(dpsi_list[1].get()[0][2])) == -0.00531168663225);
      CHECK(Approx(std::imag(psi_list[1].get()[1])) == -0.06736883166719);
      CHECK(Approx(std::imag(d2psi_list[1].get()[1])) == 0.19452200419216);
      CHECK(Approx(std::imag(dpsi_list[1].get()[1][0])) == -0.13949226666533);
      CHECK(Approx(std::imag(dpsi_list[1].get()[1][1])) == 0.05635953434634);
      CHECK(Approx(std::imag(dpsi_list[1].get()[1][2])) == -0.17227552311686);
    }


    SECTION("compare MW/SW V/VGL")
    {
      //real
      for (size_t iorb = 0; iorb < norb; iorb++)
      {
        CHECK(Approx(std::real(psi_list[0].get()[iorb])) == std::real(psiref_0[iorb]));
        CHECK(Approx(std::real(psi_list[1].get()[iorb])) == std::real(psiref_1[iorb]));
        CHECK(Approx(std::real(d2psi_list[0].get()[iorb])) == std::real(d2psiref_0[iorb]));
        CHECK(Approx(std::real(d2psi_list[1].get()[iorb])) == std::real(d2psiref_1[iorb]));
        for (size_t idim = 0; idim < SPOSet::DIM; idim++)
        {
          CHECK(Approx(std::real(dpsi_list[0].get()[iorb][idim])) == std::real(dpsiref_0[iorb][idim]));
          CHECK(Approx(std::real(dpsi_list[1].get()[iorb][idim])) == std::real(dpsiref_1[iorb][idim]));
        }
      }
      for (size_t iorb = 0; iorb < norb; iorb++)
        for (size_t iw = 0; iw < nw; iw++)
          CHECK(Approx(std::real(psi_v_list[iw].get()[iorb])) == std::real(psi_list[iw].get()[iorb]));
      //imag
      for (size_t iorb = 0; iorb < norb; iorb++)
      {
        CHECK(Approx(std::imag(psi_list[0].get()[iorb])) == std::imag(psiref_0[iorb]));
        CHECK(Approx(std::imag(psi_list[1].get()[iorb])) == std::imag(psiref_1[iorb]));
        CHECK(Approx(std::imag(d2psi_list[0].get()[iorb])) == std::imag(d2psiref_0[iorb]));
        CHECK(Approx(std::imag(d2psi_list[1].get()[iorb])) == std::imag(d2psiref_1[iorb]));
        for (size_t idim = 0; idim < SPOSet::DIM; idim++)
        {
          CHECK(Approx(std::imag(dpsi_list[0].get()[iorb][idim])) == std::imag(dpsiref_0[iorb][idim]));
          CHECK(Approx(std::imag(dpsi_list[1].get()[iorb][idim])) == std::imag(dpsiref_1[iorb][idim]));
        }
      }
      for (size_t iorb = 0; iorb < norb; iorb++)
        for (size_t iw = 0; iw < nw; iw++)
          CHECK(Approx(std::imag(psi_v_list[iw].get()[iorb])) == std::imag(psi_list[iw].get()[iorb]));
    }
    // app_log() << "vgl_refvalues: \n" << std::setprecision(14);
    // for (int iorb = 0; iorb < 2; iorb++)
    // {
    //   app_log() << "CHECK(Approx(std::real(psiref_0[" << iorb << "])) == " << std::real(psiref_0[iorb]) << ");\n";
    //   app_log() << "CHECK(Approx(std::real(d2psiref_0[" << iorb << "])) == " << std::real(d2psiref_0[iorb]) << ");\n";
    //   for (int idim = 0; idim < 3; idim++)
    //     app_log() << "CHECK(Approx(std::real(dpsiref_0[" << iorb << "][" << idim
    //               << "])) == " << std::real(dpsiref_0[iorb][idim]) << ");\n";
    // }
    // for (int iorb = 0; iorb < 2; iorb++)
    // {
    //   app_log() << "CHECK(Approx(std::real(psiref_1[" << iorb << "])) == " << std::real(psiref_1[iorb]) << ");\n";
    //   app_log() << "CHECK(Approx(std::real(d2psiref_1[" << iorb << "])) == " << std::real(d2psiref_1[iorb]) << "));\n";
    //   for (int idim = 0; idim < 3; idim++)
    //     app_log() << "CHECK(Approx(std::real(dpsiref_1[" << iorb << "][" << idim
    //               << "])) == " << std::real(dpsiref_1[iorb][idim]) << ");\n";
    // }
    // app_log() << "vgl_refvalues: \n" << std::setprecision(14);
    // for (int iorb = 0; iorb < 2; iorb++)
    // {
    //   app_log() << "CHECK(Approx(std::imag(psiref_0[" << iorb << "])) == " << std::imag(psiref_0[iorb]) << ");\n";
    //   app_log() << "CHECK(Approx(std::imag(d2psiref_0[" << iorb << "])) == " << std::imag(d2psiref_0[iorb]) << ");\n";
    //   for (int idim = 0; idim < 3; idim++)
    //     app_log() << "CHECK(Approx(std::imag(dpsiref_0[" << iorb << "][" << idim
    //               << "])) == " << std::imag(dpsiref_0[iorb][idim]) << ");\n";
    // }
    // for (int iorb = 0; iorb < 2; iorb++)
    // {
    //   app_log() << "CHECK(Approx(std::imag(psiref_1[" << iorb << "])) == " << std::imag(psiref_1[iorb]) << ");\n";
    //   app_log() << "CHECK(Approx(std::imag(d2psiref_1[" << iorb << "])) == " << std::imag(d2psiref_1[iorb]) << "));\n";
    //   for (int idim = 0; idim < 3; idim++)
    //     app_log() << "CHECK(Approx(std::imag(dpsiref_1[" << iorb << "][" << idim
    //               << "])) == " << std::imag(dpsiref_1[iorb][idim]) << ");" << std::endl;
    // }
  }
  SECTION("LCAOrbitalSet::mw_evaluateDetRatios")
  {
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
    // make virtual move of elec 1, reference ion 3
    NLPPJob<SPOSet::RealType> job_2(3, 1, ei_table_2.getDistances()[1][3], -ei_table_2.getDisplacements()[1][3]);

    // make VP refvec
    RefVectorWithLeader<VirtualParticleSet> vp_list(VP_, {VP_, VP_2});

    ResourceCollection vp_res("test_vp_res");
    VP_.createResource(vp_res);
    ResourceCollectionTeamLock<VirtualParticleSet> mw_vpset_lock(vp_res, vp_list);
    VirtualParticleSet::mw_makeMoves(vp_list, p_list, {newpos_vp_, newpos_vp_2}, {job_, job_2}, false);

    // fill invrow with dummy data for each walker
    std::vector<SPOSet::ValueType> psiMinv_data_(norb);
    std::vector<SPOSet::ValueType> psiMinv_data_2(norb);
    SPOSet::ValueVector psiMinv_ref_0(norb);
    SPOSet::ValueVector psiMinv_ref_1(norb);
    for (int i = 0; i < norb; i++)
    {
      psiMinv_data_[i]  = 0.1 * i;
      psiMinv_data_2[i] = 0.2 * i;
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
    // single-walker functions for reference
    spo->evaluateDetRatios(VP_, tmp_psi_list, psiMinv_ref_0, ratios_ref_0);
    spo_2->evaluateDetRatios(VP_2, tmp_psi_list, psiMinv_ref_1, ratios_ref_1);

    for (int ivp = 0; ivp < nvp_; ivp++)
      CHECK(Approx(std::real(ratios_list[0][ivp])) == std::real(ratios_ref_0[ivp]));
    for (int ivp = 0; ivp < nvp_2; ivp++)
      CHECK(Approx(std::real(ratios_list[1][ivp])) == std::real(ratios_ref_1[ivp]));
    for (int ivp = 0; ivp < nvp_; ivp++)
      CHECK(Approx(std::imag(ratios_list[0][ivp])) == std::imag(ratios_ref_0[ivp]));
    for (int ivp = 0; ivp < nvp_2; ivp++)
      CHECK(Approx(std::imag(ratios_list[1][ivp])) == std::imag(ratios_ref_1[ivp]));

    CHECK(Approx(std::real(ratios_ref_0[0])) == 0.0963445284);
    CHECK(Approx(std::real(ratios_ref_0[1])) == 0.0784621772);
    CHECK(Approx(std::real(ratios_ref_0[2])) == 0.0312479567);
    CHECK(Approx(std::real(ratios_ref_0[3])) == -0.0240189529);
    CHECK(Approx(std::real(ratios_ref_1[0])) == 0.1926890568);
    CHECK(Approx(std::real(ratios_ref_1[1])) == 0.1037508495);
    CHECK(Approx(std::real(ratios_ref_1[2])) == -0.0747915097);

    CHECK(Approx(std::real(ratios_list[0][0])) == 0.0963445284);
    CHECK(Approx(std::real(ratios_list[0][1])) == 0.0784621772);
    CHECK(Approx(std::real(ratios_list[0][2])) == 0.0312479567);
    CHECK(Approx(std::real(ratios_list[0][3])) == -0.0240189529);
    CHECK(Approx(std::real(ratios_list[1][0])) == 0.1926890568);
    CHECK(Approx(std::real(ratios_list[1][1])) == 0.1037508495);
    CHECK(Approx(std::real(ratios_list[1][2])) == -0.0747915097);

    CHECK(Approx(std::imag(ratios_ref_0[0])) == -0.0090812301);
    CHECK(Approx(std::imag(ratios_ref_0[1])) == -0.0385825385);
    CHECK(Approx(std::imag(ratios_ref_0[2])) == -0.0610830209);
    CHECK(Approx(std::imag(ratios_ref_0[3])) == -0.0809775403);
    CHECK(Approx(std::imag(ratios_ref_1[0])) == -0.0181624602);
    CHECK(Approx(std::imag(ratios_ref_1[1])) == -0.0856868673);
    CHECK(Approx(std::imag(ratios_ref_1[2])) == -0.1487774316);

    CHECK(Approx(std::imag(ratios_list[0][0])) == -0.0090812301);
    CHECK(Approx(std::imag(ratios_list[0][1])) == -0.0385825385);
    CHECK(Approx(std::imag(ratios_list[0][2])) == -0.0610830209);
    CHECK(Approx(std::imag(ratios_list[0][3])) == -0.0809775403);
    CHECK(Approx(std::imag(ratios_list[1][0])) == -0.0181624602);
    CHECK(Approx(std::imag(ratios_list[1][1])) == -0.0856868673);
    CHECK(Approx(std::imag(ratios_list[1][2])) == -0.1487774316);

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


TEST_CASE("LCAOrbitalSet batched PBC DiamondC", "[wavefunction]")
{
  SECTION("2x1x1 real") { test_LCAO_DiamondC_2x1x1_real(); }
#ifdef QMC_COMPLEX
  SECTION("2x1x1 cplx") { test_LCAO_DiamondC_2x1x1_cplx(); }
#endif
}
} // namespace qmcplusplus
