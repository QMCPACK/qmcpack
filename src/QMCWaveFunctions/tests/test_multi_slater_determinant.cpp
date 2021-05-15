//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2020 QMCPACK developers.
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
#include "QMCWaveFunctions/WaveFunctionFactory.h"

#include <stdio.h>
#include <string>
#include <limits>

using std::string;

namespace qmcplusplus
{
using PosType   = ParticleSet::PosType;
using ValueType = ParticleSet::ValueType;

void test_LiH_msd(const std::string& spo_xml_string,
                  const std::string& check_sponame,
                  int check_spo_size,
                  int check_basisset_size,
                  int test_batched)
{
  Communicate* c;
  c = OHMMS::Controller;

  auto ions_uptr = std::make_unique<ParticleSet>();
  auto elec_uptr = std::make_unique<ParticleSet>();
  ParticleSet& ions_(*ions_uptr);
  ParticleSet& elec_(*elec_uptr);

  ions_.setName("ion0");
  ions_.create({1, 1});
  ions_.R[0]           = {0.0, 0.0, 0.0};
  ions_.R[1]           = {0.0, 0.0, 3.0139239693};
  SpeciesSet& ispecies = ions_.getSpeciesSet();
  int LiIdx            = ispecies.addSpecies("Li");
  int HIdx             = ispecies.addSpecies("H");

  elec_.setName("elec");
  elec_.create({2, 2});
  elec_.R[0] = {0.5, 0.5, 0.5};
  elec_.R[1] = {0.1, 0.1, 1.1};
  elec_.R[2] = {-0.5, -0.5, -0.5};
  elec_.R[3] = {-0.1, -0.1, 1.5};

  SpeciesSet& tspecies       = elec_.getSpeciesSet();
  int upIdx                  = tspecies.addSpecies("u");
  int downIdx                = tspecies.addSpecies("d");
  int massIdx                = tspecies.addAttribute("mass");
  tspecies(massIdx, upIdx)   = 1.0;
  tspecies(massIdx, downIdx) = 1.0;

  // Need 1 electron and 1 proton, somehow
  //ParticleSet target = ParticleSet();
  ParticleSetPool ptcl = ParticleSetPool(c);
  ptcl.addParticleSet(std::move(elec_uptr));
  ptcl.addParticleSet(std::move(ions_uptr));

  Libxml2Document doc;
  bool okay = doc.parseFromString(spo_xml_string);
  REQUIRE(okay);

  xmlNodePtr ein_xml = doc.getRoot();

  WaveFunctionFactory wf_factory("psi0", elec_, ptcl.getPool(), c);
  wf_factory.put(ein_xml);

  SPOSet* spo_ptr(wf_factory.getSPOSet(check_sponame));
  REQUIRE(spo_ptr != nullptr);
  CHECK(spo_ptr->getOrbitalSetSize() == check_spo_size);
  CHECK(spo_ptr->getBasisSetSize() == check_basisset_size);

  ions_.update();
  elec_.update();

  auto& twf(*wf_factory.getTWF());
  twf.evaluateLog(elec_);

  std::cout << "twf.evaluateLog logpsi " << std::setprecision(16) << twf.getLogPsi() << " " << twf.getPhase()
            << std::endl;
  CHECK(std::complex<double>(twf.getLogPsi(), twf.getPhase()) ==
        LogComplexApprox(std::complex<double>(-7.646027846242066, 3.141592653589793)));

  twf.prepareGroup(elec_, 0);
  auto grad_old = twf.evalGrad(elec_, 1);
  std::cout << "twf.evalGrad grad_old " << std::setprecision(16) << grad_old << std::endl;
  CHECK(grad_old[0] == ValueApprox(0.1204183219));
  CHECK(grad_old[1] == ValueApprox(0.120821033));
  CHECK(grad_old[2] == ValueApprox(2.05904174));

  PosType delta(0.1, 0.1, 0.2);
  elec_.makeMove(1, delta);

  ParticleSet::GradType grad_new;
  auto ratio = twf.calcRatioGrad(elec_, 1, grad_new);
  std::cout << "twf.calcRatioGrad ratio " << ratio << " grad_new " << grad_new << std::endl;
  CHECK(ratio == ValueApprox(1.374307585));
  CHECK(grad_new[0] == ValueApprox(0.05732804333));
  CHECK(grad_new[1] == ValueApprox(0.05747775029));
  CHECK(grad_new[2] == ValueApprox(1.126889742));

  ratio = twf.calcRatio(elec_, 1);
  std::cout << "twf.calcRatio ratio " << ratio << std::endl;
  CHECK(ratio == ValueApprox(1.374307585));

  if (test_batched)
  {
    // set virtutal particle position
    PosType newpos(0.3, 0.2, 0.5);

    //elec_.makeVirtualMoves(newpos);
    //std::vector<ValueType> ratios(elec_.getTotalNum());
    //twf.evaluateRatiosAlltoOne(elec_, ratios);

    //CHECK(std::real(ratios[0]) == Approx());
    //CHECK(std::real(ratios[1]) == Approx());
    //CHECK(std::real(ratios[2]) == Approx());

    elec_.makeMove(0, newpos - elec_.R[0]);
    ValueType ratio_0 = twf.calcRatio(elec_, 0);
    elec_.rejectMove(0);

    CHECK(std::real(ratio_0) == Approx(2.350046921));

    VirtualParticleSet VP(elec_, 2);
    std::vector<PosType> newpos2(2);
    std::vector<ValueType> ratios2(2);
    newpos2[0] = newpos - elec_.R[1];
    newpos2[1] = PosType(0.2, 0.5, 0.3) - elec_.R[1];
    VP.makeMoves(1, elec_.R[1], newpos2);
    twf.evaluateRatios(VP, ratios2);

    CHECK(std::real(ratios2[0]) == Approx(-0.8544310407));
    CHECK(std::real(ratios2[1]) == Approx(-1.0830708458));

    //test acceptMove
    elec_.makeMove(1, newpos - elec_.R[1]);
    ValueType ratio_1 = twf.calcRatio(elec_, 1);
    twf.acceptMove(elec_, 1);
    elec_.acceptMove(1);

    CHECK(std::real(ratio_1) == Approx(-0.8544310407));
    CHECK(twf.getLogPsi() == Approx(-7.6460278462));
  }
}

TEST_CASE("LiH multi Slater dets", "[wavefunction]")
{
  app_log() << "-----------------------------------------------------------------" << std::endl;
  app_log() << "LiH_msd using the table method no precomputation" << std::endl;
  app_log() << "-----------------------------------------------------------------" << std::endl;
  const char* spo_xml_string1 = "<wavefunction name=\"psi0\" target=\"e\"> \
    <sposet_collection type=\"MolecularOrbital\" name=\"LCAOBSet\" source=\"ion0\" cuspCorrection=\"no\" href=\"LiH.orbs.h5\"> \
      <basisset name=\"LCAOBSet\" key=\"GTO\" transform=\"yes\"> \
        <grid type=\"log\" ri=\"1.e-6\" rf=\"1.e2\" npts=\"1001\"/> \
      </basisset> \
      <sposet basisset=\"LCAOBSet\" name=\"spo-up\" size=\"85\"> \
        <occupation mode=\"ground\"/> \
        <coefficient size=\"85\" spindataset=\"0\"/> \
      </sposet> \
      <sposet basisset=\"LCAOBSet\" name=\"spo-dn\" size=\"85\"> \
        <occupation mode=\"ground\"/> \
        <coefficient size=\"85\" spindataset=\"0\"/> \
      </sposet> \
    </sposet_collection> \
    <determinantset> \
      <multideterminant optimize=\"yes\" spo_up=\"spo-up\" spo_dn=\"spo-dn\" algorithm=\"table_method\"> \
        <detlist size=\"1487\" type=\"DETS\" cutoff=\"1e-20\" href=\"LiH.orbs.h5\"/> \
      </multideterminant> \
    </determinantset> \
</wavefunction> \
";
  test_LiH_msd(spo_xml_string1, "spo-up", 85, 105, true);

  app_log() << "-----------------------------------------------------------------" << std::endl;
  app_log() << "LiH_msd using the traditional slow method with all the determinants" << std::endl;
  app_log() << "-----------------------------------------------------------------" << std::endl;
  const char* spo_xml_string1_slow = "<wavefunction name=\"psi0\" target=\"e\"> \
    <sposet_collection type=\"MolecularOrbital\" name=\"LCAOBSet\" source=\"ion0\" cuspCorrection=\"no\" href=\"LiH.orbs.h5\"> \
      <basisset name=\"LCAOBSet\" key=\"GTO\" transform=\"yes\"> \
        <grid type=\"log\" ri=\"1.e-6\" rf=\"1.e2\" npts=\"1001\"/> \
      </basisset> \
      <sposet basisset=\"LCAOBSet\" name=\"spo-up\" size=\"85\"> \
        <occupation mode=\"ground\"/> \
        <coefficient size=\"85\" spindataset=\"0\"/> \
      </sposet> \
      <sposet basisset=\"LCAOBSet\" name=\"spo-dn\" size=\"85\"> \
        <occupation mode=\"ground\"/> \
        <coefficient size=\"85\" spindataset=\"0\"/> \
      </sposet> \
    </sposet_collection> \
    <determinantset> \
      <multideterminant optimize=\"yes\" spo_up=\"spo-up\" spo_dn=\"spo-dn\"  algorithm=\"all_determinants\"> \
        <detlist size=\"1487\" type=\"DETS\" cutoff=\"1e-20\" href=\"LiH.orbs.h5\"/> \
      </multideterminant> \
    </determinantset> \
</wavefunction> \
";
  test_LiH_msd(spo_xml_string1_slow, "spo-up", 85, 105, false);

  app_log() << "-----------------------------------------------------------------" << std::endl;
  app_log() << "LiH_msd using the table method with new optimization" << std::endl;
  app_log() << "-----------------------------------------------------------------" << std::endl;
  const char* spo_xml_string1_new = "<wavefunction name=\"psi0\" target=\"e\"> \
    <sposet_collection type=\"MolecularOrbital\" name=\"LCAOBSet\" source=\"ion0\" cuspCorrection=\"no\" href=\"LiH.orbs.h5\"> \
      <basisset name=\"LCAOBSet\" key=\"GTO\" transform=\"yes\"> \
        <grid type=\"log\" ri=\"1.e-6\" rf=\"1.e2\" npts=\"1001\"/> \
      </basisset> \
      <sposet basisset=\"LCAOBSet\" name=\"spo-up\" size=\"85\"> \
        <occupation mode=\"ground\"/> \
        <coefficient size=\"85\" spindataset=\"0\"/> \
      </sposet> \
      <sposet basisset=\"LCAOBSet\" name=\"spo-dn\" size=\"85\"> \
        <occupation mode=\"ground\"/> \
        <coefficient size=\"85\" spindataset=\"0\"/> \
      </sposet> \
    </sposet_collection> \
    <determinantset> \
      <multideterminant optimize=\"yes\" spo_up=\"spo-up\" spo_dn=\"spo-dn\" algorithm=\"precomputed_table_method\"> \
        <detlist size=\"1487\" type=\"DETS\" cutoff=\"1e-20\" href=\"LiH.orbs.h5\"/> \
      </multideterminant> \
    </determinantset> \
</wavefunction> \
";
  test_LiH_msd(spo_xml_string1_new, "spo-up", 85, 105, true);
}
} // namespace qmcplusplus
