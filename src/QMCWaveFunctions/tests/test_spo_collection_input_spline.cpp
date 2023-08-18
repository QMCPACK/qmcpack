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
#include "Utilities/RuntimeOptions.h"

#include <stdio.h>
#include <string>
#include <limits>

using std::string;

namespace qmcplusplus
{
void test_diamond_2x1x1_xml_input(const std::string& spo_xml_string)
{
  Communicate* c = OHMMS::Controller;

  // diamondC_2x1x1
  ParticleSet::ParticleLayout lattice;
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
  ions_.R[0] = {0.0, 0.0, 0.0};
  ions_.R[1] = {1.68658058, 1.68658058, 1.68658058};
  ions_.R[2] = {3.37316115, 3.37316115, 0.0};
  ions_.R[3] = {5.05974173, 5.05974173, 1.68658058};
  elec_.setName("elec");
  ptcl.addParticleSet(std::move(elec_uptr));
  elec_.create({2});
  elec_.R[0]                 = {0.0, 0.0, 0.0};
  elec_.R[1]                 = {0.0, 1.0, 0.0};
  SpeciesSet& tspecies       = elec_.getSpeciesSet();
  int upIdx                  = tspecies.addSpecies("u");
  int chargeIdx              = tspecies.addAttribute("charge");
  tspecies(chargeIdx, upIdx) = -1;

  Libxml2Document doc;
  bool okay = doc.parseFromString(spo_xml_string);
  REQUIRE(okay);

  xmlNodePtr ein_xml = doc.getRoot();

  WaveFunctionFactory wf_factory(elec_, ptcl.getPool(), c);
  RuntimeOptions runtime_options;
  auto twf_ptr = wf_factory.buildTWF(ein_xml, runtime_options);

  std::unique_ptr<SPOSet> spo(twf_ptr->getSPOSet("spo").makeClone());

  // for vgl
  SPOSet::ValueMatrix psiM(elec_.R.size(), spo->getOrbitalSetSize());
  SPOSet::GradMatrix dpsiM(elec_.R.size(), spo->getOrbitalSetSize());
  SPOSet::ValueMatrix d2psiM(elec_.R.size(), spo->getOrbitalSetSize());
  spo->evaluate_notranspose(elec_, 0, elec_.R.size(), psiM, dpsiM, d2psiM);

  // real part
  // due to the different ordering of bands skip the tests on CUDA+Real builds
  // checking evaluations, reference values are not independently generated.
  // value
  CHECK(std::real(psiM[1][0]) == Approx(0.9008999467));
  CHECK(std::real(psiM[1][1]) == Approx(1.2383049726));
  // grad
  CHECK(std::real(dpsiM[1][0][0]) == Approx(0.0025820041));
  CHECK(std::real(dpsiM[1][0][1]) == Approx(-0.1880052537));
  CHECK(std::real(dpsiM[1][0][2]) == Approx(-0.0025404284));
  CHECK(std::real(dpsiM[1][1][0]) == Approx(0.1069662273));
  CHECK(std::real(dpsiM[1][1][1]) == Approx(-0.4364597797));
  CHECK(std::real(dpsiM[1][1][2]) == Approx(-0.106951952));
  // lapl
  CHECK(std::real(d2psiM[1][0]) == Approx(-1.3757134676));
  CHECK(std::real(d2psiM[1][1]) == Approx(-2.4803137779));

#if defined(QMC_COMPLEX)
  // imaginary part
  // value
  CHECK(std::imag(psiM[1][0]) == Approx(0.9008999467));
  CHECK(std::imag(psiM[1][1]) == Approx(1.2383049726));
  // grad
  CHECK(std::imag(dpsiM[1][0][0]) == Approx(0.0025820041));
  CHECK(std::imag(dpsiM[1][0][1]) == Approx(-0.1880052537));
  CHECK(std::imag(dpsiM[1][0][2]) == Approx(-0.0025404284));
  CHECK(std::imag(dpsiM[1][1][0]) == Approx(0.1069453433));
  CHECK(std::imag(dpsiM[1][1][1]) == Approx(-0.43649593));
  CHECK(std::imag(dpsiM[1][1][2]) == Approx(-0.1069145575));
  // lapl
  CHECK(std::imag(d2psiM[1][0]) == Approx(-1.3757134676));
  CHECK(std::imag(d2psiM[1][1]) == Approx(-2.4919104576));
#endif
}

TEST_CASE("SPO input spline from HDF diamond_2x1x1", "[wavefunction]")
{
  // capture 3 spo input styles, the 2nd and 3rd ones will be deprecated and removed eventually.
  // the first case should be simplified using SPOSetBuilderFactory instead of WaveFunctionFactory
  app_log() << "-------------------------------------------------------------" << std::endl;
  app_log() << "diamondC_2x1x1 input style 1 using sposet_collection" << std::endl;
  app_log() << "-------------------------------------------------------------" << std::endl;
  const char* spo_xml_string1 = R"(<wavefunction name="psi0" target="elec">
<sposet_collection name="einspline_diamond_size4" type="einspline" href="diamondC_2x1x1.pwscf.h5" tilematrix="2 0 0 0 1 0 0 0 1" twistnum="0" source="ion" meshfactor="1.0" precision="float">
  <sposet name="spo" size="4" spindataset="0"/>
</sposet_collection>
</wavefunction>
)";
  test_diamond_2x1x1_xml_input(spo_xml_string1);

  app_log() << "-------------------------------------------------------------" << std::endl;
  app_log() << "diamondC_2x1x1 input style 2 sposet inside determinantset" << std::endl;
  app_log() << "-------------------------------------------------------------" << std::endl;
  const char* spo_xml_string2 = R"(<wavefunction name="psi0" target="elec">
<determinantset type="einspline" href="diamondC_2x1x1.pwscf.h5" tilematrix="2 0 0 0 1 0 0 0 1" twistnum="0" source="ion" meshfactor="1.0" precision="float">
  <sposet name="spo" size="4" spindataset="0"/>
  <slaterdeterminant>
    <determinant name="det" sposet="spo"/>
  </slaterdeterminant>
</determinantset>
</wavefunction>
)";
  test_diamond_2x1x1_xml_input(spo_xml_string2);

  app_log() << "-------------------------------------------------------------" << std::endl;
  app_log() << "diamondC_2x1x1 input style 3 sposet inside determinantset" << std::endl;
  app_log() << "-------------------------------------------------------------" << std::endl;
  const char* spo_xml_string3 = R"(<wavefunction name="psi0" target="elec">
<determinantset type="einspline" href="diamondC_2x1x1.pwscf.h5" tilematrix="2 0 0 0 1 0 0 0 1" twistnum="0" source="ion" meshfactor="1.0" precision="float">
  <slaterdeterminant>
    <determinant name="spo" size="4" spindataset="0"/>
  </slaterdeterminant>
</determinantset>
</wavefunction>
)";
  test_diamond_2x1x1_xml_input(spo_xml_string3);
}
} // namespace qmcplusplus
