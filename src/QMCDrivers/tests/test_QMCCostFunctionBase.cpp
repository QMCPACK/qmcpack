//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2022 QMCPACK developers.
//
// File developed by: Mark Dewing, mdewing@anl.gov, Argonne National Laboratory
//
// File created by: Mark Dewing, mdewing@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////

#include "catch.hpp"
#include "QMCDrivers/WFOpt/QMCCostFunctionBase.h"
#include "Utilities/RuntimeOptions.h"


namespace qmcplusplus
{

class QMCCostFunctionTest : public QMCCostFunctionBase
{
public:
  QMCCostFunctionTest(MCWalkerConfiguration& w, TrialWaveFunction& psi, QMCHamiltonian& h, Communicate* c)
      : QMCCostFunctionBase(w, psi, h, c)
  {}

  void resetPsi(bool final_reset = false) override {}
  Return_rt fillOverlapHamiltonianMatrices(Matrix<Return_rt>& Left, Matrix<Return_rt>& Right) override { return 0; }
  void getConfigurations(const std::string& aroot) override {}
  void checkConfigurations(EngineHandle& handle) override {}
  EffectiveWeight correlatedSampling(bool needGrad = true) override { return 0; }

  void callUpdateXmlNodes()
  {
    do_override_output = true;
    updateXmlNodes();
  }

  // Print the updated XML nodes. Useful for debugging.
  void printXml()
  {
    char* buf;
    int out_buflen;
    xmlDocDumpFormatMemory(m_doc_out, (xmlChar**)&buf, &out_buflen, 1);
    std::cout << "XML buffer length = " << out_buflen << std::endl;
    std::cout << "XML: " << std::endl << buf << std::endl;
    xmlFree(buf);
  }

  xmlDocPtr getDoc() { return m_doc_out; }
};


TEST_CASE("updateXmlNodes", "[drivers]")
{
  const SimulationCell simulation_cell;
  MCWalkerConfiguration w(simulation_cell);
  QMCHamiltonian h;
  RuntimeOptions runtime_options;
  TrialWaveFunction psi(runtime_options);

  Communicate* comm = OHMMS::Controller;

  QMCCostFunctionTest cost(w, psi, h, comm);

  cost.setRootName("tmp");

  const char* wf_xml = R"(
<wavefunction/>
    )";

  Libxml2Document doc;
  bool okay = doc.parseFromString(wf_xml);
  REQUIRE(okay);
  cost.setWaveFunctionNode(doc.getRoot());

  cost.callUpdateXmlNodes();
  //cost.printXml();

  xmlXPathContextPtr acontext = xmlXPathNewContext(cost.getDoc());
  OhmmsXPathObject check_elem("//override_variational_parameters", acontext);
  REQUIRE(check_elem.size() == 1);

  std::string href = getXMLAttributeValue(check_elem[0], "href");
  REQUIRE(href == "tmp.vp.h5");

  xmlXPathFreeContext(acontext);
}

TEST_CASE("updateXmlNodes with existing element", "[drivers]")
{
  const SimulationCell simulation_cell;
  MCWalkerConfiguration w(simulation_cell);
  QMCHamiltonian h;
  RuntimeOptions runtime_options;
  TrialWaveFunction psi(runtime_options);

  Communicate* comm = OHMMS::Controller;

  QMCCostFunctionTest cost(w, psi, h, comm);
  cost.setRootName("tmp2");

  const char* wf_xml = R"(
<wavefunction>
  <override_variational_parameters href="vp.tmp.h5"/>
</wavefunction>
    )";

  Libxml2Document doc;
  bool okay = doc.parseFromString(wf_xml);
  REQUIRE(okay);
  cost.setWaveFunctionNode(doc.getRoot());

  cost.callUpdateXmlNodes();
  //cost.printXml();

  xmlXPathContextPtr acontext = xmlXPathNewContext(cost.getDoc());
  OhmmsXPathObject check_elem("//override_variational_parameters", acontext);
  REQUIRE(check_elem.size() == 1);

  std::string href = getXMLAttributeValue(check_elem[0], "href");
  REQUIRE(href == "tmp2.vp.h5");

  xmlXPathFreeContext(acontext);
}


} // namespace qmcplusplus
