#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_exception.hpp>
#include "Message/Communicate.h"
#include "Message/UniformCommunicateError.h"
#include "QMCApp/QMCMain.h"
#include "OhmmsData/Libxml2Doc.h"

namespace qmcplusplus
{

class TestQMCMain : public QMCMain
{
public:
  using QMCMain::QMCMain;

  void pushDocString(const std::string& xml_str)
  {
    auto* doc = new Libxml2Document();
    REQUIRE(doc->parseFromString(xml_str));
    xml_doc_stack_.push(doc);
  }
};

TEST_CASE("QMCMain validateXML exceptions", "[qmcapp]")
{
  Communicate& c(*OHMMS::Controller);

  SECTION("include tag missing href attribute")
  {
    TestQMCMain qmc(c);
    qmc.pushDocString("<simulation><include/></simulation>");
    CHECK_THROWS_MATCHES(qmc.validateXML(), UniformCommunicateError,
                         Catch::Matchers::Message(R"(tag "include" must include an "href" attribute.)"));
  }

  SECTION("include tag invalid XML document")
  {
    TestQMCMain qmc(c);
    qmc.pushDocString(R"(<simulation><include href="nonexistent_include_file.xml"/></simulation>)");
    CHECK_THROWS_MATCHES(qmc.validateXML(), UniformCommunicateError,
                         Catch::Matchers::Message("Invalid XML document"));
  }

  SECTION("missing particleset")
  {
    TestQMCMain qmc(c);
    qmc.pushDocString("<simulation><project id=\"test\" series=\"0\"/></simulation>");
    CHECK_THROWS_MATCHES(qmc.validateXML(), UniformCommunicateError,
                         Catch::Matchers::Message("QMCMain::validateXML. Illegal input. Missing particleset."));
  }

  SECTION("missing wavefunction")
  {
    TestQMCMain qmc(c);
    const char* input = R"(<simulation>
      <project id="test" series="0"/>
      <particleset name="e" size="1">
        <group name="u">
          <parameter name="charge">-1</parameter>
        </group>
      </particleset>
    </simulation>)";
    qmc.pushDocString(input);
    CHECK_THROWS_MATCHES(qmc.validateXML(), UniformCommunicateError,
                         Catch::Matchers::Message("QMCMain::validateXML. Illegal input. Missing wavefunction."));
  }

  SECTION("missing Hamiltonian")
  {
    TestQMCMain qmc(c);
    const char* input = R"(<simulation>
      <project id="test" series="0"/>
      <particleset name="e" size="1">
        <group name="u">
          <parameter name="charge">-1</parameter>
        </group>
      </particleset>
      <wavefunction target="e" name="psi0"/>
    </simulation>)";
    qmc.pushDocString(input);
    CHECK_THROWS_MATCHES(qmc.validateXML(), UniformCommunicateError,
                         Catch::Matchers::Message("QMCMain::validateXML. Illegal input. Missing Hamiltonian."));
  }

  SECTION("multiple global estimators nodes")
  {
    TestQMCMain qmc(c);
    const char* input = R"(<simulation>
      <qmcsystem>
        <estimators/>
        <estimators/>
      </qmcsystem>
    </simulation>)";
    qmc.pushDocString(input);
    CHECK_THROWS_MATCHES(qmc.validateXML(), UniformCommunicateError,
                         Catch::Matchers::Message("QMCMain::validateXML. Illegal Input, only one global <estimators> node is permitted"));
  }
}

} // namespace qmcplusplus
