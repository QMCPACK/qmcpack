//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//////////////////////////////////////////////////////////////////////////////////////

#include <catch2/catch_test_macros.hpp>
#include "QMCTools/QMCGaussianParserBase.h"
#include "Utilities/for_testing/Catch2Approx.h"
#include "Particle/ParticleIO/XMLParticleIO.h"
#include "QMCWaveFunctions/Fermion/SlaterDetBuilder.h"
#include "QMCWaveFunctions/SPOSetBuilderFactory.h"
#include "QMCWaveFunctions/TrialWaveFunction.h"
#include "QMCWaveFunctions/tests/ConstantSPOSet.h"
#include <algorithm>
#include <bitset>
#include <cmath>
#include <complex>
#include <map>
#include <memory>
#include <string>

namespace
{
using Node = std::unique_ptr<xmlNode, decltype(&xmlFreeNode)>;

struct Parser : QMCGaussianParserBase
{
  void parse(const std::string&) override {}
};

std::string attribute(xmlNodePtr node, const char* name)
{
  auto value         = xmlGetProp(node, BAD_CAST name);
  std::string result = value ? reinterpret_cast<const char*>(value) : "";
  xmlFree(value);
  return result;
}

xmlNodePtr child(xmlNodePtr node, const char* name)
{
  for (auto current = node->children; current; current = current->next)
    if (current->type == XML_ELEMENT_NODE && xmlStrEqual(current->name, BAD_CAST name))
      return current;
  return nullptr;
}

void configureOpenShell(Parser& parser)
{
  parser.NumberOfAlpha    = 3;
  parser.NumberOfBeta     = 1;
  parser.NumberOfEls      = 4;
  parser.SpinMultiplicity = 3;
  parser.SpinRestricted   = true;
  parser.numMO            = 4;
  parser.SizeOfBasisSet   = 4;
  parser.usingCSF         = true;
  parser.ci_nstates       = 4;
  parser.ci_nea           = 3;
  parser.ci_neb           = 1;
  parser.ci_size          = 2;
  parser.ci_threshold     = 1e-5;
  parser.CSFocc           = {"2110", "2101"};
  parser.CSFalpha         = {{"1110"}, {"1101"}};
  parser.CSFbeta          = {{"1000"}, {"1000"}};
  parser.CSFexpansion     = {{1.0}, {1.0}};
  parser.coeff2csf        = {{0, 0.25}, {1, -0.75}};
}
} // namespace

TEST_CASE("CSF excitation counts distinguish reference and excited occupations", "[tools][csf]")
{
  struct Case
  {
    int alpha;
    int beta;
    std::string occupation;
    int expected;
  };
  const Case cases[] = {{2, 0, "1100", 0},   {2, 0, "1010", 1},   {2, 0, "0011", 2},   {3, 1, "211000", 0},
                        {3, 1, "210100", 1}, {3, 1, "121000", 1}, {3, 1, "111100", 1}, {2, 2, "2200", 0},
                        {2, 2, "2110", 1},   {2, 2, "1120", 2},   {2, 2, "0022", 4},   {1, 0, "0100", 1}};
  Parser parser;
  for (const auto& example : cases)
  {
    CAPTURE(example.alpha, example.beta, example.occupation);
    parser.ci_nea     = example.alpha;
    parser.ci_neb     = example.beta;
    parser.ci_nstates = example.occupation.size();
    auto occupation   = example.occupation;
    CHECK(parser.numberOfExcitationsCSF(occupation) == example.expected);
    CHECK(occupation == example.occupation);
  }
}

TEST_CASE("CSF excitation ranks match exhaustive determinant substitutions", "[tools][csf]")
{
  Parser parser;
  for (int orbitals = 1; orbitals <= 6; ++orbitals)
    for (int alpha = 0; alpha <= orbitals; ++alpha)
      for (int beta = 0; beta <= alpha; ++beta)
      {
        parser.ci_nea                  = alpha;
        parser.ci_neb                  = beta;
        parser.ci_nstates              = orbitals;
        const unsigned reference_alpha = (1u << alpha) - 1u;
        const unsigned reference_beta  = (1u << beta) - 1u;
        const unsigned limit           = 1u << orbitals;
        std::map<std::string, int> minimum_substitutions;
        // Independent spin-resolved substitutions, not the spatial-occupation formula.
        for (unsigned up = 0; up < limit; ++up)
          if (std::bitset<6>(up).count() == alpha)
            for (unsigned down = 0; down < limit; ++down)
              if (std::bitset<6>(down).count() == beta)
              {
                std::string occupation(orbitals, '0');
                for (int orbital = 0; orbital < orbitals; ++orbital)
                  occupation[orbital] += ((up >> orbital) & 1u) + ((down >> orbital) & 1u);
                const int substitutions =
                    std::bitset<6>(reference_alpha & ~up).count() + std::bitset<6>(reference_beta & ~down).count();
                auto entry          = minimum_substitutions.emplace(occupation, substitutions);
                entry.first->second = std::min(entry.first->second, substitutions);
              }
        for (const auto& entry : minimum_substitutions)
        {
          CAPTURE(orbitals, alpha, beta, entry.first);
          auto occupation = entry.first;
          CHECK(parser.numberOfExcitationsCSF(occupation) == entry.second);
        }
      }
}

TEST_CASE("CSF serialization orders excitation levels without changing coefficients", "[tools][csf]")
{
  for (bool order_by_excitation : {false, true})
  {
    CAPTURE(order_by_excitation);
    Parser parser;
    configureOpenShell(parser);
    parser.orderByExcitation = order_by_excitation;
    Node multidet(parser.createMultiDeterminantSet(), xmlFreeNode);
    auto detlist = child(multidet.get(), "detlist");
    REQUIRE(detlist != nullptr);
    auto first = child(detlist, "csf");
    REQUIRE(first != nullptr);
    CHECK(attribute(first, "occ") == (order_by_excitation ? "2110" : "2101"));
    int configurations = 0;
    for (auto csf = detlist->children; csf; csf = csf->next)
      if (csf->type == XML_ELEMENT_NODE)
      {
        const bool reference = attribute(csf, "occ") == "2110";
        CHECK(attribute(csf, "occ") == (reference ? "2110" : "2101"));
        CHECK(attribute(csf, "exctLvl") == (reference ? "0" : "1"));
        CHECK(std::stod(attribute(csf, "coeff")) == (reference ? 0.25 : -0.75));
        CHECK(std::stod(attribute(csf, "qchem_coeff")) == (reference ? 0.25 : -0.75));
        auto determinant = child(csf, "det");
        REQUIRE(determinant != nullptr);
        CHECK(attribute(determinant, "alpha") == (reference ? "1110" : "1101"));
        CHECK(attribute(determinant, "beta") == "1000");
        CHECK(attribute(determinant, "coeff") == "1");
        ++configurations;
      }
    CHECK(configurations == 2);
  }
}

TEST_CASE("CSF zero excitation cutoff retains the open-shell reference", "[tools][csf]")
{
  using Value = QMCTraits::ValueType;
  for (const std::string algorithm : {"table_method", "precomputed_table_method"})
  {
    CAPTURE(algorithm);
    Parser parser;
    configureOpenShell(parser);
    Node electron_xml(parser.createElectronSet("ion0"), xmlFreeNode);
    SimulationCell cell;
    ParticleSet electrons(cell);
    XMLParticleParser particle_reader(electrons);
    REQUIRE(particle_reader.readXML(electron_xml.get()));
    REQUIRE(electrons.groups() == 2);
    for (int electron = 0; electron < electrons.getTotalNum(); ++electron)
      electrons.R[electron] = {static_cast<double>(electron), 0.0, 0.0};
    electrons.update();
    SPOSetBuilderFactory::PSetMap particles;
    SPOSetBuilderFactory factory(*OHMMS::Controller, electrons, particles);
    auto up = std::make_unique<ConstantSPOSet<Value>>("spo-up", 3, 4);
    Matrix<Value> up_orbitals(3, 4);
    up_orbitals = Value(0);
    for (int electron = 0; electron < 3; ++electron)
    {
      up_orbitals[electron][electron] = 1;
      up_orbitals[electron][3]        = 1;
    }
    up->setRefVals(up_orbitals);
    factory.addSPOSet(std::move(up));
    auto down = std::make_unique<ConstantSPOSet<Value>>("spo-dn", 1, 4);
    Matrix<Value> down_orbitals(1, 4);
    down_orbitals = Value(2);
    down->setRefVals(down_orbitals);
    factory.addSPOSet(std::move(down));
    Node determinant_xml(xmlNewNode(nullptr, BAD_CAST "determinantset"), xmlFreeNode);
    auto multidet = parser.createMultiDeterminantSet();
    auto detlist  = child(multidet, "detlist");
    REQUIRE(detlist != nullptr);
    auto reference = child(detlist, "csf");
    while (reference && attribute(reference, "occ") != "2110")
      reference = reference->next;
    REQUIRE(reference != nullptr);
    REQUIRE(attribute(reference, "exctLvl") == "0");
    xmlSetProp(detlist, BAD_CAST "sortby", BAD_CAST "exctLvl");
    xmlSetProp(detlist, BAD_CAST "cutoff", BAD_CAST "0");
    xmlNewProp(multidet, BAD_CAST "algorithm", BAD_CAST algorithm.c_str());
    xmlAddChild(determinant_xml.get(), multidet);
    RuntimeOptions options;
    TrialWaveFunction wavefunction(options);
    SlaterDetBuilder builder(*OHMMS::Controller, factory, electrons, wavefunction, particles);
    auto component = builder.buildComponent(determinant_xml.get());
    REQUIRE(component != nullptr);
    wavefunction.addComponent(std::move(component));
    wavefunction.evaluateLog(electrons);
    const double det_reference =
        1.0 * (1.0 * 1.0 - 0.0 * 0.0) - 0.0 * (0.0 * 1.0 - 0.0 * 0.0) + 0.0 * (0.0 * 0.0 - 1.0 * 0.0);
    const double expected = 0.25 * det_reference * 2.0;
    CHECK(wavefunction.getLogPsi() == Approx(std::log(expected)));
    CHECK(std::polar(std::exp(wavefunction.getLogPsi()), wavefunction.getPhase()) ==
          ComplexApprox(std::complex<double>(expected, 0.0)));
  }
}
