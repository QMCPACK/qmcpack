//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//////////////////////////////////////////////////////////////////////////////////////

#include <catch2/catch_test_macros.hpp>
#include "QMCTools/QMCGaussianParserBase.h"
#include "Utilities/for_testing/Catch2Approx.h"
#include "Particle/ParticleIO/XMLParticleIO.h"
#include "Particle/ParticleSetPool.h"
#include "QMCWaveFunctions/Fermion/SlaterDetBuilder.h"
#include "QMCWaveFunctions/Jastrow/eeI_JastrowBuilder.h"
#include "QMCWaveFunctions/SPOSetBuilderFactory.h"
#include "QMCWaveFunctions/TrialWaveFunction.h"
#include "QMCWaveFunctions/tests/ConstantSPOSet.h"
#include <cmath>
#include <cstdio>
#include <memory>
#include <string>

namespace
{
using Node = std::unique_ptr<xmlNode, decltype(&xmlFreeNode)>;

std::string attribute(xmlNodePtr node, const char* name)
{
  auto value               = xmlGetProp(node, BAD_CAST name);
  const std::string result = value ? reinterpret_cast<const char*>(value) : "";
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

struct Converter : QMCGaussianParserBase
{
  explicit Converter(bool polarized)
  {
    NumberOfAlpha    = 2;
    NumberOfBeta     = polarized ? 0 : 1;
    NumberOfEls      = NumberOfAlpha + NumberOfBeta;
    SpinMultiplicity = NumberOfAlpha - NumberOfBeta + 1;
    SpinRestricted   = true;
    numMO            = 4;
    SizeOfBasisSet   = 4;
    EigVec.resize(32, 0.0);
    ci_nstates   = 4;
    ci_size      = 2;
    ci_nea       = 2;
    ci_neb       = NumberOfBeta;
    ci_threshold = 1e-5;
    CIcoeff      = {0.75, -0.25};
    CIalpha      = {"1100", "1010"};
    CIbeta       = polarized ? std::vector<std::string>{"0000", "0000"} : std::vector<std::string>{"1000", "0100"};
    CSFocc       = polarized ? CIalpha : std::vector<std::string>{"2100", "1110"};
    CSFalpha     = {{CIalpha[0]}, {CIalpha[1]}};
    CSFbeta      = {{CIbeta[0]}, {CIbeta[1]}};
    CSFexpansion = {{1.0}, {-1.0}};
    coeff2csf    = {{0, 0.75}, {1, -0.25}};
  }

  void parse(const std::string&) override {}
};

void check_channels(xmlNodePtr node, bool polarized)
{
  REQUIRE(node != nullptr);
  auto detlist = child(node, "detlist");
  REQUIRE(detlist != nullptr);
  CHECK(attribute(node, "spo_0") == (polarized ? "spo-up" : ""));
  CHECK(attribute(node, "spo_up") == (polarized ? "" : "spo-up"));
  CHECK(attribute(node, "spo_dn") == (polarized ? "" : "spo-dn"));
  CHECK(attribute(detlist, "nc0") == (polarized ? "0" : ""));
  CHECK(attribute(detlist, "ne0") == (polarized ? "2" : ""));
  CHECK(attribute(detlist, "nca") == (polarized ? "" : "0"));
  CHECK(attribute(detlist, "ncb") == (polarized ? "" : "0"));
  CHECK(attribute(detlist, "nea") == (polarized ? "" : "2"));
  CHECK(attribute(detlist, "neb") == (polarized ? "" : "1"));
  CHECK(attribute(detlist, "size") == "2");
  CHECK(attribute(detlist, "nstates") == "4");
}
} // namespace

TEST_CASE("convert4qmc omits empty electron groups", "[tools][converter]")
{
  for (bool polarized : {false, true})
  {
    CAPTURE(polarized);
    Converter converter(polarized);
    Node node(converter.createElectronSet("ion0"), xmlFreeNode);
    CHECK(attribute(node.get(), "spinor").empty());
    auto up = child(node.get(), "group");
    REQUIRE(up != nullptr);
    CHECK(attribute(up, "name") == "u");
    CHECK(attribute(up, "size") == "2");
    xmlNodePtr down = nullptr;
    for (auto current = up->next; current; current = current->next)
      if (current->type == XML_ELEMENT_NODE)
        down = current;
    if (polarized)
      CHECK(down == nullptr);
    else
    {
      REQUIRE(down != nullptr);
      CHECK(attribute(down, "name") == "d");
      CHECK(attribute(down, "size") == "1");
    }
    CHECK_FALSE(converter.isSpinor);
  }
}

TEST_CASE("convert4qmc XML determinant channels", "[tools][converter]")
{
  for (bool polarized : {false, true})
  {
    CAPTURE(polarized);
    Converter converter(polarized);
    Node node(converter.createMultiDeterminantSet(), xmlFreeNode);
    check_channels(node.get(), polarized);
    auto detlist = child(node.get(), "detlist");
    REQUIRE(detlist != nullptr);
    CHECK(attribute(detlist, "type") == "DETS");
    int index = 0;
    for (auto ci = detlist->children; ci; ci = ci->next)
      if (ci->type == XML_ELEMENT_NODE)
      {
        REQUIRE(index < 2);
        CHECK(std::stod(attribute(ci, "coeff")) == converter.CIcoeff[index]);
        CHECK(attribute(ci, polarized ? "occ0" : "alpha") == converter.CIalpha[index]);
        CHECK(attribute(ci, polarized ? "alpha" : "occ0").empty());
        CHECK(attribute(ci, "beta") == (polarized ? "" : converter.CIbeta[index]));
        ++index;
      }
    CHECK(index == 2);
    CHECK_FALSE(converter.isSpinor);
  }
}

TEST_CASE("convert4qmc XML CSF channels preserve expansion", "[tools][converter]")
{
  for (bool polarized : {false, true})
  {
    CAPTURE(polarized);
    Converter converter(polarized);
    converter.usingCSF = true;
    Node node(converter.createMultiDeterminantSet(), xmlFreeNode);
    check_channels(node.get(), polarized);
    auto detlist = child(node.get(), "detlist");
    REQUIRE(detlist != nullptr);
    CHECK(attribute(detlist, "type") == "CSF");
    int index = 0;
    for (auto csf = detlist->children; csf; csf = csf->next)
      if (csf->type == XML_ELEMENT_NODE)
      {
        REQUIRE(index < 2);
        CHECK(std::stod(attribute(csf, "coeff")) == converter.coeff2csf[index].second);
        CHECK(attribute(csf, "occ") == converter.CSFocc[index]);
        auto det = child(csf, "det");
        REQUIRE(det != nullptr);
        CHECK(std::stod(attribute(det, "coeff")) == converter.CSFexpansion[index][0]);
        CHECK(attribute(det, polarized ? "occ0" : "alpha") == converter.CIalpha[index]);
        CHECK(attribute(det, polarized ? "alpha" : "occ0").empty());
        CHECK(attribute(det, "beta") == (polarized ? "" : converter.CIbeta[index]));
        ++index;
      }
    CHECK(index == 2);
    CHECK_FALSE(converter.isSpinor);
  }
}

TEST_CASE("convert4qmc HDF5 input determinant channels", "[tools][converter]")
{
  for (bool polarized : {false, true})
  {
    CAPTURE(polarized);
    Converter converter(polarized);
    converter.multih5file = "orbitals.h5";
    Node node(converter.createMultiDeterminantSetFromH5(), xmlFreeNode);
    check_channels(node.get(), polarized);
    CHECK(attribute(child(node.get(), "detlist"), "href") == "orbitals.h5");
    CHECK_FALSE(converter.isSpinor);
  }
}

TEST_CASE("convert4qmc HDF5 output preserves CI payload", "[tools][converter]")
{
  const std::string filename = "test_convert4qmc_empty_spin.h5";
  struct Cleanup
  {
    const std::string& path;
    ~Cleanup() { std::remove(path.c_str()); }
  } cleanup{filename};
  for (bool polarized : {false, true})
  {
    CAPTURE(polarized);
    Converter converter(polarized);
    converter.h5file = filename;
    {
      hdf_archive file;
      REQUIRE(file.create(filename));
    }
    Node node(converter.createMultiDeterminantSetCIHDF5(), xmlFreeNode);
    check_channels(node.get(), polarized);
    hdf_archive file;
    REQUIRE(file.open(filename));
    file.push("MultiDet");
    std::vector<double> coefficients;
    REQUIRE(file.readEntry(coefficients, "Coeff"));
    CHECK(coefficients == converter.CIcoeff);
    Matrix<uint64_t> occupations;
    REQUIRE(file.readEntry(occupations, "CI_Alpha"));
    REQUIRE(occupations.rows() == 2);
    REQUIRE(occupations.cols() == 1);
    CHECK(occupations[0][0] == 3);
    CHECK(occupations[1][0] == 5);
    CHECK_FALSE(converter.isSpinor);
  }
}

TEST_CASE("convert4qmc polarized output evaluates in SlaterDetBuilder", "[tools][converter]")
{
  using Value                = QMCTraits::ValueType;
  const std::string filename = "test_convert4qmc_reader.h5";
  struct Cleanup
  {
    const std::string& path;
    ~Cleanup() { std::remove(path.c_str()); }
  } cleanup{filename};
  for (const std::string algorithm : {"table_method", "precomputed_table_method"})
    for (const std::string format : {"DETS", "CSF", "HDF5"})
    {
      CAPTURE(algorithm, format);
      Converter converter(true);
      converter.usingCSF = format == "CSF";
      converter.h5file   = filename;
      if (format == "HDF5")
      {
        hdf_archive file;
        REQUIRE(file.create(filename));
      }
      Node electron_xml(converter.createElectronSet("ion0"), xmlFreeNode);
      SimulationCell cell;
      ParticleSet electrons(cell);
      XMLParticleParser particle_reader(electrons);
      REQUIRE(particle_reader.readXML(electron_xml.get()));
      REQUIRE(electrons.groups() == 1);
      REQUIRE(electrons.getTotalNum() == 2);
      CHECK_FALSE(electrons.isSpinor());
      electrons.R[0] = {0.0, 0.0, 0.0};
      electrons.R[1] = {1.0, 0.0, 0.0};
      electrons.update();

      SPOSetBuilderFactory::PSetMap particles;
      SPOSetBuilderFactory factory(*OHMMS::Controller, electrons, particles);
      auto spo = std::make_unique<ConstantSPOSet<Value>>("spo-up", 2, 4);
      Matrix<Value> orbitals(2, 4);
      orbitals[0][0] = 1;
      orbitals[0][1] = 2;
      orbitals[0][2] = 3;
      orbitals[0][3] = 4;
      orbitals[1][0] = 2;
      orbitals[1][1] = 6;
      orbitals[1][2] = 9;
      orbitals[1][3] = 13;
      spo->setRefVals(orbitals);
      factory.addSPOSet(std::move(spo));

      Node determinant_xml(xmlNewNode(nullptr, BAD_CAST "determinantset"), xmlFreeNode);
      auto multidet =
          format == "HDF5" ? converter.createMultiDeterminantSetCIHDF5() : converter.createMultiDeterminantSet();
      xmlNewProp(multidet, BAD_CAST "algorithm", BAD_CAST algorithm.c_str());
      xmlAddChild(determinant_xml.get(), multidet);
      RuntimeOptions options;
      TrialWaveFunction wavefunction(options);
      SlaterDetBuilder builder(*OHMMS::Controller, factory, electrons, wavefunction, particles);
      auto component = builder.buildComponent(determinant_xml.get());
      REQUIRE(component != nullptr);
      wavefunction.addComponent(std::move(component));
      converter.IonSystem.setName("ion0");
      converter.IonSystem.create({1});
      converter.IonSystem.getSpeciesSet().addSpecies("H");
      converter.IonSystem.R[0] = {0.0, 1.0, 0.0};
      converter.IonSystem.update();
      Node jastrow_xml(converter.createJ3(), xmlFreeNode);
      eeI_JastrowBuilder jastrow_builder(*OHMMS::Controller, electrons, converter.IonSystem);
      auto jastrow = jastrow_builder.buildComponent(jastrow_xml.get());
      REQUIRE(jastrow != nullptr);
      wavefunction.addComponent(std::move(jastrow));
      electrons.update();
      wavefunction.evaluateLog(electrons);

      // Independent 2x2 determinants keep the conversion and reader oracle separate.
      const double det01 = 1.0 * 6.0 - 2.0 * 2.0;
      const double det02 = 1.0 * 9.0 - 3.0 * 2.0;
      const double expected =
          format == "CSF" ? 0.75 * det01 + (-0.25) * (-1.0) * det02 : 0.75 * det01 + (-0.25) * det02;
      CHECK(wavefunction.getLogPsi() == Approx(std::log(expected)));
      CHECK(std::polar(std::exp(wavefunction.getLogPsi()), wavefunction.getPhase()) ==
            ComplexApprox(std::complex<double>(expected, 0.0)));
      CHECK_FALSE(converter.isSpinor);
    }
}

TEST_CASE("convert4qmc single determinant omits empty spin", "[tools][converter]")
{
  const std::string filename = "test_convert4qmc_single.h5";
  struct Cleanup
  {
    const std::string& path;
    ~Cleanup() { std::remove(path.c_str()); }
  } cleanup{filename};
  for (bool polarized : {false, true})
    for (const std::string format : {"XML", "HDF5 input", "HDF5 output"})
    {
      CAPTURE(polarized, format);
      Converter converter(polarized);
      converter.h5file = filename;
      if (format == "HDF5 output")
      {
        hdf_archive file;
        REQUIRE(file.create(filename));
      }
      Node node(format == "XML"              ? converter.createDeterminantSet()
                    : format == "HDF5 input" ? converter.PrepareDeterminantSetFromHDF5()
                                             : converter.createDeterminantSetWithHDF5(),
                xmlFreeNode);
      auto up = child(node.get(), "determinant");
      REQUIRE(up != nullptr);
      CHECK(attribute(up, "id") == "updet");
      CHECK(attribute(up, "size") == "2");
      if (polarized)
        CHECK(up->next == nullptr);
      else
      {
        REQUIRE(up->next != nullptr);
        CHECK(attribute(up->next, "id") == "downdet");
        CHECK(attribute(up->next, "size") == "1");
      }
      CHECK_FALSE(converter.isSpinor);
    }
}

TEST_CASE("convert4qmc polarized UHF retains unoccupied HDF5 orbitals", "[tools][converter]")
{
  const std::string filename = "test_convert4qmc_uhf.h5";
  struct Cleanup
  {
    const std::string& path;
    ~Cleanup() { std::remove(path.c_str()); }
  } cleanup{filename};
  Converter converter(true);
  converter.SpinRestricted = false;
  converter.EigVal_alpha.resize(4, 0.0);
  converter.EigVal_beta.resize(4, 0.0);
  for (int index = 0; index < 32; ++index)
    converter.EigVec[index] = index + 1;
  converter.h5file = filename;
  {
    hdf_archive file;
    REQUIRE(file.create(filename));
  }
  Node node(converter.createDeterminantSetWithHDF5(), xmlFreeNode);
  auto up = child(node.get(), "determinant");
  REQUIRE(up != nullptr);
  CHECK(up->next == nullptr);
  hdf_archive file;
  REQUIRE(file.open(filename));
  file.push("Super_Twist");
  Matrix<double> down_orbitals;
  REQUIRE(file.readEntry(down_orbitals, "eigenset_1"));
  REQUIRE(down_orbitals.rows() == 4);
  REQUIRE(down_orbitals.cols() == 4);
  for (int row = 0; row < 4; ++row)
    for (int column = 0; column < 4; ++column)
      CHECK(down_orbitals[row][column] == converter.EigVec[16 + row * 4 + column]);
}

TEST_CASE("convert4qmc J3 references only occupied electron species", "[tools][converter]")
{
  for (const std::string mode : {"polarized", "two groups", "spinor"})
  {
    CAPTURE(mode);
    Converter converter(mode != "two groups");
    converter.isSpinor = mode == "spinor";
    converter.IonSystem.setName("ion0");
    converter.IonSystem.create({1});
    converter.IonSystem.getSpeciesSet().addSpecies("H");
    converter.IonSystem.R[0] = {0.0, 1.0, 0.0};
    converter.IonSystem.update();
    Node electron_xml(converter.createElectronSet("ion0"), xmlFreeNode);
    ParticleSetPool particle_pool(*OHMMS::Controller);
    xmlSetProp(electron_xml.get(), BAD_CAST "random", BAD_CAST "no");
    REQUIRE(particle_pool.put(electron_xml.get()));
    auto parsed_electrons = particle_pool.getParticleSet("e");
    REQUIRE(parsed_electrons != nullptr);
    auto& electrons = *parsed_electrons;
    for (int electron = 0; electron < electrons.getTotalNum(); ++electron)
      electrons.R[electron] = {static_cast<double>(electron), 0.0, 0.0};
    electrons.update();
    Node jastrow_xml(converter.createJ3(), xmlFreeNode);
    auto& species    = electrons.getSpeciesSet();
    int correlations = 0;
    for (auto correlation = jastrow_xml->children; correlation; correlation = correlation->next)
      if (correlation->type == XML_ELEMENT_NODE)
      {
        REQUIRE(species.findSpecies(attribute(correlation, "especies1")) < species.size());
        REQUIRE(species.findSpecies(attribute(correlation, "especies2")) < species.size());
        ++correlations;
      }
    CHECK(correlations == (mode == "two groups" ? 2 : 1));
    eeI_JastrowBuilder builder(*OHMMS::Controller, electrons, converter.IonSystem);
    auto jastrow = builder.buildComponent(jastrow_xml.get());
    REQUIRE(jastrow != nullptr);
    CHECK(electrons.isSpinor() == converter.isSpinor);
  }
}
