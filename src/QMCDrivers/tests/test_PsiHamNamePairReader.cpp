//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2026 QMCPACK developers.
//
// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////

#include "catch.hpp"
#include "PsiHamNamePairReader.h"
#include "OhmmsData/Libxml2Doc.h"

namespace qmcplusplus
{

TEST_CASE("PsiHamNamePairReader", "[drivers]")
{
  Libxml2Document doc;

  static constexpr const char* const zero_entry_input = R"(
  <qmcdriver>
  </qmcdriver>
    )";

  REQUIRE(doc.parseFromString(zero_entry_input));
  xmlNodePtr zero_entry_node = doc.getRoot();

  auto zero_entry = PsiHamNamePairReader::readOnePair(zero_entry_node);
  CHECK(!zero_entry);
  CHECK_THROWS_WITH(PsiHamNamePairReader::readMultiplePairs(zero_entry_node),
                    Catch::Matchers::Equals("Found fewer than \"two qmcsystem\" nodes in the XML file but requires two "
                                            "or more. Please fix the XML input!"));

  static constexpr const char* const one_entry_input = R"(
  <qmcdriver>
    <qmcsystem wavefunction="psi1" hamiltonian="h1"/>
  </qmcdriver>
    )";

  REQUIRE(doc.parseFromString(one_entry_input));
  xmlNodePtr one_entry_node = doc.getRoot();

  auto one_entry = PsiHamNamePairReader::readOnePair(one_entry_node);
  CHECK(one_entry);
  CHECK(one_entry->first == "psi1");
  CHECK(one_entry->second == "h1");
  CHECK_THROWS_WITH(PsiHamNamePairReader::readMultiplePairs(one_entry_node),
                    Catch::Matchers::Equals("Found fewer than \"two qmcsystem\" nodes in the XML file but requires two "
                                            "or more. Please fix the XML input!"));

  static constexpr const char* const two_entry_input = R"(
  <qmcdriver>
    <qmcsystem wavefunction="psi0"/>
    <qmcsystem wavefunction="psi1" hamiltonian="h1"/>
  </qmcdriver>
    )";

  REQUIRE(doc.parseFromString(two_entry_input));
  xmlNodePtr two_entry_node = doc.getRoot();

  CHECK_THROWS_WITH(PsiHamNamePairReader::readOnePair(two_entry_node),
                    Catch::Matchers::Equals("Found more than one \"qmcsystem\" node in the XML file but can only take "
                                            "one. Please fix the XML input!"));
  auto two_entry = PsiHamNamePairReader::readMultiplePairs(two_entry_node);
  REQUIRE(two_entry.size() == 2);
  CHECK(two_entry[0].first == "psi0");
  CHECK(two_entry[0].second == "");
  CHECK(two_entry[1].first == "psi1");
  CHECK(two_entry[1].second == "h1");
}

} // namespace qmcplusplus
