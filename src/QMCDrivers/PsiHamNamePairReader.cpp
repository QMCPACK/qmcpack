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

#include "PsiHamNamePairReader.h"
#include "OhmmsData/AttributeSet.h"
#include "Message/UniformCommunicateError.h"

namespace qmcplusplus
{
using PsiHamNamePair = typename PsiHamNamePairReader::PsiHamNamePair;

std::optional<PsiHamNamePair> PsiHamNamePairReader::readOnePair(xmlNodePtr cur)
{
  std::vector<PsiHamNamePair> name_pairs = readPairs(cur);
  if (name_pairs.size() > 1)
    throw UniformCommunicateError(
        "Found more than one \"qmcsystem\" node in the XML file but can only take one. Please fix the XML input!");
  else if (name_pairs.size() == 1)
    return name_pairs.front();
  else
    return std::nullopt;
}

std::vector<PsiHamNamePair> PsiHamNamePairReader::readMultiplePairs(xmlNodePtr cur)
{
  std::vector<PsiHamNamePair> name_pairs = readPairs(cur);
  if (name_pairs.size() < 2)
    throw UniformCommunicateError(
        "Found fewer than \"two qmcsystem\" nodes in the XML file but requires two or more. Please fix the XML input!");
  return name_pairs;
}

std::vector<PsiHamNamePair> PsiHamNamePairReader::readPairs(xmlNodePtr cur)
{
  std::string psi, ham;
  OhmmsAttributeSet attrib;
  attrib.add(psi, "wavefunction");
  attrib.add(ham, "hamiltonian");
  std::vector<PsiHamNamePair> name_pairs;
  processChildren(cur, [&](const std::string& cname, const xmlNodePtr element) {
    if (cname == "qmcsystem")
    {
      psi = "";
      ham = "";
      attrib.put(element);
      name_pairs.emplace_back(psi, ham);
    }
  });
  return name_pairs;
}

} // namespace qmcplusplus
