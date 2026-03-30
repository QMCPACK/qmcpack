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

#ifndef QMCPLUSPLUS_PSIHAMNAMEPAIRREADER_H
#define QMCPLUSPLUS_PSIHAMNAMEPAIRREADER_H

#include <string>
#include <vector>
#include <optional>
#include "OhmmsData/libxmldefs.h"

namespace qmcplusplus
{
/** Input representation for Driver base class runtime parameters
 */
class PsiHamNamePairReader
{
public:
  using PsiHamNamePair = std::pair<std::string, std::string>;

  static std::optional<PsiHamNamePair> readOnePair(xmlNodePtr cur);
  static std::vector<PsiHamNamePair> readMultiplePairs(xmlNodePtr cur);

private:
  static std::vector<PsiHamNamePair> readPairs(xmlNodePtr cur);
};

} // namespace qmcplusplus

#endif
