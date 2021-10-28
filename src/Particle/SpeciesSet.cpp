//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Raymond Clay III, j.k.rofling@gmail.com, Lawrence Livermore National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


#include "SpeciesSet.h"
#include <algorithm>
namespace qmcplusplus
{
void SpeciesSet::create(unsigned m)
{
  if (m > 0)
  {
    speciesName.insert(speciesName.end(), m, std::string("none"));
    for (auto& d : d_attrib)
    {
      d.insert(d.end(), m, 0);
    }
    TotalNum += m;
  }
}

int SpeciesSet::addSpecies(const std::string& aname)
{
  int i = findSpecies(aname); // check if the name is registered
  if (i == TotalNum)
  // if not found, add a new species
  {
    create(1);
    speciesName[i] = aname;
  }
  return i; // return an index for a species
}

int SpeciesSet::addAttribute(const std::string& aname)
{
  int i = 0;
  while (i < attribName.size())
  {
    if (attribName[i] == aname)
      return i;
    i++;
  }
  attribName.push_back(aname);
  int n = d_attrib.size();
  d_attrib.push_back(SpeciesAttrib_t(TotalNum));
  return n;
}

int SpeciesSet::getAttribute(const std::string& aname) const
{
  for (int i = 0; i < attribName.size(); i++)
  {
    if (attribName[i] == aname)
      return i;
  }
  return attribName.size();
}

bool SpeciesSet::hasAttribute(const std::string& aname) const
{
  return std::find(attribName.begin(), attribName.end(), aname) != attribName.end();
}

} // namespace qmcplusplus
