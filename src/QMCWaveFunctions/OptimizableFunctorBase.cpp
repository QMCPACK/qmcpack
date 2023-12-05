//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Luke Shulenburger, lshulen@sandia.gov, Sandia National Laboratories
//
// File created by: Luke Shulenburger, lshulen@sandia.gov, Sandia National Laboratories
//////////////////////////////////////////////////////////////////////////////////////

#include "OptimizableFunctorBase.h"
#include "OhmmsData/XMLParsingString.h"

namespace qmcplusplus
{
void print(OptimizableFunctorBase& func, std::ostream& os, double extent)
{
  using real_type = OptimizableFunctorBase::real_type;
  int n           = 1000;
  real_type d     = extent == -1.0 ? func.cutoff_radius / n : extent / n;
  real_type r     = 0;
  real_type u, du;
  for (int i = 0; i < n; ++i)
  {
    u  = func.f(r);
    du = func.df(r);
    os << std::setw(22) << r << std::setw(22) << u << std::setw(22) << du << std::endl;
    r += d;
  }
}

std::string extractCoefficientsID(xmlNodePtr cur)
{
  xmlNodePtr xmlCoefs = cur->xmlChildrenNode;
  while (xmlCoefs != NULL)
  {
    std::string cname((const char*)xmlCoefs->name);
    if (cname == "coefficients")
      return getXMLAttributeValue(xmlCoefs, "id");
    xmlCoefs = xmlCoefs->next;
  }
  return "";
}
} // namespace qmcplusplus
