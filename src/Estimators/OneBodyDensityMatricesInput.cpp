//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2021 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//
// Some code refactored from: DensityMatrices1b.cpp
//////////////////////////////////////////////////////////////////////////////////////

#include "string_utils.h"
#include "OneBodyDensityMatricesInput.h"

namespace qmcplusplus
{

OneBodyDensityMatricesInput::OneBodyDensityMatricesInput(){};
OneBodyDensityMatricesInput::OneBodyDensityMatricesInput(xmlNodePtr cur)
{
  // This results in checkParticularValidity being called on OneBodyDensityMatrixInputSection
  input_section_.readXML(cur);
}

void OneBodyDensityMatricesInput::OneBodyDensityMatrixInputSection::checkParticularValidity()
{
  if (has("scale"))
  {
    Real scale = get<Real>("scale");
    std::cout << "SCALE is :" << scale << '\n';
    if (scale > 1.0 + 1e-10)
      throw UniformCommunicateError("OneBodyDensityMatrices input: scale must be less than one");
    else if (scale < 0.0 - 1e-10)
      throw UniformCommunicateError("OneBodyDensityMatrices input: scale must be greater than zero");
  }
  std::string basis_str = get<std::string>("basis");
  auto basis_set_names  = split(basis_str);
  if (basis_set_names.size() == 0 || basis_set_names[0].size() == 0)
    throw UniformCommunicateError("OneBodyDensityMatrices input: basis must have at least one sposet");
}

} // namespace qmcplusplus
