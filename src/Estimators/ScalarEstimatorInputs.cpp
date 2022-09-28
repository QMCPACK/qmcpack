//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2022 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//
// Some code refactored from: EstimatorManagerNew.cpp
//////////////////////////////////////////////////////////////////////////////////////

#include "ScalarEstimatorInputs.h"
#include "InputSection.h"
#include "EstimatorInput.h"
#include "ModernStringUtils.hpp"

namespace qmcplusplus
{
/** make input class state consistent with legacy ignoring the distinction
 *  between the name and type of estimator for scalar estimators.
 *
 *  \param [in]      input_section_    naming here is to allow reuse of LAMBDA_setIfInInput macro
 *  \param [inout]   name              output param but assumed to be unset
 *  \param [inout]   type              output param but assumed to be unset
 */
void readNameTypeLikeLegacy(InputSection& input_section_, std::string& name, std::string& type)
{
  auto setIfInInput = LAMBDA_setIfInInput;
  setIfInInput(type, "type");
  setIfInInput(name, "name");
  if (name.empty())
    name = type;
  if (type.empty())
    type = name;
}

LocalEnergyInput::LocalEnergyInput(xmlNodePtr cur)
{
  input_section_.readXML(cur);
  readNameTypeLikeLegacy(input_section_, name_, type_);
  auto setIfInInput = LAMBDA_setIfInInput;
  setIfInInput(use_hdf5_, "hdf5");
  if (lowerCase(type_) == "elocal")
    type_ = "localenergy";
}

CSLocalEnergyInput::CSLocalEnergyInput(xmlNodePtr cur)
{
  input_section_.readXML(cur);
  readNameTypeLikeLegacy(input_section_, name_, type_);
  auto setIfInInput = LAMBDA_setIfInInput;
  setIfInInput(n_psi_, "npsi");
}

RMCLocalEnergyInput::RMCLocalEnergyInput(xmlNodePtr cur)
{
  input_section_.readXML(cur);
  readNameTypeLikeLegacy(input_section_, name_, type_);
  auto setIfInInput = LAMBDA_setIfInInput;
  setIfInInput(n_obs_, "nobs");
}

} // namespace qmcplusplus
