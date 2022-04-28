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

#include "EstimatorInput.h"
#include "ScalarEstimatorInputs.h"

namespace qmcplusplus
{

LocalEnergyInput::LocalEnergyInput(xmlNodePtr cur)
{
  input_section_.readXML(cur);
  auto setIfInInput = LAMBDA_setIfInInput;
  setIfInInput(use_hdf5_, "hdf5");
  setIfInInput(type_, "type");
  setIfInInput(name_, "name");
  if (name_.empty())
    name_ = type_;
  if (type_.empty())
    type_ = name_;
  // Support for legacy alias for localenergy
  if(lowerCase(type_) == "elocal")
    type_ = "localenergy";
}

CSLocalEnergyInput::CSLocalEnergyInput(xmlNodePtr cur)
{
  input_section_.readXML(cur);
  auto setIfInInput = LAMBDA_setIfInInput;
  setIfInInput(n_psi_, "npsi");
  setIfInInput(type_, "type");
  setIfInInput(name_, "name");
  if (name_.empty())
    name_ = type_;
  if (type_.empty())
    type_ = name_;
}

RMCLocalEnergyInput::RMCLocalEnergyInput(xmlNodePtr cur)
{
  input_section_.readXML(cur);
  auto setIfInInput = LAMBDA_setIfInInput;
  setIfInInput(n_psi_, "npsi");
  setIfInInput(type_, "type");
  setIfInInput(name_, "name");
  if (name_.empty())
    name_ = type_;
  if (type_.empty())
    type_ = name_;
}

}
