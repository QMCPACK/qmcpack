//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2022 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//
// File refactored from: EstimatorManagerNew.h
//////////////////////////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_SCALAR_ESTIMATOR_INPUTS_H
#define QMCPLUSPLUS_SCALAR_ESTIMATOR_INPUTS_H

#include "InputSection.h"

/** \file
 *  This file contains the input classes for the
 *  supported "main estimator" classes derived from ScalarEstimatorBase.
 *
 *  see: LocalEnergyEstimator.h, CSEnergyEstimator.h, RMCLocalEnergyEstimator.h
 *
 *  The legacy input format used the <estimator> name attribute instead of the type attribute
 *  to indicate the type of the scalar estimators. There would otherwise be no reason to
 *  retain the name attribute reading or native input variable.
 *
 *  The current implementation sets the native type_ or name_ equal in value to
 *  match if only name or attribute is present in input.  If both are present name is
 *  ignored and only type is used for type determination.  The name can be set to anything
 *  as is typical of name attributes in the input but at the moment that name is not
 *  used in output or anywhere else.
 *
 *  EstimatorManagerInput will delegate parsing of scalar estimator inputs to these
 *  classes based on the type of the estimator as determined by reading the type attribute in the
 *  raw input. We still fixup the name_ and type_ in the implementation to create input classes
 *  consistent with the expected "input state."
 */

namespace qmcplusplus
{

class LocalEnergyEstimator;

class LocalEnergyInput
{
  class LocalEnergyInputSection : public InputSection
  {
  public:
    LocalEnergyInputSection()
    {
      section_name = "LocalEnergy";
      section_name_alternates = {"elocal"};
      attributes   = {"name", "type", "hdf5"};
      bools        = {"hdf5"};
      strings      = {"name", "type"};
    };
  };

public:
  using Consumer = LocalEnergyEstimator;
  // required because legacy code using LocalEnergyEstimator does not use input object.
  // LocalEnergyInput must be trivially constructable.
  LocalEnergyInput() = default;
  LocalEnergyInput(xmlNodePtr cur);
  bool get_use_hdf5() const { return use_hdf5_; }
  const std::string& get_name() const { return name_; }
  const std::string& get_type() const { return type_; }

private:
  std::string name_;
  std::string type_ = "LocalEnergy";
  LocalEnergyInputSection input_section_;
  bool use_hdf5_ = true;
};

struct CSEnergyEstimator;

class CSLocalEnergyInput
{
  class CSLocalEnergyInputSection : public InputSection
  {
  public:
    CSLocalEnergyInputSection()
    {
      section_name = "CSLocalEnergy";
      attributes   = {"name", "npsi", "type"};
      integers     = {"npsi"};
      strings      = {"name", "type"};
    }
  };

public:
  using Consumer = CSEnergyEstimator;
  // required because legacy code using CSEnergyEstimator does not use input object.
  // So CSLocalEnergyInput must be trivially constructable.
  CSLocalEnergyInput() = default;
  CSLocalEnergyInput(xmlNodePtr cur);
  int get_n_psi() const { return n_psi_; }
  const std::string& get_name() const { return name_; }
  const std::string& get_type() const { return type_; }

private:
  std::string name_;
  std::string type_ = "CSLocalEnergy";
  CSLocalEnergyInputSection input_section_;
  int n_psi_ = 1;
};

class RMCLocalEnergyEstimator;

class RMCLocalEnergyInput
{
  class RMCLocalEnergyInputSection : public InputSection
  {
  public:
    RMCLocalEnergyInputSection()
    {
      section_name = "RMCLocalEnergy";
      section_name_alternates = {"rmc"};      
      attributes   = {"name", "nobs", "type"};
      integers     = {"nobs"};
      strings      = {"name", "type"};
    }
  };

public:
  using Consumer = RMCLocalEnergyEstimator;
  // required because legacy code using RMCLocalEnergyEstimator does not use input object.
  // So RMCLocalEnergyInput must be trivially constructable.

  RMCLocalEnergyInput() = default;
  RMCLocalEnergyInput(xmlNodePtr cur);
  int get_n_obs() const { return n_obs_; }
  const std::string& get_name() const { return name_; }
  const std::string& get_type() const { return type_; }

private:
  std::string name_;
  std::string type_ = "RMCLocalEnergy";
  RMCLocalEnergyInputSection input_section_;
  int n_obs_ = 1;
};

} // namespace qmcplusplus

#endif
