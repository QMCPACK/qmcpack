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
      attributes   = {"name", "type", "hdf5"};
      bools        = {"hdf5"};
      strings      = {"name","type"};
    };
  };
public:
  using Consumer = LocalEnergyEstimator;
  LocalEnergyInput() = default;
  LocalEnergyInput(xmlNodePtr cur);
  bool get_use_hdf5() const { return use_hdf5_; }
  const std::string& get_name() { return name_; }
private:
  std::string name_;
  std::string type_;
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
      attributes   = {"npsi", "type"};
      integers     = {"npsi"};
      strings      = {"type"};
    }
  };
public:
  using Consumer = CSEnergyEstimator;
  CSLocalEnergyInput() = default;
  CSLocalEnergyInput(xmlNodePtr cur);
  int get_n_psi() const { return n_psi_; }
  const std::string& get_name() { return name_; }
private:
  std::string name_;
  std::string type_;
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
      attributes   = {"npsi", "type"};
      integers     = {"npsi"};
      strings      = {"type"};
    }
  };
public:
  using Consumer = RMCLocalEnergyEstimator;
  RMCLocalEnergyInput() = default;
  RMCLocalEnergyInput(xmlNodePtr cur);
  int get_n_psi() const { return n_psi_; }
  const std::string& get_name() { return name_; }
private:
  std::string name_;
  std::string type_;
  RMCLocalEnergyInputSection input_section_;
  int n_psi_ = 1;
};


} // namespace qmcplusplus

#endif
