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

#include "Configuration.h"
#include "InputSection.h"

namespace qmcplusplus
{

class ELocalInput
{
  class ELocalInputSection : public InputSection
  {
  public:
    ELocalInputSection()
    {
      section_name = "ELocal";
      attributes   = {"use_hdf5"};
      bools        = {"use_hdf5"};
    };
  };

  ELocalInput(xmlNodePtr cur);
  bool get_use_hdf5() const { return use_hdf5_; }
private:
  ELocalInputSection input_section_;
  bool use_hdf5_ = true;
};

class CSLocalEnergyInput
{
  class CSLocalEnergyInputSection : public InputSection
  {
  public:
    CSLocalEnergyInputSection()
    {
      section_name = "CSLocalEnergy";
      attributes   = {"npsi"};
      integers     = {"nspi"};
    }
  };
  CSLocalEnergyInput(xmlNodePtr cur);
  int get_n_psi() const { return n_psi_; }

private:
  CSLocalEnergyInputSection input_section_;
  int n_psi_ = 1;
};

} // namespace qmcplusplus
