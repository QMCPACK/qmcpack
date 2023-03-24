//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2023 QMCPACK developers.
//
// File developed by: Mark Dewing, mdewing@anl.gov, Argonne National Laboratory
//
// File created by: Mark Dewing, mdewing@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////


#include "OptimizableObject.h"
#include "VariableSet.h"
#include "io/hdf/hdf_archive.h"

/**@file OptimizableObject.cpp
 *@brief OptimizableObject implementation
 */

namespace qmcplusplus
{

using opt_variables_type = optimize::VariableSet;

void OptimizableObject::openHDFToSave(const std::filesystem::path& filename, hdf_archive& hout)
{
  hout.create(filename);
  std::vector<int> vp_file_version{2, 0, 0};
  hout.write(vp_file_version, "version");
}

int OptimizableObject::openHDFToRead(const std::filesystem::path& filename, hdf_archive& hin)
{
  if (!hin.open(filename, H5F_ACC_RDONLY))
  {
    std::ostringstream err_msg;
    err_msg << "Unable to open VP file: " << filename;
    throw std::runtime_error(err_msg.str());
  }
  std::vector<int> vp_file_version(3);
  hin.read(vp_file_version, "version");
  return vp_file_version[0];
}

void OptimizableObject::saveVariationalParameters(hdf_archive& hout, opt_variables_type& opt_vars)
{
  opt_variables_type dummy;
  checkInVariablesExclusive(dummy);
  hid_t grp = hout.push(hdf_path(name_) / "name_value_lists");

  std::vector<qmcplusplus::QMCTraits::ValueType> param_values;
  std::vector<std::string> param_names;
  for (auto pair_it = dummy.begin(); pair_it != dummy.end(); pair_it++)
  {
    const std::string& param_name = pair_it->first;
    param_names.push_back(param_name);
    //auto opt_it = opt_vars.find(param_name);
    //if (opt_it == opt_vars.end())
    //  throw std::runtime_error("Inconsistent parameter " + param_name + " in " + name_);
    //param_values.push_back(opt_it->second);
    param_values.push_back(pair_it->second);
  }

  hout.write(param_names, "parameter_names");
  hout.write(param_values, "parameter_values");
  hout.pop(); // component name in name_ / name_value_lists
}

void OptimizableObject::readVariationalParameters(hdf_archive& hin, opt_variables_type& opt_vars)
{
  hid_t grp = hin.push(hdf_path(name_) / "name_value_lists", false);
  if (grp < 0)
  {
    std::ostringstream err_msg;
    err_msg << "The group name_value_lists in not present in file";
    throw std::runtime_error(err_msg.str());
  }

  std::vector<qmcplusplus::QMCTraits::ValueType> param_values;
  hin.read(param_values, "parameter_values");

  std::vector<std::string> param_names;
  hin.read(param_names, "parameter_names");

  for (int i = 0; i < param_names.size(); i++)
  {
    std::string& vp_name = param_names[i];
    // Find and set values by name.
    // Values that are not present do not get added.
    opt_vars.insert(vp_name, param_values[i]);
  }

  hin.pop(); // component name in name_ / name_value_lists
}
} // namespace qmcplusplus
