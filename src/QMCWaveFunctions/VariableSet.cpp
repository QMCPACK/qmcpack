//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


#include "VariableSet.h"
#include "io/hdf/hdf_archive.h"
#include "Host/sysutil.h"
#include <map>
#include <stdexcept>
#include <iomanip>
#include <ios>
#include <algorithm>

using std::setw;

namespace optimize
{
void VariableSet::clear()
{
  num_active_vars = 0;
  Index.clear();
  NameAndValue.clear();
  ParameterType.clear();
}

void VariableSet::insertFrom(const VariableSet& input)
{
  for (int i = 0; i < input.size(); ++i)
  {
    iterator loc = find(input.name(i));
    if (loc == NameAndValue.end())
    {
      Index.push_back(input.Index[i]);
      NameAndValue.push_back(input.NameAndValue[i]);
      ParameterType.push_back(input.ParameterType[i]);
    }
    else
      (*loc).second = input.NameAndValue[i].second;
  }
  num_active_vars = input.num_active_vars;
}

void VariableSet::removeInactive()
{
  std::vector<int> valid(Index);
  std::vector<pair_type> acopy(NameAndValue);
  std::vector<index_pair_type> ccopy(ParameterType);
  num_active_vars = 0;
  Index.clear();
  NameAndValue.clear();
  ParameterType.clear();
  for (int i = 0; i < valid.size(); ++i)
  {
    if (valid[i] > -1)
    {
      Index.push_back(num_active_vars++);
      NameAndValue.push_back(acopy[i]);
      ParameterType.push_back(ccopy[i]);
    }
  }
}

void VariableSet::resetIndex()
{
  num_active_vars = 0;
  for (int i = 0; i < Index.size(); ++i)
  {
    Index[i] = (Index[i] < 0) ? -1 : num_active_vars++;
  }
}

void VariableSet::getIndex(const VariableSet& selected)
{
  num_active_vars = 0;
  for (int i = 0; i < NameAndValue.size(); ++i)
  {
    Index[i] = selected.getIndex(NameAndValue[i].first);
    if (Index[i] >= 0)
      num_active_vars++;
  }
}

int VariableSet::getIndex(const std::string& vname) const
{
  int loc = 0;
  while (loc != NameAndValue.size())
  {
    if (NameAndValue[loc].first == vname)
      return Index[loc];
    ++loc;
  }
  return -1;
}

void VariableSet::setIndexDefault()
{
  for (int i = 0; i < Index.size(); ++i)
    Index[i] = i;
}

void VariableSet::print(std::ostream& os, int leftPadSpaces, bool printHeader) const
{
  std::string pad_str = std::string(leftPadSpaces, ' ');
  int max_name_len    = 0;
  if (NameAndValue.size() > 0)
    max_name_len =
        std::max_element(NameAndValue.begin(), NameAndValue.end(), [](const pair_type& e1, const pair_type& e2) {
          return e1.first.length() < e2.first.length();
        })->first.length();

  int max_value_len     = 28; // 6 for the precision and 7 for minus sign, leading value, period, and exponent.
  int max_type_len      = 1;
  int max_use_len       = 3;
  int max_index_len     = 1;
  if (printHeader)
  {
    max_name_len      = std::max(max_name_len, 4); // size of "Name" header
    max_type_len      = 4;
    max_index_len     = 5;
    os << pad_str << setw(max_name_len) << "Name"
       << " " << setw(max_value_len) << "Value"
       << " " << setw(max_type_len) << "Type"
       << " " << setw(max_use_len) << "Use"
       << " " << setw(max_index_len) << "Index" << std::endl;
    os << pad_str << std::setfill('-') << setw(max_name_len) << ""
       << " " << setw(max_value_len) << ""
       << " " << setw(max_type_len) << ""
       << " " << setw(max_use_len) << ""
       << " " << setw(max_index_len) << "" << std::endl;
    os << std::setfill(' ');
  }

  for (int i = 0; i < NameAndValue.size(); ++i)
  {
    os << pad_str << setw(max_name_len) << NameAndValue[i].first << " " << std::setprecision(6) << std::scientific
       << setw(max_value_len) << NameAndValue[i].second << " " << setw(max_type_len) << ParameterType[i].second << " "
       << std::defaultfloat;

    if (Index[i] < 0)
      os << setw(max_use_len) << "OFF" << std::endl;
    else
      os << setw(max_use_len) << "ON"
         << " " << setw(max_index_len) << Index[i] << std::endl;
  }
}

void VariableSet::writeToHDF(const std::string& filename, qmcplusplus::hdf_archive& hout) const
{
  hout.create(filename);

  // File Versioning
  // 1.0.0  Initial file version
  // 1.1.0  Files could have object-specific data from OptimizableObject::read/writeVariationalParameters
  std::vector<int> vp_file_version{1, 1, 0};
  hout.write(vp_file_version, "version");

  std::string timestamp(getDateAndTime("%Y-%m-%d %H:%M:%S %Z"));
  hout.write(timestamp, "timestamp");

  hout.push("name_value_lists");

  std::vector<qmcplusplus::QMCTraits::RealType> param_values;
  std::vector<std::string> param_names;
  for (auto& [name, value] : NameAndValue)
  {
    param_names.push_back(name);
    param_values.push_back(value);
  }

  hout.write(param_names, "parameter_names");
  hout.write(param_values, "parameter_values");
  hout.pop();
}

void VariableSet::readFromHDF(const std::string& filename, qmcplusplus::hdf_archive& hin)
{
  if (!hin.open(filename, H5F_ACC_RDONLY))
  {
    std::ostringstream err_msg;
    err_msg << "Unable to open VP file: " << filename;
    throw std::runtime_error(err_msg.str());
  }

  try
  {
    hin.push("name_value_lists", false);
  }
  catch (std::runtime_error&)
  {
    std::ostringstream err_msg;
    err_msg << "The group name_value_lists in not present in file: " << filename;
    throw std::runtime_error(err_msg.str());
  }

  std::vector<qmcplusplus::QMCTraits::RealType> param_values;
  hin.read(param_values, "parameter_values");

  std::vector<std::string> param_names;
  hin.read(param_names, "parameter_names");

  for (int i = 0; i < param_names.size(); i++)
  {
    std::string& vp_name = param_names[i];
    // Find and set values by name.
    // Values that are not present do not get added.
    if (find(vp_name) != end())
      (*this)[vp_name] = param_values[i];
  }

  hin.pop();
}


} // namespace optimize
