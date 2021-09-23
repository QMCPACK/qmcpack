//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2021 QMCPACK developers.
//
// File developed by: Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//
//////////////////////////////////////////////////////////////////////////////////////
#ifndef QMCPLUSPLUS_INPUTSECTION_H
#define QMCPLUSPLUS_INPUTSECTION_H

#include <typeinfo>
#include <any>
#include <unordered_set>
#include <unordered_map>

#include "Configuration.h"
#include "OhmmsData/ParameterSet.h"
#include "Containers/OhmmsPETE/TinyVector.h"

#include "Message/UniformCommunicateError.h"


namespace qmcplusplus
{


class InputSection
{
protected:
  using Real               = QMCTraits::RealType;
  using PosType            = QMCTraits::PosType;
  static constexpr int DIM = QMCTraits::DIM;

  // Internal data below comprise the input specification.
  //   Most apply attributes to input variables.
  //   Enables minimal listing of variable classification and default values in derived classes.
  //   Expand later to include allowed_values for input correctness checking
  std::string section_name; // name of the input section

  std::unordered_set<std::string> attributes; // list of attribute variables
  std::unordered_set<std::string> parameters; // list of parameter variables
  std::unordered_set<std::string> required;   // list of required variables

  std::unordered_set<std::string> strings;  // list of string variables
  std::unordered_set<std::string> bools;    // list of boolean variables
  std::unordered_set<std::string> integers; // list of integer variables
  std::unordered_set<std::string> reals;    // list of real variables

  std::unordered_map<std::string, std::any> default_values; // default values for optional variables

private:
  // Storage for variable values read from XML, etc.
  std::unordered_map<std::string, std::any> values;

public:
  // Query if a variable has been set
  bool has(const std::string& name) const { return values.find(name) != values.end(); };

  // Enable read-only access to variable values.
  //   Needs updating to allow copy-less return.
  template<typename T>
  T get(const std::string& name) const
  {
    return std::any_cast<T>(values.at(name));
  };

  // Read variable values (initialize) from XML input.
  //   Later, this should call a correctness checking function and enforce immutability.
  //   (required values are all provided, provided values fall in allowed ranges)
  void readXML(xmlNodePtr cur);

  // Initialize from unordered_map/initializer list
  void init(const std::unordered_map<std::string, std::any>& init_values);

private:
  // Query functions
  bool is_attribute(const std::string& name) const { return attributes.find(name) != attributes.end(); };
  bool is_parameter(const std::string& name) const { return parameters.find(name) != parameters.end(); };
  bool is_required(const std::string& name) const { return required.find(name) != required.end(); };

  bool is_string(const std::string& name) const { return strings.find(name) != strings.end(); };
  bool is_bool(const std::string& name) const { return bools.find(name) != bools.end(); };
  bool is_integer(const std::string& name) const { return integers.find(name) != integers.end(); };
  bool is_real(const std::string& name) const { return reals.find(name) != reals.end(); };

  bool has_default(const std::string& name) const { return default_values.find(name) != default_values.end(); };

  // Set default values for optional inputs.
  void set_defaults();

  // Perform typed read and assignment of input variables from strings
  void set_from_stream(const std::string& name, std::istringstream& svalue);

  // Perform typed assignment of input variables from intrinsic types
  template<typename T>
  void set_from_value(const std::string& name, const T& svalue);

  // Check validity of inputs
  void check_valid();


  // Simple write of contents.
  //   Developer/debugging function of limited value.
  //   May be removed at any time.
  void report() const;
};


} // namespace qmcplusplus
#endif /* INPUTSECTION_H */
