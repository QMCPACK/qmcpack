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
  std::string section_name;                    // name of the input section

  std::unordered_set<std::string> attributes;  // list of attribute variables
  std::unordered_set<std::string> parameters;  // list of parameter variables
  std::unordered_set<std::string> required;    // list of required variables

  std::unordered_set<std::string> strings;     // list of string variables
  std::unordered_set<std::string> bools;       // list of boolean variables
  std::unordered_set<std::string> integers;    // list of integer variables
  std::unordered_set<std::string> reals;       // list of real variables
  
  std::unordered_map<std::string, std::any> default_values; // default values for optional variables
  
private:
  // Storage for variable values read from XML, etc.
  std::unordered_map<std::string, std::any> values;
  
public:
  // Query if a variable has been set
  bool has(std::string name) const {return values.find(name)!=values.end();};

  // Enable read-only access to variable values.
  //   Needs updating to allow copy-less return.
  template<typename T>
  T get(std::string name) {return std::any_cast<T>(values[name]);};


  // Later replace this with initializer_list/unordered_map initializer.
  //   Enforce correctness checks and immutability/read-only access thereafter.
  //template<typename T>
  //void set(std::string name, T& value)
  //{
  //  const std::type_info& Ttype = typeid(T);
  //  bool correct_type = (Ttype == typeid(std::string)) & is_string(name);
  //  correct_type |= (Ttype == typeid(bool)) & is_bool(name);
  //  correct_type |= (Ttype == typeid(int))  & is_integer(name);
  //  correct_type |= (Ttype == typeid(Real)) & is_real(name);
  //  if (!correct_type)
  //    throw UniformCommunicateError("InputSection::set  incorrect type assignment attempted");
  //  values[name] = value;
  //};


  // Read variable values (initialize) from XML input.
  //   Later, this should call a correctness checking function and enforce immutability.
  //   (required values are all provided, provided values fall in allowed ranges)
  void readXML(xmlNodePtr cur);

  //  Simple write of contents.  Can be replaced/removed in any final implemenation.
  void report();

private:
  // Query functions
  bool is_attribute(std::string name) const {return attributes.find(name)!=attributes.end();};
  bool is_parameter(std::string name) const {return parameters.find(name)!=parameters.end();};
  bool is_required(std::string name) const {return required.find(name)!=required.end();};

  bool is_string(std::string name) const {return strings.find(name)!=strings.end();};
  bool is_bool(std::string name) const {return bools.find(name)!=bools.end();};
  bool is_integer(std::string name) const {return integers.find(name)!=integers.end();};
  bool is_real(std::string name) const {return reals.find(name)!=reals.end();};

  bool has_default(std::string name) const {return default_values.find(name)!=default_values.end();};

  // Perform typed read and assignment of input variables
  void set_from_stream(std::string name, std::istringstream& svalue);

  // Set default values for optional inputs.
  void set_defaults();

};



} // namespace qmcplusplus
#endif /* INPUTSECTION_H */
