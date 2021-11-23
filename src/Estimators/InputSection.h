//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2021 QMCPACK developers.
//
// File developed by: Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//                    Peter W. Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
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
public:
  using Real               = QMCTraits::RealType;
  using Position           = QMCTraits::PosType;
  static constexpr int DIM = QMCTraits::DIM;
protected:
  // Internal data below comprise the input specification.
  //   Most apply attributes to input variables.
  //   Enables minimal listing of variable classification and default values in derived classes.
  //   Expand later to include allowed_values for input correctness checking
  std::string section_name; // name of the input section

  std::unordered_set<std::string> attributes; // list of attribute variables
  std::unordered_set<std::string> parameters; // list of parameter variables
  std::unordered_set<std::string> required;   // list of required variables

  std::unordered_set<std::string> strings;       // list of string variables that can have one value
  std::unordered_set<std::string> multi_strings; // list of string variables that can one or more values
  std::unordered_set<std::string> bools;         // list of boolean variables
  std::unordered_set<std::string> integers;      // list of integer variables
  std::unordered_set<std::string> reals;         // list of real variables
  std::unordered_set<std::string> positions;     // list of position variables
  /** list of enum inputs which allow a finite set of strings to map to enum values
   *  The enum class types and values need only be known to IS subtypes
   */
  std::unordered_set<std::string> enums;
  std::unordered_map<std::string, std::any> default_values; // default values for optional variables

private:
  // Storage for variable values read from XML, etc.
  std::unordered_map<std::string, std::any> values;

public:
  // Query if a variable has been set
  bool has(const std::string& name) const { return values.find(name) != values.end(); }

  // Enable read-only access to variable values.
  //   Needs updating to allow copy-less return.
  template<typename T>
  T get(const std::string& name) const
  {
    if constexpr (std::is_enum<T>::value)
    {
      std::any any_enum = assignAnyEnum(name);
      return std::any_cast<T>(any_enum);
    }
    else
      return std::any_cast<T>(values.at(name));
  }

  /** set var if input section has read the tag
   *  \param[out] var  external variable to be set if tag was defined
   *  \param[in]  tag  string tag of value could be parameter or atttribute name
   *  \return          whether input section has tag
   *
   *  use this is you prefer have native c++ types for input class members
   *  as well as set default via native c++ declaration. See
   *  OneBodyDensityMatricesInput for example.
   */
  template<typename T>
  bool setIfInInput(T& var, const std::string& tag)
  {
    if (has(tag))
    {
      var = get<T>(tag);
      return true;
    }
    else
      return false;
  }

  // Read variable values (initialize) from XML input.
  //   Later, this should call a correctness checking function and enforce immutability.
  //   (required values are all provided, provided values fall in allowed ranges)
  void readXML(xmlNodePtr cur);

  // Initialize from unordered_map/initializer list
  void init(const std::unordered_map<std::string, std::any>& init_values);
  
/** Get string represtation of enum class type value from enum_val
 *  
 *  This is just a way to get around the lack of a bidirectional map type.
 */
template<typename ENUM_T>
static std::string reverseLookupInputEnumMap(ENUM_T enum_val, const std::unordered_map<std::string, std::any>& enum_map)
{
  std::string lookup_str = "not found";
  for (const auto& enum_node : enum_map)
  {
    if (enum_node.second.type() == typeid(decltype(enum_val)) &&
        enum_val == std::any_cast<decltype(enum_val)>(enum_node.second))
    {
      lookup_str = enum_node.first;
      break;
    }
  }
  return lookup_str;
}  

protected:
  /** Do validation for a particular subtype of InputSection
   *  Called by check_valid.
   *  Default implementation is noop
   */
  virtual void checkParticularValidity() {}
  /** Derived class overrides this to get proper assignment of scoped enum values.
   *
   *  In most cases all you'll need it to define the map and write:
   *    std::any DerivedInputSection::assignAnyEnum(const std::string& name) const
   *    {
   *      return lookupAnyEnum(name, get<std::string>(name), derived_input_lookup_enum);
   *    }
   *
   *  See test_InputSection.cpp and OneBodyDensityMatricesInput
   *  You really should do this if your input class has a finite set of string values for an input
   *  example: OneBodyDensityMatricesInput
   *
   *  can't be bothered then just define your enum option as a string.
   */
  virtual std::any assignAnyEnum(const std::string& tag) const
  {
    throw std::runtime_error("derived class must provide assignAnyEnum method if enum parameters are used");
    return std::any();
  }

  /** Assign any enum helper for InputSection derived class
   *  assumes enum lookup table of this form:
   *    inline static const std::unordered_map<std::string, std::any>
   *    lookup_input_enum_value{{"integrator-uniform_grid", Integrator::UNIFORM_GRID},
   *                            {"integrator-uniform", Integrator::UNIFORM},
   *                            {"integrator-density", Integrator::DENSITY},
   *                            {"evaluator-loop", Evaluator::LOOP},
   *                            {"evaluator-matrix", Evaluator::MATRIX}};
   */
  static std::any lookupAnyEnum(const std::string& enum_name, const std::string& enum_value, const std::unordered_map<std::string, std::any>& enum_map);

private:
  // Query functions
  bool is_attribute(const std::string& name) const { return attributes.find(name) != attributes.end(); }
  bool is_parameter(const std::string& name) const { return parameters.find(name) != parameters.end(); }
  bool is_required(const std::string& name) const { return required.find(name) != required.end(); }

  bool is_enum_string(const std::string& name) const { return enums.find(name) != enums.end(); }
  bool is_string(const std::string& name) const { return strings.find(name) != strings.end(); }
  bool is_multi_string(const std::string& name) const { return multi_strings.find(name) != multi_strings.end(); }
  bool is_bool(const std::string& name) const { return bools.find(name) != bools.end(); }
  bool is_integer(const std::string& name) const { return integers.find(name) != integers.end(); }
  bool is_real(const std::string& name) const { return reals.find(name) != reals.end(); }
  bool is_position(const std::string& name) const { return positions.find(name) != positions.end(); }
  bool has_default(const std::string& name) const { return default_values.find(name) != default_values.end(); }

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
