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
  using Real     = QMCTraits::RealType;
  using Position = QMCTraits::PosType;

protected:
  // Internal data below comprise the input specification.
  //   Most apply attributes to input variables.
  //   Enables minimal listing of variable classification and default values in derived classes.
  //   Expand later to include allowed_values for input correctness checking

  // Becuase it hurts to read all the trailing _ in the constructors of input section subtypes
  // NOLINTBEGIN(readability-indentifier-naming)

  std::string section_name; // name of the input section

  std::unordered_set<std::string> attributes; // list of attribute variables
  std::unordered_set<std::string> parameters; // list of parameter variables
  std::unordered_set<std::string> delegates;
  std::unordered_set<std::string> required; // list of required variables

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
  // NOLINTEND(readability-indentifier-naming)

private:
  // Storage for variable values read from XML, etc.
  std::unordered_map<std::string, std::any> values_;

public:
  // Query if a variable has been set
  bool has(const std::string& name) const { return values_.find(name) != values_.end(); }

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
      return std::any_cast<T>(values_.at(name));
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
  static std::string reverseLookupInputEnumMap(ENUM_T enum_val,
                                               const std::unordered_map<std::string, std::any>& enum_map)
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
  [[noreturn]] virtual std::any assignAnyEnum(const std::string& tag) const
  {
    throw std::runtime_error("derived class must provide assignAnyEnum method if enum parameters are used");
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
  static std::any lookupAnyEnum(const std::string& enum_name,
                                const std::string& enum_value,
                                const std::unordered_map<std::string, std::any>& enum_map);

private:
  // Query functions
  bool isAttribute(const std::string& name) const { return attributes.find(name) != attributes.end(); }
  bool isParameter(const std::string& name) const { return parameters.find(name) != parameters.end(); }
  bool isDelegate(const std::string& name) const { return delegates.find(name) != delegates.end(); }
  bool isRequired(const std::string& name) const { return required.find(name) != required.end(); }

  bool isEnumString(const std::string& name) const { return enums.find(name) != enums.end(); }
  bool isString(const std::string& name) const { return strings.find(name) != strings.end(); }
  bool isMultiString(const std::string& name) const { return multi_strings.find(name) != multi_strings.end(); }
  bool isBool(const std::string& name) const { return bools.find(name) != bools.end(); }
  bool isInteger(const std::string& name) const { return integers.find(name) != integers.end(); }
  bool isReal(const std::string& name) const { return reals.find(name) != reals.end(); }
  bool isPosition(const std::string& name) const { return positions.find(name) != positions.end(); }
  bool has_default(const std::string& name) const { return default_values.find(name) != default_values.end(); }

  // Set default values for optional inputs.
  void setDefaults();

  // Perform typed read and assignment of input variables from strings
  void setFromStream(const std::string& name, std::istringstream& svalue);

  // Perform typed assignment of input variables from intrinsic types
  template<typename T>
  void setFromValue(const std::string& name, const T& svalue);

  // Check validity of inputs
  void checkValid();


  // Simple write of contents.
  //   Developer/debugging function of limited value.
  //   May be removed at any time.
  void report() const;
};


} // namespace qmcplusplus
#endif /* INPUTSECTION_H */
