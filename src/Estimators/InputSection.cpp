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

#include "InputSection.h"
#include "Message/UniformCommunicateError.h"
#include "Utilities/string_utils.h"
namespace qmcplusplus
{

void InputSection::readXML(xmlNodePtr cur)
{
  // read attributes
  xmlAttrPtr att = cur->properties;
  while (att != NULL)
  {
    // unsafe att->name is an xmlChar, xmlChar is a UTF-8 byte
    std::string name{lowerCase(castXMLCharToChar(att->name))};
    if (!isAttribute(name))
    {
      std::stringstream error;
      error << "InputSection::readXML name " << name << " is not an attribute of " << section_name << "\n";
      throw UniformCommunicateError(error.str());
    }
    std::istringstream stream(castXMLCharToChar(att->children->content));
    setFromStream(name, stream);
    att = att->next;
  }

  // read parameters
  xmlNodePtr element = cur->xmlChildrenNode;
  while (element != NULL)
  {
    std::string ename{lowerCase(castXMLCharToChar(element->name))};
    if (ename == "parameter")
    {
      std::string name(lowerCase(getXMLAttributeValue(element, "name")));
      if (!isParameter(name))
      {
        std::stringstream error;
        error << "InputSection::readXML name " << name << " is not a parameter of " << section_name << "\n";
        throw UniformCommunicateError(error.str());
      }
      std::istringstream stream(XMLNodeString{element});
      setFromStream(name, stream);
    }
    else if (isDelegate(ename))
    {}
    element = element->next;
  }

  // assign default values for optional variables
  setDefaults();

  // check input validity
  checkValid();
  //report();
}


void InputSection::init(const std::unordered_map<std::string, std::any>& init_values)
{
  // assign inputted values
  for (auto& [name, value] : init_values)
    setFromValue(name, value);

  // assign default values for optional variables
  setDefaults();

  // check input validity
  checkValid();
  //report();
}


void InputSection::setDefaults()
{
  for (auto& [name, default_value] : default_values)
    if (!has(name))
      setFromValue(name, default_value);
}

void InputSection::setFromStream(const std::string& name, std::istringstream& svalue)
{
  if (isString(name) || isEnumString(name))
  {
    std::string value;
    svalue >> value;
    values_[name] = value;
  }
  else if (isMultiString(name))
  {
    std::vector<std::string> string_values;
    for (std::string value; svalue >> value;)
      string_values.push_back(value);
    values_[name] = string_values;
  }
  else if (isBool(name))
  {
    std::string sval;
    svalue >> sval;
    bool value    = sval == "yes" || sval == "true" || sval == "1";
    values_[name] = value;
  }
  else if (isInteger(name))
  {
    int value;
    svalue >> value;
    values_[name] = value;
  }
  else if (isReal(name))
  {
    Real value;
    svalue >> value;
    values_[name] = value;
  }
  else if (isPosition(name))
  {
    Position value;
    svalue >> value;
    values_[name] = value;
  }
  else
  {
    std::stringstream error;
    error << "InputSection::set_from_stream name " << name << " in " << section_name << " does not have a type\n";
    throw UniformCommunicateError(error.str());
  }
}

template<typename T>
void InputSection::setFromValue(const std::string& name, const T& value)
{
  if (isString(name) || isEnumString(name))
    values_[name] = std::any_cast<std::string>(value);
  else if (isMultiString(name))
    values_[name] = (std::any_cast<std::vector<std::string>>(value));
  else if (isBool(name))
    values_[name] = std::any_cast<bool>(value);
  else if (isInteger(name))
    values_[name] = std::any_cast<int>(value);
  else if (isReal(name))
    values_[name] = std::any_cast<Real>(value);
  else if (isPosition(name))
    values_[name] = std::any_cast<Position>(value);
  else
  {
    std::stringstream error;
    error << "InputSection::set_from_value name " << name << " in " << section_name << " does not have a type\n";
    throw UniformCommunicateError(error.str());
  }
}

void InputSection::checkValid()
{
  // check that all required inputs are present
  for (auto& name : required)
    if (!has(name))
    {
      std::stringstream error;
      error << "InputSection::check_valid required variable " << name << " in " << section_name
            << " has not been assigned\n";
      throw UniformCommunicateError(error.str());
    }
  checkParticularValidity();
};


void InputSection::report() const
{
  auto& out = app_log();
  out << "\n" << section_name;
  for (auto& [name, value] : values_)
  {
    out << "\n  " << name << " = ";
    if (isString(name))
      out << std::any_cast<std::string>(value);
    else if (isBool(name))
      out << std::any_cast<bool>(value);
    else if (isInteger(name))
      out << std::any_cast<int>(value);
    else if (isReal(name))
      out << std::any_cast<Real>(value);
  }
  out << "\n\n";
}

std::any InputSection::lookupAnyEnum(const std::string& enum_name,
                                     const std::string& enum_value,
                                     const std::unordered_map<std::string, std::any>& enum_map)
{
  std::string enum_value_str(lowerCase(enum_name + "-" + enum_value));
  try
  {
    return enum_map.at(enum_value_str);
  }
  catch (std::out_of_range& oor_exc)
  {
    std::throw_with_nested(std::logic_error("bad_enum_tag_value: " + enum_value_str));
  }
}

} // namespace qmcplusplus
