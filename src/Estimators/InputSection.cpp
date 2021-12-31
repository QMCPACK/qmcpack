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
    if (!is_attribute(name))
    {
      std::stringstream error;
      error << "InputSection::readXML name " << name << " is not an attribute of " << section_name << "\n";
      throw UniformCommunicateError(error.str());
    }
    std::istringstream stream(castXMLCharToChar(att->children->content));
    set_from_stream(name, stream);
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
      if (!is_parameter(name))
      {
        std::stringstream error;
        error << "InputSection::readXML name " << name << " is not a parameter of " << section_name << "\n";
        throw UniformCommunicateError(error.str());
      }
      std::istringstream stream(XMLNodeString{element});
      set_from_stream(name, stream);
    }
    element = element->next;
  }

  // assign default values for optional variables
  set_defaults();

  // check input validity
  check_valid();
  //report();
}


void InputSection::init(const std::unordered_map<std::string, std::any>& init_values)
{
  // assign inputted values
  for (auto& [name, value] : init_values)
    set_from_value(name, value);

  // assign default values for optional variables
  set_defaults();

  // check input validity
  check_valid();
  //report();
}


void InputSection::set_defaults()
{
  for (auto& [name, default_value] : default_values)
    if (!has(name))
      set_from_value(name, default_value);
}

void InputSection::set_from_stream(const std::string& name, std::istringstream& svalue)
{
  if (is_string(name) || is_enum_string(name))
  {
    std::string value;
    svalue >> value;
    values[name] = value;
  }
  else if (is_multi_string(name))
  {
    std::vector<std::string> string_values;
    for (std::string value; svalue >> value;)
      string_values.push_back(value);
    values[name] = string_values;
  }
  else if (is_bool(name))
  {
    std::string sval;
    svalue >> sval;
    bool value   = sval == "yes" || sval == "true" || sval == "1";
    values[name] = value;
  }
  else if (is_integer(name))
  {
    int value;
    svalue >> value;
    values[name] = value;
  }
  else if (is_real(name))
  {
    Real value;
    svalue >> value;
    values[name] = value;
  }
  else if (is_position(name))
  {
    Position value;
    svalue >> value;
    values[name] = value;
  }
  else
  {
    std::stringstream error;
    error << "InputSection::set_from_stream name " << name << " in " << section_name << " does not have a type\n";
    throw UniformCommunicateError(error.str());
  }
}

template<typename T>
void InputSection::set_from_value(const std::string& name, const T& value)
{
  if (is_string(name) || is_enum_string(name))
    values[name] = std::any_cast<std::string>(value);
  else if (is_multi_string(name))
    values[name] = (std::any_cast<std::vector<std::string>>(value));
  else if (is_bool(name))
    values[name] = std::any_cast<bool>(value);
  else if (is_integer(name))
    values[name] = std::any_cast<int>(value);
  else if (is_real(name))
    values[name] = std::any_cast<Real>(value);
  else if (is_position(name))
    values[name] = std::any_cast<Position>(value);
  else
  {
    std::stringstream error;
    error << "InputSection::set_from_value name " << name << " in " << section_name << " does not have a type\n";
    throw UniformCommunicateError(error.str());
  }
}

void InputSection::check_valid()
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
  for (auto& [name, value] : values)
  {
    out << "\n  " << name << " = ";
    if (is_string(name))
      out << std::any_cast<std::string>(value);
    else if (is_bool(name))
      out << std::any_cast<bool>(value);
    else if (is_integer(name))
      out << std::any_cast<int>(value);
    else if (is_real(name))
      out << std::any_cast<Real>(value);
  }
  out << "\n\n";
}

std::any InputSection::lookupAnyEnum(const std::string& enum_name, const std::string& enum_value, const std::unordered_map<std::string, std::any>& enum_map)
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
