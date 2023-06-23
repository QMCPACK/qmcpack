//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2022 QMCPACK developers.
//
// File developed by: Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//                    Peter W. Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//
//////////////////////////////////////////////////////////////////////////////////////

#include "InputSection.h"
#include "Message/UniformCommunicateError.h"
#include "ModernStringUtils.hpp"
#include "Utilities/string_utils.h"

namespace qmcplusplus
{

void InputSection::readAttributes(xmlNodePtr cur,
                                  const std::string& element_name,
                                  const std::vector<std::string>& do_not_consume)
{
  xmlAttrPtr att = cur->properties;
  while (att != NULL)
  {
    // unsafe att->name is an xmlChar, xmlChar is a UTF-8 byte
    std::string name{lowerCase(castXMLCharToChar(att->name))};
    // issue here is that we don't want to consume the name of the parameter as that has a special status in a parameter tag.
    // This is due to the <parameter name="parameter_name>  == <parameter_name> tag equivalence :(
    std::string qualified_name{((element_name.size() > 0 && element_name != "parameter") ? (element_name + "::") : "") +
                               name};

    // Since parameters don't get a qualified name this will still prevent consumption of any of the do_not_consume attributes from them
    if (std::any_of(do_not_consume.begin(), do_not_consume.end(), [&name](auto& dnc) { return name == dnc; }))
    {
      att = att->next;
      continue;
    }

    if (!isAttribute(qualified_name))
    {
      std::stringstream error;
      error << "InputSection::readXML name " << name << " is not an attribute of " << section_name << "\n";
      throw UniformCommunicateError(error.str());
    }
    std::istringstream stream(castXMLCharToChar(att->children->content));
    if (isCustom(name))
      setFromStreamCustom(element_name, qualified_name, stream);
    else
      setFromStream(qualified_name, stream);

    att = att->next;
  }
}

void InputSection::handleDelegate(const std::string& ename, const xmlNodePtr element)
{
  assert(delegate_factories_.find(ename) != delegate_factories_.end());
  std::string value_key;
  std::any value = delegate_factories_[ename](element, value_key);

  assignValue(value_key, value);
}

void InputSection::readXML(xmlNodePtr cur)
{
  // For historical reasons that actual "type" of the element/input section is expressed in a very inconsistent way.
  // It could be coded via the element name i.e. the tag, or at minimum a method, type, or name attribute.
  std::string section_ename{lowerCase(castXMLCharToChar(cur->name))};
  std::string section_method(lowerCase(getXMLAttributeValue(cur, "method")));
  std::string section_type(lowerCase(getXMLAttributeValue(cur, "type")));
  std::string section_name_actual(lowerCase(getXMLAttributeValue(cur, "name")));
  // at anyrate one of these must match the section_name.
  std::string lcase_section_name{lowerCase(section_name)};

  auto checkSectionName = [&section_name = section_name, &section_name_alternates = section_name_alternates](auto& possible_sname){
    std::string lcase_section_name{lowerCase(section_name)};
    if (possible_sname == lcase_section_name)
      return true;
    if (section_name_alternates.size() > 0)
      return std::any_of(section_name_alternates.begin(), section_name_alternates.end(),
                      [&possible_sname](auto& name_alternate) {
                        std::string lcase_alternate{lowerCase(name_alternate)};
                        return possible_sname == lcase_alternate;
                      });
    return false;
  };

  if (!(checkSectionName(section_ename) || checkSectionName(section_method) ||
        checkSectionName(section_type) || checkSectionName(section_name_actual)))
    throw UniformCommunicateError("Input is invalid  " + lcase_section_name + " does not match input node!");

  // these attributes don't get an element name passed to them because by convention we save and define them unqualified.
  readAttributes(cur, "", {"type"});
  // read parameters
  xmlNodePtr element = cur->xmlChildrenNode;
  while (element != NULL)
  {
    // ename is the "name" of the XML element i.e. <ename [attributes]>
    std::string ename{lowerCase(castXMLCharToChar(element->name))};
    // value of the elements attribute = name
    std::string name(lowerCase(getXMLAttributeValue(element, "name")));
    // we need both ename and name to figure out how to handle an element because of the <parameter name="actual_parameter"> pattern.
    // If the name attribute isn't there it should equal the element name.
    if (name.size() < 1)
      name = ename;

    if (isDelegate(ename))
    {
      handleDelegate(ename, element);
    }
    else if (isCustom(ename))
    {
      std::istringstream stream(XMLNodeString{element});
      setFromStreamCustom(ename, name, stream);
    }
    else if (ename == "parameter" || isParameter(ename))
    {
      if (ename == "parameter")
        ename = name;
      else
        // We do this because the semantics of parameters are such that name can not have a unique value because
        // if it did have one that the equivalence of <parameter_name> and <parameter name="parameter_name"> would
        // be broken. Input code being ported could depend on this an an invariant.
        name = ename;
      if (!isParameter(ename))
      {
        std::stringstream error;
        error << "InputSection::readXML name " << name << " is not a parameter of " << section_name << "\n";
        throw UniformCommunicateError(error.str());
      }

      std::istringstream stream(XMLNodeString{element});
      setFromStream(name, stream);
      readAttributes(element, name, {"name"});
    }
    else if (ename != "text")
    {
      std::stringstream error;
      error << "InputSection::readXML node name " << ename << " is not handled by InputSection subtype " << section_name
            << "\n";
      throw UniformCommunicateError(error.str());
    }
    // else can't be an error case because this is how the whitespace text nodes
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

void InputSection::registerDelegate(const std::string& tag,
                                    std::function<std::any(xmlNodePtr cur, std::string& value_key)> factory)
{
  delegate_factories_[tag] = factory;
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
    assignValue(name, value);
  }
  else if (isMultiString(name))
  {
    std::vector<std::string> string_values;
    for (std::string value; svalue >> value;)
      string_values.push_back(value);
    assignValue(name, string_values);
  }
  else if (isMultiReal(name))
  {
    std::vector<Real> real_values;
    for (FullPrecReal value; svalue >> value;)
      real_values.push_back(static_cast<Real>(value));
    assignValue(name, real_values);
  }
  else if (isBool(name))
  {
    std::string sval;
    svalue >> sval;
    bool value = sval == "yes" || sval == "true" || sval == "1";
    assignValue(name, value);
  }
  else if (isInteger(name))
  {
    int value;
    svalue >> value;
    assignValue(name, value);
  }
  else if (isReal(name))
  {
    FullPrecReal value;
    svalue >> value;
    assignValue(name, Real(value));
  }
  else if (isPosition(name))
  {
    Position value;
    svalue >> value;
    assignValue(name, value);
  }
  else
  {
    std::stringstream error;
    error << "InputSection::set_from_stream name " << name << " in " << section_name << " does not have a type\n";
    throw UniformCommunicateError(error.str());
  }
}

template<typename T>
void InputSection::assignValue(const std::string& name, const T& value)
{
  if (has(name) && !isMultiple(name))
    throw UniformCommunicateError("Input is invalid  " + section_name + " contains " + name +
                                  " node with duplicate name!");

  if (!isMultiple(name))
    values_[name] = value;
  else
  {
    if (has(name))
      std::any_cast<std::vector<T>>(values_[name]).push_back(value);
    else
      values_[name] = std::vector<T>{value};
  }
}

void InputSection::setFromValue(const std::string& name, const std::any& value)
{
  try
  {
    if (isString(name) || isEnumString(name))
      assignValue(name, std::any_cast<std::string>(value));
    else if (isMultiString(name))
      assignValue(name, std::any_cast<std::vector<std::string>>(value));
    else if (isMultiReal(name))
      assignValue(name, std::any_cast<std::vector<std::string>>(value));
    else if (isBool(name))
      assignValue(name, std::any_cast<bool>(value));
    else if (isInteger(name))
      assignValue(name, std::any_cast<int>(value));
    else if (isReal(name))
      assignValue(name, std::any_cast<Real>(value));
    else if (isPosition(name))
      assignValue(name, std::any_cast<Position>(value));
    else
    {
      std::stringstream error;
      error << "InputSection::set_from_value name " << name << " in " << section_name << " does not have a type\n";
      throw UniformCommunicateError(error.str());
    }
  }
  catch (const std::bad_cast& exc)
  {
    std::throw_with_nested(UniformCommunicateError("std::any_cast failed in setFromValue for name:" + name));
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

void InputSection::report(std::ostream& out) const
{
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
    std::throw_with_nested(UniformCommunicateError("bad_enum_tag_value: " + enum_value_str));
  }
}

} // namespace qmcplusplus
