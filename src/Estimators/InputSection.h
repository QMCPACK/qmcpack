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

public:
  using Real               = QMCTraits::RealType;
  using PosType            = QMCTraits::PosType;
  static constexpr int DIM = QMCTraits::DIM;

protected:
  std::string section_name;

  std::unordered_set<std::string> attributes;
  std::unordered_set<std::string> parameters;
  std::unordered_set<std::string> required;
  std::unordered_set<std::string> derived;

  std::unordered_set<std::string> strings;
  std::unordered_set<std::string> bools;
  std::unordered_set<std::string> integers;
  std::unordered_set<std::string> reals;
  
  std::unordered_map<std::string, std::any> default_values;
  std::unordered_map<std::string, std::any> values;
  
public:
  bool is_attribute(std::string name){return attributes.find(name)!=attributes.end();};
  bool is_parameter(std::string name){return parameters.find(name)!=parameters.end();};
  bool is_required(std::string name){return required.find(name)!=required.end();};
  bool is_derived(std::string name){return derived.find(name)!=derived.end();};

  bool is_string(std::string name){return strings.find(name)!=strings.end();};
  bool is_bool(std::string name){return bools.find(name)!=bools.end();};
  bool is_integer(std::string name){return integers.find(name)!=integers.end();};
  bool is_real(std::string name){return reals.find(name)!=reals.end();};

  bool has_default(std::string name){return default_values.find(name)!=default_values.end();};
  bool has(std::string name){return values.find(name)!=values.end();};


  template<typename T>
  T get(std::string name){return std::any_cast<T>(values[name]);};


  template<typename T>
  void set(std::string name, T& value)
  {
    const std::type_info& Ttype = typeid(T);
    bool correct_type = (Ttype == typeid(std::string)) & is_string(name);
    correct_type |= (Ttype == typeid(int)) & is_integer(name);
    correct_type |= (Ttype == typeid(Real)) & is_real(name);
    if (!correct_type)
      throw UniformCommunicateError("InputSection::set  incorrect type assignment attempted");
    values[name] = value;
  };


  void set_defaults()
  {
    for (auto &[name, default_value] : default_values)
      if (!has(name))
        values[name] = default_value;
  };


  void readXML(xmlNodePtr cur)
  {
    // read attributes
    xmlAttrPtr att = cur->properties;
    while (att != NULL)
    {
      std::string name = (const char*)(att->name);
      for(auto& c : name)
        c = tolower(c);
      if (!is_attribute(name))
      {
        std::stringstream error;
        error << "InputSection::readXML name "<<name<<" is not an attribute of "<<section_name<<"\n";
        throw UniformCommunicateError(error.str());
      }
      std::istringstream stream((const char*)(att->children->content));
      set_from_stream(name,stream);
      att = att->next;
    }

    // read parameters
    xmlNodePtr element = cur->xmlChildrenNode;
    while (element != NULL)
    {
      std::string ename((const char*)element->name);
      for(auto& c : ename)
        c = tolower(c);
      if (ename == "parameter")
      {
        XMLAttrString name(element, "name");
        for(auto& c : name)
          c = tolower(c);
        if (!is_parameter(name))
        {
          std::stringstream error;
          error << "InputSection::readXML name "<<name<<" is not a parameter of "<<section_name<<"\n";
          throw UniformCommunicateError(error.str());
        }
        std::istringstream stream(XMLNodeString{element});
        set_from_stream(name,stream);
      }
      element = element->next;
    }

    set_defaults();

    //report();
  };


  void set_from_stream(std::string name, std::istringstream& svalue)
  {
    if (is_string(name))
    {
      std::string value;
      svalue >> value;
      values[name] = value;
    }
    else if(is_bool(name))
    {
      std::string sval;
      svalue >> sval;
      bool value = sval=="yes" || sval=="true" || sval=="1";
      values[name] = value;
    }
    else if(is_integer(name))
    {
      int value;
      svalue >> value;
      values[name] = value;
    }
    else if(is_real(name))
    {
      Real value;
      svalue >> value;
      values[name] = value;
    }
    else
    {
      std::stringstream error;
      error << "InputSection::set_from_stream name "<<name<<" in "<<section_name<<" does not have a type\n";
      throw UniformCommunicateError(error.str());
    }
  };


  void report()
  {
    auto& out = app_log();
    out<<"\n"<<section_name;
    for (auto &[name, value] : values)
    {
      out<<"\n  "<<name<<" = ";
      if (is_string(name))
        out<<std::any_cast<std::string>(value);
      else if (is_bool(name))
        out<<std::any_cast<bool>(value);
      else if (is_integer(name))
        out<<std::any_cast<int>(value);
      else if (is_real(name))
        out<<std::any_cast<Real>(value);
    }
    out<<"\n\n";
  }

};



} // namespace qmcplusplus
#endif /* INPUTSECTION_H */
