//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2021 QMCPACK developers.
//
// File developed by: Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//
//////////////////////////////////////////////////////////////////////////////////////

#include "InputSection.h"
#include "Message/UniformCommunicateError.h"

namespace qmcplusplus
{

  void InputSection::readXML(xmlNodePtr cur)
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

    // assign default values for optional variables
    set_defaults();

    //report();
  };


  void InputSection::report() const
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
  };


  void InputSection::set_from_stream(const std::string& name, std::istringstream& svalue)
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


  void InputSection::set_defaults()
  {
    for (auto &[name, default_value] : default_values)
      if (!has(name))
        values[name] = default_value;
  };

}
