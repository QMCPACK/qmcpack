//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


#ifndef OHMMS_OHMMSPARAMETER_H
#define OHMMS_OHMMSPARAMETER_H

#include <stdexcept>
#include "OhmmsData/OhmmsElementBase.h"
#include "Host/OutputManager.h"
#include "ModernStringUtils.hpp"

/** generic class for parameter xmlNode
 *
 * <parameter/> node is used to generically add a named parameter whose
 *value is the content of an xmlNode. The definition confirms docbook::parameter
 *\htmlonly
 Usage is:
 &lt;parameter name="<b>aname</b>" &gt; <b>value</b> &lt;/parameter&gt;
 <br>
  Example is:
 &lt;parameter name="temperature" &gt; 200 &lt;/parameter&gt;
 <br>
 <ul>
 <li> <b>name</b> is the name of the parameter.
 <li> <b>value</b> is the content of the parameter of type T.
 </ul>
 <ul> Two kinds of template parameter T are valid:
 <li> intrinsic C/C variables, such as int,double
 <li> Ohmms basic data types in OhmmsPETE, such as  TinyVector&lt;T,D&gt; and  Tensor&lt;T,D&gt;
 </ul>
 Using other types are valid, as far as the operators << and >> can handle the xmlChar*.

 \endhtmlonly
 */

enum class TagStatus
{
  OPTIONAL,   // there is a default value
  REQUIRED,   // no default value, input required.
  DEPRECATED, // deprecated input tag. Input becomes optional
  DELETED,    // deleted input tag. Hard stop if detected.
  UNSUPPORTED // not supported input tags due to incompatible changes. Hard stop if detected.
};

inline void checkTagStatus(const std::string& tagname, TagStatus status)
{
  if (status == TagStatus::DELETED)
  {
    std::ostringstream msg;
    msg << "Input tag \"" << tagname << "\" has been deleted. Please remove it from the input file!" << std::endl;
    throw std::runtime_error(msg.str());
  }
  else if (status == TagStatus::UNSUPPORTED)
  {
    std::ostringstream msg;
    msg << "Input tag \"" << tagname << "\" is not supported. Please remove it from the input file!" << std::endl;
    throw std::runtime_error(msg.str());
  }
  else if (status == TagStatus::DEPRECATED)
    qmcplusplus::app_warning() << "Input tag \"" << tagname
                               << "\" has been deprecated and will be deleted in the next release!" << std::endl;
}


template<class T>
class OhmmsParameter : public OhmmsElementBase
{
  //@{
  ///reference to a value of type T
  T& ref_;
  ///valid candidate values, if not empty, the first one is the default.
  std::vector<T> candidate_values_;
  ///pointer to the corresponding xmlNode
  xmlNodePtr node_;
  /// The status of the tag
  TagStatus tag_staus_;
  //@}

  void checkValues()
  {
    // when the size is 1, values are unbounded and the given value is the default value
    if (candidate_values_.size() > 1 &&
        std::find(candidate_values_.begin(), candidate_values_.end(), ref_) == candidate_values_.end())
    {
      std::ostringstream msg;
      msg << "Input tag \"" << myName << "\" value \"" << ref_ << "\" is not valid. Candidate values are : ";
      for (const auto& value : candidate_values_)
        msg << " \"" << value << "\"";
      msg << std::endl;
      throw std::runtime_error(msg.str());
    }
  }

public:
  /*!\fn OhmmsParameter(...) constructor
   *\param a the value to be referenced
   *\param aname the name of this object
   *\param candidate_values candidate values to be checked against, the first element is the default value. It can be left empty for backward compatibility or unbounded input.
   *\param status Tag status, See TagStatus enum.
   */
  OhmmsParameter(T& a,
                 const std::string& aname,
                 std::vector<T>&& candidate_values = {},
                 TagStatus status                  = TagStatus::OPTIONAL)
      : OhmmsElementBase(aname.c_str()),
        ref_(a),
        candidate_values_(std::move(candidate_values)),
        node_(NULL),
        tag_staus_(status)
  {
    // set default value
    if (status != TagStatus::REQUIRED && !candidate_values_.empty())
      ref_ = candidate_values_[0];
  }

  ///print to an std::ostream
  inline bool get(std::ostream& os) const override
  {
    os << "<parameter name=\"" << myName << "\">" << ref_ << "</parameter>" << std::endl;
    return true;
  }

  /*!inline bool put(xmlNodePtr cur), used by ParameterSet
   *\param cur the current xmlNode whose content is assigned to ref_
   */
  inline bool put(xmlNodePtr cur) override
  {
    checkTagStatus(myName, tag_staus_);
    node_ = cur;
    putContent(ref_, cur);
    checkValues();
    return true;
  }


  ///read from std::istream, used by OhmmsAttributeSet
  inline bool put(std::istream& is) override
  {
    checkTagStatus(myName, tag_staus_);
    is >> ref_;
    checkValues();
    return true;
  }

  /*!\fn bool add(xmlNodePtr parent)
   *\param parent the parent node to which a xmlNode for this object is appended.
   *\brief This function is used by the users to add a xmlNode, when the
   *input file does not contain the corresponding <parameter/>. The content
   *of the new xmlNode is the current value of ref_.
   */
  bool add(xmlNodePtr parent) override
  {
    if (!node_)
    {
      node_ = xmlNewChild(parent, parent->ns, (const xmlChar*)"parameter", NULL);
      xmlNewProp(node_, (const xmlChar*)"name", (const xmlChar*)(myName.c_str()));
      getContent(ref_, node_);
    }
    return true;
  }

  inline void setValue(T x) { ref_ = x; }

  ///reset member data
  inline void reset() override { getContent(ref_, node_); }
};

/*!\class OhmmsParameter<bool>
 *\brief A specialization of OhmmsParameter<T> for T = boolean.
 */
template<>
class OhmmsParameter<bool> : public OhmmsElementBase
{
  //@{
  ///reference to a value of type T
  bool& ref_;
  ///valid candidate values, if not empty, the first one is the default.
  std::vector<bool> candidate_values_;
  ///pointer to the corresponding xmlNode
  xmlNodePtr node_;
  /// The status of the tag
  TagStatus tag_staus_;
  //@}

public:
  /*!\fn OhmmsParameter(bool& a, const char* aname, const char* uname)
   *\param a the boolean to be referenced.
   *\param aname the name of this object
   *\param candidate_values valid input values
   */
  OhmmsParameter(bool& a,
                 const std::string& aname,
                 std::vector<bool>&& candidate_values = {},
                 TagStatus status                     = TagStatus::OPTIONAL)
      : OhmmsElementBase(aname.c_str()),
        ref_(a),
        candidate_values_(std::move(candidate_values)),
        node_(NULL),
        tag_staus_(status)
  {
    // set default value.
    if (!candidate_values_.empty())
      ref_ = candidate_values_[0];
  }

  ///print to an std::ostream
  inline bool get(std::ostream& os) const override
  {
    os << "<parameter name=\"" << myName << "\"> " << (ref_ ? "yes" : "no") << " </parameter>" << std::endl;
    return true;
  }

  /*!inline bool put(xmlNodePtr cur)
   *\param cur the current xmlNode whose content is assigned to ref_
   *\brief If the content is empty, the negation of the current value is taken.
   *Example is <parameter name="force"/> to turn on the force-evaluation flag
   *of NoPropagator.
   */
  inline bool put(xmlNodePtr cur) override
  {
    checkTagStatus(myName, tag_staus_);
    node_ = cur;
    const XMLNodeString ac(cur);
    std::istringstream stream(ac);
    return put(stream);
  }

  inline void setValue(bool x) { ref_ = x; }

  ///read from std::istream
  inline bool put(std::istream& is) override
  {
    checkTagStatus(myName, tag_staus_);
    std::string input_value;
    is >> input_value;
    std::string input_value_lower_case = qmcplusplus::lowerCase(input_value);

    if (input_value_lower_case == "yes" || input_value_lower_case == "true")
      ref_ = true;
    else if (input_value_lower_case == "no" || input_value_lower_case == "false")
      ref_ = false;
    else if (input_value.empty())
      throw std::runtime_error(myName + " requires a single value input.");
    else
      throw std::runtime_error(myName + " only accepts 'yes'/'no'/'true'/'false' but the input value is '" +
                               input_value + "'.");

    std::string dummy;
    is >> dummy;
    if (!dummy.empty())
      throw std::runtime_error(myName + " only accepts a single value input.");
    return true;
  }

  /*!\fn bool add(xmlNodePtr parent)
   *\param parent the parent node to which a xmlNode for this object is appended.
   *\brief This function is used by the users to add a xmlNode, when the
   *input file does not contain the corresponding <parameter/>. The content
   *of the new xmlNode is the current value of ref_.
   */
  bool add(xmlNodePtr parent) override
  {
    if (!node_)
    {
      node_ = xmlNewChild(parent, parent->ns, (const xmlChar*)"parameter", NULL);
      xmlNewProp(node_, (const xmlChar*)"name", (const xmlChar*)(myName.c_str()));
      getContent(ref_, node_);
    }
    return true;
  }

  ///reset member data
  inline void reset() override { getContent(ref_, node_); }
};
#endif /*OHMMS_OHMMSPARAMETER_H*/
