//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


#include "VariableSet.h"
#include <map>
#include <stdexcept>
#include <iomanip>
#include <ios>
#include <algorithm>
#include "io/OhmmsData/Libxml2Doc.h"

using std::setw;

using qmcplusplus::app_log;

namespace optimize
{
//   VariableSet::VariableSet(variable_map_type& input):num_active_vars(0)
//   {
//     Index.resize(input.size(),-1);
//     NameAndValue.resize(input.size());
//     copy(input.begin(),input.end(),NameAndValue.begin());
//     for(int i=0; i<Index.size(); ++i) Index[i]=i;
//
//     ParameterType.resize(0); Recompute.resize(0);
//     for(int i=0; i<Index.size(); ++i) ParameterType.push_back(index_pair_type(NameAndValue[i].first,0));
//     for(int i=0; i<Index.size(); ++i) Recompute.push_back(index_pair_type(NameAndValue[i].first,1));
//   }

void VariableSet::clear()
{
  num_active_vars = 0;
  Index.clear();
  NameAndValue.clear();
  Recompute.clear();
  ParameterType.clear();
}

/** insert name-value pairs  of this object to output
 * @param output parameters to be added
 */
//   void VariableSet::insertTo(variable_map_type& output) const
//   {
//     for(int i=0; i<Index.size(); ++i)
//     {
//       if(Index[i]>-1)
//       {
//         output[NameAndValue[i].first]=NameAndValue[i].second;
//         output[Recompute[i].first]=Recompute[i].second;
//         output[ParameterType[i].first]=ParameterType[i].second;
//       }
//     }
//   }

void VariableSet::insertFrom(const VariableSet& input)
{
  for (int i = 0; i < input.size(); ++i)
  {
    iterator loc = find(input.name(i));
    if (loc == NameAndValue.end())
    {
      Index.push_back(input.Index[i]);
      NameAndValue.push_back(input.NameAndValue[i]);
      ParameterType.push_back(input.ParameterType[i]);
      Recompute.push_back(input.Recompute[i]);
    }
    else
      (*loc).second = input.NameAndValue[i].second;
  }
  num_active_vars = input.num_active_vars;
}

void VariableSet::insertFromSum(const VariableSet& input_1, const VariableSet& input_2)
{
  value_type sum_val;
  std::string vname;

  // Check that objects to be summed together have the same number of active
  // variables.
  if (input_1.num_active_vars != input_2.num_active_vars)
    throw std::runtime_error("Inconsistent number of parameters in two provided "
                             "variable sets.");

  for (int i = 0; i < input_1.size(); ++i)
  {
    // Check that each of the equivalent variables in both VariableSet objects
    // have the same name - otherwise we certainly shouldn't be adding them.
    if (input_1.NameAndValue[i].first != input_2.NameAndValue[i].first)
      throw std::runtime_error("Inconsistent parameters exist in the two provided "
                               "variable sets.");

    sum_val = input_1.NameAndValue[i].second + input_2.NameAndValue[i].second;

    iterator loc = find(input_1.name(i));
    if (loc == NameAndValue.end())
    {
      Index.push_back(input_1.Index[i]);
      ParameterType.push_back(input_1.ParameterType[i]);
      Recompute.push_back(input_1.Recompute[i]);

      // We can reuse the above values, which aren't summed between the
      // objects, but the parameter values themselves need to use the summed
      // values.
      vname = input_1.NameAndValue[i].first;
      NameAndValue.push_back(pair_type(vname, sum_val));
    }
    else
      (*loc).second = sum_val;
  }
  num_active_vars = input_1.num_active_vars;
}

void VariableSet::insertFromDiff(const VariableSet& input_1, const VariableSet& input_2)
{
  value_type diff_val;
  std::string vname;

  // Check that objects to be subtracted have the same number of active
  // variables.
  if (input_1.num_active_vars != input_2.num_active_vars)
    throw std::runtime_error("Inconsistent number of parameters in two provided "
                             "variable sets.");

  for (int i = 0; i < input_1.size(); ++i)
  {
    // Check that each of the equivalent variables in both VariableSet objects
    // have the same name - otherwise we certainly shouldn't be subtracting them.
    if (input_1.NameAndValue[i].first != input_2.NameAndValue[i].first)
      throw std::runtime_error("Inconsistent parameters exist in the two provided "
                               "variable sets.");

    diff_val = input_1.NameAndValue[i].second - input_2.NameAndValue[i].second;

    iterator loc = find(input_1.name(i));
    if (loc == NameAndValue.end())
    {
      Index.push_back(input_1.Index[i]);
      ParameterType.push_back(input_1.ParameterType[i]);
      Recompute.push_back(input_1.Recompute[i]);

      // We can reuse the above values, which aren't subtracted between the
      // objects, but the parameter values themselves need to use the
      // subtracted values.
      vname = input_1.NameAndValue[i].first;
      NameAndValue.push_back(pair_type(vname, diff_val));
    }
    else
      (*loc).second = diff_val;
  }
  num_active_vars = input_1.num_active_vars;
}

void VariableSet::activate(const variable_map_type& selected)
{
  //activate the variables
  variable_map_type::const_iterator it(selected.begin()), it_end(selected.end());
  while (it != it_end)
  {
    iterator loc = find((*it++).first);
    if (loc != NameAndValue.end())
    {
      int i = loc - NameAndValue.begin();
      if (Index[i] < 0)
        Index[i] = num_active_vars++;
    }
  }
}

void VariableSet::disable(const variable_map_type& selected)
{
  variable_map_type::const_iterator it(selected.begin()), it_end(selected.end());
  while (it != it_end)
  {
    int loc = find((*it++).first) - NameAndValue.begin();
    if (loc < NameAndValue.size())
      Index[loc] = -1;
  }
}


void VariableSet::removeInactive()
{
  std::vector<int> valid(Index);
  std::vector<pair_type> acopy(NameAndValue);
  std::vector<index_pair_type> bcopy(Recompute), ccopy(ParameterType);
  num_active_vars = 0;
  Index.clear();
  NameAndValue.clear();
  Recompute.clear();
  ParameterType.clear();
  for (int i = 0; i < valid.size(); ++i)
  {
    if (valid[i] > -1)
    {
      Index.push_back(num_active_vars++);
      NameAndValue.push_back(acopy[i]);
      Recompute.push_back(bcopy[i]);
      ParameterType.push_back(ccopy[i]);
    }
  }
}

void VariableSet::resetIndex()
{
  num_active_vars = 0;
  for (int i = 0; i < Index.size(); ++i)
  {
    Index[i] = (Index[i] < 0) ? -1 : num_active_vars++;
  }
}

void VariableSet::getIndex(const VariableSet& selected)
{
  for (int i = 0; i < NameAndValue.size(); ++i)
  {
    Index[i] = selected.getIndex(NameAndValue[i].first);
    if (Index[i] >= 0)
      num_active_vars++;
  }
}

void VariableSet::setDefaults(bool optimize_all)
{
  for (int i = 0; i < Index.size(); ++i)
    Index[i] = optimize_all ? i : -1;
}

void VariableSet::print(std::ostream& os, int leftPadSpaces, bool printHeader) const
{
  std::string pad_str = std::string(leftPadSpaces, ' ');
  int max_name_len    = 0;
  if (NameAndValue.size() > 0)
    max_name_len =
        std::max_element(NameAndValue.begin(), NameAndValue.end(), [](const pair_type& e1, const pair_type& e2) {
          return e1.first.length() < e2.first.length();
        })->first.length();

  int max_value_len     = 28; // 6 for the precision and 7 for minus sign, leading value, period, and exponent.
  int max_type_len      = 1;
  int max_recompute_len = 1;
  int max_use_len       = 3;
  int max_index_len     = 1;
  if (printHeader)
  {
    max_name_len      = std::max(max_name_len, 4); // size of "Name" header
    max_type_len      = 4;
    max_recompute_len = 9;
    max_index_len     = 5;
    os << pad_str << setw(max_name_len) << "Name"
       << " " << setw(max_value_len) << "Value"
       << " " << setw(max_type_len) << "Type"
       << " " << setw(max_recompute_len) << "Recompute"
       << " " << setw(max_use_len) << "Use"
       << " " << setw(max_index_len) << "Index" << std::endl;
    os << pad_str << std::setfill('-') << setw(max_name_len) << ""
       << " " << setw(max_value_len) << ""
       << " " << setw(max_type_len) << ""
       << " " << setw(max_recompute_len) << ""
       << " " << setw(max_use_len) << ""
       << " " << setw(max_index_len) << "" << std::endl;
    os << std::setfill(' ');
  }

  for (int i = 0; i < NameAndValue.size(); ++i)
  {
    os << pad_str << setw(max_name_len) << NameAndValue[i].first << " " << std::setprecision(6) << std::scientific
       << setw(max_value_len) << NameAndValue[i].second << " " << setw(max_type_len) << ParameterType[i].second << " "
       << setw(max_recompute_len) << Recompute[i].second << " ";

    os << std::defaultfloat;

    if (Index[i] < 0)
      os << setw(max_use_len) << "OFF" << std::endl;
    else
      os << setw(max_use_len) << "ON"
         << " " << setw(max_index_len) << Index[i] << std::endl;
  }
}

void VariableSet::saveAsXML(const std::string& filename) const
{
  Libxml2Document doc;
  doc.newDoc("variational_parameters");
  xmlNodePtr pairs = doc.addChild(doc.getRoot(), "name_value_pairs");
  for (auto pair_it : NameAndValue)
  {
    xmlNodePtr pair = doc.addChild(pairs, "vp");
    doc.addChild(pair, "name", pair_it.first);
    doc.addChild(pair, "value", pair_it.second);
  }
  doc.dump(filename);
}

bool VariableSet::readFromXML(const std::string& filename)
{
  Libxml2Document doc;
  bool parse_okay = doc.parse(filename);
  if (!parse_okay)
  {
    app_log() << "Failed to load or parse variational parameter file: " << filename << std::endl;
    return false;
  }

  OhmmsXPathObject vp_base("//variational_parameters", doc.getXPathContext());
  if (vp_base.empty())
  {
    app_log() << "No variational_parameters tag in file:" << filename << std::endl;
    return false;
  }
  if (vp_base.size() > 1)
  {
    app_log() << "More than one variational_parameters tag in file:" << filename << std::endl;
    app_log() << "Using only the first one" << std::endl;
  }

  processChildren(vp_base[0], [&](const std::string& cname, const xmlNodePtr element) {
    if (cname == "name_value_pairs")
    {
      processChildren(element, [&](const std::string cname2, const xmlNodePtr element2) {
        if (cname2 == "vp")
        {
          std::string vp_name;
          value_type vp_value;
          processChildren(element2, [&](const std::string cname3, const xmlNodePtr element3) {
            if (cname3 == "name")
              putContent(vp_name, element3);
            else if (cname3 == "value")
              putContent(vp_value, element3);
          });
          if (find(vp_name) != end())
            (*this)[vp_name] = vp_value;
        }
      });
    }
  });

  return true;
}


} // namespace optimize
