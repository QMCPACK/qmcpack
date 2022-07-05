//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2019 QMCPACK developers.
//
// File developed by: Leon Otis, leon_otis@berkeley.edu University, University of California Berkeley
//                    Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Leon Otis, leon_otis@berkeley.edu University, University of California Berkeley
//////////////////////////////////////////////////////////////////////////////////////

#include <vector>
#include <string>
#include <numeric>
#include "HybridEngine.h"
#include "OhmmsData/ParameterSet.h"
#include "Message/CommOperators.h"
#include "OhmmsData/XMLParsingString.h"


namespace qmcplusplus
{
HybridEngine::HybridEngine(Communicate* comm, const xmlNodePtr cur) : myComm(comm)
{
  step_num_ = -1;
  processXML(cur);
}

//Hybrid Engine processes xml input to determine which optimizations are being used and for how many iterations
bool HybridEngine::processXML(const xmlNodePtr opt_xml)
{
  opt_methods_.clear();
  saved_xml_opt_methods_.clear();
  num_updates_opt_methods_.clear();

  xmlNodePtr cur = opt_xml->children;
  while (cur != NULL)
  {
    std::string cname((const char*)(cur->name));
    if (cname == "optimizer")
    {
      std::string children_MinMethod;
      ParameterSet m_param;
      m_param.add(children_MinMethod, "MinMethod");
      m_param.put(cur);

      if (children_MinMethod.empty())
        throw std::runtime_error("MinMethod must be given!\n");
      std::string updates_string(getXMLAttributeValue(cur, "num_updates"));
      app_log() << "HybridEngine saved MinMethod " << children_MinMethod << " num_updates = " << updates_string << '\n';
      auto iter = OptimizerNames.find(children_MinMethod);
      if (iter == OptimizerNames.end())
        throw std::runtime_error("Unknown MinMethod!\n");
      opt_methods_.push_back(iter->second);
      saved_xml_opt_methods_.push_back(cur);
      num_updates_opt_methods_.push_back(std::stoi(updates_string));
    }
    cur = cur->next;
  }

  if (saved_xml_opt_methods_.size() != 2)
    throw std::runtime_error("MinMethod hybrid needs two optimizer input blocks!\n");

  return true;
}

//Retrieves the appropriate XML input for the current optimization method
xmlNodePtr HybridEngine::getSelectedXML()
{
  return saved_xml_opt_methods_[identifyMethodIndex()];
}

//Queries whether information from one optimization method(currently needs to be descent) should be stored.
bool HybridEngine::queryStore(int store_num, const OptimizerType method) const
{
  bool store = false;

  int idx = 0;

  auto iter = std::find(opt_methods_.begin(), opt_methods_.end(), method);
  if (iter == opt_methods_.end())
    throw std::runtime_error("Unknown MinMethod!\n");
  else
    idx = std::distance(opt_methods_.begin(), iter);


  const int tot_micro_it = std::accumulate(num_updates_opt_methods_.begin(), num_updates_opt_methods_.end(), 0);
  const int pos          = step_num_ % tot_micro_it;
  //interval is the number of steps between stores based on the number of stores desired
  int interval = num_updates_opt_methods_[idx] / store_num;

  if (interval == 0)
  {
    app_log()
        << "Requested Number of Stored Vectors greater than number of descent steps. Storing a vector on each step."
        << std::endl;
    interval = 1;
  }
  if ((pos + 1) % interval == 0)
  {
    store = true;
  }

  return store;
}

//Determines what the current optimization method is based on the current step number and returns the associated index
int HybridEngine::identifyMethodIndex() const
{
  const int tot_micro_it = std::accumulate(num_updates_opt_methods_.begin(), num_updates_opt_methods_.end(), 0);
  const int pos          = step_num_ % tot_micro_it;

  int run_sum    = 0;
  int select_idx = 0;

  //Compare pos to running sum of microiterations of different methods to determine which method is being used
  for (int i = 0; i < num_updates_opt_methods_.size(); i++)
  {
    run_sum += num_updates_opt_methods_[i];
    if (run_sum > pos)
    {
      select_idx = i;
      break;
    }
  }

  return select_idx;
}

} // namespace qmcplusplus
