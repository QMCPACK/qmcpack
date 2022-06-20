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

//Code for an engine to handle hybrid optimization


#ifndef QMCPLUSPLUS_HYBRID_ENGINE_HEADER
#define QMCPLUSPLUS_HYBRID_ENGINE_HEADER

#include <vector>
#include <libxml/tree.h>
#include "Message/Communicate.h"
#include "VariableSet.h"
#include "QMCDrivers/Optimizers/OptimizerTypes.h"


namespace qmcplusplus
{
class HybridEngine
{
  using FullPrecValueType = qmcplusplus::QMCTraits::FullPrecValueType;
  using ValueType         = qmcplusplus::QMCTraits::ValueType;

private:
  Communicate* myComm;


  ///number of optimization steps taken
  int step_num_;

  ///process xml node
  bool processXML(const xmlNodePtr cur);

  ///list of methods known by hybrid
  std::vector<OptimizerType> opt_methods_;

  ///xml saved node
  std::vector<xmlNodePtr> saved_xml_opt_methods_;

  //number of updates in each individual method of hybrid method
  std::vector<int> num_updates_opt_methods_;

  //Helper method for determining the index in opt_methods_ that corresponds to a method based on current step
  int identifyMethodIndex() const;

public:
  //Constructor for engine
  HybridEngine(Communicate* comm, const xmlNodePtr cur);

  //Fake constructor only used in unit test for hybrid engine.
  HybridEngine() { step_num_ = 0; }

  //Returns the appropriate XML for a optimization method to be handled inside QMCFixedSampleLinearOptimize
  xmlNodePtr getSelectedXML();

  //Determines whether to store a vector based on how many have been requested and what the current step is
  bool queryStore(int store_num, OptimizerType method_type) const;

  //adds an optimization method to list. Used only in unit test for hybrid engine.
  void addMethod(OptimizerType method) { opt_methods_.push_back(method); }

  //adds a number of update steps to list. Used only in unit test for hybrid engine.
  void addUpdates(int num_steps) { num_updates_opt_methods_.push_back(num_steps); }
};

} // namespace qmcplusplus
#endif
