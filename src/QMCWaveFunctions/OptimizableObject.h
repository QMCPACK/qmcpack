//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2022 QMCPACK developers.
//
// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////


#ifndef QMCPLUSPLUS_OPTIMIZABLEOBJECT_H
#define QMCPLUSPLUS_OPTIMIZABLEOBJECT_H

#include "VariableSet.h"

/**@file OptimizableObject.h
 *@brief Declaration of OptimizableObject
 */
namespace qmcplusplus
{
using opt_variables_type = optimize::VariableSet;

class OptimizableObject
{
public:
   OptimizableObject(const std::string& name, bool optimizable)
     : name_(name), optimizable_(optimizable) {}

   bool isOptimized() const { return is_optimized_; }

   const std::string& getName() const { return name_; }
protected:
  /** Name of the optimizable object
   */
  const std::string name_;
  /** If true, this object is actively modified during WFOpt
   */
  bool is_optimized_ = false;
  /** boolean to set optimization
   *
   * If true, this object is actively modified during optimization
   */
  const bool optimizable_;
public:
  /** check in optimizable parameters
   * @param active a super set of optimizable variables
   *
   * Add the paramemters this WaveFunctionComponent manage to active.
   */
  virtual void checkInVariables(opt_variables_type& active) {}

  /** check out optimizable variables
   *
   * Update myVars index map
   */
  virtual void checkOutVariables(const opt_variables_type& active) {}

  /** reset the parameters during optimizations
   */
  virtual void resetParameters(const opt_variables_type& active) {}

  /** print the state, e.g., optimizables */
  virtual void reportStatus(std::ostream& os) {}
};
} // namespace qmcplusplus
#endif
