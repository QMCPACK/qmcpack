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
#include "type_traits/template_types.hpp"

/**@file OptimizableObject.h
 *@brief Declaration of OptimizableObject
 */
namespace qmcplusplus
{
using opt_variables_type = optimize::VariableSet;

class OptimizableObject
{
public:
  OptimizableObject(const std::string& name) : name_(name) {}

  const std::string& getName() const { return name_; }
  bool isOptimized() const { return is_optimized_; }

private:
  /** Name of the optimizable object
   */
  const std::string name_;
  /** If true, this object is actively modified during WFOpt
   */
  bool is_optimized_ = false;

public:
  /** check in variational parameters to the global list of parameters used by the optimizer.
   * @param active a super set of optimizable variables
   *
   * This is a query function and should never be implemented as a feature blocker.
   * If an OptimizableObject derived class doesn't support optimization, use the base class fallback.
   */
  virtual void checkInVariables(opt_variables_type& active)
  {
    throw std::logic_error("BUG!! OptimizableObject::checkInVariables should not be called!");
  }

  /** check out variational optimizable variables
   * @param active a super set of optimizable variables
   */
  virtual void checkOutVariables(const opt_variables_type& active)
  {
    throw std::logic_error("BUG!! OptimizableObject::checkOutVariables should not be called!");
  }

  /** reset the parameters during optimizations
   */
  virtual void resetParameters(const opt_variables_type& active)
  {
    throw std::logic_error("BUG!! OptimizableObject::resetParameters should not be called!");
  }

  /** print the state, e.g., optimizables */
  virtual void reportStatus(std::ostream& os) {}
};

class UniqueOptObjRefs : public RefVector<OptimizableObject>
{
public:
  OptimizableObject& operator[](size_t i) const { return RefVector<OptimizableObject>::operator[](i); }

  void push_back(OptimizableObject& obj)
  {
    auto result =
        std::find_if(begin(), end(), [&](OptimizableObject& element) { return element.getName() == obj.getName(); });
    if (result == end())
      RefVector<OptimizableObject>::push_back(obj);
  }
};

} // namespace qmcplusplus
#endif
