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
   * The existing checkInVariables implementation in WFC/SPO/.. are inclusive and it calls checkInVariables of its members
   * class A: public SPOSet {}
   * class B: public WFC
   * {
   *   A objA;
   *   checkInVariables() { objA.checkInVariables(); }
   * };
   *
   * With OptimizableObject,
   * class A: public OptimizableObject {}
   * class B: public OptimizableObject
   * {
   *   A objA;
   *   checkInVariablesExclusive() { // should not call objA.checkInVariablesExclusive() if objA has been extracted; }
   * };
   * A vector of OptimizableObject, will be created by calling extractOptimizableObjects().
   * All the checkInVariablesExclusive() will be called through this vector and thus
   * checkInVariablesExclusive implementation should only handle non-OptimizableObject members.
   */
  virtual void checkInVariablesExclusive(opt_variables_type& active) = 0;

  /** reset the parameters during optimizations. Exclusive, see checkInVariablesExclusive
   */
  virtual void resetParametersExclusive(const opt_variables_type& active) = 0;

  /** print the state, e.g., optimizables */
  virtual void reportStatus(std::ostream& os) {}

  void setOptimization(bool state) { is_optimized_ = state; }

  /** Write the variational parameters for this object to the VP HDF file
   *
   * The hout parameter should come from VariableSet::writeToHDF
   *
   * Objects can use this function to store additional information to the file.
   *
   * By default the parameters are saved in VariableSet::writeToHDF, and objects
   * do not need to implement this function (yet).
   *
   */
  virtual void writeVariationalParameters(hdf_archive& hout){};

  /** Read the variational parameters for this object from the VP HDF file
   *
   * The hin parameter should come from VariableSet::readFromHDF
   *
   * By default the parameters are read in VariableSet::readFromHDF, and objects
   * do not need to implement this function (yet).
   */
  virtual void readVariationalParameters(hdf_archive& hin){};
};

class UniqueOptObjRefs : public RefVector<OptimizableObject>
{
public:
  OptimizableObject& operator[](size_t i) const { return RefVector<OptimizableObject>::operator[](i); }

  void push_back(OptimizableObject& obj)
  {
    if (obj.getName().empty())
      throw std::logic_error("BUG!! Only named OptimizableObject object can be added to UniqueOptObjRefs!");
    auto result =
        std::find_if(begin(), end(), [&](OptimizableObject& element) { return element.getName() == obj.getName(); });
    if (result == end())
      RefVector<OptimizableObject>::push_back(obj);
  }
};

} // namespace qmcplusplus
#endif
