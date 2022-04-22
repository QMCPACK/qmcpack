//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2022 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//
// File refactored from: EstimatorManagerNew.h
//////////////////////////////////////////////////////////////////////////////////////

/** \file
 *  declares the list  supported estimator input types and declares the input type for 
 *  EstimatorManagerNew.
 */
#ifndef QMCPLUSPLUS_ESIMATORMANAGERINPUT_H
#define QMCPLUSPLUS_ESIMATORMANAGERINPUT_H


#include "type_traits/template_types.hpp"
#include <functional>
#include <vector>
#include <variant>
#include <libxml/tree.h>

namespace qmcplusplus
{

namespace testing
{
/** class to allow private method testing in unit tests
 */
class EstimatorManagerInputTests;
} // namespace testing

/** These are the estimator input types EstimatorManagerInput delegates to.
 *  We of course know all the estimator types at compile time and it is useful to have type safety for their usage.
 */
class SpinDensityInput;
class MomentumDistributionInput;
class OneBodyDensityMatricesInput;
using EstimatorInput  = std::variant<MomentumDistributionInput, SpinDensityInput, OneBodyDensityMatricesInput>;
using EstimatorInputs = std::vector<EstimatorInput>;

/** The scalar esimtator inputs
 */
class LocalEnergyInput;
class CSLocalEnergyInput;
using ScalarEstimatorInput  = std::variant<LocalEnergyInput, CSLocalEnergyInput>;
using ScalarEstimatorInputs = std::vector<ScalarEstimatorInput>;

/** Input type for EstimatorManagerNew
 *  Parses Estimators level of input and and delegates child estimator nodes
 *  to appropriate estimator input class
 *  Will later provide access to estimator input instances for esimator construction.
 */
class EstimatorManagerInput
{
public:
  EstimatorManagerInput()                            = default;
  EstimatorManagerInput(EstimatorManagerInput&& emi) = default;
  EstimatorManagerInput(xmlNodePtr cur);
  EstimatorInputs& get_estimator_inputs() { return estimator_inputs_; }
  ScalarEstimatorInputs& get_scalar_estimator_inputs() { return scalar_estimator_inputs_; }
  /** read <Estimators> node
   *  Note that this can be done multiple times to combine global and section local estimator inputs.
   */
  void readXML(xmlNodePtr cur);

private:
  /// this is a vector of variants for typesafe access to the estimator inputs
  EstimatorInputs estimator_inputs_;
  ScalarEstimatorInputs scalar_estimator_inputs_;

  template<typename T, typename... Args>
  void appendEstimatorInput(Args&&... args)
  {
    estimator_inputs_.emplace_back(std::in_place_type<T>, std::forward<Args>(args)...);
    ;
  }

  template<typename T, typename... Args>
  void appendScalarEstimatorInput(Args&&... args)
  {
    scalar_estimator_inputs_.emplace_back(std::in_place_type<T>, std::forward<Args>(args)...);
    ;
  }

  friend class testing::EstimatorManagerInputTests;
};

} // namespace qmcplusplus

#endif
