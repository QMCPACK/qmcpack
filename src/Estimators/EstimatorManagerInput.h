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

#ifndef QMCPLUSPLUS_ESIMATORMANAGERINPUT_H
#define QMCPLUSPLUS_ESIMATORMANAGERINPUT_H

#include "type_traits/template_types.hpp"
#include <functional>
#include <vector>
#include <variant>
#include <libxml/tree.h>
#include "io/InputNode.hpp"

namespace qmcplusplus
{

class SpinDensityInput;
class MomentumDistributionInput;
class OneBodyDensityMatricesInput;

using EstimatorInputs = std::vector<
    std::variant<RefW<SpinDensityInput>, RefW<MomentumDistributionInput>, RefW<OneBodyDensityMatricesInput>>>;

namespace testing
{
class EstimatorManagerInputTests;
}

/** Input representation for Driver base class runtime parameters
 */
class EstimatorManagerInput
{
public:
  EstimatorManagerInput() = default;
  EstimatorManagerInput(EstimatorManagerInput&& emi) = default;
  EstimatorManagerInput(xmlNodePtr cur);
  const EstimatorInputs& get_estimator_inputs() const { return estimator_inputs; };
protected:
  /** read <Estimators> node
   */
  void readXML(xmlNodePtr cur);

  EstimatorInputs estimator_inputs;
  UPtrVector<InputNode> estimator_input_storage;

  template<typename T, typename... Args>
  void appendEstimatorInput(Args&&... args)
  {
    estimator_input_storage.push_back(std::make_unique<T>(std::forward<T>(args)...));
    estimator_inputs.push_back(static_cast<T&>(*estimator_input_storage.back()));
  }

  friend class testing::EstimatorManagerInputTests;

};

} // namespace qmcplusplus

#endif
