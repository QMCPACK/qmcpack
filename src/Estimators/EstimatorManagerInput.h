//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2025 QMCPACK developers.
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
 *  We know all the estimator types at compile time and it is useful to have type safety for their usage.
 *  All input clasess must satisfy std::is_trivially_copyable..
 */
class EnergyDensityInput;
class SpinDensityInput;
class MomentumDistributionInput;
class OneBodyDensityMatricesInput;
class SelfHealingOverlapInput;
class MagnetizationDensityInput;
class PerParticleHamiltonianLoggerInput;
class StructureFactorInput;
using EstimatorInput = std::variant<std::monostate,
                                    MomentumDistributionInput,
                                    SpinDensityInput,
                                    OneBodyDensityMatricesInput,
                                    SelfHealingOverlapInput,
                                    MagnetizationDensityInput,
                                    PerParticleHamiltonianLoggerInput,
                                    EnergyDensityInput,
                                    StructureFactorInput>;

using EstimatorInputs = std::vector<EstimatorInput>;

/** The scalar esimtator inputs
 */
class LocalEnergyInput;
class CSLocalEnergyInput;
class RMCLocalEnergyInput;
using ScalarEstimatorInput  = std::variant<std::monostate, LocalEnergyInput, CSLocalEnergyInput, RMCLocalEnergyInput>;
using ScalarEstimatorInputs = std::vector<ScalarEstimatorInput>;

/** Input type for EstimatorManagerNew
 *  Parses Estimators level of input and and delegates child estimator nodes
 *  to appropriate estimator input class
 *  Will later provide access to estimator input instances for estimator construction.
 */
class EstimatorManagerInput
{
public:
  EstimatorManagerInput()                                            = default;
  EstimatorManagerInput(const EstimatorManagerInput& emi)            = default;
  EstimatorManagerInput(EstimatorManagerInput&& emi)                 = default;
  EstimatorManagerInput& operator=(const EstimatorManagerInput& emi) = default;
  EstimatorManagerInput& operator=(EstimatorManagerInput&& emi)      = default;
  // This constructor is for merging the global and local EstimatorManagerInput.
  // but you could also merge more instances
  EstimatorManagerInput(std::initializer_list<EstimatorManagerInput> emil);
  EstimatorManagerInput(xmlNodePtr cur);
  EstimatorInputs& get_estimator_inputs() { return estimator_inputs_; }
  ScalarEstimatorInputs& get_scalar_estimator_inputs() { return scalar_estimator_inputs_; }

  /** read <estimators> node or (<estimator> node for legacy support)
   *  This can be done multiple times with <estimators> nodes
   *  or with <estimator> nodes to support deprecated bare <estimator> definitions
   */
  void readXML(xmlNodePtr cur);

  /** typed appending of already parsed inputs.
   *  only used in testing.
   */
  void append(const EstimatorInput& ei);
  void append(const ScalarEstimatorInput& sei);

  /** Get Indexes of inputs of type T
   *
   *  Right now This is only used to prevent magic numbers in tests,
   *  magic numbers that will result in test maintence issues as the
   *  EstimatorManagerInputTest::create.. cases are
   *  expanded.
   */
  template<typename T>
  std::vector<int> getEstimatorTypeIndexes() const;

private:
  /// this is a vector of variants for typesafe access to the estimator inputs
  EstimatorInputs estimator_inputs_;
  ScalarEstimatorInputs scalar_estimator_inputs_;

  template<typename T>
  void appendEstimatorInput(xmlNodePtr node)
  {
    estimator_inputs_.emplace_back(std::in_place_type<T>, node);
  }

  template<typename T>
  void appendScalarEstimatorInput(xmlNodePtr node)
  {
    scalar_estimator_inputs_.emplace_back(std::in_place_type<T>, node);
  }

  friend class testing::EstimatorManagerInputTests;
};

// Explicity instantiate these to prevent type bleed into the header.
// There is one of these for each forward declaration of a type in the
// EstimaterInput variant
extern template std::vector<int> EstimatorManagerInput::getEstimatorTypeIndexes<EnergyDensityInput>() const;
extern template std::vector<int> EstimatorManagerInput::getEstimatorTypeIndexes<SpinDensityInput>() const;
extern template std::vector<int> EstimatorManagerInput::getEstimatorTypeIndexes<MomentumDistributionInput>() const;
extern template std::vector<int> EstimatorManagerInput::getEstimatorTypeIndexes<OneBodyDensityMatricesInput>() const;
extern template std::vector<int> EstimatorManagerInput::getEstimatorTypeIndexes<SelfHealingOverlapInput>() const;
extern template std::vector<int> EstimatorManagerInput::getEstimatorTypeIndexes<MagnetizationDensityInput>() const;
extern template std::vector<int> EstimatorManagerInput::getEstimatorTypeIndexes<PerParticleHamiltonianLoggerInput>()
    const;

} // namespace qmcplusplus

#endif
