//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2019 QMCPACK developers.
//
// File developed by: Leon Otis, leon_otis@berkeley.edu University, University of California Berkeley
//		      Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Leon Otis, leon_otis@berkeley.edu University, University of California Berkeley
//////////////////////////////////////////////////////////////////////////////////////

#include "catch.hpp"

#include "OhmmsData/Libxml2Doc.h"
#include "QMCDrivers/Optimizers/DescentEngine.h"
#include "VariableSet.h"
#include "Configuration.h"
#include "Message/Communicate.h"

namespace qmcplusplus
{

using FullPrecValueType = qmcplusplus::QMCTraits::FullPrecValueType;
using ValueType         = qmcplusplus::QMCTraits::ValueType;

#if !defined(MIXED_PRECISION)
///This provides a basic test of the descent engine's parameter update algorithm
TEST_CASE("DescentEngine RMSprop update", "[drivers][descent]")
{
  Communicate* c = OHMMS::Controller;


  const std::string engine_input("<tmp> </tmp>");

  Libxml2Document doc;
  bool okay = doc.parseFromString(engine_input);
  REQUIRE(okay);

  xmlNodePtr fakeXML = doc.getRoot();

  std::unique_ptr<DescentEngine> descentEngineObj = std::make_unique<DescentEngine>(c, fakeXML);

  optimize::VariableSet myVars;

  //Two fake parameters are specified
  optimize::VariableSet::value_type first_param(1.0);
  optimize::VariableSet::value_type second_param(-2.0);

  myVars.insert("first", first_param);
  myVars.insert("second", second_param);

  std::vector<ValueType> LDerivs;

  //Corresponding fake derivatives are specified and given to the engine
  ValueType first_deriv  = 5;
  ValueType second_deriv = 1;

  LDerivs.push_back(first_deriv);
  LDerivs.push_back(second_deriv);

  descentEngineObj->setDerivs(LDerivs);

  descentEngineObj->setupUpdate(myVars);

  descentEngineObj->storeDerivRecord();
  descentEngineObj->updateParameters();

  std::vector<ValueType> results = descentEngineObj->retrieveNewParams();

  app_log() << "Descent engine test of parameter update" << std::endl;
  app_log() << "First parameter: " << results[0] << std::endl;
  app_log() << "Second parameter: " << results[1] << std::endl;

  //The engine should update the parameters using the generic default step size of .001 and obtain these values.
  REQUIRE(std::real(results[0]) == Approx(.995));
  REQUIRE(std::real(results[1]) == Approx(-2.001));

  //Provide fake data to test mpi_unbiased_ratio_of_means
  int n              = 2;
  ValueType mean     = 0;
  ValueType variance = 0;
  ValueType stdErr   = 0;

  std::vector<ValueType> weights;
  weights.push_back(1.0);
  weights.push_back(1.0);
  std::vector<ValueType> numerSamples;
  numerSamples.push_back(-2.0);
  numerSamples.push_back(-2.0);
  std::vector<ValueType> denomSamples;
  denomSamples.push_back(1.0);
  denomSamples.push_back(1.0);

  descentEngineObj->mpi_unbiased_ratio_of_means(n, weights, numerSamples, denomSamples, mean, variance, stdErr);
  app_log() << "Descent engine test of mpi_unbiased_ratio_of_means" << std::endl;
  app_log() << "Mean: " << mean << std::endl;
  app_log() << "Variance: " << variance << std::endl;
  app_log() << "Standard Error: " << stdErr << std::endl;

  //mpi_unbiased_ratio_of_means should calculate the mean, variance, and standard error and obtain the values below
  REQUIRE(std::real(mean) == Approx(-2.0));
  REQUIRE(std::real(variance) == Approx(0.0));
  REQUIRE(std::real(stdErr) == Approx(0.0));
}
#endif
} // namespace qmcplusplus
