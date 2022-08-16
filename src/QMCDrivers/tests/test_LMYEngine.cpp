//////////////////////////////////////////////////////////////////////////////////////
//// This file is distributed under the University of Illinois/NCSA Open Source License.
//// See LICENSE file in top directory for details.
////
//// Copyright (c) 2022 QMCPACK developers.
////
//// File developed by: Leon Otis, leon_otis@berkeley.edu, University of California Berkeley
////
//// File created by: Leon Otis, leon_otis@berkeley.edu, University of California Berkeley
////////////////////////////////////////////////////////////////////////////////////////
//

#ifdef HAVE_LMY_ENGINE

#include "catch.hpp"

#include "VariableSet.h"
#include "formic/utils/lmyengine/engine.h"
#include "formic/utils/lmyengine/var_dependencies.h"

namespace qmcplusplus
{
using FullPrecValueType = qmcplusplus::QMCTraits::FullPrecValueType;
using ValueType         = qmcplusplus::QMCTraits::ValueType;

TEST_CASE("LMYEngine Sample Storage", "[drivers][lm]")
{
  //Construct LM engine as in QMCFixedSampleLinearOptimize for testing
  formic::VarDeps vdeps(1, std::vector<double>());
  std::vector<double> shift_scales(3, 1.0);
  cqmc::engine::LMYEngine<ValueType>* EngineObj =
      new cqmc::engine::LMYEngine<ValueType>(&vdeps,
                                             false, // exact sampling
                                             true,  // ground state?
                                             false, // variance correct,
                                             true,
                                             true,  // print matrices,
                                             true,  // build matrices
                                             false, // spam
                                             false, // use var deps?
                                             true,  // chase lowest
                                             false, // chase closest
                                             false, // eom
                                             false,
                                             false,  // eom related
                                             false,  // eom related
                                             false,  // use block?
                                             120000, // number of samples
                                             0,      // number of parameters
                                             60,     // max krylov iter
                                             0,      // max spam inner iter
                                             1,      // spam appro degree
                                             0,      // eom related
                                             0,      // eom related
                                             0,      // eom related
                                             0.0,    // omega
                                             0.0,    // var weight
                                             1.0e-6, // convergence threshold
                                             0.99,   // minimum S singular val
                                             0.0, 0.0,
                                             10.0, // max change allowed
                                             1.00, // identity shift
                                             1.00, // overlap shift
                                             0.3,  // max parameter change
                                             shift_scales, app_log());


  EngineObj->setFiltering(true);
  EngineObj->setStoringSamples(true);
  EngineObj->setThreshold(1.0);
  EngineObj->setFilterInfo(true);

  app_log() << "Creating fake data to test LM sample storage and parameter filtration" << std::endl;
  int fakeParamNum   = 1;
  int fakeNumSamples = 2;

  EngineObj->setUpStorage(fakeParamNum, fakeNumSamples);
  std::vector<FullPrecValueType> der_rat_samp;
  std::vector<FullPrecValueType> le_der_samp;

  der_rat_samp.resize(fakeParamNum + 1, 0.0);
  le_der_samp.resize(fakeParamNum + 1, 0.0);
  der_rat_samp[0] = 1.0;
  le_der_samp[0] - 1.0;

  der_rat_samp[1] = 0.5;
  le_der_samp[1]  = -2.0;

  EngineObj->store_sample(der_rat_samp, le_der_samp, le_der_samp, 1.0, 1.0, 0);
  le_der_samp[0]  = -1.5;
  der_rat_samp[1] = 1.5;
  le_der_samp[1]  = -0.5;
  EngineObj->store_sample(der_rat_samp, le_der_samp, le_der_samp, 1.0, 1.0, 1);
  EngineObj->selectParameters();

  bool paramOn = EngineObj->getParameterSetting(0);
  //Parameter should be left off
  REQUIRE(paramOn == false);

  delete EngineObj;
}

} // namespace qmcplusplus
#endif
