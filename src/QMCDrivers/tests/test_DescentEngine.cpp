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

#include <libxml/tree.h>
#include "QMCDrivers/Optimizers/DescentEngine.h"
#include "Optimize/VariableSet.h"
#include "Configuration.h"
#include "Message/Communicate.h"

namespace qmcplusplus
{

typedef qmcplusplus::QMCTraits::FullPrecValueType FullPrecValueType;
typedef qmcplusplus::QMCTraits::ValueType ValueType;

///This provides a basic test of the descent engine's parameter update algorithm
TEST_CASE("DescentEngine RMSprop update","[drivers][descent]")
{

Communicate* myComm;

xmlNodePtr fakeXML;

std::unique_ptr<DescentEngine> descentEngineObj = std::make_unique<DescentEngine>(myComm, fakeXML);

optimize::VariableSet myVars;

//Two fake parameters are specified
optimize::VariableSet::value_type first_param(1.0);
optimize::VariableSet::value_type second_param(-2.0);

myVars.insert("first",first_param);
myVars.insert("second",second_param);

std::vector<ValueType> LDerivs;

//Corresponding fake derivatives are specified and given to the engine
ValueType first_deriv = 5;
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

//The engine should update the parameters using the generic defualt step size of .001 and obtain these values.
REQUIRE(std::real(results[0]) == Approx(.995));
REQUIRE(std::real(results[1]) == Approx(-2.001));

}

}


