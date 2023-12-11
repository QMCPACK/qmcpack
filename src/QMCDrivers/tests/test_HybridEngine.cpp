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

#include "QMCDrivers/Optimizers/HybridEngine.h"
#include "QMCDrivers/Optimizers/OptimizerTypes.h"

namespace qmcplusplus
{

/// This provides a basic test of the hybrid engine's check on whether vectors need to be stored
TEST_CASE("Hybrid Engine query store","[drivers][hybrid]")
{

    //A fake constructor is used that sets the hybrid engine's step_num_ variable to 0 instead of the usual -1
    std::unique_ptr<HybridEngine> hybridEngineObj = std::make_unique<HybridEngine>();

    //Descent and adaptive methods with associated numbers of update steps are added to the engine
    hybridEngineObj->addMethod(OptimizerType::DESCENT);
    hybridEngineObj->addUpdates(100);
    
    hybridEngineObj->addMethod(OptimizerType::ADAPTIVE);
    hybridEngineObj->addUpdates(3);

    //A typical number of vectors that might be requested over the course of a descent optimization in the hybrid method
    int test_store_num = 5;

    //Engine checks whether a vector should be stored. With step_num_ set to 0 this should be false.
    int result = hybridEngineObj->queryStore(test_store_num,OptimizerType::DESCENT);

    app_log() << "Hybrid Test Result: " << result << std::endl;

    REQUIRE(result == false);


}


}

