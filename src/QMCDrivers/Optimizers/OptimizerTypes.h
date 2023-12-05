//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2019 QMCPACK developers.
//
// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////

//Defines names for different types of available optimization methods
#ifndef QMCPLUSPLUS_OPTIMIZER_TYPE_HEADER
#define QMCPLUSPLUS_OPTIMIZER_TYPE_HEADER
namespace qmcplusplus
{
enum class OptimizerType
{
  NONE,
  QUARTIC,
  RESCALE,
  LINEMIN,
  ONESHIFTONLY,
  ADAPTIVE,
  DESCENT,
  HYBRID,
  GRADIENT_TEST
};

const std::map<std::string, OptimizerType> OptimizerNames = {{"quartic", OptimizerType::QUARTIC},
                                                             {"rescale", OptimizerType::RESCALE},
                                                             {"linemin", OptimizerType::LINEMIN},
                                                             {"OneShiftOnly", OptimizerType::ONESHIFTONLY},
                                                             {"adaptive", OptimizerType::ADAPTIVE},
                                                             {"descent", OptimizerType::DESCENT},
                                                             {"hybrid", OptimizerType::HYBRID},
                                                             {"gradient_test", OptimizerType::GRADIENT_TEST}};

} // namespace qmcplusplus
#endif
