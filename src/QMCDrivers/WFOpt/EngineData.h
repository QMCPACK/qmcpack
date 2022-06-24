//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2022 QMCPACK developers.
//
// File developed by: Leon Otis, leon_otis@berkeley.edu, UC Berkeley
//
// File created by: Leon Otis, leon_otis@berkeley.edu, UC Berkeley
//////////////////////////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_ENGINE_DATA_HEADER
#define QMCPLUSPLUS_ENGINE_DATA_HEADER  

#include "QMCDrivers/Optimizers/DescentEngine.h"
namespace qmcplusplus
{
  struct engineData
  {
#ifdef HAVE_LMY_ENGINE 
using ValueType = QMCTraits::ValueType; 
    cqmc::engine::LMYEngine<ValueType>* lmEngine;
  
#endif

    DescentEngine* descentEngine;
    std::string method;
  };
} // namespace qmcplusplus
#endif

 
