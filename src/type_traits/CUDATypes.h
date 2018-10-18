 //////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2018 QMCPACK developers
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Lab
//
// File created by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Lab
//////////////////////////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_CUDA_TYPES_H
#define QMCPLUSPLUS_CUDA_TYPES_H

#include "config.h"
#include "type_traits/QMCTypes.h"

namespace qmcplusplus
{

// Right now we build only one CUDA precision at a time.
// This was set once in QMCTraits then inherited everywhere
// The CUDAGlobalTypeAliases reproduces this old behavior for now
using CUDAGlobalTypeAliases = QMCTypes<CUDA_PRECISION, OHMMS_DIM>;
 
}
#endif
