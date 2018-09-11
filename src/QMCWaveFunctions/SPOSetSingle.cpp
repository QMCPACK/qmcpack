//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2018 QMCPACK developers.
//
// File developed by: Peter Doak, Oak Ridge National Laboratory
//                      refactored from SPOSet.cpp
//
// File created by: Peter Doak, Oak Ridge National Laboratory
//////////////////////////////////////////////////////////////////////////////////////

#include "SPOSetSingle.h"

namespace qmcplusplus
{

SPOSetSingle* SPOSetSingle::makeClone() const
{
  APP_ABORT("Missing  SPOSet::makeClone for "+className);
  return nullptr;
}

}
