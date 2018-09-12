//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2018 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//////////////////////////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_SPOSET_BATCHED_H
#define QMCPLUSPLUS_SPOSET_BATCHED_H

#include "QMCWaveFunctions/BsplineFactory/temp_batch_type.h"
#include "QMCWaveFunctions/SPOSetTypeAliases.h"
//! SPOSet evaluation interface depends on walker batching strategy
/*!
  SPOSet inherits this so the correct evaluation function signatures and translations
  for the walker batching strategy are specified.  The actual evaluations
  are implement by SplineAdoptor descendents.
*/
namespace qmcplusplus
{

class SPOSetBatched : public SPOSet,
		      public SPOSetEvaluation<Batching::BATCHED>
{
public:
  using SPOSetPtr = SPOSetBatched*;
};

}

#endif
