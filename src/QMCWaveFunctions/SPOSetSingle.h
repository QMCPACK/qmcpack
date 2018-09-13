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

#ifndef QMCPLUSPLUS_SPOSET_SINGLE_H
#define QMCPLUSPLUS_SPOSET_SINGLE_H

#include "QMCWaveFunctions/Batching.h"
#include "QMCWaveFunctions/SPOSetTypeAliases.h"
#include "QMCWaveFunctions/SPOSet.h"
#include "QMCWaveFunctions/SPOSetEvaluation.h"
/*! SPOSet evaluation interface depends on walker batching strategy
  SPOSet inherits this so the correct evaluation function signatures and translations
  for the walker batching strategy are specified.  The actual evaluations
  are implemented by SplineAdoptor children.
*/
namespace qmcplusplus
{

class SPOSetSingle : public SPOSet,
		     public SPOSetEvaluation<Batching::SINGLE>
{
  /*! Implements SPOSet for single walker evaluation
    uses templated evaluation interface
  */
public:
  using SPOSetPtr = SPOSetSingle*;

  virtual SPOSetSingle* makeClone() const;
};

  

}

#endif
