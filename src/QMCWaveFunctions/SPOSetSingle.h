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

#include "QMCWaveFunctions/BsplineFactory/temp_batch_type.h"
#include "QMCWaveFunctions/SPOSetTypeAliases.h"
#include "QMCWaveFunctions/SPOSet.h"
//! SPOSet evaluation interface depends on walker batching strategy
/*!
  SPOSet inherits this so the correct evaluation function signatures and translations
  for the walker batching strategy are specified.  The actual evaluations
  are implement by SplineAdoptor descendents.
*/
namespace qmcplusplus
{

class SPOSetSingle : public SPOSet,
		     public SPOSetEvaluation<Batching::SINGLE>
{
public:
  //Because so many SPOSetSingle children exist and
  //expect to inherit typedefs
  using SSTA = SPOSetTypeAliases;
  using ValueType = QMCTraits::ValueType;
  using IndexVector_t = SSTA::IndexVector_t;
  using ValueVector_t = SSTA::ValueVector_t;
  using ValueAlignedVector_t = SSTA::ValueAlignedVector_t;
  using ValueMatrix_t = SSTA::ValueMatrix_t;
  using GradVector_t = SSTA::GradVector_t;
  using GradMatrix_t = SSTA::GradMatrix_t;
  using HessVector_t = SSTA::HessVector_t;
  using HessMatrix_t = SSTA::HessMatrix_t;
  using HessType = SSTA::HessType;
  using HessArray_t = SSTA::HessArray_t;
  using GradHessType = SSTA::GGGType;
  //using GradHessVector_t = SSTA::GGGVector_t;
  using GGGVector_t = SSTA::GGGVector_t;
  //using GradHessMatrix_t = SSTA::GGGMatrix_t;
  using GGGMatrix_t = SSTA::GGGMatrix_t;
  using VGLVector_t = SSTA::VGLVector_t;
  using Walker_t = SSTA::Walker_t;
  using SPOSetPtr = SPOSetSingle*;
  virtual SPOSetSingle* makeClone() const;
  using GGGType = SSTA::GGGType;
};

  

}

#endif
