//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2019 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//
// File created: consult git log
//////////////////////////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_ESTIMATORMANAGERINTERFACE_H
#define QMCPLUSPLUS_ESTIMATORMANAGERINTERFACE_H

#include "Configuration.h"

namespace qmcplusplus
{
  
class EstimatorManagerInterface
{
public:
  using RealType = QMCTraits::FullPrecRealType;
  using EstimatorType =  ScalarEstimatorBase;
  using BufferType = std::vector<RealType>;

  virtual bool is_manager() const = 0;
  ///return the number of ScalarEstimators
  virtual int size() const = 0;

  virtual void start(int blocks, bool record) = 0;
  /** stop a qmc run
   *
   * Replace finalize();
   */
  virtual void stop() = 0;

  /** start  a block
   * @param steps number of steps in a block
   */
  virtual void startBlock(int steps) = 0;

    /** stop a block
   * @param accept acceptance rate of this block
   */


};
}
#endif
