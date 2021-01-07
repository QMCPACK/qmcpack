//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2020 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//
// File refactored from: OperatorEstBase.cpp
//////////////////////////////////////////////////////////////////////////////////////

/**@file
 */
#include "Message/Communicate.h"
#include "OperatorEstBase.h"
#include "QMCHamiltonians/QMCHamiltonian.h"

namespace qmcplusplus
{
OperatorEstBase::OperatorEstBase() : walkers_weight_(0)
{
  data_locality_ = DataLocality::crowd;
}

OperatorEstBase::OperatorEstBase(const OperatorEstBase& oth) : walkers_weight_(0)
{
  data_locality_ = oth.data_locality_;     
}

// I suspect this can be a pure function outside of the class.
// In this case at least we don't care to copy the data_ as we are going to reduce these later and don't want
// to end up with a multiplicative factor if we already have data.
OperatorEstBase::Data OperatorEstBase::createLocalData(size_t size, DataLocality data_locality)
{
  Data new_data;
  if (data_locality == DataLocality::crowd)
  {
    new_data = std::make_unique<std::vector<QMCT::RealType>>(size, 0);
  }
  else
  {
    throw std::runtime_error("currently SpinDensityNew only supports crowd level datalocality");
  }
  return new_data;
}

} // namespace qmcplusplus
