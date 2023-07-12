//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Miguel Morales, moralessilva2@llnl.gov, Lawrence Livermore National Laboratory
//
// File created by: Miguel Morales, moralessilva2@llnl.gov, Lawrence Livermore National Laboratory
//////////////////////////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_AFQMC_SERIALWALKERSET_HPP
#define QMCPLUSPLUS_AFQMC_SERIALWALKERSET_HPP

#include <random>
#include <type_traits>
#include <memory>

#include "Configuration.h"
#include "OhmmsData/libxmldefs.h"
#include "Utilities/RandomGenerator.h"

#include "AFQMC/config.h"
#include "AFQMC/Utilities/taskgroup.h"

#include "AFQMC/Walkers/Walkers.hpp"
#include "AFQMC/Walkers/WalkerControl.hpp"
#include "AFQMC/Walkers/WalkerConfig.hpp"
#include "AFQMC/Walkers/WalkerSetBase.h"

namespace qmcplusplus
{
namespace afqmc
{
/*
 * Class that contains and handles walkers.
 * Implements communication, load balancing, and I/O operations.   
 * Walkers are always accessed through the handler.
 */
class SerialWalkerSet : public WalkerSetBase<device_allocator<ComplexType>, device_ptr<ComplexType>>

{
  using Base = WalkerSetBase<device_allocator<ComplexType>, device_ptr<ComplexType>>;

public:
  /// constructor
  SerialWalkerSet(afqmc::TaskGroup_& tg_, xmlNodePtr cur, AFQMCInfo& info, RandomBase<RealType>& r)
      : Base(tg_, cur, info, r, device_allocator<ComplexType>{}, shared_allocator<ComplexType>{tg_.TG_local()})
  {}

  /// destructor
  ~SerialWalkerSet() {}

  SerialWalkerSet(SerialWalkerSet const& other) = delete;
  SerialWalkerSet(SerialWalkerSet&& other)      = default;
  SerialWalkerSet& operator=(SerialWalkerSet const& other) = delete;
  SerialWalkerSet& operator=(SerialWalkerSet&& other) = delete;
};

} // namespace afqmc

} // namespace qmcplusplus

//#include "AFQMC/Walkers/SerialWalkerSet.icc"

#endif
