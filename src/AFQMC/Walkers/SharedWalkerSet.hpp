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

#ifndef QMCPLUSPLUS_AFQMC_SHAREDWALKERSET_H
#define QMCPLUSPLUS_AFQMC_SHAREDWALKERSET_H

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
class SharedWalkerSet : public WalkerSetBase<shared_allocator<ComplexType>, ComplexType*>
{
public:
  using Base = WalkerSetBase<shared_allocator<ComplexType>, ComplexType*>;

  /// constructor
  SharedWalkerSet(afqmc::TaskGroup_& tg_, xmlNodePtr cur, AFQMCInfo& info, RandomBase<RealType>& r)
      : Base(tg_,
             cur,
             info,
             r,
             shared_allocator<ComplexType>{tg_.TG_local()},
             shared_allocator<ComplexType>{tg_.TG_local()})
  {}

  /// destructor
  ~SharedWalkerSet() {}

  SharedWalkerSet(SharedWalkerSet const& other) = delete;
  SharedWalkerSet(SharedWalkerSet&& other)      = default;
  SharedWalkerSet& operator=(SharedWalkerSet const& other) = delete;
  SharedWalkerSet& operator=(SharedWalkerSet&& other) = delete;
};

} // namespace afqmc

} // namespace qmcplusplus

//#include "AFQMC/Walkers/SharedWalkerSet.icc"

#endif
