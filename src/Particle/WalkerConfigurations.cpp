//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2020 QMCPACK developers.
//
// File developed by: Jordan E. Vincent, University of Illinois at Urbana-Champaign
//                    Bryan Clark, bclark@Princeton.edu, Princeton University
//                    Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Cynthia Gu, zg1@ornl.gov, Oak Ridge National Laboratory
//                    Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////


#include "WalkerConfigurations.h"
#include <map>
#include "Utilities/IteratorUtility.h"

namespace qmcplusplus
{
WalkerConfigurations::WalkerConfigurations() = default;

///default destructor
WalkerConfigurations::~WalkerConfigurations() { destroyWalkers(walker_list_.begin(), walker_list_.end()); }


void WalkerConfigurations::createWalkers(int n, size_t numPtcls)
{
  if (walker_list_.empty())
  {
    while (n)
    {
      walker_list_.push_back(std::make_unique<Walker_t>(numPtcls));
      --n;
    }
  }
  else
  {
    if (walker_list_.size() >= n)
    {
      int iw = walker_list_.size(); //copy from the back
      for (int i = 0; i < n; ++i)
      {
        walker_list_.push_back(std::make_unique<Walker_t>(*walker_list_[--iw]));
      }
    }
    else
    {
      int nc  = n / walker_list_.size();
      int nw0 = walker_list_.size();
      for (int iw = 0; iw < nw0; ++iw)
      {
        for (int ic = 0; ic < nc; ++ic)
          walker_list_.push_back(std::make_unique<Walker_t>(*walker_list_[iw]));
      }
      n -= nc * nw0;
      while (n > 0)
      {
        walker_list_.push_back(std::make_unique<Walker_t>(*walker_list_[--nw0]));
        --n;
      }
    }
  }
}


void WalkerConfigurations::resize(int numWalkers, size_t numPtcls)
{
  int dn = numWalkers - walker_list_.size();
  if (dn > 0)
    createWalkers(dn, numPtcls);
  if (dn < 0)
  {
    int nw = -dn;
    if (nw < walker_list_.size())
    {
      walker_list_.erase(walker_list_.begin(), walker_list_.begin() - dn);
    }
  }
}

///returns the next valid iterator
WalkerConfigurations::iterator WalkerConfigurations::destroyWalkers(iterator first, iterator last)
{
  return walker_list_.erase(first, last);
}

void WalkerConfigurations::createWalkers(iterator first, iterator last)
{
  destroyWalkers(walker_list_.begin(), walker_list_.end());
  while (first != last)
  {
    walker_list_.push_back(std::make_unique<Walker_t>(**first));
    ++first;
  }
}

void WalkerConfigurations::destroyWalkers(int nw)
{
  if (nw > walker_list_.size())
  {
    app_warning() << "  Cannot remove walkers. Current Walkers = " << walker_list_.size() << std::endl;
    return;
  }
  nw     = walker_list_.size() - nw;
  int iw = nw;
  walker_list_.erase(walker_list_.begin() + nw, walker_list_.end());
}

void WalkerConfigurations::copyWalkers(iterator first, iterator last, iterator it)
{
  while (first != last)
  {
    (*it++)->makeCopy(**first++);
  }
}

/** Make Metropolis move to the walkers and save in a temporary array.
 * @param it the iterator of the first walker to work on
 * @param tauinv  inverse of the time step
 *
 * R + D + X
 */
void WalkerConfigurations::reset()
{
  for (auto& walker : walker_list_)
  {
    walker->Weight       = 1.0;
    walker->Multiplicity = 1.0;
  }
}

void WalkerConfigurations::putConfigurations(Walker_t::RealType* target, QMCTraits::FullPrecRealType* weights) const
{
  for (const auto& walker : walker_list_)
  {
    std::copy(get_first_address(walker->R), get_last_address(walker->R), target);
    target += get_last_address(walker->R) - get_first_address(walker->R);
    *weights = walker->Weight;
    ++weights;
  }
}

} // namespace qmcplusplus
