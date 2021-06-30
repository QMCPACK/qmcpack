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
WalkerConfigurations::WalkerConfigurations() : LocalNumWalkers(0), GlobalNumWalkers(0) {}

///default destructor
WalkerConfigurations::~WalkerConfigurations() { destroyWalkers(WalkerList.begin(), WalkerList.end()); }


void WalkerConfigurations::createWalkers(int n, size_t numPtcls)
{
  if (WalkerList.empty())
  {
    while (n)
    {
      WalkerList.push_back(std::make_unique<Walker_t>(numPtcls));
      --n;
    }
  }
  else
  {
    if (WalkerList.size() >= n)
    {
      int iw = WalkerList.size(); //copy from the back
      for (int i = 0; i < n; ++i)
      {
        WalkerList.push_back(std::make_unique<Walker_t>(*WalkerList[--iw]));
      }
    }
    else
    {
      int nc  = n / WalkerList.size();
      int nw0 = WalkerList.size();
      for (int iw = 0; iw < nw0; ++iw)
      {
        for (int ic = 0; ic < nc; ++ic)
          WalkerList.push_back(std::make_unique<Walker_t>(*WalkerList[iw]));
      }
      n -= nc * nw0;
      while (n > 0)
      {
        WalkerList.push_back(std::make_unique<Walker_t>(*WalkerList[--nw0]));
        --n;
      }
    }
  }
}


void WalkerConfigurations::resize(int numWalkers, size_t numPtcls)
{
  int dn = numWalkers - WalkerList.size();
  if (dn > 0)
    createWalkers(dn, numPtcls);
  if (dn < 0)
  {
    int nw = -dn;
    if (nw < WalkerList.size())
    {
      WalkerList.erase(WalkerList.begin(), WalkerList.begin() - dn);
    }
  }
}

///returns the next valid iterator
WalkerConfigurations::iterator WalkerConfigurations::destroyWalkers(iterator first, iterator last)
{
  return WalkerList.erase(first, last);
}

void WalkerConfigurations::createWalkers(iterator first, iterator last)
{
  destroyWalkers(WalkerList.begin(), WalkerList.end());
  while (first != last)
  {
    WalkerList.push_back(std::make_unique<Walker_t>(**first));
    ++first;
  }
}

void WalkerConfigurations::destroyWalkers(int nw)
{
  if (nw > WalkerList.size())
  {
    app_warning() << "  Cannot remove walkers. Current Walkers = " << WalkerList.size() << std::endl;
    return;
  }
  nw     = WalkerList.size() - nw;
  int iw = nw;
  WalkerList.erase(WalkerList.begin() + nw, WalkerList.end());
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
  iterator it(WalkerList.begin()), it_end(WalkerList.end());
  while (it != it_end)
  //(*it)->reset();++it;}
  {
    (*it)->Weight       = 1.0;
    (*it)->Multiplicity = 1.0;
    ++it;
  }
}

} // namespace qmcplusplus
