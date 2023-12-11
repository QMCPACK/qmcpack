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

#ifndef QMCPLUSPLUS_AFQMC_WALKERCONTROL_HPP
#define QMCPLUSPLUS_AFQMC_WALKERCONTROL_HPP


#include <tuple>
#include <cassert>
#include <memory>
#include <stack>
#include <mpi.h>
#include "AFQMC/config.h"
#include "Utilities/FairDivide.h"

#include "AFQMC/Walkers/WalkerConfig.hpp"
#include "AFQMC/Walkers/WalkerUtilities.hpp"

#include "mpi3/communicator.hpp"
#include "mpi3/request.hpp"

namespace qmcplusplus
{
namespace afqmc
{
/** swap Walkers with Recv/Send
 *
 * The algorithm ensures that the load per node can differ only by one walker.
 * The communication is one-dimensional.
 * Wexcess is an object with multi::array concept which contains walkers beyond the expected
 * pupolation target.
 */
template<class WlkBucket, class Mat, class IVec = std::vector<int>>
inline int swapWalkersSimple(WlkBucket& wset,
                             Mat&& Wexcess,
                             IVec& CurrNumPerNode,
                             IVec& NewNumPerNode,
                             communicator& comm)
{
  int wlk_size = wset.single_walker_size() + wset.single_walker_bp_size();
  int NumContexts, MyContext;
  NumContexts = comm.size();
  MyContext   = comm.rank();
  static_assert(std::decay<Mat>::type::dimensionality == 2, "Wrong dimensionality");
  if (wlk_size != std::get<1>(Wexcess.sizes()))
    throw std::runtime_error("Array dimension error in swapWalkersSimple().");
  if (1 != Wexcess.stride(1))
    throw std::runtime_error("Array shape error in swapWalkersSimple().");
  if (CurrNumPerNode.size() < NumContexts || NewNumPerNode.size() < NumContexts)
    throw std::runtime_error("Array dimension error in swapWalkersSimple().");
  if (wset.capacity() < NewNumPerNode[MyContext])
    throw std::runtime_error("Insufficient capacity in swapWalkersSimple().");
  std::vector<int> minus, plus;
  int deltaN = 0;
  for (int ip = 0; ip < NumContexts; ip++)
  {
    int dn = CurrNumPerNode[ip] - NewNumPerNode[ip];
    if (ip == MyContext)
      deltaN = dn;
    if (dn > 0)
    {
      plus.insert(plus.end(), dn, ip);
    }
    else if (dn < 0)
    {
      minus.insert(minus.end(), -dn, ip);
    }
  }
  int nswap = std::min(plus.size(), minus.size());
  int nsend = 0;
  if (deltaN <= 0 && wset.size() != CurrNumPerNode[MyContext])
    throw std::runtime_error("error in swapWalkersSimple().");
  if (deltaN > 0 && (wset.size() != NewNumPerNode[MyContext] || int(std::get<0>(Wexcess.sizes())) != deltaN))
    throw std::runtime_error("error in swapWalkersSimple().");
  std::vector<ComplexType> buff;
  if (deltaN < 0)
    buff.resize(wlk_size);
  for (int ic = 0; ic < nswap; ic++)
  {
    if (plus[ic] == MyContext)
    {
      comm.send_n(Wexcess[nsend].origin(), Wexcess[nsend].size(), minus[ic], plus[ic] + 999);
      ++nsend;
    }
    if (minus[ic] == MyContext)
    {
      comm.receive_n(buff.data(), buff.size(), plus[ic], plus[ic] + 999);
      wset.push_walkers(boost::multi::array_ref<ComplexType, 2>(buff.data(), {1, wlk_size}));
    }
  }
  return nswap;
}

/** swap Walkers with Irecv/Send
 *
 * The algorithm ensures that the load per node can differ only by one walker.
 * The communication is one-dimensional.
 */
template<class WlkBucket, class Mat, class IVec = std::vector<int>>
// eventually generalize MPI_Comm to a MPI wrapper
inline int swapWalkersAsync(WlkBucket& wset,
                            Mat&& Wexcess,
                            IVec& CurrNumPerNode,
                            IVec& NewNumPerNode,
                            communicator& comm)
{
  int wlk_size = wset.single_walker_size() + wset.single_walker_bp_size();
  int NumContexts, MyContext;
  NumContexts = comm.size();
  MyContext   = comm.rank();
  static_assert(std::decay<Mat>::type::dimensionality == 2, "Wrong dimensionality");
  if (wlk_size != std::get<1>(Wexcess.sizes()))
    throw std::runtime_error("Array dimension error in swapWalkersAsync().");
  if (1 != Wexcess.stride(1) || (std::get<0>(Wexcess.sizes()) > 0 && std::get<1>(Wexcess.sizes()) != Wexcess.stride(0)))
    throw std::runtime_error("Array shape error in swapWalkersAsync().");
  if (CurrNumPerNode.size() < NumContexts || NewNumPerNode.size() < NumContexts)
    throw std::runtime_error("Array dimension error in swapWalkersAsync().");
  if (wset.capacity() < NewNumPerNode[MyContext])
    throw std::runtime_error("Insufficient capacity in swapWalkersAsync().");
  std::vector<int> minus, plus;
  int deltaN = 0;
  for (int ip = 0; ip < NumContexts; ip++)
  {
    int dn = CurrNumPerNode[ip] - NewNumPerNode[ip];
    if (ip == MyContext)
      deltaN = dn;
    if (dn > 0)
    {
      plus.insert(plus.end(), dn, ip);
    }
    else if (dn < 0)
    {
      minus.insert(minus.end(), -dn, ip);
    }
  }
  int nswap     = std::min(plus.size(), minus.size());
  int nsend     = 0;
  int countSend = 1;
  if (deltaN <= 0 && wset.size() != CurrNumPerNode[MyContext])
    throw std::runtime_error("error(1) in swapWalkersAsync().");
  if (deltaN > 0 && (wset.size() != NewNumPerNode[MyContext] || int(std::get<0>(Wexcess.sizes())) != deltaN))
    throw std::runtime_error("error(2) in swapWalkersAsync().");
  std::vector<ComplexType*> buffers;
  std::vector<boost::mpi3::request> requests;
  std::vector<int> recvCounts;
  for (int ic = 0; ic < nswap; ic++)
  {
    if (plus[ic] == MyContext)
    {
      if ((ic < nswap - 1) && (plus[ic] == plus[ic + 1]) && (minus[ic] == minus[ic + 1]))
      {
        countSend++;
      }
      else
      {
        requests.emplace_back(comm.isend(Wexcess[nsend].origin(), Wexcess[nsend].origin() + countSend * std::get<1>(Wexcess.sizes()),
                                         minus[ic], plus[ic] + 1999));
        nsend += countSend;
        countSend = 1;
      }
    }
    if (minus[ic] == MyContext)
    {
      if ((ic < nswap - 1) && (plus[ic] == plus[ic + 1]) && (minus[ic] == minus[ic + 1]))
      {
        countSend++;
      }
      else
      {
        ComplexType* bf = new ComplexType[countSend * wlk_size];
        buffers.push_back(bf);
        recvCounts.push_back(countSend);
        requests.emplace_back(comm.ireceive_n(bf, countSend * wlk_size, plus[ic], plus[ic] + 1999));
        countSend = 1;
      }
    }
  }
  if (deltaN < 0)
  {
    // receiving nodes
    for (int ip = 0; ip < requests.size(); ++ip)
    {
      requests[ip].wait();
      wset.push_walkers(boost::multi::array_ref<ComplexType, 2>(buffers[ip], {recvCounts[ip], wlk_size}));
      delete[] buffers[ip];
    }
  }
  else
  {
    // sending nodes
    for (int ip = 0; ip < requests.size(); ++ip)
      requests[ip].wait();
  }
  return nswap;
}


/**
 * Implements Cafarrel's minimum branching algorithm.
 *   - buff: array of walker info (weight,num).
 */
inline void min_branch(std::vector<std::pair<double, int>>& buff,
                       RandomBase<QMCTraits::FullPrecRealType>& rng,
                       double max_c,
                       double min_c)
{
  APP_ABORT(" Error: min_branch not implemented yet. \n\n\n");
}

/**
 * Implements Cafarrel's minimum branching algorithm.
 *   - buff: array of walker info (weight,num).
 */
inline void serial_comb(std::vector<std::pair<double, int>>& buff, RandomBase<QMCTraits::FullPrecRealType>& rng)
{
  APP_ABORT(" Error: serial_comb not implemented yet. \n\n\n");
}

/**
 * Implements the paired branching algorithm on a popultion of walkers,
 * given a list of walker weights. For each walker in the list, returns the weight
 * and number of times the walker should appear in the new list.
 *   - buff: array of walker info (weight,num).
 */
inline void pair_branch(std::vector<std::pair<double, int>>& buff,
                        RandomBase<QMCTraits::FullPrecRealType>& rng,
                        double max_c,
                        double min_c)
{
  using tp    = std::tuple<double, int, int>;
  using tp_it = std::vector<tp>::iterator;
  // slow for now, not efficient!!!
  int nw = buff.size();
  std::vector<tp> wlks(nw);
  for (int i = 0; i < nw; i++)
    wlks[i] = tp{buff[i].first, 1, i};

  std::sort(wlks.begin(), wlks.end(), [](const tp& a, const tp& b) { return std::get<0>(a) < std::get<0>(b); });

  tp_it it_s = wlks.begin();
  tp_it it_l = wlks.end() - 1;

  while (it_s < it_l)
  {
    if (std::abs(std::get<0>(*it_s)) < min_c || std::abs(std::get<0>(*it_l)) > max_c)
    {
      double w12 = std::get<0>(*it_s) + std::get<0>(*it_l);
      if (rng() < std::get<0>(*it_l) / w12)
      {
        std::get<0>(*it_l) = 0.5 * w12;
        std::get<0>(*it_s) = 0.0;
        std::get<1>(*it_l) = 2;
        std::get<1>(*it_s) = 0;
      }
      else
      {
        std::get<0>(*it_s) = 0.5 * w12;
        std::get<0>(*it_l) = 0.0;
        std::get<1>(*it_s) = 2;
        std::get<1>(*it_l) = 0;
      }
      it_s++;
      it_l--;
    }
    else
      break;
  }

  int nnew  = 0;
  int nzero = 0;
  for (auto& w : wlks)
  {
    buff[std::get<2>(w)] = {std::get<0>(w), std::get<1>(w)};
    nnew += std::get<1>(w);
    if (std::get<1>(w) > 0 && std::abs(std::get<0>(w)) < 1e-7)
      nzero++;
  }
  if (nzero > 0)
  {
    app_error() << " Error in pair_branch: nzero>0: " << nzero << std::endl;
    app_error() << " Found walkers with zero weight after branch.\n"
                << " Try reducing subSteps or reducing the time step." << std::endl;
    APP_ABORT("Error in pair_branch.");
  }
  if (nw != nnew)
    APP_ABORT("Error: Problems with pair_branching.\n");
}

/**
 * Implements the serial branching algorithm on the set of walkers.
 * Serial branch involves gathering the list of weights on the root node
 * and making the decisions locally. The new list of walker weights is then bcasted.
 * This implementation requires contiguous walkers and fixed population walker sets.
 */
template<class WalkerSet,
         class Mat,
         typename = typename std::enable_if<(WalkerSet::contiguous_walker)>::type,
         typename = typename std::enable_if<(WalkerSet::fixed_population)>::type>
inline void SerialBranching(WalkerSet& wset,
                            BRANCHING_ALGORITHM type,
                            double min_,
                            double max_,
                            std::vector<int>& wlk_counts,
                            Mat& Wexcess,
                            RandomBase<QMCTraits::FullPrecRealType>& rng,
                            communicator& comm)
{
  std::vector<std::pair<double, int>> buffer(wset.get_global_target_population());

  // assemble list of weights
  getGlobalListOfWalkerWeights(wset, buffer, comm);

  // using global weight list, use pair branching algorithm
  if (comm.root())
  {
    if (type == PAIR)
      pair_branch(buffer, rng, max_, min_);
    else if (type == MIN_BRANCH)
      min_branch(buffer, rng, max_, min_);
    else if (type == SERIAL_COMB)
      serial_comb(buffer, rng);
    else
      APP_ABORT("Error: Unknown branching type in SerialBranching. \n");
  }

  // bcast walker information and calculate new walker counts locally
  comm.broadcast_n(buffer.data(), buffer.size());

  int target = wset.get_TG_target_population();
  wlk_counts.resize(comm.size());
  for (int i = 0, p = 0; i < comm.size(); i++)
  {
    int cnt = 0;
    for (int k = 0; k < target; k++, p++)
      cnt += buffer[p].second;
    wlk_counts[i] = cnt;
  }
  if (wset.get_global_target_population() != std::accumulate(wlk_counts.begin(), wlk_counts.end(), 0))
  {
    app_error() << " Error: targetN != nwold: " << target << " "
                << std::accumulate(wlk_counts.begin(), wlk_counts.end(), 0) << std::endl;
    APP_ABORT(" Error: targetN != nwold.");
  }

  // reserve space for extra walkers
  if (wlk_counts[comm.rank()] > target)
    Wexcess.reextent(
        {std::max(0, wlk_counts[comm.rank()] - target), wset.single_walker_size() + wset.single_walker_bp_size()});

  // perform local branching
  // walkers beyond target go in Wexcess
  wset.branch(buffer.begin() + target * comm.rank(), buffer.begin() + target * (comm.rank() + 1), Wexcess);
}

/**
 * Implements the distributed comb branching algorithm.
 */
template<class WalkerSet,
         class Mat,
         typename = typename std::enable_if<(WalkerSet::contiguous_walker)>::type,
         typename = typename std::enable_if<(WalkerSet::fixed_population)>::type>
inline void CombBranching(WalkerSet& wset,
                          BRANCHING_ALGORITHM type,
                          std::vector<int>& wlk_counts,
                          Mat& Wexcess,
                          RandomBase<QMCTraits::FullPrecRealType>& rng,
                          communicator& comm)
{
  APP_ABORT("Error: comb not implemented yet. \n");
}

} // namespace afqmc

} // namespace qmcplusplus

#endif
