//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2020 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Lab
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


#include <array>
#include <cmath>
#include <sstream>

#include "WalkerControlMPI.h"
#include "Utilities/IteratorUtility.h"
#include "Utilities/FairDivide.h"
#include "QMCWaveFunctions/TrialWaveFunction.h"

namespace qmcplusplus
{
//#define MCWALKERSET_MPI_DEBUG

enum DMC_MPI_Timers
{
  DMC_MPI_branch,
  DMC_MPI_imbalance,
  DMC_MPI_prebalance,
  DMC_MPI_copyWalkers,
  DMC_MPI_allreduce,
  DMC_MPI_loadbalance,
  DMC_MPI_send,
  DMC_MPI_recv,
};

TimerNameList_t<DMC_MPI_Timers> DMCMPITimerNames = {{DMC_MPI_branch, "WalkerControlMPI::branch"},
                                                    {DMC_MPI_imbalance, "WalkerControlMPI::imbalance"},
                                                    {DMC_MPI_prebalance, "WalkerControlMPI::pre-loadbalance"},
                                                    {DMC_MPI_copyWalkers, "WalkerControlMPI::copyWalkers"},
                                                    {DMC_MPI_allreduce, "WalkerControlMPI::allreduce"},
                                                    {DMC_MPI_loadbalance, "WalkerControlMPI::loadbalance"},
                                                    {DMC_MPI_send, "WalkerControlMPI::send"},
                                                    {DMC_MPI_recv, "WalkerControlMPI::recv"}};

/** default constructor
 *
 * set SwapMode? SwapMode is set to 1 but what does that mean?
 * This object persists inside the SFNB which also persists
 * The zeroing here will not happen in later QMC sections...
 * This seems problematic in that NumWalkersSent will start at a 
 * value of no concern to the current section.
 *
 * In the new drivers SFNB should throw an except if there is attempted 
 * reuse of WalkerController
 */
WalkerControlMPI::WalkerControlMPI(Communicate* c)
    : WalkerControlBase(c), myTimers(getGlobalTimerManager(), DMCMPITimerNames, timer_level_medium)
{
  NumWalkersSent = 0;
  SwapMode       = 1;
  Cur_min        = 0;
  Cur_max        = 0;
}

/** Perform branch and swap walkers as required
 *
 *  It takes 5 steps:
 *    1. sortWalkers marks good and bad walkers.
 *    2. allreduce collects the number of good walkers + copies on every rank.
 *    3. applyNmaxNmin avoids too large or too small global population.
 *    4. swapWalkersSimple makes a decision of load balancing and send/recv walkers.
 *       Receiving side recycles bad walkers' memory first.
 *    5. copyWalkers generates copies of good walkers.
 *  In order to minimize the memory footprint fluctuation
 *  the walker copying is placed as the last step.
 *  In order to reduce the time for allocating walker memory,
 *  this algorithm does not destroy the bad walkers in step 1.
 *  All the bad walkers are recycled as much as possible in step 3/4.
 */
int WalkerControlMPI::branch(int iter, MCWalkerConfiguration& W, FullPrecRealType trigger)
{
  ScopedTimer local_timer(myTimers[DMC_MPI_branch]);
  {
    ScopedTimer local_timer(myTimers[DMC_MPI_prebalance]);
    std::fill(curData.begin(), curData.end(), 0.0);
    sortWalkers(W);
    //use NumWalkersSent from the previous exchange
    curData[SENTWALKERS_INDEX] = NumWalkersSent;
    //update the number of walkers for this rank
    //Causes implicit conversion to FullPrecRealType
    curData[LE_MAX + MyContext] = NumWalkers;
    //{ ScopedTimer local_timer(myTimers[DMC_MPI_imbalance]);
    //}
    {
      ScopedTimer local_timer(myTimers[DMC_MPI_allreduce]);
      myComm->allreduce(curData);
    }
    measureProperties(iter);
    W.EnsembleProperty = ensemble_property_;
    for (int i = 0, j = LE_MAX; i < num_contexts_; i++, j++)
      NumPerRank[i] = static_cast<int>(curData[j]);
    int current_population = std::accumulate(NumPerRank.begin(), NumPerRank.end(), 0);

    Cur_pop = applyNmaxNmin(current_population);
  }
  {
    ScopedTimer local_timer(myTimers[DMC_MPI_loadbalance]);
    swapWalkersSimple(W);
  }
  {
    ScopedTimer local_timer(myTimers[DMC_MPI_copyWalkers]);
    copyWalkers(W);
  }
  //set Weight and Multiplicity to default values
  for (auto& walker : W)
  {
    walker->Weight       = 1.0;
    walker->Multiplicity = 1.0;
  }
  //update the walkers offsets
  W.setWalkerOffsets(FairOffSet);

  return Cur_pop;
}

// determine new walker population on each MPI rank
void WalkerControlMPI::determineNewWalkerPopulation(int cur_pop,
                                                    int num_contexts,
                                                    int my_context,
                                                    std::vector<int>& num_per_rank,
                                                    std::vector<int>& fair_offset,
                                                    std::vector<int>& minus,
                                                    std::vector<int>& plus)
{
  FairDivideLow(cur_pop, num_contexts, fair_offset);
  for (int ip = 0; ip < num_contexts; ip++)
  {
    // (FairOffSet[ip + 1] - FairOffSet[ip]) gives the partiion ip walker pop
    int dn = num_per_rank[ip] - (fair_offset[ip + 1] - fair_offset[ip]);
    num_per_rank[ip] -= dn;
    if (dn > 0)
    {
      plus.insert(plus.end(), dn, ip);
    }
    else if (dn < 0)
    {
      minus.insert(minus.end(), -dn, ip);
    }
  }
#ifndef NDEBUG
  if (plus.size() != minus.size())
  {
    app_error() << "Walker send/recv pattern doesn't match. "
                << "The send size " << plus.size() << " is not equal to the recv size " << minus.size() << " ."
                << std::endl;
    throw std::runtime_error("Trying to swap in WalkerControlMPI::swapWalkersSimple with mismatched queues");
  }
#endif
}

/** swap Walkers with Recv/Send or Irecv/Isend
 *
 * The algorithm ensures that the load per rank can differ only by one walker.
 * Each MPI rank can only send or receive or be silent.
 * The communication is one-dimensional and very local.
 * If multiple copies of a walker need to be sent to the target rank, only send one.
 * The number of copies is communicated ahead via blocking send/recv.
 * Then the walkers are transferred via blocking or non-blocking send/recv.
 * The blocking send/recv may become serialized and worsen load imbalance.
 * Non blocking send/recv algorithm avoids serialization completely.
 */
void WalkerControlMPI::swapWalkersSimple(MCWalkerConfiguration& W)
{
  std::vector<int> minus, plus;
  //legacy code does not modify NumPerRank in this call so we copy NumPerRank
  std::vector<int> num_per_rank(NumPerRank);
  determineNewWalkerPopulation(Cur_pop, num_contexts_, MyContext, num_per_rank, FairOffSet, minus, plus);

  if (good_w.empty() && bad_w.empty())
  {
    app_error() << "It should never happen that no walkers, "
                << "neither good nor bad, exist on a rank. "
                << "Please report to developers. " << std::endl;
    APP_ABORT("WalkerControlMPI::swapWalkersSimple no existing walker");
  }

  Walker_t& wRef(*(good_w.empty() ? bad_w[0] : good_w[0]));
  std::vector<std::unique_ptr<Walker_t>> newW;
  std::vector<int> ncopy_newW;
#ifdef MCWALKERSET_MPI_DEBUG
  std::array<char, 128> fname;
  if (std::snprintf(fname.data(), fname.size() "test.%d", MyContext) < 0)
    throw std::runtime_error("Error generating filename");
  std::ofstream fout(fname.data(), std::ios::app);
  //fout << NumSwaps << " " << Cur_pop << " ";
  //for(int ic=0; ic<NumContexts; ic++) fout << NumPerRank[ic] << " ";
  //fout << " | ";
  //for(int ic=0; ic<NumContexts; ic++) fout << FairOffSet[ic+1]-FairOffSet[ic] << " ";
  //fout << " | ";
  for (int ic = 0; ic < plus.size(); ic++)
  {
    fout << plus[ic] << " ";
  }
  fout << " | ";
  for (int ic = 0; ic < minus.size(); ic++)
  {
    fout << minus[ic] << " ";
  }
  fout << std::endl;
#endif
  int nswap = plus.size();
  // sort good walkers by the number of copies
  assert(good_w.size() == ncopy_w.size());
  std::vector<std::pair<int, int>> ncopy_pairs;
  for (int iw = 0; iw < ncopy_w.size(); iw++)
    ncopy_pairs.push_back(std::make_pair(ncopy_w[iw], iw));
  std::sort(ncopy_pairs.begin(), ncopy_pairs.end());

  int nsend = 0;
  struct job
  {
    const int walkerID;
    const int target;
    job(int wid, int target_in) : walkerID(wid), target(target_in){};
  };
  std::vector<job> job_list;
  for (int ic = 0; ic < nswap; ic++)
  {
    if (plus[ic] == MyContext)
    {
      // always send the last good walker
      auto& awalker = good_w[ncopy_pairs.back().second];
      // count the possible copies in one send
      int nsentcopy = 0;

      for (int id = ic + 1; id < nswap; id++)
        if (plus[ic] == plus[id] && minus[ic] == minus[id] && ncopy_pairs.back().first > 0)
        { // increment copy counter
          ncopy_pairs.back().first--;
          nsentcopy++;
        }
        else
        { // not enough copies to send or not the same send/recv pair
          break;
        }

      // send the number of copies to the target
      myComm->comm.send_value(nsentcopy, minus[ic]);
      job_list.push_back(job(ncopy_pairs.back().second, minus[ic]));
#ifdef MCWALKERSET_MPI_DEBUG
      fout << "rank " << plus[ic] << " sends a walker with " << nsentcopy << " copies to rank " << minus[ic]
           << std::endl;
#endif

      // update counter and cursor
      ++nsend;
      ic += nsentcopy;

      // update copy counter
      if (ncopy_pairs.back().first > 0)
      {
        ncopy_pairs.back().first--;
        std::sort(ncopy_pairs.begin(), ncopy_pairs.end());
      }
      else
      {
        ncopy_pairs.pop_back();
        bad_w.push_back(std::make_unique<Walker_t>(*awalker));
      }
    }
    if (minus[ic] == MyContext)
    {
      std::unique_ptr<Walker_t> awalker;
      if (!bad_w.empty())
      {
        awalker = std::move(bad_w.back());
        bad_w.pop_back();
      }

      int nsentcopy = 0;
      // recv the number of copies from the target
      myComm->comm.receive_n(&nsentcopy, 1, plus[ic]);
      job_list.push_back(job(newW.size(), plus[ic]));
      if (plus[ic] != plus[ic + nsentcopy] || minus[ic] != minus[ic + nsentcopy])
        APP_ABORT("WalkerControlMPI::swapWalkersSimple send/recv pair checking failed!");
#ifdef MCWALKERSET_MPI_DEBUG
      fout << "rank " << minus[ic] << " recvs a walker with " << nsentcopy << " copies from rank " << plus[ic]
           << std::endl;
#endif

      // save the new walker
      if (awalker)
      {
        newW.push_back(std::make_unique<Walker_t>(*awalker));
      }
      else
      {
        newW.push_back(nullptr);
      }
      ncopy_newW.push_back(nsentcopy);
      // update cursor
      ic += nsentcopy;
    }
  }

  if (nsend > 0)
  {
    std::vector<mpi3::request> requests;
    // mark all walkers not in send
    for (auto jobit = job_list.begin(); jobit != job_list.end(); jobit++)
      good_w[jobit->walkerID]->SendInProgress = false;
    for (auto jobit = job_list.begin(); jobit != job_list.end(); jobit++)
    {
      // pack data and send
      auto& awalker   = good_w[jobit->walkerID];
      size_t byteSize = awalker->byteSize();
      if (!awalker->SendInProgress)
      {
        awalker->updateBuffer();
        awalker->SendInProgress = true;
      }
      if (use_nonblocking)
        requests.push_back(myComm->comm.isend_n(awalker->DataSet.data(), byteSize, jobit->target));
      else
      {
        ScopedTimer local_timer(myTimers[DMC_MPI_send]);
        myComm->comm.send_n(awalker->DataSet.data(), byteSize, jobit->target);
      }
    }
    if (use_nonblocking)
    {
      // wait all the isend
      for (int im = 0; im < requests.size(); im++)
      {
        ScopedTimer local_timer(myTimers[DMC_MPI_send]);
        requests[im].wait();
      }
      requests.clear();
    }
  }
  else
  {
    std::vector<mpi3::request> requests;
    for (auto jobit = job_list.begin(); jobit != job_list.end(); jobit++)
    {
      // recv and unpack data
      auto& awalker = newW[jobit->walkerID];
      if (!awalker)
        awalker = std::make_unique<Walker_t>(wRef);
      size_t byteSize = awalker->byteSize();
      if (use_nonblocking)
        requests.push_back(myComm->comm.ireceive_n(awalker->DataSet.data(), byteSize, jobit->target));
      else
      {
        ScopedTimer local_timer(myTimers[DMC_MPI_recv]);
        myComm->comm.receive_n(awalker->DataSet.data(), byteSize, jobit->target);
        awalker->copyFromBuffer();
      }
    }
    if (use_nonblocking)
    {
      std::vector<bool> not_completed(requests.size(), true);
      bool completed = false;
      while (!completed)
      {
        completed = true;
        for (int im = 0; im < requests.size(); im++)
          if (not_completed[im])
          {
            if (requests[im].completed())
            {
              newW[job_list[im].walkerID]->copyFromBuffer();
              not_completed[im] = false;
            }
            else
              completed = false;
          }
      }
      requests.clear();
    }
  }
  //save the number of walkers sent
  NumWalkersSent = nsend;
  // rebuild good_w and ncopy_w
  std::vector<std::unique_ptr<Walker_t>> good_w_temp(std::move(good_w));
  good_w.resize(ncopy_pairs.size());
  ncopy_w.resize(ncopy_pairs.size());
  for (int iw = 0; iw < ncopy_pairs.size(); iw++)
  {
    good_w[iw]  = std::move(good_w_temp[ncopy_pairs[iw].second]);
    ncopy_w[iw] = ncopy_pairs[iw].first;
  }
  //add walkers from other rank
  if (newW.size())
  {
    good_w.insert(good_w.end(), std::make_move_iterator(newW.begin()), std::make_move_iterator(newW.end()));
    ncopy_w.insert(ncopy_w.end(), ncopy_newW.begin(), ncopy_newW.end());
  }

  assert(std::accumulate(ncopy_w.begin(), ncopy_w.end(), ncopy_w.size()) == num_per_rank[MyContext]);
}

} // namespace qmcplusplus
