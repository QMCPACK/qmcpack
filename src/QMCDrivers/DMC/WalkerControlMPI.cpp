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
WalkerControlMPI::WalkerControlMPI(Communicate* c) : WalkerControlBase(c)
{
  NumWalkersSent = 0;
  SwapMode       = 1;
  Cur_min        = 0;
  Cur_max        = 0;
  setup_timers(myTimers, DMCMPITimerNames, timer_level_medium);
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
    std::fill(curData.begin(), curData.end(), 0);
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
  MCWalkerConfiguration::iterator it(W.begin()), it_end(W.end());
  while (it != it_end)
  {
    (*it)->Weight       = 1.0;
    (*it)->Multiplicity = 1.0;
    ++it;
  }
  //update the global number of walkers and offsets
  W.setGlobalNumWalkers(Cur_pop);
  W.setWalkerOffsets(FairOffSet);

  return Cur_pop;
}

/** Unified Driver version
 *
 *  It takes 5 steps:
 *    1. calcPopulationAdjustment produces a PopulationAdjustment
 *    2. allreduce collects the number of good walkers + copies on every rank.
 *    3. properties are measured which updates the ensemble properties.
 *    4. adjustPopulation applies nMinNmax if ranks have exceeded min or max limits

 *    5. swapWalkersSimple enacts the population adjustment, sending and recieving walkers.
 *       Receiving side recycles bad walkers' memory first.
 *    6. onRankspawnkill kill's remaining bad walkers, spawns any walkers who have not yet been copied.
 *
 *  In order to reduce the time for allocating walker memory,
 *  this algorithm does not destroy the bad walkers in step 1.
 *  All the dead walkers are reused in step 5 & 6. None are ever GC'd
 */
QMCTraits::FullPrecRealType WalkerControlMPI::branch(int iter, MCPopulation& pop)
{
  ScopedTimer local_timer(myTimers[DMC_MPI_branch]);
  // This has the same ridiculous side effect as SortWalkers
  // i.e. it updates most of curData
  PopulationAdjustment adjust(calcPopulationAdjustment(pop));

  {
    ScopedTimer local_timer(myTimers[DMC_MPI_prebalance]);

    //use NumWalkersSent from the previous exchange
    //You need another copy because curData is zeroed out defensively.
    curData[SENTWALKERS_INDEX] = NumWalkersSent;

    //This should not be used by the new driver code
    //curData[LE_MAX + MyContext] = -1000;
    {
      ScopedTimer local_timer(myTimers[DMC_MPI_allreduce]);
      // You might think we are just reducing LE and sent walkers but
      // see calcPopulationAdjustments massive side effects.
      myComm->allreduce(curData);
    }
    measureProperties(iter);

    pop.set_ensemble_property(ensemble_property_);

    limitPopulation(adjust);
  }

  auto num_per_rank = WalkerControlBase::syncFutureWalkersPerRank(myComm, adjust.num_walkers);

  {
    ScopedTimer local_timer(myTimers[DMC_MPI_loadbalance]);
    NumWalkersSent = swapWalkersSimple(pop, adjust, num_per_rank);
  }

  WalkerControlBase::onRankKill(pop, adjust);
  WalkerControlBase::onRankSpawn(pop, adjust);

  if (adjust.num_walkers != num_per_rank[MyContext])
  {
    std::ostringstream error_message;
    error_message << "failure MPI population control pop.get_num_local_walkers() " << pop.get_num_local_walkers()
                  << " != "
                  << "num_per_rank[" << num_per_rank[MyContext] << "]\n";
    throw std::runtime_error(error_message.str());
  }

  // Update to the current population
  pop.syncWalkersPerRank(myComm);

  for (UPtr<MCPWalker>& walker : pop.get_walkers())
  {
    walker->Weight       = 1.0;
    walker->Multiplicity = 1.0;
  }

  return pop.get_num_global_walkers();
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
  std::vector<Walker_t*> newW;
  std::vector<int> ncopy_newW;
#ifdef MCWALKERSET_MPI_DEBUG
  char fname[128];
  sprintf(fname, "test.%d", MyContext);
  std::ofstream fout(fname, std::ios::app);
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
      Walker_t*& awalker = good_w[ncopy_pairs.back().second];

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
        bad_w.push_back(awalker);
      }
    }
    if (minus[ic] == MyContext)
    {
      Walker_t* awalker(nullptr);
      if (!bad_w.empty())
      {
        awalker = bad_w.back();
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
      newW.push_back(awalker);
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
      Walker_t*& awalker = good_w[jobit->walkerID];
      size_t byteSize    = awalker->byteSize();
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
      Walker_t*& awalker = newW[jobit->walkerID];
      if (!awalker)
        awalker = new Walker_t(wRef);
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
  std::vector<Walker_t*> good_w_temp(good_w);
  good_w.resize(ncopy_pairs.size());
  ncopy_w.resize(ncopy_pairs.size());
  for (int iw = 0; iw < ncopy_pairs.size(); iw++)
  {
    good_w[iw]  = good_w_temp[ncopy_pairs[iw].second];
    ncopy_w[iw] = ncopy_pairs[iw].first;
  }
  //add walkers from other rank
  if (newW.size())
  {
    good_w.insert(good_w.end(), newW.begin(), newW.end());
    ncopy_w.insert(ncopy_w.end(), ncopy_newW.begin(), ncopy_newW.end());
  }

  assert(std::accumulate(ncopy_w.begin(), ncopy_w.end(), ncopy_w.size()) == num_per_rank[MyContext]);
}

/** swap Walkers between rank MCPopulations
 *
 *  MCPopulation is sufficiently different from MCWalkerConfiguration that this is 
 *  basically a rewrite.
 *  \param[inout] pop
 *  \param[inout] adjust
 *
 *  This method should not be dependent on legacy state variables of
 *  Cur_pop, NumPerRank
 *
 */
int WalkerControlMPI::swapWalkersSimple(MCPopulation& pop,
                                        PopulationAdjustment& adjust,
                                        std::vector<IndexType>& num_per_rank)
{
  int expanded_population = std::accumulate(num_per_rank.begin(), num_per_rank.end(), 0);
  std::vector<int> minus, plus;
  determineNewWalkerPopulation(expanded_population, num_contexts_, MyContext, num_per_rank, FairOffSet, plus, minus);

  // local struct for sort vector
  // std::sort requires reference_wrapper for w_elem
  struct CopiesAndWE
  {
    int copies;
    std::reference_wrapper<WalkerElementsRef> w_elem;
  };
  // sort good walkers by the number of copies
  std::vector<CopiesAndWE> sorted_good_walkers;
  for (int iw = 0; iw < adjust.copies_to_make.size(); iw++)
    sorted_good_walkers.push_back(CopiesAndWE{adjust.copies_to_make[iw], adjust.good_walkers[iw]});

  // Sort only on the number of copies
  std::sort(sorted_good_walkers.begin(), sorted_good_walkers.end(),
            [](CopiesAndWE& a, CopiesAndWE& b) { return a.copies < b.copies; });

  int nswap = plus.size();

  std::vector<WalkerMessage> send_message_list;
  // overallocated by number of duplicate messages but message isn't big.
  std::vector<WalkerMessage> recv_message_list;
  std::vector<WalkerElementsRef> new_walkers;
  // Their data needs to not get written over until we are done.
  RefVector<WalkerElementsRef> zombies;
  for (int ic = 0; ic < nswap; ic++)
  {
    if (minus[ic] == MyContext)
    {
      // always send the last good walker or if we're out send the first zombie
      if (!sorted_good_walkers.empty())
      {
        send_message_list.push_back(WalkerMessage{sorted_good_walkers.back().w_elem, minus[ic], plus[ic]});
        --(sorted_good_walkers.back().copies);
        if (sorted_good_walkers.back().copies < 0)
        {
          // Danger possible race condition if this dead walker ends up back in the pool
          // so temporary refvector for those to be killed.
          zombies.push_back(sorted_good_walkers.back().w_elem);
          sorted_good_walkers.pop_back();
        }
      }
      else
      {
        send_message_list.push_back(WalkerMessage{zombies.front(), minus[ic], plus[ic]});
        app_warning() << "Rank " << myComm->rank() << "Had to send best zombie for population control.\n";
      }
    }
    else if (plus[ic] == MyContext)
    {
      if (adjust.bad_walkers.size() > 0)
      {
        pop.killWalker(adjust.bad_walkers.back().walker);
        adjust.bad_walkers.pop_back();
      }
      new_walkers.push_back(pop.spawnWalker());
      recv_message_list.push_back(WalkerMessage{new_walkers.back(), minus[ic], plus[ic]});
    }
  }

  //create send requests
  std::vector<mpi3::request> send_requests;

  if (send_message_list.size() > 0)
  {
    std::for_each(send_message_list.begin(), send_message_list.end(), [&send_requests, this](WalkerMessage& message) {
      MCPWalker& this_walker = message.walker_elements.walker;
      ParticleSet& this_pset = message.walker_elements.pset;
      // Most of these calls are unecessary,
      // evaluateLog definitely is but is invaluable for checking for the state of the walker before and after transfer
      // \todo narrow these down to a minimum and manage to reason out the state of a valid fat walker.
      this_pset.saveWalker(this_walker);
      this_walker.updateBuffer();
      send_requests.emplace_back(myComm->comm.isend_n(message.walker_elements.walker.DataSet.data(),
                                                      message.walker_elements.walker.DataSet.size(),
                                                      message.target_rank));
    });
  }

  //create recv requests
  std::vector<mpi3::request> recv_requests;
  if (recv_message_list.size() > 0)
  {
    std::for_each(recv_message_list.begin(), recv_message_list.end(), [&recv_requests, this](WalkerMessage& message) {
      recv_requests.emplace_back(myComm->comm.ireceive_n(message.walker_elements.walker.DataSet.data(),
                                                         message.walker_elements.walker.DataSet.size(),
                                                         message.source_rank));
      size_t dsize = message.walker_elements.walker.DataSet.size();
    });
  }

  if (recv_message_list.size() > 0)
  {
    ScopedTimer local_timer(myTimers[DMC_MPI_recv]);
    std::vector<int> recv_completed(recv_message_list.size(), 0);
    std::vector<int> recv_waited(recv_message_list.size(), 0);

    while (std::any_of(recv_completed.begin(), recv_completed.end(), [](int i) { return i == 0; }))
    {
      for (int im = 0; im < recv_requests.size(); ++im)
      {
        if (!recv_waited[im])
        {
          //recv_requests[im].wait();
          recv_waited[im] = 1;
        }

        if (!recv_completed[im] && recv_requests[im].completed())
        {
          MCPWalker& this_walker = recv_message_list[im].walker_elements.walker;
          // This sequence of calls is our best effort to go from the wire to a working fat walker.
          // \todo narrow these down to a minimum and manage to reason out the state of a valid fat walker.
          this_walker.copyFromBuffer();
#ifndef NDEBUG
          this_walker.set_has_been_on_wire(true);
#endif
          ParticleSet& this_pset = recv_message_list[im].walker_elements.pset;
          this_pset.loadWalker(this_walker, true);
          // If this update isn't called then the Jastrow's will not match those in the sent walker.
          // The update call is required to update the internal state of pset used by TWF to do the
          // Jastrow evaluations in evaluateLog.
          this_pset.update();
          TrialWaveFunction& this_twf = recv_message_list[im].walker_elements.twf;
          this_twf.evaluateLog(this_pset);
          recv_completed[im] = 1;
        }
      }
    }
  }

  if (send_message_list.size() > 0)
  {
    ScopedTimer local_timer(myTimers[DMC_MPI_send]);
    std::vector<int> send_completed(send_message_list.size(), 0);
    std::vector<int> send_waited(send_message_list.size(), 0);

    // After we've got all our receives wait if we're not done sending.
    while (std::any_of(send_completed.begin(), send_completed.end(), [](int i) { return i == 0; }))
    {
      for (int im = 0; im < send_requests.size(); im++)
      {
        if (!send_waited[im])
        {
          //send_requests[im].wait();
          send_waited[im] = 1;
        }

        send_requests[im].wait();
        if (!send_completed[im] && send_requests[im].completed())
          send_completed[im] = 1;
      }
    }
  }


  std::for_each(zombies.begin(), zombies.end(), [&pop](WalkerElementsRef& zombie) { pop.killWalker(zombie.walker); });
  adjust.good_walkers.clear();
  adjust.copies_to_make.clear();
  for (int iw = 0; iw < sorted_good_walkers.size(); ++iw)
  {
    assert(sorted_good_walkers[iw].copies >= 0);
    adjust.good_walkers.push_back(sorted_good_walkers[iw].w_elem);
    adjust.copies_to_make.push_back(sorted_good_walkers[iw].copies);
  }
  for (int iw = 0; iw < new_walkers.size(); ++iw)
  {
    assert(new_walkers[iw].walker.get_has_been_on_wire());
    adjust.good_walkers.push_back(new_walkers[iw]);
    adjust.copies_to_make.push_back(0);
  }

  adjust.num_walkers =
      std::accumulate(adjust.copies_to_make.begin(), adjust.copies_to_make.end(), adjust.copies_to_make.size());

  //myComm->barrier();

  assert(adjust.num_walkers == num_per_rank[MyContext]);

  // if ( send_message_list.empty() )
  //   return 0;
  // else
  return send_message_list.size();
}


} // namespace qmcplusplus
