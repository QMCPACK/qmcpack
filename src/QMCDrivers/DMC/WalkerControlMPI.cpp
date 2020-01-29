//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2019 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Lab
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


#include <cmath>
#include <queue>
#include <QMCDrivers/DMC/WalkerControlMPI.h>
#include <Utilities/IteratorUtility.h>
#include <Utilities/FairDivide.h>
#include <Utilities/NewTimer.h>
#include <Utilities/Timer.h>

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
 * The zeroing here will not happen in late QMC sections...
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
  myTimers[DMC_MPI_branch]->start();
  myTimers[DMC_MPI_prebalance]->start();
  std::fill(curData.begin(), curData.end(), 0);
  sortWalkers(W);
  //use NumWalkersSent from the previous exchange
  curData[SENTWALKERS_INDEX] = NumWalkersSent;
  //update the number of walkers for this node
  //Causes implicit conversion to FullPrecRealType
  curData[LE_MAX + MyContext] = NumWalkers;
  //myTimers[DMC_MPI_imbalance]->start();
  //myComm->barrier();
  //myTimers[DMC_MPI_imbalance]->stop();
  myTimers[DMC_MPI_allreduce]->start();
  myComm->allreduce(curData);
  myTimers[DMC_MPI_allreduce]->stop();
  measureProperties(iter);
  W.EnsembleProperty = ensemble_property_;
  for (int i = 0, j = LE_MAX; i < num_contexts_; i++, j++)
    NumPerNode[i] = static_cast<int>(curData[j]);
  int current_population = std::accumulate(NumPerNode.begin(), NumPerNode.end(), 0);

  Cur_pop = applyNmaxNmin(current_population);
  myTimers[DMC_MPI_prebalance]->stop();
  myTimers[DMC_MPI_loadbalance]->start();
  swapWalkersSimple(W);
  myTimers[DMC_MPI_loadbalance]->stop();
  myTimers[DMC_MPI_copyWalkers]->start();
  copyWalkers(W);
  myTimers[DMC_MPI_copyWalkers]->stop();
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

  myTimers[DMC_MPI_branch]->stop();
  return Cur_pop;
}

/** Unified Driver version
 *
 *  It takes 5 steps:
 *    1. calcPopulationAdjustment produces a PopulationAdjustment
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
int WalkerControlMPI::branch(int iter, MCPopulation& pop, FullPrecRealType trigger)
{
  myTimers[DMC_MPI_branch]->start();
  myTimers[DMC_MPI_prebalance]->start();
  std::fill(curData.begin(), curData.end(), 0);
  // This has the same ridiculous side effect as SortWalkers
  // i.e. it updates most of curData
  PopulationAdjustment adjust(calcPopulationAdjustment(pop));

  //use NumWalkersSent from the previous exchange
  //You need another copy because curData is zeroed out defensively.
  curData[SENTWALKERS_INDEX] = NumWalkersSent;

  //This should not be used by the new driver code
  //curData[LE_MAX + MyContext] = -1000;
  myTimers[DMC_MPI_allreduce]->start();
  // You might think we are just reducing LE and sent walkers but
  // see calcPopulationAdjustments massive side effects.
  myComm->allreduce(curData);
  myTimers[DMC_MPI_allreduce]->stop();
  measureProperties(iter);
  pop.set_ensemble_property(ensemble_property_);

  //All of this should really just accomplish what onRankSpawnKill does for a nonmpi job.
  auto num_per_node = WalkerControlBase::syncFutureWalkersPerRank(this->getCommunicator(), adjust.num_walkers);

  myTimers[DMC_MPI_prebalance]->stop();
  myTimers[DMC_MPI_loadbalance]->start();
  swapWalkersSimple(pop, adjust, num_per_node);
  myTimers[DMC_MPI_loadbalance]->stop();

  adjustPopulation(pop, adjust);

  onRankSpawnKill(pop, adjust);
  
  for (UPtr<MCPWalker>& walker : pop.get_walkers())
  {
    walker->Weight       = 1.0;
    walker->Multiplicity = 1.0;
  }
  // //update the global number of walkers and offsets

  myComm->allreduce(curData);

  // Discover the population now
  pop.syncWalkersPerNode(getCommunicator());

  myTimers[DMC_MPI_branch]->stop();

  return pop.get_num_global_walkers();
}

// determine new walker population on each node
void WalkerControlMPI::determineNewWalkerPopulation(int Cur_pop,
                                                    int NumContexts,
                                                    int MyContext,
                                                    const std::vector<int>& NumPerNode,
                                                    std::vector<int>& FairOffSet,
                                                    std::vector<int>& minus,
                                                    std::vector<int>& plus)
{
  // Cur_pop - in - current population
  // NumContexts -in - number of MPI processes
  // MyContext - in - my MPI rank
  // NumPerNode - in - current walkers per node
  // FairOffSet - out - running population count at each partition boundary
  // minus - out - list of partition indexes one occurance for each walker removed
  // plus -  out - list of partition indexes one occurance for each walker added

  FairDivideLow(Cur_pop, NumContexts, FairOffSet);
  for (int ip = 0; ip < NumContexts; ip++)
  {
    // (FairOffSet[ip + 1] - FairOffSet[ip]) gives the partiion ip walker pop
    int dn = NumPerNode[ip] - (FairOffSet[ip + 1] - FairOffSet[ip]);
    if (dn > 0)
    {
      plus.insert(plus.end(), dn, ip);
    }
    else if (dn < 0)
    {
      minus.insert(minus.end(), -dn, ip);
    }
  }
}

/** swap Walkers with Recv/Send or Irecv/Isend
 *
 * The algorithm ensures that the load per node can differ only by one walker.
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
  determineNewWalkerPopulation(Cur_pop, num_contexts_, MyContext, NumPerNode, FairOffSet, minus, plus);

  if (good_w.empty() && bad_w.empty())
  {
    app_error() << "It should never happen that no walkers, "
                << "neither good nor bad, exist on a node. "
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
  //for(int ic=0; ic<NumContexts; ic++) fout << NumPerNode[ic] << " ";
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
  if (plus.size() != minus.size())
  {
    app_error() << "Walker send/recv pattern doesn't match. "
                << "The send size " << plus.size() << " is not equal to the receive size " << minus.size() << " ."
                << std::endl;
    APP_ABORT("WalkerControlMPI::swapWalkersSimple");
  }
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
        myTimers[DMC_MPI_send]->start();
        myComm->comm.send_n(awalker->DataSet.data(), byteSize, jobit->target);
        myTimers[DMC_MPI_send]->stop();
      }
    }
    if (use_nonblocking)
    {
      // wait all the isend
      for (int im = 0; im < requests.size(); im++)
      {
        myTimers[DMC_MPI_send]->start();
        requests[im].wait();
        myTimers[DMC_MPI_send]->stop();
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
        myTimers[DMC_MPI_recv]->start();
        myComm->comm.receive_n(awalker->DataSet.data(), byteSize, jobit->target);
        awalker->copyFromBuffer();
        myTimers[DMC_MPI_recv]->stop();
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
  //add walkers from other node
  if (newW.size())
  {
    good_w.insert(good_w.end(), newW.begin(), newW.end());
    ncopy_w.insert(ncopy_w.end(), ncopy_newW.begin(), ncopy_newW.end());
  }
}

/** swap Walkers between rank MCPopulations
 *
 *  MCPopulation is sufficiently different from MCWalkerConfiguration that this is 
 *  basically a rewrite.
 *  \param[inout] pop
 *  \param[inout] adjust
 *
 *  This method should not be dependent on legacy state variables of
 *  Cur_pop, NumPernode
 *
 */
void WalkerControlMPI::swapWalkersSimple(MCPopulation& pop,
                                         PopulationAdjustment& adjust,
                                         std::vector<IndexType> num_per_node)
{
  int expanded_population = std::accumulate(num_per_node.begin(), num_per_node.end(), 0);
  std::vector<int> minus, plus;
  determineNewWalkerPopulation(expanded_population, num_contexts_, MyContext, num_per_node, FairOffSet, plus, minus);

  if (adjust.good_walkers.empty() && adjust.bad_walkers.empty())
  {
    app_error() << "It should never happen that no walkers, "
                << "neither good nor bad, exist on a node. "
                << "Please report to developers. " << std::endl;
    APP_ABORT("WalkerControlMPI::swapWalkersSimple no existing walker");
  }

  // Looks like a just in case update should be justified.
  for (MCPWalker& walker : adjust.good_walkers)
    walker.updateBuffer();

  if (plus.size() != minus.size())
  {
    app_error() << "Walker send/recv pattern doesn't match. "
                << "The send size " << plus.size() << " is not equal to the recv size " << minus.size() << " ."
                << std::endl;
    throw std::runtime_error("Trying to swap in WalkerControlMPI::swapWalkersSimple with mismatched queues");
  }

  // sort good walkers by the number of copies
  std::vector<std::pair<int, MCPWalker&>> sorted_good_walkers;
  for (int iw = 0; iw < adjust.copies_to_make.size(); iw++)
    sorted_good_walkers.push_back(std::make_pair(adjust.copies_to_make[iw], adjust.good_walkers[iw]));

  // Sort only on the number of copies
  std::sort(sorted_good_walkers.begin(), sorted_good_walkers.end(), [](auto& a, auto& b) { return a.first < b.first; });

  //useful counts
  int nswap       = plus.size();
  int local_sends = 0;
  int local_recvs = 0;

  for (int ic = 0; ic < nswap; ic++)
  {
    if (minus[ic] == MyContext)
      ++local_sends;
    else if (plus[ic] == MyContext)
      ++local_recvs;
  }

  std::vector<WalkerMessage> send_message_list;
  // overallocated by number of duplicate messages but message isn't big.
  send_message_list.reserve(local_sends);
  std::vector<WalkerMessage> recv_message_list;
  recv_message_list.reserve(local_recvs);
  RefVector<MCPWalker> new_walkers;
  new_walkers.reserve(local_recvs);
  // Their data needs to not get written over until we are done.
  RefVector<MCPWalker> zombies;
  for (int ic = 0; ic < nswap; ic++)
  {
    if (minus[ic] == MyContext)
    {
      // always send the last good walker or if we're out send the first zombie
      if (!sorted_good_walkers.empty())
      {
        send_message_list.push_back(WalkerMessage{sorted_good_walkers.back().second, minus[ic], plus[ic]});
        --(sorted_good_walkers.back().first);
        if (sorted_good_walkers.back().first < 0)
        {
          // Danger possible race condition if this dead walker ends up back in the pool
          // so temporary refvector for those to be killed.
          zombies.push_back(sorted_good_walkers.back().second);
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
        pop.killWalker(adjust.bad_walkers.back());
        adjust.bad_walkers.pop_back();
      }
      MCPWalker& spawned_walker = *(pop.spawnWalker());
      new_walkers.push_back(spawned_walker);

      recv_message_list.push_back(WalkerMessage{new_walkers.back(), minus[ic], plus[ic]});
    }
  }

  //create send requests
  std::vector<mpi3::request> send_requests;

  if (send_message_list.size() > 0)
  {
    std::for_each(send_message_list.begin(), send_message_list.end(), [&send_requests, this](WalkerMessage& message) {
      MCPWalker& this_walker = message.walker;
      send_requests.emplace_back(myComm->comm.isend_n(message.walker.DataSet.data(),
                                                      message.walker.DataSet.size(), message.target_rank));
    });
  }

  //create recv requests
  std::vector<mpi3::request> recv_requests;
  if (recv_message_list.size() > 0)
  {
    std::for_each(recv_message_list.begin(), recv_message_list.end(), [&recv_requests, this](WalkerMessage& message) {
      MCPWalker& walker = message.walker;
      recv_requests.emplace_back(myComm->comm.ireceive_n(message.walker.DataSet.data(),
                                                         message.walker.DataSet.size(),
                                                         message.source_rank));
    });
  }

  RefVector<MCPWalker> recv_walkers;
  if (local_recvs > 0)
  {
    for (int im = 0; im < recv_requests.size(); ++im)
    {
      recv_requests[im].wait();
      MCPWalker& walker_to_check = recv_message_list[im].walker;
      recv_message_list[im].walker.copyFromBuffer();
      // for (auto property : walker_to_check.Properties)
      // {
      //   if (std::isnan(property))
      //     throw std::runtime_error("received property is nan!");
      // }
      assert( walker_to_check.Multiplicity > 0 );
      assert( walker_to_check.Weight > 0 );
      
      recv_walkers.push_back(walker_to_check);
    }
  }

  if (local_sends > 0)
  {
    std::vector<int> send_completed(local_sends, 0);
    // After we've got all our receives wait if we're not done sending.
    myTimers[DMC_MPI_send]->start();
    while(std::any_of(send_completed.begin(), send_completed.end(), [](int i){ return i == 0; }))
    {
      for (int im = 0; im < send_requests.size(); im++)
      {
        send_requests[im].wait();
        if( ! send_completed[im] && send_requests[im].completed() )
          send_completed[im] = 1;        
      }
    }
    myTimers[DMC_MPI_send]->stop();
  }

  std::for_each(zombies.begin(), zombies.end(), [&pop](MCPWalker& zombie) { pop.killWalker(zombie); });

  adjust.good_walkers.clear();
  adjust.copies_to_make.clear();
  for(int iw = 0; iw < sorted_good_walkers.size(); ++iw)
  {
    adjust.good_walkers.push_back(sorted_good_walkers[iw].second);
    adjust.copies_to_make.push_back(sorted_good_walkers[iw].first);
  }
  for(int iw = 0; iw < recv_walkers.size(); ++iw)
  {
    adjust.good_walkers.push_back(recv_walkers[iw]);
    adjust.copies_to_make.push_back(0);
  }
  adjust.num_walkers = std::accumulate(adjust.copies_to_make.begin(),adjust.copies_to_make.end(),adjust.copies_to_make.size());

  NumWalkersSent = local_sends;
}


} // namespace qmcplusplus
