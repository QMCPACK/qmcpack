//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2020 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//                    Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////

#include <cassert>
#include <stdexcept>
#include <numeric>
#include <sstream>

#include "WalkerControl.h"
#include "QMCDrivers/WalkerProperties.h"
#include "Particle/HDFWalkerIO.h"
#include "OhmmsData/ParameterSet.h"
#include "type_traits/template_types.hpp"

namespace qmcplusplus
{
using WP = WalkerProperties::Indexes;

enum WC_Timers
{
  WC_branch,
  WC_imbalance,
  WC_prebalance,
  WC_copyWalkers,
  WC_allreduce,
  WC_loadbalance,
  WC_send,
  WC_recv,
};

TimerNameList_t<WC_Timers> WalkerControlTimerNames = {{WC_branch, "WalkerControlMPI::branch"},
                                                           {WC_imbalance, "WalkerControlMPI::imbalance"},
                                                           {WC_prebalance, "WalkerControlMPI::pre-loadbalance"},
                                                           {WC_copyWalkers, "WalkerControlMPI::copyWalkers"},
                                                           {WC_allreduce, "WalkerControlMPI::allreduce"},
                                                           {WC_loadbalance, "WalkerControlMPI::loadbalance"},
                                                           {WC_send, "WalkerControlMPI::send"},
                                                           {WC_recv, "WalkerControlMPI::recv"}};

WalkerControl::WalkerControl(Communicate* c, RandomGenerator_t& rng, bool use_fixed_pop)
    : MPIObjectBase(c),
      rng_(rng),
      use_fixed_pop_(use_fixed_pop),
      n_min_(1),
      n_max_(10),
      MaxCopy(2),
      target_sigma_(10),
      dmcStream(0),
      NumWalkersCreated(0),
      SwapMode(0),
      NumWalkersSent(0)
{
  num_contexts_ = myComm->size();
  MyContext     = myComm->rank();
  curData.resize(LE_MAX + num_contexts_);
  NumPerNode.resize(num_contexts_);
  OffSet.resize(num_contexts_ + 1);
  FairOffSet.resize(num_contexts_ + 1);

  setup_timers(myTimers, WalkerControlTimerNames, timer_level_medium);
}

WalkerControl::~WalkerControl()
{
  if (dmcStream)
  {
    // without this its possible to end up without all data flushed to dmc.dat.
    (*dmcStream) << std::endl;
    delete dmcStream;
  }
}

void WalkerControl::start()
{
  if (MyContext == 0)
  {
    std::string hname(myComm->getName());
    hname.append(".dmc.dat");
    if (hname != dmcFname)
    {
      if (dmcStream)
      {
        *dmcStream << std::endl;
        delete dmcStream;
      }
      dmcStream = new std::ofstream(hname.c_str());
      //oa = new boost::archive::binary_oarchive (*dmcStream);
      dmcStream->setf(std::ios::scientific, std::ios::floatfield);
      dmcStream->precision(10);
      (*dmcStream) << "# Index " << std::setw(20) << "LocalEnergy" << std::setw(20) << "Variance" << std::setw(20)
                   << "Weight" << std::setw(20) << "NumOfWalkers" << std::setw(20)
                   << "AvgSentWalkers"; //add the number of walkers
      (*dmcStream) << std::setw(20) << "TrialEnergy" << std::setw(20) << "DiffEff";
      (*dmcStream) << std::setw(20) << "LivingFraction";
      (*dmcStream) << std::endl;
      dmcFname = hname;
    }
  }
}

void WalkerControl::setWalkerID(MCPopulation& population)
{
  start(); //do the normal start
}

/** Depends on alot of state
 *
 *  Does not depend on state refactored to  PopulationAdjustment directly
 */
void WalkerControl::measureProperties(int iter)
{
  //taking average over the walkers
  FullPrecRealType wgtInv(1.0 / curData[WEIGHT_INDEX]);
  FullPrecRealType eavg         = curData[ENERGY_INDEX] * wgtInv;
  ensemble_property_.Energy     = eavg;
  ensemble_property_.Weight     = curData[WEIGHT_INDEX];
  ensemble_property_.Variance   = (curData[ENERGY_SQ_INDEX] * wgtInv - eavg * eavg);
  ensemble_property_.NumSamples = curData[WALKERSIZE_INDEX];
  ensemble_property_.R2Accepted = curData[R2ACCEPTED_INDEX];
  ensemble_property_.R2Proposed = curData[R2PROPOSED_INDEX];
  ensemble_property_.LivingFraction =
      static_cast<FullPrecRealType>(curData[FNSIZE_INDEX]) / static_cast<FullPrecRealType>(curData[WALKERSIZE_INDEX]);
  ensemble_property_.AlternateEnergy = FullPrecRealType(0);
  ensemble_property_.RNSamples       = FullPrecRealType(0);
  // \\todo If WalkerControl is not exclusively for dmc then this shouldn't be here.
  // If it is it shouldn't be in QMDrivers but QMCDrivers/DMC
  if (dmcStream)
  {
    //boost::archive::text_oarchive oa(*dmcStream);
    //(*oa) & iter  & eavg_cur & wgt_cur & Etrial  & pop_old;
    (*dmcStream) << std::setw(10) << iter << std::setw(20) << ensemble_property_.Energy << std::setw(20)
                 << ensemble_property_.Variance << std::setw(20) << ensemble_property_.Weight << std::setw(20)
                 << ensemble_property_.NumSamples << std::setw(20)
                 << curData[SENTWALKERS_INDEX] / static_cast<double>(num_contexts_);
    (*dmcStream) << std::setw(20) << trialEnergy << std::setw(20)
                 << ensemble_property_.R2Accepted / ensemble_property_.R2Proposed;
    (*dmcStream) << std::setw(20) << ensemble_property_.LivingFraction;
    // Work around for bug with deterministic scalar trace test on select compiler/architectures.
    // While WalkerControl appears to have exclusive ownership of the dmcStream pointer,
    // this is not actually true. Apparently it doesn't actually and can loose ownership then it is
    // either leaked or not flushed before it is destroyed.
    // \todo fix this, you don't want to flush every step since you really hope that could be very rapid.
    (*dmcStream)
        << std::endl; //'\n'; // this is definitely not a place to put an endl as that is also a signal for a flush.
  }
}

void WalkerControl::reset() { std::fill(curData.begin(), curData.end(), 0.0); }

QMCTraits::FullPrecRealType WalkerControl::branch(int iter, MCPopulation& pop)
{
  /* dynamic population
    1. compute multiplicity. If iter 0, multiplicity = 1
    2. compute curData, collect multiplicity on every rank

     fix population
    1. compute curData, collect weight on every rank
    2. compute multiplicity by comb method

    3. figure out final distribution
    4. collect good, bad walkers
    5. communicate walkers
    6. unpack received walkers, apply walker count floor and ceiling.
   */

  auto& walkers = pop.get_walkers();

  if (use_fixed_pop_)
  {
    computeCurData(walkers);
    // convert  node local num of walkers after combing
    // curData[LE_MAX + MyContext] = wsum to num_total_copies
    // calculate walker->Multiplicity;
  }
  else
  {
    // no branching at the first iteration to avoid large population change.
    if (iter == 0)
      for (auto& walker : walkers)
        walker->Multiplicity = 1.0;
    else
      for (auto& walker : walkers)
        walker->Multiplicity = static_cast<int>(walker->Weight + rng_());
    computeCurData(walkers);
    for (int i = 0, j = LE_MAX; i < num_contexts_; i++, j++)
      NumPerNode[i] = static_cast<int>(curData[j]);
  }
  // at this point, curData[LE_MAX + MyContext] and walker->Multiplicity are ready.

  //const int current_population = std::accumulate(NumPerNode.begin(), NumPerNode.end(), 0);
  //std::cout << "debug " << iter << " current_population " << current_population << std::endl;

  // kill walkers, actually put them in deadlist
  RefVector<MCPWalker> bad_walkers;
  bad_walkers.reserve(walkers.size());
  for (auto& walker : walkers)
    if (static_cast<int>(walker->Multiplicity) == 0)
      bad_walkers.push_back(*walker);
  for (MCPWalker& bad_walker : bad_walkers)
    pop.killWalker(bad_walker);
  bad_walkers.clear();

  // copy good walkers
  const size_t good_walkers = walkers.size();
  for (size_t iw = 0; iw < good_walkers; iw++)
  {
    size_t num_copies = static_cast<int>(walkers[iw]->Multiplicity);
    while (num_copies > 1)
    {
      auto walker_elements = pop.spawnWalker();
      pop.checkIntegrity();
      walker_elements.walker = *walkers[iw];
      walker_elements.pset.loadWalker(walker_elements.walker, true);
      walker_elements.pset.update();
      walker_elements.twf.evaluateLog(walker_elements.pset);
      num_copies--;
    }
  }

  measureProperties(iter);
  pop.set_ensemble_property(ensemble_property_);

  pop.syncWalkersPerNode(myComm);

  for (UPtr<MCPWalker>& walker : pop.get_walkers())
  {
    walker->Weight       = 1.0;
    walker->Multiplicity = 1.0;
  }

  // At this point Weight == global_walkers
  return pop.get_num_global_walkers();
}

void WalkerControl::computeCurData(const UPtrVector<MCPWalker>& walkers)
{
  FullPrecRealType esum = 0.0, e2sum = 0.0, wsum = 0.0;
  FullPrecRealType r2_accepted = 0.0, r2_proposed = 0.0;
  int num_good_walkers(0), num_total_copies(0);
  for (const auto& walker : walkers)
  {
    const int num_copies = static_cast<int>(walker->Multiplicity);
    if (num_copies > 0)
    {
      num_good_walkers++;
      num_total_copies += num_copies;
      // Ye : not sure about these r2
      r2_accepted += walker->Properties(WP::R2ACCEPTED);
      r2_proposed += walker->Properties(WP::R2PROPOSED);
      FullPrecRealType e   = walker->Properties(WP::LOCALENERGY);
      FullPrecRealType wgt = walker->Weight;
      esum += wgt * e;
      e2sum += wgt * e * e;
      wsum += wgt;
    }
  }
  //temp is an array to perform reduction operations
  std::fill(curData.begin(), curData.end(), 0);
  curData[ENERGY_INDEX]      = esum;
  curData[ENERGY_SQ_INDEX]   = e2sum;
  curData[WALKERSIZE_INDEX]  = walkers.size(); // num of all the current walkers (good+bad)
  curData[WEIGHT_INDEX]      = wsum;
  curData[R2ACCEPTED_INDEX]  = r2_accepted;
  curData[R2PROPOSED_INDEX]  = r2_proposed;
  curData[FNSIZE_INDEX]      = num_good_walkers; // num of good walkers before branching
  curData[SENTWALKERS_INDEX] = NumWalkersSent;
  if (use_fixed_pop_)
    curData[LE_MAX + MyContext] = wsum; // node sum of walker weights
  else
    curData[LE_MAX + MyContext] = num_total_copies; // node num of walkers after local branching

  myComm->allreduce(curData);
}

void WalkerControl::Write2XYZ(MCWalkerConfiguration& W)
{
  std::ofstream fout("bad.xyz");
  MCWalkerConfiguration::iterator it(W.begin());
  MCWalkerConfiguration::iterator it_end(W.end());
  int nptcls(W.getTotalNum());
  while (it != it_end)
  {
    fout << nptcls << std::endl
         << "# E = " << (*it)->Properties(WP::LOCALENERGY) << " Wgt= " << (*it)->Weight << std::endl;
    for (int i = 0; i < nptcls; i++)
      fout << "H " << (*it)->R[i] << std::endl;
    ++it;
  }
}

/** legacy population limiting
 */
int WalkerControl::applyNmaxNmin(int current_population)
{
  // limit Nmax
  int current_max = (current_population + num_contexts_ - 1) / num_contexts_;
  if (current_max > n_max_)
  {
    app_warning() << "Exceeding Max Walkers per MPI rank : " << n_max_ << ". Ceiling is applied" << std::endl;
    int nsub = current_population - n_max_ * num_contexts_;
    for (int inode = 0; inode < num_contexts_; inode++)
      if (NumPerNode[inode] > n_max_)
      {
        int n_remove = std::min(nsub, NumPerNode[inode] - n_max_);
        NumPerNode[inode] -= n_remove;
        nsub -= n_remove;

        if (inode == MyContext)
        {
          for (int iw = 0; iw < ncopy_w.size(); iw++)
          {
            int n_remove_walker = std::min(ncopy_w[iw], n_remove);
            ncopy_w[iw] -= n_remove_walker;
            n_remove -= n_remove_walker;
            if (n_remove == 0)
              break;
          }

          if (n_remove > 0)
          {
            app_warning() << "Removing copies of good walkers is not enough. "
                          << "Removing good walkers." << std::endl;
            do
            {
              bad_w.push_back(good_w.back());
              good_w.pop_back();
              ncopy_w.pop_back();
              --n_remove;
            } while (n_remove > 0 && !good_w.empty());
          }

          if (n_remove)
            APP_ABORT("WalkerControl::applyNmaxNmin not able to remove sufficient walkers on a node!");
          if (std::accumulate(ncopy_w.begin(), ncopy_w.end(), ncopy_w.size()) != NumPerNode[inode])
            APP_ABORT("WalkerControl::applyNmaxNmin walker removal mismatch!");
        }

        if (nsub == 0)
          break;
      }

    if (nsub)
      APP_ABORT("WalkerControl::applyNmaxNmin not able to remove sufficient walkers overall!");
  }

  // limit Nmin
  if (current_population / num_contexts_ < n_min_)
  {
    app_warning() << "The number of walkers is running lower than Min Walkers per MPI rank : " << n_min_
                  << ". Floor is applied" << std::endl;
    int nadd = n_min_ * num_contexts_ - current_population;
    for (int inode = 0; inode < num_contexts_; inode++)
      if (NumPerNode[inode] > 0 && NumPerNode[inode] < n_min_)
      {
        int n_insert = std::min(nadd, n_min_ - NumPerNode[inode]);
        NumPerNode[inode] += n_insert;
        nadd -= n_insert;

        if (inode == MyContext)
        {
          int n_avg_insert_per_walker = (n_insert + ncopy_w.size() - 1) / ncopy_w.size();
          for (int iw = 0; iw < ncopy_w.size(); iw++)
          {
            int n_insert_walker = std::min(n_avg_insert_per_walker, n_insert);
            ncopy_w[iw] += n_insert_walker;
            n_insert -= n_insert_walker;
            if (n_insert == 0)
              break;
          }

          if (std::accumulate(ncopy_w.begin(), ncopy_w.end(), ncopy_w.size()) != NumPerNode[inode])
            APP_ABORT("WalkerControl::applyNmaxNmin walker insertion mismatch!");
        }

        if (nadd == 0)
          break;
      }

    if (nadd)
      app_warning() << "WalkerControl::applyNmaxNmin not able to add sufficient walkers overall!" << std::endl;
  }

  // check current population
  current_population = std::accumulate(NumPerNode.begin(), NumPerNode.end(), 0);
  // at least one walker after load-balancing
  if (current_population / num_contexts_ == 0)
  {
    app_error() << "Some MPI ranks have no walkers after load balancing. This should not happen."
                << "Improve the trial wavefunction or adjust the simulation parameters." << std::endl;
    APP_ABORT("WalkerControl::sortWalkers");
  }

  return current_population;
}

// determine new walker population on each node
void WalkerControl::determineNewWalkerPopulation(const std::vector<int>& num_per_node,
                                                 std::vector<int>& fair_offset,
                                                 std::vector<int>& minus,
                                                 std::vector<int>& plus)
{
  const int num_contexts       = num_per_node.size();
  const int current_population = std::accumulate(num_per_node.begin(), num_per_node.end(), 0);
  FairDivideLow(current_population, num_contexts, fair_offset);
  for (int ip = 0; ip < num_contexts; ip++)
  {
    // (FairOffSet[ip + 1] - FairOffSet[ip]) gives the partiion ip walker pop
    int dn = num_per_node[ip] - (fair_offset[ip + 1] - fair_offset[ip]);
    if (dn > 0)
      plus.insert(plus.end(), dn, ip);
    else if (dn < 0)
      minus.insert(minus.end(), -dn, ip);
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
 * The algorithm ensures that the load per node can differ only by one walker.
 * Each MPI rank can only send or receive or be silent.
 * The communication is one-dimensional and very local.
 * If multiple copies of a walker need to be sent to the target rank, only send one.
 * The number of copies is communicated ahead via blocking send/recv.
 * Then the walkers are transferred via blocking or non-blocking send/recv.
 * The blocking send/recv may become serialized and worsen load imbalance.
 * Non blocking send/recv algorithm avoids serialization completely.
 */
void WalkerControl::swapWalkersSimple(MCPopulation& pop)
{
  std::vector<int> minus, plus;
  determineNewWalkerPopulation(NumPerNode, FairOffSet, minus, plus);

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
        myTimers[WC_send]->start();
        myComm->comm.send_n(awalker->DataSet.data(), byteSize, jobit->target);
        myTimers[WC_send]->stop();
      }
    }
    if (use_nonblocking)
    {
      // wait all the isend
      for (int im = 0; im < requests.size(); im++)
      {
        myTimers[WC_send]->start();
        requests[im].wait();
        myTimers[WC_send]->stop();
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
        myTimers[WC_recv]->start();
        myComm->comm.receive_n(awalker->DataSet.data(), byteSize, jobit->target);
        awalker->copyFromBuffer();
        myTimers[WC_recv]->stop();
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

  assert(std::accumulate(ncopy_w.begin(), ncopy_w.end(), ncopy_w.size()) == num_per_node[MyContext]);
}

std::vector<WalkerControl::IndexType> WalkerControl::syncFutureWalkersPerRank(Communicate* comm, IndexType n_walkers)
{
  int ncontexts = comm->size();
  std::vector<IndexType> future_walkers(ncontexts, 0);
  future_walkers[comm->rank()] = n_walkers;
  comm->allreduce(future_walkers);
  return future_walkers;
}

bool WalkerControl::put(xmlNodePtr cur)
{
  int nw_target = 0, nw_max = 0;
  std::string nonblocking = "yes";
  ParameterSet params;
  params.add(target_sigma_, "sigmaBound", "double");
  params.add(MaxCopy, "maxCopy", "int");
  params.add(nw_target, "targetwalkers", "int");
  params.add(nw_max, "max_walkers", "int");
  params.add(nonblocking, "use_nonblocking", "string");

  bool success = params.put(cur);

  // validating input
  if (nonblocking == "yes")
  {
    use_nonblocking = true;
  }
  else if (nonblocking == "no")
  {
    use_nonblocking = false;
  }
  else
  {
    APP_ABORT("WalkerControl::put unknown use_nonblocking option " + nonblocking);
  }

  setMinMax(nw_target, nw_max);

  app_log() << "  WalkerControl parameters " << std::endl;
  //app_log() << "    energyBound = " << targetEnergyBound << std::endl;
  //app_log() << "    sigmaBound = " << targetSigma << std::endl;
  app_log() << "    maxCopy = " << MaxCopy << std::endl;
  app_log() << "    Max Walkers per MPI rank " << n_max_ << std::endl;
  app_log() << "    Min Walkers per MPI rank " << n_min_ << std::endl;
  app_log() << "    Using " << (use_nonblocking ? "non-" : "") << "blocking send/recv" << std::endl;
  return true;
}

void WalkerControl::setMinMax(int nw_in, int nmax_in)
{
  if (nw_in > 0)
  {
    int npernode = nw_in / num_contexts_;
    if (use_fixed_pop_)
    {
      n_max_ = npernode;
      n_min_ = npernode;
    }
    else
    {
      n_max_ = MaxCopy * npernode + 1;
      n_min_ = npernode / 5 + 1;
      if (nmax_in > 0)
        n_max_ = nmax_in;
    }
  }
}

} // namespace qmcplusplus
