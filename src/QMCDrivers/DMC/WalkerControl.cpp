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


#include <array>
#include <cassert>
#include <stdexcept>
#include <numeric>
#include <sstream>

#include "WalkerControl.h"
#include "QMCDrivers/WalkerProperties.h"
#include "OhmmsData/ParameterSet.h"
#include "type_traits/template_types.hpp"
#include "QMCWaveFunctions/TrialWaveFunction.h"

namespace qmcplusplus
{
using WP = WalkerProperties::Indexes;

enum WC_Timers
{
  WC_branch,
  WC_imbalance,
  WC_prebalance,
  WC_copyWalkers,
  WC_recomputing,
  WC_allreduce,
  WC_loadbalance,
  WC_send,
  WC_recv,
};

TimerNameList_t<WC_Timers> WalkerControlTimerNames = {{WC_branch, "WalkerControl::branch"},
                                                      {WC_imbalance, "WalkerControl::imbalance"},
                                                      {WC_prebalance, "WalkerControl::pre-loadbalance"},
                                                      {WC_copyWalkers, "WalkerControl::copyWalkers"},
                                                      {WC_recomputing, "WalkerControl::recomputing"},
                                                      {WC_allreduce, "WalkerControl::allreduce"},
                                                      {WC_loadbalance, "WalkerControl::loadbalance"},
                                                      {WC_send, "WalkerControl::send"},
                                                      {WC_recv, "WalkerControl::recv"}};

WalkerControl::WalkerControl(Communicate* c, RandomBase<FullPrecRealType>& rng, bool use_fixed_pop)
    : MPIObjectBase(c),
      rng_(rng),
      use_fixed_pop_(use_fixed_pop),
      n_min_(1),
      n_max_(10),
      max_copy_(2),
      rank_num_(c->rank()),
      num_ranks_(c->size()),
      SwapMode(0),
      use_nonblocking_(true),
      debug_disable_branching_(false),
      my_timers_(getGlobalTimerManager(), WalkerControlTimerNames, timer_level_medium),
      saved_num_walkers_sent_(0)
{
  num_per_rank_.resize(num_ranks_);
  fair_offset_.resize(num_ranks_ + 1);
}

WalkerControl::~WalkerControl() = default;

void WalkerControl::start()
{
  if (rank_num_ == 0)
  {
    std::filesystem::path hname(myComm->getName());
    hname.concat(".dmc.dat");
    if (hname != dmcFname)
    {
      dmcStream = std::make_unique<std::ofstream>(hname);
      dmcStream->setf(std::ios::scientific, std::ios::floatfield);
      dmcStream->precision(10);
      (*dmcStream) << "# Index " << std::setw(20) << "LocalEnergy" << std::setw(20) << "Variance" << std::setw(20)
                   << "Weight" << std::setw(20) << "NumOfWalkers" << std::setw(20)
                   << "AvgSentWalkers"; //add the number of walkers
      (*dmcStream) << std::setw(20) << "TrialEnergy" << std::setw(20) << "DiffEff";
      (*dmcStream) << std::setw(20) << "LivingFraction";
      (*dmcStream) << std::endl;
      dmcFname = std::move(hname);
    }
  }
}

void WalkerControl::writeDMCdat(int iter, const std::vector<FullPrecRealType>& curData)
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
  // \\todo If WalkerControl is not exclusively for dmc then this shouldn't be here.
  // If it is it shouldn't be in QMDrivers but QMCDrivers/DMC
  if (dmcStream)
  {
    //boost::archive::text_oarchive oa(*dmcStream);
    //(*oa) & iter  & eavg_cur & wgt_cur & Etrial  & pop_old;
    (*dmcStream) << std::setw(10) << iter << std::setw(20) << ensemble_property_.Energy << std::setw(20)
                 << ensemble_property_.Variance << std::setw(20) << ensemble_property_.Weight << std::setw(20)
                 << ensemble_property_.NumSamples << std::setw(20)
                 << curData[SENTWALKERS_INDEX] / static_cast<double>(num_ranks_);
    (*dmcStream) << std::setw(20) << trial_energy_ << std::setw(20)
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

void WalkerControl::branch(int iter, MCPopulation& pop, bool do_not_branch)
{
  if (debug_disable_branching_)
    do_not_branch = true;
  /* dynamic population
    1. compute multiplicity. If iter 0, multiplicity = 1
    2. compute curData, collect multiplicity on every rank

     fix population
    1. compute curData, collect weight on every rank
    2. compute multiplicity by comb method

    3. figure out final distribution, apply walker count ceiling
    4. collect good, bad walkers
    5. communicate walkers
    6. unpack received walkers, apply walker count floor
   */

  ScopedTimer branch_timer(my_timers_[WC_branch]);
  auto& walkers = pop.get_walkers();

  {
    ScopedTimer prebalance_timer(my_timers_[WC_prebalance]);
    ///any temporary data includes many ridiculous conversions of integral types to and from fp
    std::vector<FullPrecRealType> curData(LE_MAX + num_ranks_, 0.0);

    if (use_fixed_pop_)
    {
      computeCurData(walkers, curData);
      // convert  node local num of walkers after combing
      // curData[LE_MAX + rank_num_] = wsum to num_total_copies
      // calculate walker->Multiplicity;
    }
    else
    {
      // no branching at the first iteration to avoid large population change.
      if (do_not_branch)
        for (auto& walker : walkers)
          walker->Multiplicity = 1.0;
      else
        for (auto& walker : walkers)
          walker->Multiplicity = static_cast<int>(walker->Weight + rng_());
      computeCurData(walkers, curData);
      for (int i = 0, j = LE_MAX; i < num_ranks_; i++, j++)
        num_per_rank_[i] = static_cast<int>(curData[j]);
    }
    // at this point, curData[LE_MAX + rank_num_] and walker->Multiplicity are ready.

    writeDMCdat(iter, curData);
    pop.set_ensemble_property(ensemble_property_);
  }

  auto untouched_walkers = walkers.size();
#if defined(HAVE_MPI)
  {
    ScopedTimer loadbalance_timer(my_timers_[WC_loadbalance]);
    // kill walkers, actually put them in deadlist for be recycled for receiving walkers
    killDeadWalkersOnRank(pop);
    // ranks receiving walkers from other ranks have the lowest walker count now.
    untouched_walkers = std::min(untouched_walkers, walkers.size());

    // load balancing over MPI
    swapWalkersSimple(pop);
  }
#endif

  // kill dead walker to be recycled by the following copy
  killDeadWalkersOnRank(pop);
  // ranks sending walkers from other ranks have the lowest walker count now.
  untouched_walkers = std::min(untouched_walkers, walkers.size());

  { // copy good walkers
    ScopedTimer copywalkers_timer(my_timers_[WC_copyWalkers]);
    const size_t good_walkers = walkers.size();
    for (size_t iw = 0; iw < good_walkers; iw++)
    {
      size_t num_copies = static_cast<int>(walkers[iw]->Multiplicity);
      while (num_copies > 1)
      {
        auto walker_elements = pop.spawnWalker();
        // save this walkers ID
        // \todo revisit Walker assignment operator after legacy drivers removed.
        // but in the modern scheme walker IDs are permanent after creation, what walker they
        // were copied from is in ParentID.
        long save_id                    = walker_elements.walker.ID;
        walker_elements.walker          = *walkers[iw];
        walker_elements.walker.ParentID = walker_elements.walker.ID;
        walker_elements.walker.ID       = save_id;
        num_copies--;
      }
    }
  }

  const int current_num_global_walkers = std::accumulate(num_per_rank_.begin(), num_per_rank_.end(), 0);
  pop.set_num_global_walkers(current_num_global_walkers);
#ifndef NDEBUG
  pop.checkIntegrity();
  pop.syncWalkersPerRank(myComm);
  if (current_num_global_walkers != pop.get_num_global_walkers())
    throw std::runtime_error("Potential bug! Population num_global_walkers mismatched!");
#endif

  if (!do_not_branch)
    for (UPtr<MCPWalker>& walker : pop.get_walkers())
    {
      walker->Weight       = 1.0;
      walker->Multiplicity = 1.0;
    }

  for (int iw = 0; iw < untouched_walkers; iw++)
    pop.get_walkers()[iw]->wasTouched = false;

  for (int iw = untouched_walkers; iw < pop.get_num_local_walkers(); iw++)
    pop.get_walkers()[iw]->wasTouched = true;
}

void WalkerControl::computeCurData(const UPtrVector<MCPWalker>& walkers, std::vector<FullPrecRealType>& curData)
{
  FullPrecRealType esum = 0.0, e2sum = 0.0, wsum = 0.0;
  FullPrecRealType r2_accepted = 0.0, r2_proposed = 0.0;
  int num_good_walkers(0), num_total_copies(0);
  for (const auto& walker : walkers)
  {
    const int num_copies = static_cast<int>(walker->Multiplicity);
    num_good_walkers += num_copies > 0 ? 1 : 0;
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
  //temp is an array to perform reduction operations
  std::fill(curData.begin(), curData.end(), 0.0);
  curData[ENERGY_INDEX]      = esum;
  curData[ENERGY_SQ_INDEX]   = e2sum;
  curData[WALKERSIZE_INDEX]  = walkers.size(); // num of all the current walkers (good+bad)
  curData[WEIGHT_INDEX]      = wsum;
  curData[R2ACCEPTED_INDEX]  = r2_accepted;
  curData[R2PROPOSED_INDEX]  = r2_proposed;
  curData[FNSIZE_INDEX]      = num_good_walkers; // num of good walkers before branching
  curData[SENTWALKERS_INDEX] = saved_num_walkers_sent_;
  if (use_fixed_pop_)
    curData[LE_MAX + rank_num_] = wsum; // node sum of walker weights
  else
    curData[LE_MAX + rank_num_] = num_total_copies; // node num of walkers after local branching

  {
    ScopedTimer allreduce_timer(my_timers_[WC_allreduce]);
    myComm->allreduce(curData);
  }
}

// determine new walker population on each node
void WalkerControl::determineNewWalkerPopulation(const std::vector<int>& num_per_rank,
                                                 std::vector<int>& fair_offset,
                                                 std::vector<int>& minus,
                                                 std::vector<int>& plus)
{
  const int num_contexts       = num_per_rank.size();
  const int current_population = std::accumulate(num_per_rank.begin(), num_per_rank.end(), 0);
  FairDivideLow(current_population, num_contexts, fair_offset);
  for (int ip = 0; ip < num_contexts; ip++)
  {
    int dn = num_per_rank[ip] - (fair_offset[ip + 1] - fair_offset[ip]);
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
    throw std::runtime_error("Trying to swap in WalkerControl::swapWalkersSimple with mismatched queues");
  }
#endif
}

#if defined(HAVE_MPI)
void WalkerControl::swapWalkersSimple(MCPopulation& pop)
{
  std::vector<int> minus, plus;
  determineNewWalkerPopulation(num_per_rank_, fair_offset_, minus, plus);

#ifdef MCWALKERSET_MPI_DEBUG
  std::array<char, 128> fname;
  if (std::snprintf(fname.data(), fname.size(), "test.%d", rank_num_) < 0)
    throw std::runtime_error("Error generating filename");
  std::ofstream fout(fname.data(), std::ios::app);

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

  auto& good_walkers = pop.get_walkers();
  const int nswap    = plus.size();
  // sort good walkers by the number of copies
  std::vector<std::pair<int, int>> ncopy_pairs;
  for (int iw = 0; iw < good_walkers.size(); iw++)
    ncopy_pairs.push_back(std::make_pair(static_cast<int>(good_walkers[iw]->Multiplicity), iw));
  std::sort(ncopy_pairs.begin(), ncopy_pairs.end());

  struct job
  {
    const int walkerID;
    const int target;
    job(int wid, int target_in) : walkerID(wid), target(target_in){};
  };

  int nsend = 0;
  std::vector<job> job_list;
  std::vector<WalkerElementsRef> newW;
  std::vector<int> ncopy_newW;

  for (int ic = 0; ic < nswap; ic++)
  {
    int nsentcopy = 0;
    if (plus[ic] == rank_num_)
    {
      // always send the last good walker with most copies
      // count the possible copies in one send
      for (int id = ic + 1; id < nswap; id++)
        if (plus[ic] == plus[id] && minus[ic] == minus[id] && ncopy_pairs.back().first > 1)
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

      // update copy counter
      if (ncopy_pairs.back().first > 1)
      {
        ncopy_pairs.back().first--;
        std::sort(ncopy_pairs.begin(), ncopy_pairs.end());
      }
      else
      {
        good_walkers[ncopy_pairs.back().second]->Multiplicity = 0.0;
        ncopy_pairs.pop_back();
      }
    }

    if (minus[ic] == rank_num_)
    {
      newW.push_back(pop.spawnWalker());

      // recv the number of copies from the target
      myComm->comm.receive_n(&nsentcopy, 1, plus[ic]);
      job_list.push_back(job(newW.size() - 1, plus[ic]));
      if (plus[ic] != plus[ic + nsentcopy] || minus[ic] != minus[ic + nsentcopy])
        throw std::runtime_error("WalkerControl::swapWalkersSimple send/recv pair checking failed!");
#ifdef MCWALKERSET_MPI_DEBUG
      fout << "rank " << minus[ic] << " recvs a walker with " << nsentcopy << " copies from rank " << plus[ic]
           << std::endl;
#endif

      ncopy_newW.push_back(nsentcopy);
    }

    // update cursor
    ic += nsentcopy;
  }

  if (nsend > 0)
  {
    std::vector<mpi3::request> requests;
    // mark all walkers not in send
    for (auto jobit = job_list.begin(); jobit != job_list.end(); jobit++)
      good_walkers[jobit->walkerID]->SendInProgress = false;
    for (auto jobit = job_list.begin(); jobit != job_list.end(); jobit++)
    {
      // pack data and send
      auto& awalker   = good_walkers[jobit->walkerID];
      size_t byteSize = awalker->byteSize();
      if (!awalker->SendInProgress)
      {
        awalker->updateBuffer();
        awalker->SendInProgress = true;
      }
      if (use_nonblocking_)
        requests.push_back(myComm->comm.isend_n(awalker->DataSet.data(), byteSize, jobit->target));
      else
      {
        ScopedTimer local_timer(my_timers_[WC_send]);
        myComm->comm.send_n(awalker->DataSet.data(), byteSize, jobit->target);
      }
    }
    if (use_nonblocking_)
    {
      // wait all the isend
      for (int im = 0; im < requests.size(); im++)
      {
        ScopedTimer local_timer(my_timers_[WC_send]);
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
      auto& walker_elements = newW[jobit->walkerID];
      auto& awalker         = walker_elements.walker;
      size_t byteSize       = awalker.byteSize();
      if (use_nonblocking_)
        requests.push_back(myComm->comm.ireceive_n(awalker.DataSet.data(), byteSize, jobit->target));
      else
      {
        ScopedTimer local_timer(my_timers_[WC_recv]);
        myComm->comm.receive_n(awalker.DataSet.data(), byteSize, jobit->target);
        awalker.copyFromBuffer();
      }
    }
    if (use_nonblocking_)
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
              auto& walker_elements = newW[job_list[im].walkerID];
              walker_elements.walker.copyFromBuffer();
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
  saved_num_walkers_sent_ = nsend;

  // rebuild Multiplicity
  for (int iw = 0; iw < ncopy_pairs.size(); iw++)
    good_walkers[ncopy_pairs[iw].second]->Multiplicity = ncopy_pairs[iw].first;

  for (int iw = 0; iw < newW.size(); iw++)
    newW[iw].walker.Multiplicity = ncopy_newW[iw] + 1;

#ifndef NDEBUG
  FullPrecRealType TotalMultiplicity = 0;
  for (int iw = 0; iw < good_walkers.size(); iw++)
    TotalMultiplicity += good_walkers[iw]->Multiplicity;
  if (static_cast<int>(TotalMultiplicity) != fair_offset_[rank_num_ + 1] - fair_offset_[rank_num_])
    throw std::runtime_error("Multiplicity check failed in WalkerControl::swapWalkersSimple!");
#endif
}
#endif

void WalkerControl::killDeadWalkersOnRank(MCPopulation& pop)
{
  // kill walkers, actually put them in deadlist
  RefVector<MCPWalker> bad_walkers;
  auto& walkers = pop.get_walkers();
  bad_walkers.reserve(walkers.size());
  for (auto& walker : walkers)
    if (static_cast<int>(walker->Multiplicity) == 0)
      bad_walkers.push_back(*walker);
  for (MCPWalker& bad_walker : bad_walkers)
    pop.killWalker(bad_walker);
#ifndef NDEBUG
  pop.checkIntegrity();
#endif
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
  ParameterSet params;
  params.add(max_copy_, "maxCopy");
  params.add(nw_target, "targetwalkers");
  params.add(nw_max, "max_walkers");
  params.add(use_nonblocking_, "use_nonblocking", {true});
  params.add(debug_disable_branching_, "debug_disable_branching", {false});

  try
  {
    bool success = params.put(cur);
  }
  catch (const std::runtime_error& re)
  {
    myComm->barrier_and_abort("WalkerControl::put parsing error. " + std::string(re.what()));
  }

  setMinMax(nw_target, nw_max);

  app_log() << "  WalkerControl parameters " << std::endl;
  //app_log() << "    energyBound = " << targetEnergyBound << std::endl;
  //app_log() << "    sigmaBound = " << targetSigma << std::endl;
  app_log() << "    maxCopy = " << max_copy_ << std::endl;
  app_log() << "    Max Walkers per MPI rank " << n_max_ << std::endl;
  app_log() << "    Min Walkers per MPI rank " << n_min_ << std::endl;
  app_log() << "    Using " << (use_nonblocking_ ? "non-" : "") << "blocking send/recv" << std::endl;
  if (debug_disable_branching_)
    app_log() << "    Disable branching for debugging as the user input request." << std::endl;
  return true;
}

void WalkerControl::setMinMax(int nw_in, int nmax_in)
{
  if (nw_in > 0)
  {
    int npernode = nw_in / num_ranks_;
    if (use_fixed_pop_)
    {
      n_max_ = npernode;
      n_min_ = npernode;
    }
    else
    {
      n_max_ = max_copy_ * npernode + 1;
      n_min_ = npernode / 5 + 1;
      if (nmax_in > 0)
        n_max_ = nmax_in;
    }
  }
}

} // namespace qmcplusplus
