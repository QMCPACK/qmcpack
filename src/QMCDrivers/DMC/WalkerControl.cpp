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
#include <cmath>
#include <stdexcept>
#include <numeric>
#include <sstream>

#include "WalkerControl.h"
#include "QMCDrivers/WalkerProperties.h"
#include "OhmmsData/ParameterSet.h"
#include "type_traits/template_types.hpp"
#include "QMCWaveFunctions/TrialWaveFunction.h"
#include "Message/UniformCommunicateError.h"
#include "CPU/math.hpp"

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

WalkerControl::WalkerControl(Communicate& c, RandomBase<FullPrecRealType>& rng, bool use_fixed_pop)
    : MPIObjectBase(c),
      rng_(rng),
      use_fixed_pop_(use_fixed_pop),
      rank_num_(c.rank()),
      num_ranks_(c.size()),
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

void WalkerControl::start(const std::string& root_name)
{
  if (rank_num_ == 0)
  {
    std::filesystem::path hname(root_name);
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
  FullPrecRealType wgtInv(use_fixed_pop_ ? 0.0 : 1.0 / curData[WEIGHT_INDEX]);
  FullPrecRealType eavg =
      use_fixed_pop_ ? curData[ENERGY_INDEX] / curData[WEIGHT_INDEX] : curData[ENERGY_INDEX] * wgtInv;
  const FullPrecRealType second_moment =
      use_fixed_pop_ ? curData[ENERGY_SQ_INDEX] / curData[WEIGHT_INDEX] : curData[ENERGY_SQ_INDEX] * wgtInv;
  ensemble_property_.Energy     = eavg;
  ensemble_property_.Weight     = curData[WEIGHT_INDEX];
  ensemble_property_.Variance   = second_moment - eavg * eavg;
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

/** unified: perform branch and swap (balance) walkers as required
   *  **This has many side effects**
   *  ## For:
   *  ### dynamic population
   *  1. compute multiplicity. If iter 0 and in warmup -> multiplicity = 1
   *     Multiplicity in normal branching is walker->Weight + rng()
   *  2. compute curData, collect multiplicity on every rank
   *  ### fixed population
   *  1. compute curData, collect weight on every rank
   *  2. compute multiplicity by comb method
   *  ---
   *  3. figure out final distribution, apply walker count ceiling
   *  4. collect good, bad walkers
   *  5. communicate walkers
   *  6. unpack received walkers, apply walker count floor
   *  7. call MCPopulation to amplify walkers with Multiplicity > 1
   */
void WalkerControl::branch(int iter, MCPopulation& pop, bool do_not_branch)
{
  if (debug_disable_branching_)
    do_not_branch = true;

  ScopedTimer branch_timer(my_timers_[WC_branch]);
  auto& walkers = pop.get_walkers();
  FullPrecRealType offspring_weight = 1.0;

  {
    ScopedTimer prebalance_timer(my_timers_[WC_prebalance]);
    ///any temporary data includes many ridiculous conversions of integral types to and from fp
    // Fixed population mode packs its common offset and invalid-input flag into the single reduction.
    std::vector<FullPrecRealType> curData(LE_MAX + num_ranks_ + (use_fixed_pop_ ? 2 : 0), 0.0);

    if (use_fixed_pop_)
    {
      computeCurData(walkers, curData, do_not_branch);
      const int population = pop.get_num_global_walkers();
      if (population <= 0 || curData[WALKERSIZE_INDEX] != population)
        throw UniformCommunicateError("WalkerControl: unexpected inconsistent fixed population. BUG.");
      if (curData[WEIGHT_INDEX] / population == 0)
        throw UniformCommunicateError("WalkerControl: fixed population descendant weight underflow.");
      if (!do_not_branch)
      {
        const std::vector<FullPrecRealType> rank_weights(curData.begin() + LE_MAX,
                                                         curData.begin() + LE_MAX + num_ranks_);
        assignStochasticComb(walkers, rank_weights, rank_num_, population, curData[LE_MAX + num_ranks_], num_per_rank_);
        offspring_weight = curData[WEIGHT_INDEX] / population;
      }
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
      computeCurData(walkers, curData, do_not_branch);
      for (int i = 0, j = LE_MAX; i < num_ranks_; i++, j++)
        num_per_rank_[i] = static_cast<int>(curData[j]);

      // A rank may have no survivors; only global extinction is fatal. The reduced counts
      // make this error uniform across ranks, before diagnostics or walker deletion.
      if (std::accumulate(num_per_rank_.begin(), num_per_rank_.end(), 0) == 0)
        throw UniformCommunicateError("WalkerControl::branch: All walkers were killed in DMC at iteration " +
                                      std::to_string(iter) + ".");
    }
    // at this point, curData[LE_MAX + rank_num_] and walker->Multiplicity are ready.

    writeDMCdat(iter, curData);
    pop.set_ensemble_property(ensemble_property_);
  }

  // Do not migrate walkers with fixed population and a skipped comb.
  if (use_fixed_pop_ && do_not_branch)
    return;

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

  {
    ScopedTimer copywalkers_timer(my_timers_[WC_copyWalkers]);
    pop.fissionHighMultiplicityWalkers();
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
      walker->Weight       = offspring_weight;
      walker->Multiplicity = 1.0;
    }

  for (int iw = 0; iw < untouched_walkers; iw++)
    pop.get_walkers()[iw]->wasTouched = false;

  for (int iw = untouched_walkers; iw < pop.get_num_local_walkers(); iw++)
    pop.get_walkers()[iw]->wasTouched = true;
}

void WalkerControl::computeCurData(const UPtrVector<MCPWalker>& walkers,
                                   std::vector<FullPrecRealType>& curData,
                                   bool do_not_branch)
{
  FullPrecRealType esum = 0.0, e2sum = 0.0, wsum = 0.0;
  FullPrecRealType r2_accepted = 0.0, r2_proposed = 0.0;
  int num_good_walkers(0), num_total_copies(0);
  bool invalid_weights = false;
  for (const auto& walker : walkers)
  {
    const int num_copies = use_fixed_pop_ ? 0 : static_cast<int>(walker->Multiplicity);
    // In fixed mode LivingFraction describes the positive-weight input
    num_good_walkers += (use_fixed_pop_ ? walker->Weight > 0 : num_copies > 0) ? 1 : 0;
    num_total_copies += num_copies;
    // Ye : not sure about these r2
    r2_accepted += walker->Properties(WP::R2ACCEPTED);
    r2_proposed += walker->Properties(WP::R2PROPOSED);
    FullPrecRealType e   = walker->Properties(WP::LOCALENERGY);
    FullPrecRealType wgt = walker->Weight;
    if (use_fixed_pop_ && (!qmcplusplus::isfinite(wgt) || wgt < 0))
    {
      invalid_weights = true;
      continue;
    }
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
  {
    curData[LE_MAX + rank_num_] = wsum; // node sum of walker weights
    curData[LE_MAX + num_ranks_]     = rank_num_ == 0 && !do_not_branch ? rng_() : 0.0; // random offset for stochastic comb computed only on rank 0
    curData[LE_MAX + num_ranks_ + 1] = invalid_weights ? 1.0 : 0.0;
  }
  else
    curData[LE_MAX + rank_num_] = num_total_copies; // node num of walkers after local branching

  {
    ScopedTimer allreduce_timer(my_timers_[WC_allreduce]);
    myComm.allreduce(curData);
  }

  if (use_fixed_pop_)
  {
    curData[WEIGHT_INDEX] =
        std::accumulate(curData.begin() + LE_MAX, curData.begin() + LE_MAX + num_ranks_, FullPrecRealType(0));
    if (curData[LE_MAX + num_ranks_ + 1] != 0 || !qmcplusplus::isfinite(curData[WEIGHT_INDEX]) ||
        curData[WEIGHT_INDEX] <= 0 || !qmcplusplus::isfinite(curData[ENERGY_INDEX]) ||
        !qmcplusplus::isfinite(curData[ENERGY_SQ_INDEX]))
      throw UniformCommunicateError("WalkerControl: fixed population found erroneous weights.");
  }
}

void WalkerControl::assignStochasticComb(const UPtrVector<MCPWalker>& walkers,
                                         const std::vector<FullPrecRealType>& rank_weights,
                                         int rank,
                                         int population,
                                         FullPrecRealType offset,
                                         std::vector<int>& num_per_rank)
{
  // Stochastic comb implementation for assigning walker multiplicities based on rank weights and a shared offset.
  // No communication is performed. All ranks independently compute the entire comb using the shared rank_weights and offset.
  // For N=population and W=total, comb teeth are x_k=(k+offset)*W/N, k=0,...,N-1. Each weight w walker owns the interval [C, C+w).
  // branch() later balances/copies the selected walkers and given every descendant the appropriate weight W/N.

  // 1. Compute W from the shared rank_weights, using the same rank order on every MPI rank.
  // Validate N, the normalized offset, and W before evaluating the comb spacing or tooth indices.
  const FullPrecRealType total = std::accumulate(rank_weights.begin(), rank_weights.end(), FullPrecRealType(0));
  if (population <= 0 || !qmcplusplus::isfinite(total) || total <= 0 || total / population == 0 ||
      !qmcplusplus::isfinite(offset) || offset < 0 || offset >= 1 ||
      std::any_of(
          rank_weights.begin(), rank_weights.end(),
          [](FullPrecRealType weight) { return !qmcplusplus::isfinite(weight) || weight < 0; }))
    throw UniformCommunicateError("WalkerControl: invalid stochastic comb input.");

  // 2. Express the comb in power-of-two-scaled weight coordinates. frexp writes exponent=e such that
  // W=mantissa*2^e; spacing is (W*2^-e)/N, not the unscaled W/N.
  // Scaling both C and W equally leaves C/(W/N) unchanged; binary exponent scaling is exact unless it underflows.
  // Keeping spacing near 1/N avoids reciprocal overflow/underflow if the compiler rewrites x/y as x*(1/y).
  int exponent;
  const FullPrecRealType spacing = std::frexp(total, &exponent) / population;

  // toothBoundary(C) counts teeth strictly below the cumulative global weight C.
  // Subtracting boundary counts at a walker's endpoints implements the half-open interval convention,
  // so a tooth exactly at an endpoint belongs to the interval *beginning* there, not ending there.
  const auto toothBoundary = [spacing, exponent, population, offset](FullPrecRealType cumulative) {
    // scaled is C/(W/N), the endpoint in units of comb spacing.
    const FullPrecRealType scaled = std::ldexp(cumulative, -exponent) / spacing;
    if (scaled == 0 && cumulative > 0 && offset == 0)
      return 1; // The tooth at zero belongs to the first positive interval even if C/delta underflows.
    // whole is the integral part of scaled. The fractional part determines whether one more tooth
    // precedes C. This is ceil(scaled-offset), clamped to [0,N], without subtracting offset directly:
    // direct subtraction can round down to an integer when offset is immediately below one.
    const FullPrecRealType whole = std::floor(scaled);
    return static_cast<int>(
        std::clamp(whole + (scaled - whole > offset ? 1 : 0), FullPrecRealType(0), FullPrecRealType(population)));
  };

  // 3. boundaries[r] is the global tooth index at the beginning of source rank r.
  // prefix accumulates the global weight through each rank ir; local_start records the weight
  // preceding this rank's walkers, needed to convert local cumulative weights to global endpoints.
  std::vector<int> boundaries(rank_weights.size() + 1, 0);
  FullPrecRealType prefix      = 0;
  FullPrecRealType local_start = 0;
  for (int ir = 0; ir < rank_weights.size(); ++ir)
  {
    if (ir == rank)
      local_start = prefix;
    prefix += rank_weights[ir];
    boundaries[ir + 1] = toothBoundary(prefix);
  }
  // The final endpoint must include all N teeth, independently of floating-point boundary arithmetic.
  boundaries.back() = population;

  // Adjacent rank boundaries determine how many offspring each source rank supplies, with no further
  // MPI communication. These num_per_rank counts are not the balanced destination-rank walker counts.
  num_per_rank.resize(rank_weights.size());
  for (int ir = 0; ir < rank_weights.size(); ++ir)
    num_per_rank[ir] = boundaries[ir + 1] - boundaries[ir];

  // 4. Find the last positive-weight local walker. last_positive stays at walkers.end() for an all-zero
  // rank; otherwise its interval will consume exactly the remaining teeth allocated to this source rank.
  // Pinning this endpoint avoids disagreement between local and shared rank sums, without assigning
  // offspring to any trailing zero-weight walkers.
  auto last_positive = walkers.end();
  for (auto it = walkers.begin(); it != walkers.end(); ++it)
    if ((*it)->Weight > 0)
      last_positive = it;

  // 5. Walk the local intervals in order. cumulative is the local weight through the current walker it;
  // previous is the tooth index at its left endpoint, initially this rank's first tooth index.
  FullPrecRealType cumulative = 0;
  int previous                = boundaries[rank];
  for (auto it = walkers.begin(); it != walkers.end(); ++it)
  {
    // weight is this walker's w; local_start+cumulative becomes its global right endpoint C+w.
    const auto weight = (*it)->Weight;
    cumulative += weight;

    // next is the tooth index at the right endpoint. Zero weight leaves the boundary unchanged.
    // The final positive walker uses the shared rank endpoint; other walkers use the computed boundary,
    // clamped between previous and the rank endpoint to keep multiplicities nonnegative and on this rank.
    const int next = weight == 0 ? previous
        : it == last_positive    ? boundaries[rank + 1]
                                 : std::clamp(toothBoundary(local_start + cumulative), previous, boundaries[rank + 1]);
    // The boundary difference counts this walker's teeth/offspring. Carry next forward as the left
    // endpoint for the following walker. Local multiplicities sum to num_per_rank[rank], globally to N.
    (*it)->Multiplicity = next - previous;
    previous            = next;
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
  // first --> multiplicity
  // second -->  walker index in good_walkers
  std::vector<std::pair<int, int>> ncopy_pairs;
  for (int iw = 0; iw < good_walkers.size(); iw++)
  {
    ncopy_pairs.push_back(std::make_pair(static_cast<int>(good_walkers[iw]->Multiplicity), iw));
  }
  // sort good walkers by the number of copies
  std::sort(ncopy_pairs.begin(), ncopy_pairs.end());

  struct job
  {
    // Walker_index is just its index not its "walker_id"
    const int walker_index;
    const int target;
    job(int wid, int target_in) : walker_index(wid), target(target_in) {};
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
      myComm.getCommMPI3().send_value(nsentcopy, minus[ic]);
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
      myComm.getCommMPI3().receive_n(&nsentcopy, 1, plus[ic]);
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
      good_walkers[jobit->walker_index]->SendInProgress = false;
    for (auto jobit = job_list.begin(); jobit != job_list.end(); jobit++)
    {
      // pack data and send
      auto& awalker   = good_walkers[jobit->walker_index];
      size_t byteSize = awalker->byteSize();
      if (!awalker->SendInProgress)
      {
        awalker->updateBuffer();
        awalker->SendInProgress = true;
      }
      if (use_nonblocking_)
        requests.push_back(myComm.getCommMPI3().isend_n(awalker->DataSet.data(), byteSize, jobit->target));
      else
      {
        ScopedTimer local_timer(my_timers_[WC_send]);
        myComm.getCommMPI3().send_n(awalker->DataSet.data(), byteSize, jobit->target);
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
    auto unpackWalker = [](auto& awalker) {
      auto walker_id = awalker.getWalkerID();
      // Walker::copyFromBuffer overwrites the walker_id
      awalker.copyFromBuffer();
      auto parent_id = awalker.getWalkerID();
      // The logic of walker_id's is that they do not change once a walker becomes living so we write it back
      awalker.setWalkerID(walker_id);
      // And set the parent_id to this walkers parent.
      awalker.setParentID(parent_id);
    };

    std::vector<mpi3::request> requests;
    for (auto jobit = job_list.begin(); jobit != job_list.end(); jobit++)
    {
      // recv and unpack data
      auto& walker_elements = newW[jobit->walker_index];
      auto& awalker         = walker_elements.walker;
      size_t byteSize       = awalker.byteSize();
      if (use_nonblocking_)
        requests.push_back(myComm.getCommMPI3().ireceive_n(awalker.DataSet.data(), byteSize, jobit->target));
      else
      {
        ScopedTimer local_timer(my_timers_[WC_recv]);
        myComm.getCommMPI3().receive_n(awalker.DataSet.data(), byteSize, jobit->target);
        unpackWalker(awalker);
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
              auto& walker_elements = newW[job_list[im].walker_index];
              // Here the walker ID from the spawn is overwritten with the copy.
              auto& awalker = walker_elements.walker;
              unpackWalker(walker_elements.walker);
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
  {
    std::ostringstream error_msg;
    error_msg << "Multiplicity check failed in WalkerControl::swapWalkersSimple!\n"
              << "for rank: " << rank_num_ << " total_multiplicity: " << TotalMultiplicity
              << "  fair_offset_[rank_num_ + 1] - fair_offset_[rank_num_]: " << fair_offset_[rank_num_ + 1] << " - "
              << fair_offset_[rank_num_] << '\n';
    throw std::runtime_error(error_msg.str());
  }
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

std::vector<WalkerControl::IndexType> WalkerControl::syncFutureWalkersPerRank(Communicate& comm, IndexType n_walkers)
{
  int ncontexts = comm.size();
  std::vector<IndexType> future_walkers(ncontexts, 0);
  future_walkers[comm.rank()] = n_walkers;
  comm.allreduce(future_walkers);
  return future_walkers;
}

bool WalkerControl::put(xmlNodePtr cur)
{
  ParameterSet params;
  params.add(use_nonblocking_, "use_nonblocking", {true});
  params.add(debug_disable_branching_, "debug_disable_branching", {false});

  try
  {
    bool success = params.put(cur);
  }
  catch (const std::runtime_error& re)
  {
    throw UniformCommunicateError("WalkerControl::put parsing error. " + std::string(re.what()));
  }

  app_log() << "  WalkerControl parameters " << std::endl;
  app_log() << "    Using " << (use_nonblocking_ ? "non-" : "") << "blocking send/recv" << std::endl;
  if (debug_disable_branching_)
    app_log() << "    Disable branching for debugging as the user input request." << std::endl;
  return true;
}

} // namespace qmcplusplus
