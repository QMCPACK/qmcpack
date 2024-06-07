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


#ifndef QMCPLUSPLUS_WALKER_CONTROL_BASE_H
#define QMCPLUSPLUS_WALKER_CONTROL_BASE_H

#include "Configuration.h"
#include "Particle/MCWalkerConfiguration.h"
#include "QMCDrivers/MCPopulation.h"
#include "Message/MPIObjectBase.h"
#include "Message/CommOperators.h"
#include "Utilities/RandomGenerator.h"

#include <filesystem>

namespace qmcplusplus
{
namespace testing
{
class UnifiedDriverWalkerControlMPITest;
}

/** Class for controlling the walkers for DMC simulations.
 * w and w/o MPI. Fixed and dynamic population in one place.
 */
class WalkerControl : public MPIObjectBase
{
public:
  ///typedef of Walker_t
  using Walker_t = MCWalkerConfiguration::Walker_t;
  /// distinct type for "new" walker, currently same as Walker_t
  using MCPWalker = MCPopulation::MCPWalker;
  ///typedef of FullPrecRealType
  using FullPrecRealType = QMCTraits::FullPrecRealType;
  ///typedef of IndexType
  using IndexType = QMCTraits::IndexType;

  /** default constructor
   *
   * Set the SwapMode to zero so that instantiation can be done
   */
  WalkerControl(Communicate* c, RandomBase<FullPrecRealType>& rng, bool use_fixed_pop = false);

  /** empty destructor to clean up the derived classes */
  ~WalkerControl();

  /** start a block */
  void start();

  /** take averages and writes to a file */
  void writeDMCdat(int iter, const std::vector<FullPrecRealType>& curData);

  /** set the trial energy for writing to dmc.dat
   */
  inline void setTrialEnergy(FullPrecRealType et) { trial_energy_ = et; }

  /** unified: perform branch and swap (balance) walkers as required
   *  **This has many side effects**
   *  ## For:
   *  ### dynamic population
   *  1. compute multiplicity. If iter 0 and all of warmup -> multiplicity = 1
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
  void branch(int iter, MCPopulation& pop, bool do_not_branch);

  bool put(xmlNodePtr cur);

  void setMinMax(int nw_in, int nmax_in);

  int get_n_max() const { return n_max_; }
  int get_n_min() const { return n_min_; }
  MCDataType<FullPrecRealType>& get_ensemble_property() { return ensemble_property_; }
  void set_ensemble_property(MCDataType<FullPrecRealType>& ensemble_property)
  {
    ensemble_property_ = ensemble_property;
  }
  IndexType get_num_contexts() const { return num_ranks_; }
  const std::vector<int>& getNumPerRank() { return num_per_rank_; }
private:
  // Can't really unit test without this.
  void setNumPerRank(const std::vector<int>& num_per_rank) { num_per_rank_ = num_per_rank; }

  /// kill dead walkers in the population
  static void killDeadWalkersOnRank(MCPopulation& pop);

  static std::vector<IndexType> syncFutureWalkersPerRank(Communicate* comm, IndexType n_walkers);

  /** update the curData state buffer.
   *  weighted sum over walker properties then all reduce.  see discussion #curData
   */
  void computeCurData(const UPtrVector<MCPWalker>& walkers, std::vector<FullPrecRealType>& curData);

  /** creates the distribution plan
   *
   *  populates the minus and plus vectors they contain 1 copy of a partition index 
   *  for each adjustment in population to the context.
   *  \param[in] #num_per_rank
   *  \param[out] fair_offset running population count at each partition boundary
   *  \param[out] minus list of partition indexes one occurrence for each walker removed
   *  \param[out] plus list of partition indexes one occurrence for each walker added
   */
  static void determineNewWalkerPopulation(const std::vector<int>& num_per_rank,
                                           std::vector<int>& fair_offset,
                                           std::vector<int>& minus,
                                           std::vector<int>& plus);

#if defined(HAVE_MPI)
  /** swap Walkers with Recv/Send or Irecv/Isend
   *
   * The algorithm ensures that the load per node can differ only by one walker.
   * Each MPI rank can only send or receive or be silent.

   * If multiple copies of a walker need to be sent to the target rank, only send one.
   * The number of copies is communicated ahead via blocking send/recv.

   * The communication is one-dimensional and depending on your node layout and network topology
   * could be very local.
   * Then the walkers are transferred via blocking or non-blocking send/recv.
   * The blocking send/recv may become serialized and worsen load imbalance.
   * Non blocking send/recv algorithm avoids serialization completely.
   */
  void swapWalkersSimple(MCPopulation& pop);
#endif

  /** An enum to access/document curData's elements, this is just a subset of curData's indexes
   * curData[LE_MAX:LE_MAX+num_ranks-1] are sums of walker multiplicty by rank
   */
  enum
  {
    ENERGY_INDEX = 0,
    ENERGY_SQ_INDEX,
    WALKERSIZE_INDEX,
    WEIGHT_INDEX,
    R2ACCEPTED_INDEX,
    R2PROPOSED_INDEX,
    FNSIZE_INDEX,
    SENTWALKERS_INDEX,
    LE_MAX
  };

  ///random number generator
  RandomBase<FullPrecRealType>& rng_;
  ///if true, use fixed population
  bool use_fixed_pop_;
  ///minimum number of walkers
  IndexType n_min_;
  ///maximum number of walkers
  IndexType n_max_;
  ///maximum copy per walker
  IndexType max_copy_;
  ///trial energy energy
  FullPrecRealType trial_energy_;
  /** number of walkers on each MPI wrt to the branching algorithm before load balancing
   *  not the same as WalkerControl::branch.
   *  For dynamic population: std::floor(sum over walkers of  std::floor(walker weight + rng))
   *  For fixed population: std::floor(sum over walkers weight)

   */
  std::vector<int> num_per_rank_;
  ///offset of the particle index for a fair distribution
  std::vector<int> fair_offset_;
  ///filename for dmc.dat
  std::filesystem::path dmcFname;
  ///file to save energy histogram
  std::unique_ptr<std::ofstream> dmcStream;
  ///context id
  const IndexType rank_num_;
  ///number of contexts
  const IndexType num_ranks_;
  ///0 is default
  IndexType SwapMode;
  /** large amount of state data MPI_AllReduce once per step size is LE_MAX + number of ranks.
   *  partial enum to access elements is the enum above dumped into WalkerControl's namespace
   *  includes many integral types converted to and from fp
   *  \todo MPI in the 21st centure allows user defined data types that would allow all this to benefit from actual
   *        C++ type safety and still can have performant collective operation efficiency.
   *        The desired operation isn't actually an all reduce for potentially the largest number of elements.
   *        From LE_MAX:LEMAX_+num_ranks-1 is an "integer" allgather for dynamic populations
   *        and a floating point all gather for fixed population.  
   */
  std::vector<FullPrecRealType> curData;
  ///Use non-blocking isend/irecv
  bool use_nonblocking_;
  ///disable branching for debugging
  bool debug_disable_branching_;
  ///ensemble properties
  MCDataType<FullPrecRealType> ensemble_property_;
  ///timers
  TimerList_t my_timers_;
  ///Number of walkers sent during the exchange
  IndexType saved_num_walkers_sent_;

  friend testing::UnifiedDriverWalkerControlMPITest;
};

} // namespace qmcplusplus
#endif
