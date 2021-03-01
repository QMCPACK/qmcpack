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
#include "QMCDrivers/WalkerElementsRef.h"
#include "QMCDrivers/MCPopulation.h"
#include "Message/MPIObjectBase.h"
#include "Message/CommOperators.h"
// #include "QMCDrivers/ForwardWalking/ForwardWalkingStructure.h"

//#include <boost/archive/binary_oarchive.hpp>

namespace qmcplusplus
{
namespace testing
{
class UnifiedDriverWalkerControlMPITest;
}


/** Base class to control the walkers for DMC simulations.
 *
 * The virtual functions are implemented for a serial execution with a usual birth/death
 * process. Inherited classes implement other WalkerControl algorithms by implementing
 * branch function.
 * - 2006/10/18
 * -- curData and accumData are added.
 * -- branch function handles averaging of the local energies to apply the weight over the
 *   current walker sample correctly.
 * -- removes an extra global reduction.
 *   --- curData should be used for global reductions
 * - 2007/11/07
 * -- added functions to dump dmc histograms
 */
class WalkerControlBase : public MPIObjectBase
{
public:
  ///typedef of Walker_t
  typedef MCWalkerConfiguration::Walker_t Walker_t;
  /// distinct type for "new" walker, currently same as Walker_t
  using MCPWalker = MCPopulation::MCPWalker;
  ///typedef of FullPrecRealType
  using FullPrecRealType = QMCTraits::FullPrecRealType;
  ///typedef of IndexType
  typedef QMCTraits::IndexType IndexType;

  /** An enum to access curData and accumData for reduction
   *
   * curData is larger than this //LE_MAX + n_node * T
   */
  enum
  {
    ENERGY_INDEX = 0,
    ENERGY_SQ_INDEX,
    WALKERSIZE_INDEX,
    WEIGHT_INDEX,
    EREF_INDEX,
    R2ACCEPTED_INDEX,
    R2PROPOSED_INDEX,
    FNSIZE_INDEX,
    RNONESIZE_INDEX,
    RNSIZE_INDEX,
    B_ENERGY_INDEX,
    B_WGT_INDEX,
    SENTWALKERS_INDEX,
    LE_MAX
  };

  struct PopulationAdjustment
  {
    int num_walkers{0}; // This is the number of walkers we are adjusting to
    std::vector<WalkerElementsRef> good_walkers;
    std::vector<int> copies_to_make;
    std::vector<WalkerElementsRef> bad_walkers;
  };

  struct WalkerAdjustmentCriteria
  {
    FullPrecRealType esum{0.0};
    FullPrecRealType e2sum{0.0};
    FullPrecRealType wsum{0.0};
    FullPrecRealType ecum{0.0};
    FullPrecRealType w2sum{0.0};
    FullPrecRealType besum{0.0};
    FullPrecRealType bwgtsum{0.0};
    FullPrecRealType r2_accepted{0.0};
    FullPrecRealType r2_proposed{0.0};
    int nfn{0};
    int nrn{0};
    int ngoodfn{0};
    int ncr{0};
    int nc{0};
  };


  /** default constructor
   *
   * Set the SwapMode to zero so that instantiation can be done
   */
  WalkerControlBase(Communicate* c, bool rn = false);

  /** empty destructor to clean up the derived classes */
  virtual ~WalkerControlBase();

  /** start a block */
  void start();

  /** legacy: start controller  and initialize the IDs of walkers*/
  void setWalkerID(MCWalkerConfiguration& walkers);

  /** unified driver: start controller
   *
   *  WalkerID's are initialized by MCPopulation, SOC
   */
  void setWalkerID(MCPopulation& population);

  /** take averages and writes to a file */
  void measureProperties(int iter);

  /** set the trial energy
   */
  inline void setTrialEnergy(FullPrecRealType et) { trialEnergy = et; }

  /** return a value accumulated during a block
   * @param i index of the data
   *
   * use enum for i, see DMCEnergyEstimator
   */
  inline FullPrecRealType getValue(int i) { return accumData[i]; }

  /** return a current value
   * @param i index of the data
   *
   * use enum for i, see DMCEnergyEstimator
   */
  inline FullPrecRealType getCurrentValue(int i) { return curData[i]; }

  /** legacy: return global population
   *  update properties without branching 
   *
   */
  int doNotBranch(int iter, MCWalkerConfiguration& W);

  /** unified driver: return global population
   *  update properties without branching
   */
  FullPrecRealType doNotBranch(int iter, MCPopulation& pop);

  /** legacy: sort Walkers between good and bad and prepare branching
   *
   *  not a sort changes internal state of walkers to copy and how many of each copy
   */
  int sortWalkers(MCWalkerConfiguration& W);

  static std::vector<IndexType> syncFutureWalkersPerRank(Communicate* comm, IndexType n_walkers);

  /** unified: create data structure needed to do population adjustment
   *
   *  refactored sortWalkers
   *  This data structure contains what was updated in the state.
   */
  PopulationAdjustment calcPopulationAdjustment(MCPopulation& pop);

  /** unified: do the actual adjustment
   *
   *  unfortunately right now this requires knowledge of the global context, seems unnecessary
   *  but this is why MCPopulation is handed in.
   *  This does was applyNmaxNmin used to.
   */
  void limitPopulation(PopulationAdjustment& adjust);

  /** legacy: apply per node limit Nmax and Nmin
   */
  int applyNmaxNmin(int current_population);

  /** legacy: copy good walkers to W
   */
  int copyWalkers(MCWalkerConfiguration& W);

  /** reset to accumulate data */
  virtual void reset();

  /** legacy: perform branch and swap walkers as required 
   *
   *  \return global population
   */
  virtual int branch(int iter, MCWalkerConfiguration& W, FullPrecRealType trigger);

  /** unified: perform branch and swap walkers as required 
   *
   *  \return global population
   */
  virtual FullPrecRealType branch(int iter, MCPopulation& pop);

  bool put(xmlNodePtr cur);

  void setMinMax(int nw_in, int nmax_in);

  int get_n_max() const { return n_max_; }
  int get_n_min() const { return n_min_; }
  FullPrecRealType get_target_sigma() const { return target_sigma_; }
  MCDataType<FullPrecRealType>& get_ensemble_property() { return ensemble_property_; }
  void set_ensemble_property(MCDataType<FullPrecRealType>& ensemble_property)
  {
    ensemble_property_ = ensemble_property;
  }
  IndexType get_num_contexts() const { return num_contexts_; }
  void set_write_release_nodes(bool write_release_nodes) { write_release_nodes_ = write_release_nodes; }
  IndexType get_method() const { return method_; }
  void set_method(IndexType method) { method_ = method; }

protected:
  /** makes adjustments to local population based on adjust
   *
   * \param[inout] pop the local population
   * \param[inout] the population adjustment 
   */
  static void onRankKill(MCPopulation& pop, PopulationAdjustment& adjust);
  /** makes adjustments to local population based on adjust
   *
   * \param[inout] pop the local population
   * \param[in]    the population adjustment, it is not updated to reflect local state and is now invalid.
   */
  static void onRankSpawn(MCPopulation& pop, PopulationAdjustment& adjust);


  ///id for the method
  IndexType method_;
  ///minimum number of walkers
  IndexType n_min_;
  ///maximum number of walkers
  IndexType n_max_;
  ///maximum copy per walker
  IndexType MaxCopy;
  ///current number of walkers per processor
  IndexType NumWalkers;
  ///trial energy energy
  FullPrecRealType trialEnergy;
  ///target sigma to limit fluctuations of the trial energy
  FullPrecRealType target_sigma_;
  ///number of particle per node
  std::vector<int> NumPerNode;
  ///offset of the particle index
  std::vector<int> OffSet;
  ///offset of the particle index for a fair distribution
  std::vector<int> FairOffSet;

  ///filename for dmc.dat
  std::string dmcFname;
  ///file to save energy histogram
  std::ofstream* dmcStream;
  ///Number of walkers created by this node
  IndexType NumWalkersCreated;
  ///context id
  IndexType MyContext;
  ///number of contexts
  IndexType num_contexts_;
  ///0 is default
  IndexType SwapMode;
  ///any accumulated data over a block
  std::vector<FullPrecRealType> accumData;
  ///any temporary data includes many ridiculous conversions of integral types to and from fp
  std::vector<FullPrecRealType> curData;
  ///temporary storage for good and bad walkers
  std::vector<Walker_t*> good_w, bad_w;
  ///temporary storage for copy counters
  std::vector<int> ncopy_w;
  ///Add released-node fields to .dmc.dat file
  bool write_release_nodes_;
  ///Use non-blocking isend/irecv
  bool use_nonblocking;

  ///ensemble properties
  MCDataType<FullPrecRealType> ensemble_property_;

  friend class qmcplusplus::testing::UnifiedDriverWalkerControlMPITest;
private:
  /** unified: Refactoring possibly dead releaseNodesCode out
   * @{
   */
  static auto rn_walkerCalcAdjust(MCPWalker& walker, WalkerAdjustmentCriteria wac);

  static auto addReleaseNodeWalkers(PopulationAdjustment& adjust,
                                    WalkerAdjustmentCriteria& wac,
                                    std::vector<WalkerElementsRef>& good_walkers_rn,
                                    std::vector<int>& copies_to_make_rn);
  /**}@*/

  /** unified: CalcAdjust segmenting
   * @{
   */
  static auto walkerCalcAdjust(MCPWalker& walker, WalkerAdjustmentCriteria wac);

  static void updateCurDataWithCalcAdjust(std::vector<FullPrecRealType>& data,
                                          WalkerAdjustmentCriteria wac,
                                          PopulationAdjustment& adjust,
                                          MCPopulation& pop);
  /**}@*/
};

} // namespace qmcplusplus
#endif
