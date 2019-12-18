//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
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
// #include "QMCDrivers/ForwardWalking/ForwardWalkingStructure.h"

//#include <boost/archive/binary_oarchive.hpp>

namespace qmcplusplus
{
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

  /** default constructor
   *
   * Set the SwapMode to zero so that instantiation can be done
   */
  WalkerControlBase(Communicate* c, bool rn = false);

  /** empty destructor to clean up the derived classes */
  virtual ~WalkerControlBase();

  /** start a block */
  void start();

  /** start controller  and initialize the IDs of walkers*/
  void setWalkerID(MCWalkerConfiguration& walkers);

  /** start controller
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

  /** update properties without branching */
  int doNotBranch(int iter, MCWalkerConfiguration& W);

  /** update properties without branching */
  int doNotBranch(int iter, MCPopulation& pop);

  /** sort Walkers between good and bad and prepare branching
   *
   *  not a sort changes internal state of walkers to copy and how many of each copy
   */
  int sortWalkers(MCWalkerConfiguration& W);

  struct PopulationAdjustment {
    int num_walkers; // This is the number of walkers we are adjusting to
    RefVector<MCPWalker> good_walkers;
    std::vector<int> copies_to_make;
    RefVector<MCPWalker> bad_walkers;
  };

  static std::vector<IndexType> syncFutureWalkersPerRank(Communicate* comm, IndexType n_walkers );

  /** create data structure needed to do population adjustment
   *
   *  refactored sortWalkers
   *  This data structure contains what was updated in the state.
   */
  PopulationAdjustment calcPopulationAdjustment(MCPopulation& pop);

  /** do the actual adjustment
   *
   *  unfortunately right now this requires knowledge of the global context, seems unecessary
   *  but this is why MCPopulation is handed in.
   */
  int adjustPopulation(MCPopulation& pop, PopulationAdjustment& adjust);

  /** apply per node limit Nmax and Nmin
   */
  int applyNmaxNmin(int current_population);

  /** copy good walkers to W
   */
  int copyWalkers(MCWalkerConfiguration& W);

  void Write2XYZ(MCWalkerConfiguration& W);

  /** reset to accumulate data */
  virtual void reset();

  /** perform branch and swap walkers as required */
  virtual int branch(int iter, MCWalkerConfiguration& W, FullPrecRealType trigger);

  /** perform branch and swap walkers as required */
  virtual int branch(int iter, MCPopulation& pop, FullPrecRealType trigger);

  virtual FullPrecRealType getFeedBackParameter(int ngen, FullPrecRealType tau) { return 1.0 / (static_cast<FullPrecRealType>(ngen) * tau); }

  bool put(xmlNodePtr cur);

  void setMinMax(int nw_in, int nmax_in);

  int get_n_max() const { return n_max_; }
  int get_n_min() const { return n_min_; }
  FullPrecRealType get_target_sigma() const { return target_sigma_; }
  MCDataType<FullPrecRealType>& get_ensemble_property() { return ensemble_property_; }
  void set_ensemble_property(MCDataType<FullPrecRealType>& ensemble_property) { ensemble_property_ = ensemble_property; }
  IndexType get_num_contexts() const { return num_contexts_; }
  void set_write_release_nodes(bool write_release_nodes) { write_release_nodes_ = write_release_nodes; }
  IndexType get_method() const { return method_; }
  void set_method(IndexType method) { method_ = method; }

protected:
  static void onRankSpawnKill(MCPopulation& pop, PopulationAdjustment& adjust);
  
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

};

} // namespace qmcplusplus
#endif
