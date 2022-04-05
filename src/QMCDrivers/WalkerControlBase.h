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
#include "Message/MPIObjectBase.h"
#include "Message/CommOperators.h"
// #include "QMCDrivers/ForwardWalking/ForwardWalkingStructure.h"

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
  using Walker_t = MCWalkerConfiguration::Walker_t;
  ///typedef of FullPrecRealType
  using FullPrecRealType = QMCTraits::FullPrecRealType;
  ///typedef of IndexType
  using IndexType = QMCTraits::IndexType;

  /** An enum to access curData and accumData for reduction
   *
   * curData is larger than this //LE_MAX + n_rank * T
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

  /** legacy: start controller  and initialize the IDs of walkers*/
  void setWalkerID(MCWalkerConfiguration& walkers);

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

  /** legacy: sort Walkers between good and bad and prepare branching
   *
   *  not a sort changes internal state of walkers to copy and how many of each copy
   */
  int sortWalkers(MCWalkerConfiguration& W);

  static std::vector<IndexType> syncFutureWalkersPerRank(Communicate* comm, IndexType n_walkers);

  /** legacy: apply per rank limit Nmax and Nmin
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
  ///number of particle per rank
  std::vector<int> NumPerRank;
  ///offset of the particle index
  std::vector<int> OffSet;
  ///offset of the particle index for a fair distribution
  std::vector<int> FairOffSet;

  ///filename for dmc.dat
  std::string dmcFname;
  ///file to save energy histogram
  std::unique_ptr<std::ofstream> dmcStream;
  ///Number of walkers created by this rank
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
  std::vector<std::unique_ptr<Walker_t>> good_w, bad_w;
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
