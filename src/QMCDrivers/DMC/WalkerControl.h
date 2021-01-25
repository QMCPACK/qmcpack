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
 * -- curData is added.
 * -- branch function handles averaging of the local energies to apply the weight over the
 *   current walker sample correctly.
 * -- removes an extra global reduction.
 *   --- curData should be used for global reductions
 * - 2007/11/07
 * -- added functions to dump dmc histograms
 */
class WalkerControl : public MPIObjectBase
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

  /** An enum to access curData for reduction
   *
   * curData is larger than this //LE_MAX + n_node * T
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
    FullPrecRealType r2_accepted{0.0};
    FullPrecRealType r2_proposed{0.0};
    int nfn{0};
    int ngoodfn{0};
    int nc{0};
  };


  /** default constructor
   *
   * Set the SwapMode to zero so that instantiation can be done
   */
  WalkerControl(Communicate* c, RandomGenerator_t& rng, bool use_fixed_pop = false);

  /** empty destructor to clean up the derived classes */
  virtual ~WalkerControl();

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

  static std::vector<IndexType> syncFutureWalkersPerRank(Communicate* comm, IndexType n_walkers);

  /** legacy: apply per node limit Nmax and Nmin
   */
  int applyNmaxNmin(int current_population);

  void Write2XYZ(MCWalkerConfiguration& W);

  /** reset to accumulate data */
  virtual void reset();

  /** unified: perform branch and swap walkers as required 
   *
   *  \return global population
   */
  virtual FullPrecRealType branch(int iter, MCPopulation& pop);

  virtual FullPrecRealType getFeedBackParameter(int ngen, FullPrecRealType tau)
  {
    return 1.0 / (static_cast<FullPrecRealType>(ngen) * tau);
  }

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

protected:
  ///random number generator
  RandomGenerator_t& rng_;
  ///if true, use fixed population
  bool use_fixed_pop_;
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
  ///number of walkers on each node after branching before load balancing
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
  ///any temporary data includes many ridiculous conversions of integral types to and from fp
  std::vector<FullPrecRealType> curData;
  ///temporary storage for good and bad walkers
  std::vector<Walker_t*> good_w, bad_w;
  ///temporary storage for copy counters
  std::vector<int> ncopy_w;
  ///Use non-blocking isend/irecv
  bool use_nonblocking;

  ///ensemble properties
  MCDataType<FullPrecRealType> ensemble_property_;

  friend class qmcplusplus::testing::UnifiedDriverWalkerControlMPITest;

private:
  void computeCurData(const UPtrVector<MCPWalker>& walkers);
};

} // namespace qmcplusplus
#endif
