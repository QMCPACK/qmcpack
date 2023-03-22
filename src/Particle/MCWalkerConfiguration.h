//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jordan E. Vincent, University of Illinois at Urbana-Champaign
//                    Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Cynthia Gu, zg1@ornl.gov, Oak Ridge National Laboratory
//                    Raymond Clay III, j.k.rofling@gmail.com, Lawrence Livermore National Laboratory
//                    Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


/** @file MCWalkerConfiguration.h
 * @brief Declaration of a MCWalkerConfiguration
 */
#ifndef QMCPLUSPLUS_MCWALKERCONFIGURATION_H
#define QMCPLUSPLUS_MCWALKERCONFIGURATION_H
#include "Particle/ParticleSet.h"
#include "Particle/WalkerConfigurations.h"
#include "Particle/Walker.h"
#include "Particle/SampleStack.h"
#include "Utilities/IteratorUtility.h"

namespace qmcplusplus
{
//Forward declaration
class MultiChain;
class HDFWalkerOutput;
class Reptile;

/** A set of walkers that are to be advanced by Metropolis Monte Carlo.
 *
 *As a derived class from ParticleSet, MCWalkerConfiguration interacts with
 *QMCHamiltonian and TrialWaveFunction as a ParticleSet, while QMCDrivers
 *use it as multiple walkers whose configurations are advanced according
 to MC algorithms.
 *
 Each walker is represented by Walker<PosVector_t> and
 *MCWalkerConfiguration contains a list of
 *the walkers.  This class enables two possible moves:
 *<ul>
 *<li> move the entire active walkers, similarly to molecu. Suitable for
 *small and big systems with a small time step.
 *<li> move a particle for each walker. Suitable for large systems.

 *</ul>
 */
class MCWalkerConfiguration : public ParticleSet, public WalkerConfigurations
{
public:
  /**enumeration for update*/
  enum
  {
    Update_All = 0, ///move all the active walkers
    Update_Walker,  ///move a walker by walker
    Update_Particle ///move a particle by particle
  };

  using Walker_t = WalkerConfigurations::Walker_t;
  ///container type of the Properties of a Walker
  using PropertyContainer_t = Walker_t::PropertyContainer_t;
  ///container type of Walkers
  using WalkerList_t = std::vector<std::unique_ptr<Walker_t>>;
  /// FIX: a type alias of iterator for an object should not be for just one of many objects it holds.
  using iterator = WalkerList_t::iterator;
  ///const_iterator of Walker container
  using const_iterator = WalkerList_t::const_iterator;

  using ReptileList_t = UPtrVector<Reptile>;

  ///default constructor
  MCWalkerConfiguration(const SimulationCell& simulation_cell,
                        const DynamicCoordinateKind kind = DynamicCoordinateKind::DC_POS);

  ///default constructor: copy only ParticleSet
  MCWalkerConfiguration(const MCWalkerConfiguration& mcw);
  ~MCWalkerConfiguration();
  /** create numWalkers Walkers
   *
   * Append Walkers to WalkerList.
   */
  void createWalkers(int numWalkers);
  ///clean up the walker list and make a new list
  void resize(int numWalkers, int numPtcls);

  ///clean up the walker list
  using WalkerConfigurations::clear;
  ///resize Walker::PropertyHistory and Walker::PHindex:
  void resizeWalkerHistories();

  ///make random moves for all the walkers
  //void sample(iterator first, iterator last, value_type tauinv);
  ///make a random move for a walker
  void sample(iterator it, RealType tauinv);

  ///return the number of particles per walker
  inline int getParticleNum() const { return R.size(); }
  /**@}*/

  /** set LocalEnergy
   * @param e current average Local Energy
   */
  inline void setLocalEnergy(RealType e) { LocalEnergy = e; }

  /** return LocalEnergy
   */
  inline RealType getLocalEnergy() const { return LocalEnergy; }

  inline MultiChain* getPolymer() { return Polymer; }

  inline void setPolymer(MultiChain* chain) { Polymer = chain; }

  void resetWalkerProperty(int ncopy = 1);

  inline bool updatePbyP() const { return ReadyForPbyP; }

  //@{save/load/clear function for optimization
  //
  int numSamples() const { return samples.getNumSamples(); }
  ///set the number of max samples
  void setNumSamples(int n);
  ///save the position of current walkers to SampleStack
  void saveEnsemble();
  ///save the position of current walkers
  void saveEnsemble(iterator first, iterator last);
  /// load a single sample from SampleStack
  void loadSample(ParticleSet& pset, size_t iw) const;
  /// load SampleStack data to the current list of walker configurations
  void loadEnsemble();
  /// load the SampleStacks of others to the current list of walker configurations
  void loadEnsemble(std::vector<MCWalkerConfiguration*>& others, bool doclean = true);
  /** dump Samples to a file
   * @param others MCWalkerConfigurations whose samples will be collected
   * @param out engine to write the samples to state_0/walkers
   * @param np number of processors
   * @return true with non-zero samples
   *
   * CAUTION: The current implementation assumes the same amount of active walkers on all the MPI ranks.
   */
  static bool dumpEnsemble(std::vector<MCWalkerConfiguration*>& others, HDFWalkerOutput& out, int np, int nBlock);
  ///clear the ensemble
  void clearEnsemble();

  const SampleStack& getSampleStack() const { return samples; }
  SampleStack& getSampleStack() { return samples; }

  /// Transitional forwarding methods
  int getMaxSamples() const;
  //@}

protected:
  ///true if the buffer is ready for particle-by-particle updates
  bool ReadyForPbyP;
  ///update-mode index
  int UpdateMode;

  RealType LocalEnergy;

public:
  ///a collection of reptiles contained in MCWalkerConfiguration.
  ReptileList_t ReptileList;
  Reptile* reptile;

  friend class MCPopulation;

private:
  MultiChain* Polymer;

  SampleStack samples;
};
} // namespace qmcplusplus
#endif
