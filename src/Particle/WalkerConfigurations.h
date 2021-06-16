//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2020 QMCPACK developers.
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
// File created by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////


/** @file WalkerConfigurations.h
 * @brief Declaration of a WalkerConfigurations
 */
#ifndef QMCPLUSPLUS_WALKERCONFIGURATIONS_H
#define QMCPLUSPLUS_WALKERCONFIGURATIONS_H
#include "Configuration.h"
#include "Particle/Walker.h"
#include "Utilities/IteratorUtility.h"

namespace qmcplusplus
{
/** Monte Carlo Data of an ensemble
 *
 * The quantities are shared by all the nodes in a group
 * - NumSamples number of samples
 * - Weight     total weight of a sample
 * - Energy     average energy of a sample
 * - Variance   variance
 * - LivingFraction fraction of walkers alive each step.
 */
template<typename T>
struct MCDataType
{
  T NumSamples;
  T RNSamples;
  T Weight;
  T Energy;
  T AlternateEnergy;
  T Variance;
  T R2Accepted;
  T R2Proposed;
  T LivingFraction;
};

/** A set of light weight walkers that are carried between driver sections and restart
 */
class WalkerConfigurations
{
public:
  /// walker type
  using Walker_t         = Walker<QMCTraits, PtclOnLatticeTraits>;
  using FullPrecRealType = QMCTraits::FullPrecRealType;
  ///container type of Walkers
  using WalkerList_t = std::vector<std::unique_ptr<Walker_t>>;
  /// FIX: a type alias of iterator for an object should not be for just one of many objects it holds.
  using iterator = WalkerList_t::iterator;
  ///const_iterator of Walker container
  using const_iterator = WalkerList_t::const_iterator;

  /** starting index of the walkers in a processor group
   *
   * WalkerOffsets[0]=0 and WalkerOffsets[WalkerOffsets.size()-1]=total number of walkers in a group
   * WalkerOffsets[processorid+1]-WalkerOffsets[processorid] is equal to the number of walkers on a processor,
   * i.e., W.getActiveWalkers().
   * WalkerOffsets is added to handle parallel I/O with hdf5
   */
  std::vector<int> WalkerOffsets;

  MCDataType<FullPrecRealType> EnsembleProperty;

  WalkerConfigurations();
  ~WalkerConfigurations();
  WalkerConfigurations(const WalkerConfigurations&) = delete;
  WalkerConfigurations& operator=(const WalkerConfigurations&) = delete;
  WalkerConfigurations(WalkerConfigurations&&)                 = default;
  WalkerConfigurations& operator=(WalkerConfigurations&&) = default;

  /** create numWalkers Walkers
   *
   * Append Walkers to WalkerList.
   */
  void createWalkers(int numWalkers, size_t numPtcls);
  /** create walkers
   * @param first walker iterator
   * @param last walker iterator
   */
  void createWalkers(iterator first, iterator last);
  /** copy walkers
   * @param first input walker iterator
   * @param last input walker iterator
   * @param start first target iterator
   *
   * No memory allocation is allowed.
   */
  void copyWalkers(iterator first, iterator last, iterator start);

  /** destroy Walkers from itstart to itend
   *@param first starting iterator of the walkers
   *@param last ending iterator of the walkers
   */
  iterator destroyWalkers(iterator first, iterator last);

  /** destroy Walkers
   *@param nw number of walkers to be destroyed
   */
  void destroyWalkers(int nw);

  ///clean up the walker list and make a new list
  void resize(int numWalkers, size_t numPtcls);

  ///return the number of active walkers
  inline size_t getActiveWalkers() const { return WalkerList.size(); }
  ///return the total number of active walkers among a MPI group
  inline size_t getGlobalNumWalkers() const { return GlobalNumWalkers; }
  ///return the total number of active walkers among a MPI group
  inline void setGlobalNumWalkers(size_t nw)
  {
    GlobalNumWalkers            = nw;
    EnsembleProperty.NumSamples = nw;
    EnsembleProperty.Weight     = nw;
  }

  inline void setWalkerOffsets(const std::vector<int>& o) { WalkerOffsets = o; }

  /// return the first iterator
  inline iterator begin() { return WalkerList.begin(); }
  /// return the last iterator, [begin(), end())
  inline iterator end() { return WalkerList.end(); }

  /// return the first const_iterator
  inline const_iterator begin() const { return WalkerList.begin(); }

  /// return the last const_iterator  [begin(), end())
  inline const_iterator end() const { return WalkerList.end(); }
  /**@}*/

  /** clear the WalkerList without destroying them
   *
   * Provide std::vector::clear interface
   */
  inline void clear() { WalkerList.clear(); }

  /** insert elements
   * @param it locator where the inserting begins
   * @param first starting iterator
   * @param last ending iterator
   *
   * Provide std::vector::insert interface
   */
  template<class INPUT_ITER>
  inline void insert(iterator it, INPUT_ITER first, INPUT_ITER last)
  {
    WalkerList.insert(it, first, last);
  }

  /** add Walker_t* at the end
   * @param awalker pointer to a walker
   *
   * Provide std::vector::push_back interface
   */
  inline void push_back(std::unique_ptr<Walker_t> awalker) { WalkerList.push_back(std::move(awalker)); }

  /** delete the last Walker_t*
   *
   * Provide std::vector::pop_back interface
   */
  inline void pop_back()
  {
    WalkerList.pop_back();
  }

  inline Walker_t* operator[](int i) { return WalkerList[i].get(); }

  inline const Walker_t* operator[](int i) const { return WalkerList[i].get(); }

  /** reset the Walkers
   */
  void reset();

protected:
  ///number of walkers on a node
  int LocalNumWalkers;
  ///number of walkers shared by a MPI group
  size_t GlobalNumWalkers;

public:
  ///a collection of walkers
  WalkerList_t WalkerList;
};
} // namespace qmcplusplus
#endif
