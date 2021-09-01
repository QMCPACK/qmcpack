//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2020 QMCPACK developers.
//
// File developed by: Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//
// File refactored from: MomentumEstimator.h
//////////////////////////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_MOMENTUMDISTRIBUTION_H
#define QMCPLUSPLUS_MOMENTUMDISTRIBUTION_H
#include <vector>

#include "Configuration.h"
#include "OperatorEstBase.h"
#include "Containers/OhmmsPETE/TinyVector.h"

#include "MomentumDistributionInput.h"


namespace qmcplusplus
{
/** Class that collects momentum distribution of electrons
 *  
 */
class MomentumDistribution : public OperatorEstBase
{
public:
  using POLT    = PtclOnLatticeTraits;
  using Lattice = POLT::ParticleLayout_t;
  using QMCT    = QMCTraits;

  //data members
  MomentumDistributionInput input_;

  /** @ingroup MomentumDistribution mutable parameters
   */

  /** Constructor for MomentumDistributionInput 
   */
  MomentumDistribution(MomentumDistributionInput&& mdi, DataLocality dl = DataLocality::crowd);

  MomentumDistribution(const MomentumDistribution& md);

  /** This allows us to allocate the necessary data for the DataLocality::queue 
   */
  void startBlock(int steps) override;

  /** standard interface
   */
  MomentumDistribution* clone() override;

  /** accumulate 1 or more walkers of MomentumDistribution samples
   */
  void accumulate(const RefVector<MCPWalker>& walkers, const RefVector<ParticleSet>& psets) override;

  /** this allows the EstimatorManagerNew to reduce without needing to know the details
   *  of MomentumDistribution's data.
   *
   *  can use base class default until crowd level MomentumDistribution
   *  estimators don't have a copy of the density grid.
   */
  void collect(const RefVector<OperatorEstBase>& operator_estimators) override;

  /** this allows the EstimatorManagerNew to reduce without needing to know the details
   *  of MomentumDistribution's data.
   *
   *  can use base class default until crowd level MomentumDistribution estimators don't have a copy of the density grid.
   */
  //void collect(const OperatorEstBase&  oeb);

  /** this gets us into the hdf5 file
   *
   *  Just parroting for now don't fully understand.
   *, needs to be unraveled and simplified the hdf5 output is another 
   *  big state big coupling design.
   */
  void registerOperatorEstimator(hid_t gid) override;

};

} // namespace qmcplusplus

#endif /* QMCPLUSPLUS_MOMENTUMDISTRIBUTION_H */
