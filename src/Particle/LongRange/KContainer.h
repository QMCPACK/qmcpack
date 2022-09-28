//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


#ifndef QMCPLUSPLUS_KCONTAINER_H
#define QMCPLUSPLUS_KCONTAINER_H

#include "Configuration.h"
#include "OhmmsSoA/VectorSoaContainer.h"
#include "OMPTarget/OffloadAlignedAllocators.hpp"

namespace qmcplusplus
{
/** Container for k-points
 *
 * It generates a set of k-points that are unit-translations of the reciprocal-space
 * cell. K-points are generated within a spherical cutoff set by the supercell
 */
class KContainer : public QMCTraits
{
private:
  /// The cutoff up to which k-vectors are generated.
  RealType kcutoff;

public:
  //Typedef for the lattice-type
  using ParticleLayout = PtclOnLatticeTraits::ParticleLayout;

  ///number of k-points
  int numk;

  /** maximum integer translations of reciprocal cell within kc.
   *
   * Last index is max. of first dimension+1
   */
  TinyVector<int, DIM + 1> mmax;

  /** K-vector in reduced coordinates
   */
  std::vector<TinyVector<int, DIM>> kpts;
  /** K-vector in Cartesian coordinates
   */
  std::vector<PosType> kpts_cart;
  /** squre of kpts in Cartesian coordniates
   */
  std::vector<RealType> ksq;
  /** Given a k index, return index to -k
   */
  std::vector<int> minusk;
  /** kpts which belong to the ith-shell [kshell[i], kshell[i+1]) */
  std::vector<int> kshell;

  /** k points sorted by the |k|  excluding |k|=0
   *
   * The first for |k|
   * The second for a map to the full index. The size of the second is the degeneracy.
   */
  //std::map<int,std::vector<int>*>  kpts_sorted;

  /** update k-vectors 
   * @param sc supercell
   * @param kc cutoff radius in the K
   * @param twist shifts the center of the grid of k-vectors
   * @param useSphere if true, use the |K|
   */
  void updateKLists(const ParticleLayout& lattice,
                    RealType kc,
                    unsigned ndim,
                    const PosType& twist = PosType(),
                    bool useSphere       = true);

  const auto& get_kpts_cart_soa() const { return kpts_cart_soa_; }

private:
  /** compute approximate parallelpiped that surrounds kc
   * @param lattice supercell
   */
  void findApproxMMax(const ParticleLayout& lattice, unsigned ndim);
  /** construct the container for k-vectors */
  void BuildKLists(const ParticleLayout& lattice, const PosType& twist, bool useSphere);

  /** K-vector in Cartesian coordinates in SoA layout
   */
  VectorSoaContainer<RealType, DIM, OffloadAllocator<RealType>> kpts_cart_soa_;
};

} // namespace qmcplusplus

#endif
