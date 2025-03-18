//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2025 QMCPACK developers.
//
// File developed by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//                    Peter W. Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
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
template<typename REAL>
class KContainerT
{
public:
  using Real               = REAL;
  using FullPrecReal       = QMCTraits::FullPrecRealType;
  using PositionFull       = QMCTraits::QTFull::PosType;
  static constexpr int DIM = OHMMS_DIM;
  using AppPosition        = typename QMCTraits::PosType;
  using Position           = typename QMCTypes<Real, DIM>::PosType;

private:
  /// The cutoff up to which k-vectors are generated.
  FullPrecReal kcutoff;

public:
  //Typedef for the lattice-type
  using ParticleLayout = CrystalLattice<Real, OHMMS_DIM>;

  const auto& get_kpts_cart_soa() const { return kpts_cart_soa_; }

  const std::vector<TinyVector<int, DIM>>& getKpts() const { return kpts_; }
  const std::vector<AppPosition>& getKptsCartWorking() const;

  const std::vector<int>& getKShell() const { return kshell; };

  const std::vector<Real>& getKSQWorking() const;

  int getMinusK(int k) const;

  int getNumK() const { return numk; }

  /** update k-vectors 
   * @param sc supercell
   * @param kc cutoff radius in the K
   * @param twist shifts the center of the grid of k-vectors
   * @param useSphere if true, use the |K|
   */
  template<typename TT>
  void updateKLists(const CrystalLattice<TT, OHMMS_DIM>& lattice,
                    FullPrecReal kc,
                    unsigned ndim,
                    const Position& twist = Position(),
                    bool useSphere        = true);

private:
  ///number of k-points
  int numk;

  /** maximum integer translations of reciprocal cell within kc.
   *
   * Last index is max. of first dimension+1
   */
  TinyVector<int, DIM + 1> mmax;

  /** K-vector in reduced coordinates
   */
  std::vector<TinyVector<int, DIM>> kpts_;
  /** K-vector in Cartesian coordinates
   */
  std::vector<Position> kpts_cart_;
  std::vector<AppPosition> kpts_cart_working_;
  /** squre of kpts in Cartesian coordniates
   */
  std::vector<FullPrecReal> ksq_;
  std::vector<Real> ksq_working_;

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

private:
  /** compute approximate parallelpiped that surrounds kc
   * @param lattice supercell
   */
  template<typename TT>
  void findApproxMMax(const CrystalLattice<TT, OHMMS_DIM>& lattice, unsigned ndim);
  /** construct the container for k-vectors */
  template<typename TT>
  void BuildKLists(const CrystalLattice<TT, OHMMS_DIM>& lattice, const Position& twist, bool useSphere);

  /** K-vector in Cartesian coordinates in SoA layout
   */
  VectorSoaContainer<FullPrecReal, DIM, OffloadAllocator<FullPrecReal>> kpts_cart_soa_;
};

using KContainer = KContainerT<QMCTraits::RealType>;

#ifdef MIXED_PRECISION
extern template class KContainerT<float>;
extern template void KContainerT<float>::updateKLists<float>(const CrystalLattice<float, OHMMS_DIM>& lattice,
                                                             FullPrecReal kc,
                                                             unsigned ndim,
                                                             const Position& twist = Position(),
                                                             bool useSphere        = true);
extern template void KContainerT<float>::updateKLists<double>(const CrystalLattice<double, OHMMS_DIM>& lattice,
                                                              FullPrecReal kc,
                                                              unsigned ndim,
                                                              const Position& twist = Position(),
                                                              bool useSphere        = true);
extern template void KContainerT<float>::findApproxMMax<float>(const CrystalLattice<float, OHMMS_DIM>& lattice,
                                                               unsigned ndim);
extern template void KContainerT<float>::findApproxMMax<double>(const CrystalLattice<double, OHMMS_DIM>& lattice,
                                                                unsigned ndim);

extern template class KContainerT<double>;
extern template void KContainerT<double>::updateKLists<float>(const CrystalLattice<float, OHMMS_DIM>& lattice,
                                                              FullPrecReal kc,
                                                              unsigned ndim,
                                                              const Position& twist = Position(),
                                                              bool useSphere        = true);
extern template void KContainerT<double>::updateKLists<double>(const CrystalLattice<double, OHMMS_DIM>& lattice,
                                                               FullPrecReal kc,
                                                               unsigned ndim,
                                                               const Position& twist = Position(),
                                                               bool useSphere        = true);
extern template void KContainerT<double>::findApproxMMax<float>(const CrystalLattice<float, OHMMS_DIM>& lattice,
                                                                unsigned ndim);
extern template void KContainerT<double>::findApproxMMax<double>(const CrystalLattice<double, OHMMS_DIM>& lattice,
                                                                 unsigned ndim);
#else
extern template class KContainerT<double>;
extern template void KContainerT<double>::updateKLists<float>(const CrystalLattice<float, OHMMS_DIM>& lattice,
                                                              FullPrecReal kc,
                                                              unsigned ndim,
                                                              const Position& twist = Position(),
                                                              bool useSphere        = true);
extern template void KContainerT<double>::updateKLists<double>(const CrystalLattice<double, OHMMS_DIM>& lattice,
                                                               FullPrecReal kc,
                                                               unsigned ndim,
                                                               const Position& twist = Position(),
                                                               bool useSphere        = true);
extern template void KContainerT<double>::findApproxMMax<float>(const CrystalLattice<float, OHMMS_DIM>& lattice,
                                                                unsigned ndim);
extern template void KContainerT<double>::findApproxMMax<double>(const CrystalLattice<double, OHMMS_DIM>& lattice,
                                                                 unsigned ndim);
#endif

} // namespace qmcplusplus

#endif
