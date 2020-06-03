//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2019 QMCPACK developers.
//
// File developed by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@intel.com, University of Illinois at Urbana-Champaign
//                    Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//                    Anouar Benali, benali@anl.gov, Argonne National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


/** @file
 *
 * class to handle complex splines to real orbitals with splines of arbitrary precision
 */
#ifndef QMCPLUSPLUS_SPLINE_C2R_H
#define QMCPLUSPLUS_SPLINE_C2R_H

#include <memory>
#include <QMCWaveFunctions/BsplineFactory/BsplineSet.h>
#include <OhmmsSoA/Container.h>
#include <spline2/MultiBspline.hpp>
#include "Utilities/FairDivide.h"

namespace qmcplusplus
{
/** class to match std::complex<ST> spline with BsplineSet::ValueType (real) SPOs
 * @tparam ST precision of spline
 *
 * Requires temporage storage and multiplication of phase vectors
 * Internal storage use double sized arrays of ST type, aligned and padded.
 */
template<typename ST>
class SplineC2R : public BsplineSet
{
public:
  using SplineType       = typename bspline_traits<ST, 3>::SplineType;
  using BCType           = typename bspline_traits<ST, 3>::BCType;
  using DataType         = ST;
  using PointType        = TinyVector<ST, 3>;
  using SingleSplineType = UBspline_3d_d;
  // types for evaluation results
  using TT = typename BsplineSet::ValueType;
  using BsplineSet::GGGVector_t;
  using BsplineSet::GradVector_t;
  using BsplineSet::HessVector_t;
  using BsplineSet::ValueVector_t;

  using vContainer_type  = Vector<ST, aligned_allocator<ST>>;
  using gContainer_type  = VectorSoaContainer<ST, 3>;
  using hContainer_type  = VectorSoaContainer<ST, 6>;
  using ghContainer_type = VectorSoaContainer<ST, 10>;

private:
  ///primitive cell
  CrystalLattice<ST, 3> PrimLattice;
  ///\f$GGt=G^t G \f$, transformation for tensor in LatticeUnit to CartesianUnit, e.g. Hessian
  Tensor<ST, 3> GGt;
  ///number of complex bands
  int nComplexBands;
  ///multi bspline set
  std::shared_ptr<MultiBspline<ST>> SplineInst;

  vContainer_type mKK;
  VectorSoaContainer<ST, 3> myKcart;

  ///thread private ratios for reduction when using nested threading, numVP x numThread
  Matrix<TT> ratios_private;

protected:
  /// intermediate result vectors
  vContainer_type myV;
  vContainer_type myL;
  gContainer_type myG;
  hContainer_type myH;
  ghContainer_type mygH;

public:
  SplineC2R() : nComplexBands(0)
  {
    is_complex = true;
    className  = "SplineC2R";
    KeyWord    = "SplineC2R";
  }

  virtual SPOSet* makeClone() const override { return new SplineC2R(*this); }

  inline void resizeStorage(size_t n, size_t nvals)
  {
    init_base(n);
    size_t npad = getAlignedSize<ST>(2 * n);
    myV.resize(npad);
    myG.resize(npad);
    myL.resize(npad);
    myH.resize(npad);
    mygH.resize(npad);
  }

  void bcast_tables(Communicate* comm) { chunked_bcast(comm, SplineInst->getSplinePtr()); }

  void gather_tables(Communicate* comm)
  {
    if (comm->size() == 1)
      return;
    const int Nbands      = kPoints.size();
    const int Nbandgroups = comm->size();
    offset.resize(Nbandgroups + 1, 0);
    FairDivideLow(Nbands, Nbandgroups, offset);

    for (size_t ib = 0; ib < offset.size(); ib++)
      offset[ib] = offset[ib] * 2;
    gatherv(comm, SplineInst->getSplinePtr(), SplineInst->getSplinePtr()->z_stride, offset);
  }

  template<typename GT, typename BCT>
  void create_spline(GT& xyz_g, BCT& xyz_bc)
  {
    resize_kpoints();
    SplineInst = std::make_shared<MultiBspline<ST>>();
    SplineInst->create(xyz_g, xyz_bc, myV.size());

    app_log() << "MEMORY " << SplineInst->sizeInByte() / (1 << 20) << " MB allocated "
              << "for the coefficients in 3D spline orbital representation" << std::endl;
  }

  inline void flush_zero() { SplineInst->flush_zero(); }

  /** remap kPoints to pack the double copy */
  inline void resize_kpoints()
  {
#ifndef QMC_CUDA
    // GPU CUDA code doesn't allow a change of the ordering
    nComplexBands = this->remap_kpoints();
#endif
    int nk = kPoints.size();
    mKK.resize(nk);
    myKcart.resize(nk);
    for (size_t i = 0; i < nk; ++i)
    {
      mKK[i]     = -dot(kPoints[i], kPoints[i]);
      myKcart(i) = kPoints[i];
    }
  }

  void set_spline(SingleSplineType* spline_r, SingleSplineType* spline_i, int twist, int ispline, int level);

  bool read_splines(hdf_archive& h5f);

  bool write_splines(hdf_archive& h5f);

  void assign_v(const PointType& r, const vContainer_type& myV, ValueVector_t& psi, int first, int last) const;

  virtual void evaluateValue(const ParticleSet& P, const int iat, ValueVector_t& psi) override;

  virtual void evaluateDetRatios(const VirtualParticleSet& VP,
                                 ValueVector_t& psi,
                                 const ValueVector_t& psiinv,
                                 std::vector<TT>& ratios) override;

  /** assign_vgl
   */
  void assign_vgl(const PointType& r, ValueVector_t& psi, GradVector_t& dpsi, ValueVector_t& d2psi, int first, int last)
      const;

  /** assign_vgl_from_l can be used when myL is precomputed and myV,myG,myL in cartesian
   */
  void assign_vgl_from_l(const PointType& r, ValueVector_t& psi, GradVector_t& dpsi, ValueVector_t& d2psi);

  virtual void evaluateVGL(const ParticleSet& P,
                           const int iat,
                           ValueVector_t& psi,
                           GradVector_t& dpsi,
                           ValueVector_t& d2psi) override;

  void assign_vgh(const PointType& r,
                  ValueVector_t& psi,
                  GradVector_t& dpsi,
                  HessVector_t& grad_grad_psi,
                  int first,
                  int last) const;

  virtual void evaluateVGH(const ParticleSet& P,
                           const int iat,
                           ValueVector_t& psi,
                           GradVector_t& dpsi,
                           HessVector_t& grad_grad_psi) override;

  void assign_vghgh(const PointType& r,
                    ValueVector_t& psi,
                    GradVector_t& dpsi,
                    HessVector_t& grad_grad_psi,
                    GGGVector_t& grad_grad_grad_psi,
                    int first = 0,
                    int last  = -1) const;

  virtual void evaluateVGHGH(const ParticleSet& P,
                             const int iat,
                             ValueVector_t& psi,
                             GradVector_t& dpsi,
                             HessVector_t& grad_grad_psi,
                             GGGVector_t& grad_grad_grad_psi) override;

  template<class BSPLINESPO>
  friend class SplineSetReader;
  friend class BsplineReaderBase;
};

extern template class SplineC2R<float>;
extern template class SplineC2R<double>;

} // namespace qmcplusplus
#endif
