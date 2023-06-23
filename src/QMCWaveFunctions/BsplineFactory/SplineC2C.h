//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2019 QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@intel.com, Intel Corp.
//                    Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@intel.com, Intel Corp.
//////////////////////////////////////////////////////////////////////////////////////


/** @file
 *
 * class to handle complex splines to complex orbitals with splines of arbitrary precision
 */
#ifndef QMCPLUSPLUS_SPLINE_C2C_H
#define QMCPLUSPLUS_SPLINE_C2C_H

#include <memory>
#include "QMCWaveFunctions/BsplineFactory/BsplineSet.h"
#include "OhmmsSoA/VectorSoaContainer.h"
#include "spline2/MultiBspline.hpp"
#include "Utilities/FairDivide.h"

namespace qmcplusplus
{
/** class to match std::complex<ST> spline with BsplineSet::ValueType (complex) SPOs
 * @tparam ST precision of spline
 *
 * Requires temporage storage and multiplication of phase vectors
 * The internal storage of complex spline coefficients uses double sized real arrays of ST type, aligned and padded.
 * All the output orbitals are complex.
 */
template<typename ST>
class SplineC2C : public BsplineSet
{
public:
  using SplineType       = typename bspline_traits<ST, 3>::SplineType;
  using BCType           = typename bspline_traits<ST, 3>::BCType;
  using DataType         = ST;
  using PointType        = TinyVector<ST, 3>;
  using SingleSplineType = UBspline_3d_d;

  // types for evaluation results
  using ComplexT = typename BsplineSet::ValueType;
  using BsplineSet::GGGVector;
  using BsplineSet::GradVector;
  using BsplineSet::HessVector;
  using BsplineSet::ValueVector;

  using vContainer_type  = Vector<ST, aligned_allocator<ST>>;
  using gContainer_type  = VectorSoaContainer<ST, 3>;
  using hContainer_type  = VectorSoaContainer<ST, 6>;
  using ghContainer_type = VectorSoaContainer<ST, 10>;

private:
  ///primitive cell
  CrystalLattice<ST, 3> PrimLattice;
  ///\f$GGt=G^t G \f$, transformation for tensor in LatticeUnit to CartesianUnit, e.g. Hessian
  Tensor<ST, 3> GGt;
  ///multi bspline set
  std::shared_ptr<MultiBspline<ST>> SplineInst;

  ///Copy of original splines for orbital rotation
  std::shared_ptr<std::vector<RealType>> coef_copy_;

  vContainer_type mKK;
  VectorSoaContainer<ST, 3> myKcart;

  ///thread private ratios for reduction when using nested threading, numVP x numThread
  Matrix<ComplexT> ratios_private;

protected:
  /// intermediate result vectors
  vContainer_type myV;
  vContainer_type myL;
  gContainer_type myG;
  hContainer_type myH;
  ghContainer_type mygH;

public:
  SplineC2C(const std::string& my_name) : BsplineSet(my_name) {}

  SplineC2C(const SplineC2C& in);
  virtual std::string getClassName() const override { return "SplineC2C"; }
  virtual std::string getKeyword() const override { return "SplineC2C"; }
  bool isComplex() const override { return true; };


  std::unique_ptr<SPOSet> makeClone() const override { return std::make_unique<SplineC2C>(*this); }

  bool isRotationSupported() const override { return true; }

  /// Store an original copy of the spline coefficients for orbital rotation
  void storeParamsBeforeRotation() override;

  /*
    Implements orbital rotations via [1,2].
    Should be called by RotatedSPOs::apply_rotation()
    This implementation requires that NSPOs > Nelec. In other words,
    if you want to run a orbopt wfn, you must include some virtual orbitals!
    Some results (using older Berkeley branch) were published in [3].
    [1] Filippi & Fahy, JCP 112, (2000)
    [2] Toulouse & Umrigar, JCP 126, (2007)
    [3] Townsend et al., PRB 102, (2020)
  */
  void applyRotation(const ValueMatrix& rot_mat, bool use_stored_copy) override;

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
      offset[ib] *= 2;
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
    const size_t nk = kPoints.size();
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

  void assign_v(const PointType& r, const vContainer_type& myV, ValueVector& psi, int first, int last) const;

  void evaluateValue(const ParticleSet& P, const int iat, ValueVector& psi) override;

  void evaluateDetRatios(const VirtualParticleSet& VP,
                         ValueVector& psi,
                         const ValueVector& psiinv,
                         std::vector<ValueType>& ratios) override;

  /** assign_vgl
   */
  void assign_vgl(const PointType& r, ValueVector& psi, GradVector& dpsi, ValueVector& d2psi, int first, int last)
      const;

  /** assign_vgl_from_l can be used when myL is precomputed and myV,myG,myL in cartesian
   */
  void assign_vgl_from_l(const PointType& r, ValueVector& psi, GradVector& dpsi, ValueVector& d2psi);

  void evaluateVGL(const ParticleSet& P,
                   const int iat,
                   ValueVector& psi,
                   GradVector& dpsi,
                   ValueVector& d2psi) override;

  void assign_vgh(const PointType& r,
                  ValueVector& psi,
                  GradVector& dpsi,
                  HessVector& grad_grad_psi,
                  int first,
                  int last) const;

  void evaluateVGH(const ParticleSet& P,
                   const int iat,
                   ValueVector& psi,
                   GradVector& dpsi,
                   HessVector& grad_grad_psi) override;

  void assign_vghgh(const PointType& r,
                    ValueVector& psi,
                    GradVector& dpsi,
                    HessVector& grad_grad_psi,
                    GGGVector& grad_grad_grad_psi,
                    int first = 0,
                    int last  = -1) const;

  void evaluateVGHGH(const ParticleSet& P,
                     const int iat,
                     ValueVector& psi,
                     GradVector& dpsi,
                     HessVector& grad_grad_psi,
                     GGGVector& grad_grad_grad_psi) override;

  template<class BSPLINESPO>
  friend struct SplineSetReader;
  friend struct BsplineReaderBase;
};

extern template class SplineC2C<float>;
extern template class SplineC2C<double>;

} // namespace qmcplusplus
#endif
