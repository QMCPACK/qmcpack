//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2019 QMCPACK developers.
//
// File developed by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//                    Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


#ifndef QMCPLUSPLUS_SPLINE_R2R_H
#define QMCPLUSPLUS_SPLINE_R2R_H

#include <memory>
#include "QMCWaveFunctions/BsplineFactory/BsplineSet.h"
#include "OhmmsSoA/VectorSoaContainer.h"
#include "spline2/MultiBspline.hpp"
#include "Utilities/FairDivide.h"

namespace qmcplusplus
{
/** class to match ST real spline with BsplineSet::ValueType (real) SPOs
 * @tparam ST precision of spline
 *
 * Requires temporage storage and multiplication of the sign of the real part of the phase
 * Internal storage ST type arrays are aligned and padded.
 */
template<typename ST>
class SplineR2R : public BsplineSet
{
public:
  using SplineType       = typename bspline_traits<ST, 3>::SplineType;
  using BCType           = typename bspline_traits<ST, 3>::BCType;
  using DataType         = ST;
  using PointType        = TinyVector<ST, 3>;
  using SingleSplineType = UBspline_3d_d;
  // types for evaluation results
  using TT = typename BsplineSet::ValueType;
  using BsplineSet::GGGVector;
  using BsplineSet::GradVector;
  using BsplineSet::HessVector;
  using BsplineSet::ValueVector;

  using vContainer_type  = Vector<ST, aligned_allocator<ST>>;
  using gContainer_type  = VectorSoaContainer<ST, 3>;
  using hContainer_type  = VectorSoaContainer<ST, 6>;
  using ghContainer_type = VectorSoaContainer<ST, 10>;

private:
  bool IsGamma;
  ///\f$GGt=G^t G \f$, transformation for tensor in LatticeUnit to CartesianUnit, e.g. Hessian
  Tensor<ST, 3> GGt;
  ///multi bspline set
  std::shared_ptr<MultiBspline<ST>> SplineInst;

  ///thread private ratios for reduction when using nested threading, numVP x numThread
  Matrix<TT> ratios_private;

protected:
  ///primitive cell
  CrystalLattice<ST, 3> PrimLattice;
  /// intermediate result vectors
  vContainer_type myV;
  vContainer_type myL;
  gContainer_type myG;
  hContainer_type myH;
  ghContainer_type mygH;

public:
  SplineR2R()
  {
    is_complex = false;
    className  = "SplineR2R";
    KeyWord    = "SplineR2R";
  }

  std::unique_ptr<SPOSet> makeClone() const override { return std::make_unique<SplineR2R>(*this); }

  inline void resizeStorage(size_t n, size_t nvals)
  {
    init_base(n);
    const size_t npad = getAlignedSize<ST>(n);
    myV.resize(npad);
    myG.resize(npad);
    myL.resize(npad);
    myH.resize(npad);
    mygH.resize(npad);

    IsGamma = ((HalfG[0] == 0) && (HalfG[1] == 0) && (HalfG[2] == 0));
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
    gatherv(comm, SplineInst->getSplinePtr(), SplineInst->getSplinePtr()->z_stride, offset);
  }

  template<typename GT, typename BCT>
  void create_spline(GT& xyz_g, BCT& xyz_bc)
  {
    GGt        = dot(transpose(PrimLattice.G), PrimLattice.G);
    SplineInst = std::make_shared<MultiBspline<ST>>();
    SplineInst->create(xyz_g, xyz_bc, myV.size());

    app_log() << "MEMORY " << SplineInst->sizeInByte() / (1 << 20) << " MB allocated "
              << "for the coefficients in 3D spline orbital representation" << std::endl;
  }

  inline void flush_zero() { SplineInst->flush_zero(); }

  void set_spline(SingleSplineType* spline_r, SingleSplineType* spline_i, int twist, int ispline, int level);

  bool read_splines(hdf_archive& h5f);

  bool write_splines(hdf_archive& h5f);

  /** convert position in PrimLattice unit and return sign */
  inline int convertPos(const PointType& r, PointType& ru)
  {
    ru          = PrimLattice.toUnit(r);
    int bc_sign = 0;
    for (int i = 0; i < D; i++)
      if (-std::numeric_limits<ST>::epsilon() < ru[i] && ru[i] < 0)
        ru[i] = ST(0.0);
      else
      {
        ST img = std::floor(ru[i]);
        ru[i] -= img;
        bc_sign += HalfG[i] * (int)img;
      }
    return bc_sign;
  }

  void assign_v(int bc_sign, const vContainer_type& myV, ValueVector& psi, int first, int last) const;

  void evaluateValue(const ParticleSet& P, const int iat, ValueVector& psi) override;

  void evaluateDetRatios(const VirtualParticleSet& VP,
                         ValueVector& psi,
                         const ValueVector& psiinv,
                         std::vector<TT>& ratios) override;

  void assign_vgl(int bc_sign, ValueVector& psi, GradVector& dpsi, ValueVector& d2psi, int first, int last) const;

  /** assign_vgl_from_l can be used when myL is precomputed and myV,myG,myL in cartesian
   */
  void assign_vgl_from_l(int bc_sign, ValueVector& psi, GradVector& dpsi, ValueVector& d2psi);

  void evaluateVGL(const ParticleSet& P,
                   const int iat,
                   ValueVector& psi,
                   GradVector& dpsi,
                   ValueVector& d2psi) override;

  void assign_vgh(int bc_sign, ValueVector& psi, GradVector& dpsi, HessVector& grad_grad_psi, int first, int last)
      const;

  void evaluateVGH(const ParticleSet& P,
                   const int iat,
                   ValueVector& psi,
                   GradVector& dpsi,
                   HessVector& grad_grad_psi) override;

  void assign_vghgh(int bc_sign,
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

extern template class SplineR2R<float>;
extern template class SplineR2R<double>;

} // namespace qmcplusplus
#endif
