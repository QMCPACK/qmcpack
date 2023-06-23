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


#include "Concurrency/OpenMP.h"
#include "SplineR2R.h"
#include "spline2/MultiBsplineEval.hpp"
#include "QMCWaveFunctions/BsplineFactory/contraction_helper.hpp"

namespace qmcplusplus
{
template<typename ST>
SplineR2R<ST>::SplineR2R(const SplineR2R& in) = default;

template<typename ST>
inline void SplineR2R<ST>::set_spline(SingleSplineType* spline_r,
                                      SingleSplineType* spline_i,
                                      int twist,
                                      int ispline,
                                      int level)
{
  SplineInst->copy_spline(spline_r, ispline);
}

template<typename ST>
bool SplineR2R<ST>::read_splines(hdf_archive& h5f)
{
  std::ostringstream o;
  o << "spline_" << MyIndex;
  einspline_engine<SplineType> bigtable(SplineInst->getSplinePtr());
  return h5f.readEntry(bigtable, o.str().c_str()); //"spline_0");
}

template<typename ST>
bool SplineR2R<ST>::write_splines(hdf_archive& h5f)
{
  std::ostringstream o;
  o << "spline_" << MyIndex;
  einspline_engine<SplineType> bigtable(SplineInst->getSplinePtr());
  return h5f.writeEntry(bigtable, o.str().c_str()); //"spline_0");
}

template<typename ST>
void SplineR2R<ST>::storeParamsBeforeRotation()
{
  const auto spline_ptr     = SplineInst->getSplinePtr();
  const auto coefs_tot_size = spline_ptr->coefs_size;
  coef_copy_                = std::make_shared<std::vector<RealType>>(coefs_tot_size);

  std::copy_n(spline_ptr->coefs, coefs_tot_size, coef_copy_->begin());
}

/*
  ~~ Notes for rotation ~~
  spl_coefs      = Raw pointer to spline coefficients
  basis_set_size = Number of spline coefs per orbital
  OrbitalSetSize = Number of orbitals (excluding padding)

  spl_coefs has a complicated layout depending on dimensionality of splines.
  Luckily, for our purposes, we can think of spl_coefs as pointing to a
  matrix of size BasisSetSize x (OrbitalSetSize + padding), with the spline
  index adjacent in memory. The orbital index is SIMD aligned and therefore
  may include padding.

  As a result, due to SIMD alignment, Nsplines may be larger than the
  actual number of splined orbitals. This means that in practice rot_mat
  may be smaller than the number of 'columns' in the coefs array!

      SplineR2R spl_coef layout:
             ^         | sp1 | ... | spN | pad |
             |         |=====|=====|=====|=====|
             |         | c11 | ... | c1N | 0   |
      basis_set_size   | c21 | ... | c2N | 0   |
             |         | ... | ... | ... | 0   |
             |         | cM1 | ... | cMN | 0   |
             v         |=====|=====|=====|=====|
                       <------ Nsplines ------>

      SplineC2C spl_coef layout:
             ^         | sp1_r | sp1_i |  ...  | spN_r | spN_i |  pad  |
             |         |=======|=======|=======|=======|=======|=======|
             |         | c11_r | c11_i |  ...  | c1N_r | c1N_i |   0   |
      basis_set_size   | c21_r | c21_i |  ...  | c2N_r | c2N_i |   0   |
             |         |  ...  |  ...  |  ...  |  ...  |  ...  |  ...  |
             |         | cM1_r | cM1_i |  ...  | cMN_r | cMN_i |   0   |
             v         |=======|=======|=======|=======|=======|=======|
                       <------------------ Nsplines ------------------>

  NB: For splines (typically) BasisSetSize >> OrbitalSetSize, so the spl_coefs
  "matrix" is very tall and skinny.
*/
template<typename ST>
void SplineR2R<ST>::applyRotation(const ValueMatrix& rot_mat, bool use_stored_copy)
{
  // SplineInst is a MultiBspline. See src/spline2/MultiBspline.hpp
  const auto spline_ptr = SplineInst->getSplinePtr();
  assert(spline_ptr != nullptr);
  const auto spl_coefs      = spline_ptr->coefs;
  const auto Nsplines       = spline_ptr->num_splines; // May include padding
  const auto coefs_tot_size = spline_ptr->coefs_size;
  const auto BasisSetSize   = coefs_tot_size / Nsplines;
  const auto TrueNOrbs      = rot_mat.size1(); // == Nsplines - padding
  assert(OrbitalSetSize == rot_mat.rows());
  assert(OrbitalSetSize == rot_mat.cols());

  if (!use_stored_copy)
  {
    assert(coef_copy_ != nullptr);
    std::copy_n(spl_coefs, coefs_tot_size, coef_copy_->begin());
  }

  // Apply rotation the dumb way b/c I can't get BLAS::gemm to work...
  for (auto i = 0; i < BasisSetSize; i++)
  {
    for (auto j = 0; j < OrbitalSetSize; j++)
    {
      const auto cur_elem = Nsplines * i + j;
      auto newval{0.};
      for (auto k = 0; k < OrbitalSetSize; k++)
      {
        const auto index = i * Nsplines + k;
        newval += (*coef_copy_)[index] * rot_mat[k][j];
      }
      spl_coefs[cur_elem] = newval;
    }
  }
}


template<typename ST>
inline void SplineR2R<ST>::assign_v(int bc_sign, const vContainer_type& myV, ValueVector& psi, int first, int last)
    const
{
  // protect last
  last = last > kPoints.size() ? kPoints.size() : last;

  const ST signed_one = (bc_sign & 1) ? -1 : 1;
#pragma omp simd
  for (size_t j = first; j < last; ++j)
    psi[first_spo + j] = signed_one * myV[j];
}

template<typename ST>
void SplineR2R<ST>::evaluateValue(const ParticleSet& P, const int iat, ValueVector& psi)
{
  const PointType& r = P.activeR(iat);
  PointType ru;
  int bc_sign = convertPos(r, ru);

#pragma omp parallel
  {
    int first, last;
    FairDivideAligned(psi.size(), getAlignment<ST>(), omp_get_num_threads(), omp_get_thread_num(), first, last);

    spline2::evaluate3d(SplineInst->getSplinePtr(), ru, myV, first, last);
    assign_v(bc_sign, myV, psi, first, last);
  }
}

template<typename ST>
void SplineR2R<ST>::evaluateDetRatios(const VirtualParticleSet& VP,
                                      ValueVector& psi,
                                      const ValueVector& psiinv,
                                      std::vector<TT>& ratios)
{
  const bool need_resize = ratios_private.rows() < VP.getTotalNum();

#pragma omp parallel
  {
    int tid = omp_get_thread_num();
    // initialize thread private ratios
    if (need_resize)
    {
      if (tid == 0) // just like #pragma omp master, but one fewer call to the runtime
        ratios_private.resize(VP.getTotalNum(), omp_get_num_threads());
#pragma omp barrier
    }
    int first, last;
    FairDivideAligned(psi.size(), getAlignment<ST>(), omp_get_num_threads(), tid, first, last);
    const int last_real = kPoints.size() < last ? kPoints.size() : last;

    for (int iat = 0; iat < VP.getTotalNum(); ++iat)
    {
      const PointType& r = VP.activeR(iat);
      PointType ru;
      int bc_sign = convertPos(r, ru);

      spline2::evaluate3d(SplineInst->getSplinePtr(), ru, myV, first, last);
      assign_v(bc_sign, myV, psi, first, last_real);
      ratios_private[iat][tid] = simd::dot(psi.data() + first, psiinv.data() + first, last_real - first);
    }
  }

  // do the reduction manually
  for (int iat = 0; iat < VP.getTotalNum(); ++iat)
  {
    ratios[iat] = TT(0);
    for (int tid = 0; tid < ratios_private.cols(); tid++)
      ratios[iat] += ratios_private[iat][tid];
  }
}

template<typename ST>
inline void SplineR2R<ST>::assign_vgl(int bc_sign,
                                      ValueVector& psi,
                                      GradVector& dpsi,
                                      ValueVector& d2psi,
                                      int first,
                                      int last) const
{
  // protect last
  last = last > kPoints.size() ? kPoints.size() : last;

  const ST signed_one = (bc_sign & 1) ? -1 : 1;
  const ST g00 = PrimLattice.G(0), g01 = PrimLattice.G(1), g02 = PrimLattice.G(2), g10 = PrimLattice.G(3),
           g11 = PrimLattice.G(4), g12 = PrimLattice.G(5), g20 = PrimLattice.G(6), g21 = PrimLattice.G(7),
           g22      = PrimLattice.G(8);
  const ST symGG[6] = {GGt[0], GGt[1] + GGt[3], GGt[2] + GGt[6], GGt[4], GGt[5] + GGt[7], GGt[8]};

  const ST* restrict g0  = myG.data(0);
  const ST* restrict g1  = myG.data(1);
  const ST* restrict g2  = myG.data(2);
  const ST* restrict h00 = myH.data(0);
  const ST* restrict h01 = myH.data(1);
  const ST* restrict h02 = myH.data(2);
  const ST* restrict h11 = myH.data(3);
  const ST* restrict h12 = myH.data(4);
  const ST* restrict h22 = myH.data(5);

#pragma omp simd
  for (size_t j = first; j < last; ++j)
  {
    const size_t psiIndex = first_spo + j;
    psi[psiIndex]         = signed_one * myV[j];
    dpsi[psiIndex][0]     = signed_one * (g00 * g0[j] + g01 * g1[j] + g02 * g2[j]);
    dpsi[psiIndex][1]     = signed_one * (g10 * g0[j] + g11 * g1[j] + g12 * g2[j]);
    dpsi[psiIndex][2]     = signed_one * (g20 * g0[j] + g21 * g1[j] + g22 * g2[j]);
    d2psi[psiIndex]       = signed_one * SymTrace(h00[j], h01[j], h02[j], h11[j], h12[j], h22[j], symGG);
  }
}

/** assign_vgl_from_l can be used when myL is precomputed and myV,myG,myL in cartesian
   */
template<typename ST>
inline void SplineR2R<ST>::assign_vgl_from_l(int bc_sign, ValueVector& psi, GradVector& dpsi, ValueVector& d2psi)
{
  const ST signed_one   = (bc_sign & 1) ? -1 : 1;
  const ST* restrict g0 = myG.data(0);
  const ST* restrict g1 = myG.data(1);
  const ST* restrict g2 = myG.data(2);

#pragma omp simd
  for (int psiIndex = first_spo; psiIndex < last_spo; ++psiIndex)
  {
    const size_t j    = psiIndex - first_spo;
    psi[psiIndex]     = signed_one * myV[j];
    dpsi[psiIndex][0] = signed_one * g0[j];
    dpsi[psiIndex][1] = signed_one * g1[j];
    dpsi[psiIndex][2] = signed_one * g2[j];
    d2psi[psiIndex]   = signed_one * myL[j];
  }
}

template<typename ST>
void SplineR2R<ST>::evaluateVGL(const ParticleSet& P,
                                const int iat,
                                ValueVector& psi,
                                GradVector& dpsi,
                                ValueVector& d2psi)
{
  const PointType& r = P.activeR(iat);
  PointType ru;
  int bc_sign = convertPos(r, ru);

#pragma omp parallel
  {
    int first, last;
    FairDivideAligned(psi.size(), getAlignment<ST>(), omp_get_num_threads(), omp_get_thread_num(), first, last);

    spline2::evaluate3d_vgh(SplineInst->getSplinePtr(), ru, myV, myG, myH, first, last);
    assign_vgl(bc_sign, psi, dpsi, d2psi, first, last);
  }
}

template<typename ST>
void SplineR2R<ST>::assign_vgh(int bc_sign,
                               ValueVector& psi,
                               GradVector& dpsi,
                               HessVector& grad_grad_psi,
                               int first,
                               int last) const
{
  // protect last
  last = last > kPoints.size() ? kPoints.size() : last;

  const ST signed_one = (bc_sign & 1) ? -1 : 1;
  const ST g00 = PrimLattice.G(0), g01 = PrimLattice.G(1), g02 = PrimLattice.G(2), g10 = PrimLattice.G(3),
           g11 = PrimLattice.G(4), g12 = PrimLattice.G(5), g20 = PrimLattice.G(6), g21 = PrimLattice.G(7),
           g22 = PrimLattice.G(8);

  const ST* restrict g0  = myG.data(0);
  const ST* restrict g1  = myG.data(1);
  const ST* restrict g2  = myG.data(2);
  const ST* restrict h00 = myH.data(0);
  const ST* restrict h01 = myH.data(1);
  const ST* restrict h02 = myH.data(2);
  const ST* restrict h11 = myH.data(3);
  const ST* restrict h12 = myH.data(4);
  const ST* restrict h22 = myH.data(5);

#pragma omp simd
  for (size_t j = first; j < last; ++j)
  {
    //dot(PrimLattice.G,myG[j])
    const ST dX_r = g00 * g0[j] + g01 * g1[j] + g02 * g2[j];
    const ST dY_r = g10 * g0[j] + g11 * g1[j] + g12 * g2[j];
    const ST dZ_r = g20 * g0[j] + g21 * g1[j] + g22 * g2[j];

    const size_t psiIndex = j + first_spo;
    psi[psiIndex]         = signed_one * myV[j];
    dpsi[psiIndex][0]     = signed_one * dX_r;
    dpsi[psiIndex][1]     = signed_one * dY_r;
    dpsi[psiIndex][2]     = signed_one * dZ_r;

    const ST h_xx_r = v_m_v(h00[j], h01[j], h02[j], h11[j], h12[j], h22[j], g00, g01, g02, g00, g01, g02);
    const ST h_xy_r = v_m_v(h00[j], h01[j], h02[j], h11[j], h12[j], h22[j], g00, g01, g02, g10, g11, g12);
    const ST h_xz_r = v_m_v(h00[j], h01[j], h02[j], h11[j], h12[j], h22[j], g00, g01, g02, g20, g21, g22);
    const ST h_yx_r = v_m_v(h00[j], h01[j], h02[j], h11[j], h12[j], h22[j], g10, g11, g12, g00, g01, g02);
    const ST h_yy_r = v_m_v(h00[j], h01[j], h02[j], h11[j], h12[j], h22[j], g10, g11, g12, g10, g11, g12);
    const ST h_yz_r = v_m_v(h00[j], h01[j], h02[j], h11[j], h12[j], h22[j], g10, g11, g12, g20, g21, g22);
    const ST h_zx_r = v_m_v(h00[j], h01[j], h02[j], h11[j], h12[j], h22[j], g20, g21, g22, g00, g01, g02);
    const ST h_zy_r = v_m_v(h00[j], h01[j], h02[j], h11[j], h12[j], h22[j], g20, g21, g22, g10, g11, g12);
    const ST h_zz_r = v_m_v(h00[j], h01[j], h02[j], h11[j], h12[j], h22[j], g20, g21, g22, g20, g21, g22);

    grad_grad_psi[psiIndex][0] = signed_one * h_xx_r;
    grad_grad_psi[psiIndex][1] = signed_one * h_xy_r;
    grad_grad_psi[psiIndex][2] = signed_one * h_xz_r;
    grad_grad_psi[psiIndex][3] = signed_one * h_yx_r;
    grad_grad_psi[psiIndex][4] = signed_one * h_yy_r;
    grad_grad_psi[psiIndex][5] = signed_one * h_yz_r;
    grad_grad_psi[psiIndex][6] = signed_one * h_zx_r;
    grad_grad_psi[psiIndex][7] = signed_one * h_zy_r;
    grad_grad_psi[psiIndex][8] = signed_one * h_zz_r;
  }
}

template<typename ST>
void SplineR2R<ST>::evaluateVGH(const ParticleSet& P,
                                const int iat,
                                ValueVector& psi,
                                GradVector& dpsi,
                                HessVector& grad_grad_psi)
{
  const PointType& r = P.activeR(iat);
  PointType ru;
  int bc_sign = convertPos(r, ru);

#pragma omp parallel
  {
    int first, last;
    FairDivideAligned(psi.size(), getAlignment<ST>(), omp_get_num_threads(), omp_get_thread_num(), first, last);

    spline2::evaluate3d_vgh(SplineInst->getSplinePtr(), ru, myV, myG, myH, first, last);
    assign_vgh(bc_sign, psi, dpsi, grad_grad_psi, first, last);
  }
}

template<typename ST>
void SplineR2R<ST>::assign_vghgh(int bc_sign,
                                 ValueVector& psi,
                                 GradVector& dpsi,
                                 HessVector& grad_grad_psi,
                                 GGGVector& grad_grad_grad_psi,
                                 int first,
                                 int last) const
{
  // protect last
  last = last < 0 ? kPoints.size() : (last > kPoints.size() ? kPoints.size() : last);

  const ST signed_one = (bc_sign & 1) ? -1 : 1;
  const ST g00 = PrimLattice.G(0), g01 = PrimLattice.G(1), g02 = PrimLattice.G(2), g10 = PrimLattice.G(3),
           g11 = PrimLattice.G(4), g12 = PrimLattice.G(5), g20 = PrimLattice.G(6), g21 = PrimLattice.G(7),
           g22 = PrimLattice.G(8);

  const ST* restrict g0  = myG.data(0);
  const ST* restrict g1  = myG.data(1);
  const ST* restrict g2  = myG.data(2);
  const ST* restrict h00 = myH.data(0);
  const ST* restrict h01 = myH.data(1);
  const ST* restrict h02 = myH.data(2);
  const ST* restrict h11 = myH.data(3);
  const ST* restrict h12 = myH.data(4);
  const ST* restrict h22 = myH.data(5);

  const ST* restrict gh000 = mygH.data(0);
  const ST* restrict gh001 = mygH.data(1);
  const ST* restrict gh002 = mygH.data(2);
  const ST* restrict gh011 = mygH.data(3);
  const ST* restrict gh012 = mygH.data(4);
  const ST* restrict gh022 = mygH.data(5);
  const ST* restrict gh111 = mygH.data(6);
  const ST* restrict gh112 = mygH.data(7);
  const ST* restrict gh122 = mygH.data(8);
  const ST* restrict gh222 = mygH.data(9);

  //SIMD doesn't work quite right yet.  Comment out until further debugging.
  //#pragma omp simd
  for (size_t j = first; j < last; ++j)
  {
    const ST val_r = myV[j];


    //dot(PrimLattice.G,myG[j])
    const ST dX_r = g00 * g0[j] + g01 * g1[j] + g02 * g2[j];
    const ST dY_r = g10 * g0[j] + g11 * g1[j] + g12 * g2[j];
    const ST dZ_r = g20 * g0[j] + g21 * g1[j] + g22 * g2[j];

    const size_t psiIndex = j + first_spo;
    psi[psiIndex]         = signed_one * val_r;
    dpsi[psiIndex][0]     = signed_one * dX_r;
    dpsi[psiIndex][1]     = signed_one * dY_r;
    dpsi[psiIndex][2]     = signed_one * dZ_r;

    //intermediates for computation of hessian. \partial_i \partial_j phi in cartesian coordinates.
    const ST f_xx_r = v_m_v(h00[j], h01[j], h02[j], h11[j], h12[j], h22[j], g00, g01, g02, g00, g01, g02);
    const ST f_xy_r = v_m_v(h00[j], h01[j], h02[j], h11[j], h12[j], h22[j], g00, g01, g02, g10, g11, g12);
    const ST f_xz_r = v_m_v(h00[j], h01[j], h02[j], h11[j], h12[j], h22[j], g00, g01, g02, g20, g21, g22);
    const ST f_yy_r = v_m_v(h00[j], h01[j], h02[j], h11[j], h12[j], h22[j], g10, g11, g12, g10, g11, g12);
    const ST f_yz_r = v_m_v(h00[j], h01[j], h02[j], h11[j], h12[j], h22[j], g10, g11, g12, g20, g21, g22);
    const ST f_zz_r = v_m_v(h00[j], h01[j], h02[j], h11[j], h12[j], h22[j], g20, g21, g22, g20, g21, g22);

    /*    const ST h_xx_r=f_xx_r;
      const ST h_xy_r=f_xy_r+(kX*dY_i+kY*dX_i)-kX*kY*val_r;
      const ST h_xz_r=f_xz_r+(kX*dZ_i+kZ*dX_i)-kX*kZ*val_r;
      const ST h_yy_r=f_yy_r+2*kY*dY_i-kY*kY*val_r;
      const ST h_yz_r=f_yz_r+(kY*dZ_i+kZ*dY_i)-kY*kZ*val_r;
      const ST h_zz_r=f_zz_r+2*kZ*dZ_i-kZ*kZ*val_r; */

    grad_grad_psi[psiIndex][0] = f_xx_r * signed_one;
    grad_grad_psi[psiIndex][1] = f_xy_r * signed_one;
    grad_grad_psi[psiIndex][2] = f_xz_r * signed_one;
    grad_grad_psi[psiIndex][4] = f_yy_r * signed_one;
    grad_grad_psi[psiIndex][5] = f_yz_r * signed_one;
    grad_grad_psi[psiIndex][8] = f_zz_r * signed_one;

    //symmetry:
    grad_grad_psi[psiIndex][3] = grad_grad_psi[psiIndex][1];
    grad_grad_psi[psiIndex][6] = grad_grad_psi[psiIndex][2];
    grad_grad_psi[psiIndex][7] = grad_grad_psi[psiIndex][5];
    //These are the real and imaginary components of the third SPO derivative.  _xxx denotes
    // third derivative w.r.t. x, _xyz, a derivative with resepect to x,y, and z, and so on.

    const ST f3_xxx_r = t3_contract(gh000[j], gh001[j], gh002[j], gh011[j], gh012[j], gh022[j], gh111[j], gh112[j],
                                    gh122[j], gh222[j], g00, g01, g02, g00, g01, g02, g00, g01, g02);
    const ST f3_xxy_r = t3_contract(gh000[j], gh001[j], gh002[j], gh011[j], gh012[j], gh022[j], gh111[j], gh112[j],
                                    gh122[j], gh222[j], g00, g01, g02, g00, g01, g02, g10, g11, g12);
    const ST f3_xxz_r = t3_contract(gh000[j], gh001[j], gh002[j], gh011[j], gh012[j], gh022[j], gh111[j], gh112[j],
                                    gh122[j], gh222[j], g00, g01, g02, g00, g01, g02, g20, g21, g22);
    const ST f3_xyy_r = t3_contract(gh000[j], gh001[j], gh002[j], gh011[j], gh012[j], gh022[j], gh111[j], gh112[j],
                                    gh122[j], gh222[j], g00, g01, g02, g10, g11, g12, g10, g11, g12);
    const ST f3_xyz_r = t3_contract(gh000[j], gh001[j], gh002[j], gh011[j], gh012[j], gh022[j], gh111[j], gh112[j],
                                    gh122[j], gh222[j], g00, g01, g02, g10, g11, g12, g20, g21, g22);
    const ST f3_xzz_r = t3_contract(gh000[j], gh001[j], gh002[j], gh011[j], gh012[j], gh022[j], gh111[j], gh112[j],
                                    gh122[j], gh222[j], g00, g01, g02, g20, g21, g22, g20, g21, g22);
    const ST f3_yyy_r = t3_contract(gh000[j], gh001[j], gh002[j], gh011[j], gh012[j], gh022[j], gh111[j], gh112[j],
                                    gh122[j], gh222[j], g10, g11, g12, g10, g11, g12, g10, g11, g12);
    const ST f3_yyz_r = t3_contract(gh000[j], gh001[j], gh002[j], gh011[j], gh012[j], gh022[j], gh111[j], gh112[j],
                                    gh122[j], gh222[j], g10, g11, g12, g10, g11, g12, g20, g21, g22);
    const ST f3_yzz_r = t3_contract(gh000[j], gh001[j], gh002[j], gh011[j], gh012[j], gh022[j], gh111[j], gh112[j],
                                    gh122[j], gh222[j], g10, g11, g12, g20, g21, g22, g20, g21, g22);
    const ST f3_zzz_r = t3_contract(gh000[j], gh001[j], gh002[j], gh011[j], gh012[j], gh022[j], gh111[j], gh112[j],
                                    gh122[j], gh222[j], g20, g21, g22, g20, g21, g22, g20, g21, g22);

    //Here is where we build up the components of the physical hessian gradient, namely, d^3/dx^3(e^{-ik*r}\phi(r)
    /*     const ST gh_xxx_r= f3_xxx_r + 3*kX*f_xx_i - 3*kX*kX*dX_r - kX*kX*kX*val_i;
      const ST gh_xxy_r= f3_xxy_r +(kY*f_xx_i+2*kX*f_xy_i) - (kX*kX*dY_r+2*kX*kY*dX_r)-kX*kX*kY*val_i; 
      const ST gh_xxz_r= f3_xxz_r +(kZ*f_xx_i+2*kX*f_xz_i) - (kX*kX*dZ_r+2*kX*kZ*dX_r)-kX*kX*kZ*val_i; 
      const ST gh_xyy_r= f3_xyy_r +(2*kY*f_xy_i+kX*f_yy_i) - (2*kX*kY*dY_r+kY*kY*dX_r)-kX*kY*kY*val_i;
      const ST gh_xyz_r= f3_xyz_r +(kX*f_yz_i+kY*f_xz_i+kZ*f_xy_i)-(kX*kY*dZ_r+kY*kZ*dX_r+kZ*kX*dY_r) - kX*kY*kZ*val_i;
      const ST gh_xzz_r= f3_xzz_r +(2*kZ*f_xz_i+kX*f_zz_i) - (2*kX*kZ*dZ_r+kZ*kZ*dX_r)-kX*kZ*kZ*val_i;
      const ST gh_yyy_r= f3_yyy_r + 3*kY*f_yy_i - 3*kY*kY*dY_r - kY*kY*kY*val_i;
      const ST gh_yyz_r= f3_yyz_r +(kZ*f_yy_i+2*kY*f_yz_i) - (kY*kY*dZ_r+2*kY*kZ*dY_r)-kY*kY*kZ*val_i; 
      const ST gh_yzz_r= f3_yzz_r +(2*kZ*f_yz_i+kY*f_zz_i) - (2*kY*kZ*dZ_r+kZ*kZ*dY_r)-kY*kZ*kZ*val_i;
      const ST gh_zzz_r= f3_zzz_r + 3*kZ*f_zz_i - 3*kZ*kZ*dZ_r - kZ*kZ*kZ*val_i;*/
    //[x][xx] //These are the unique entries
    grad_grad_grad_psi[psiIndex][0][0] = signed_one * f3_xxx_r;
    grad_grad_grad_psi[psiIndex][0][1] = signed_one * f3_xxy_r;
    grad_grad_grad_psi[psiIndex][0][2] = signed_one * f3_xxz_r;
    grad_grad_grad_psi[psiIndex][0][4] = signed_one * f3_xyy_r;
    grad_grad_grad_psi[psiIndex][0][5] = signed_one * f3_xyz_r;
    grad_grad_grad_psi[psiIndex][0][8] = signed_one * f3_xzz_r;

    //filling in the symmetric terms.  Filling out the xij terms
    grad_grad_grad_psi[psiIndex][0][3] = grad_grad_grad_psi[psiIndex][0][1];
    grad_grad_grad_psi[psiIndex][0][6] = grad_grad_grad_psi[psiIndex][0][2];
    grad_grad_grad_psi[psiIndex][0][7] = grad_grad_grad_psi[psiIndex][0][5];

    //Now for everything that's a permutation of the above:
    grad_grad_grad_psi[psiIndex][1][0] = grad_grad_grad_psi[psiIndex][0][1];
    grad_grad_grad_psi[psiIndex][1][1] = grad_grad_grad_psi[psiIndex][0][4];
    grad_grad_grad_psi[psiIndex][1][2] = grad_grad_grad_psi[psiIndex][0][5];
    grad_grad_grad_psi[psiIndex][1][3] = grad_grad_grad_psi[psiIndex][0][4];
    grad_grad_grad_psi[psiIndex][1][6] = grad_grad_grad_psi[psiIndex][0][5];

    grad_grad_grad_psi[psiIndex][2][0] = grad_grad_grad_psi[psiIndex][0][2];
    grad_grad_grad_psi[psiIndex][2][1] = grad_grad_grad_psi[psiIndex][0][5];
    grad_grad_grad_psi[psiIndex][2][2] = grad_grad_grad_psi[psiIndex][0][8];
    grad_grad_grad_psi[psiIndex][2][3] = grad_grad_grad_psi[psiIndex][0][5];
    grad_grad_grad_psi[psiIndex][2][6] = grad_grad_grad_psi[psiIndex][0][8];

    grad_grad_grad_psi[psiIndex][1][4] = signed_one * f3_yyy_r;
    grad_grad_grad_psi[psiIndex][1][5] = signed_one * f3_yyz_r;
    grad_grad_grad_psi[psiIndex][1][8] = signed_one * f3_yzz_r;

    grad_grad_grad_psi[psiIndex][1][7] = grad_grad_grad_psi[psiIndex][1][5];
    grad_grad_grad_psi[psiIndex][2][4] = grad_grad_grad_psi[psiIndex][1][5];
    grad_grad_grad_psi[psiIndex][2][5] = grad_grad_grad_psi[psiIndex][1][8];
    grad_grad_grad_psi[psiIndex][2][7] = grad_grad_grad_psi[psiIndex][1][8];

    grad_grad_grad_psi[psiIndex][2][8] = signed_one * f3_zzz_r;
  }
}

template<typename ST>
void SplineR2R<ST>::evaluateVGHGH(const ParticleSet& P,
                                  const int iat,
                                  ValueVector& psi,
                                  GradVector& dpsi,
                                  HessVector& grad_grad_psi,
                                  GGGVector& grad_grad_grad_psi)
{
  const PointType& r = P.activeR(iat);
  PointType ru;
  int bc_sign = convertPos(r, ru);

#pragma omp parallel
  {
    int first, last;
    FairDivideAligned(psi.size(), getAlignment<ST>(), omp_get_num_threads(), omp_get_thread_num(), first, last);

    spline2::evaluate3d_vghgh(SplineInst->getSplinePtr(), ru, myV, myG, myH, mygH, first, last);
    assign_vghgh(bc_sign, psi, dpsi, grad_grad_psi, grad_grad_grad_psi, first, last);
  }
}

template class SplineR2R<float>;
template class SplineR2R<double>;

} // namespace qmcplusplus
