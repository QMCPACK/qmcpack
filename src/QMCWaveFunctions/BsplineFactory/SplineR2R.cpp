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


#include "Message/OpenMP.h"
#include "SplineR2R.h"
#include "spline2/MultiBsplineEval.hpp"
#include "QMCWaveFunctions/BsplineFactory/contraction_helper.hpp"
#include "Numerics/MatrixOperators.h"
#include "CPU/BLAS.hpp"

namespace qmcplusplus
{

  /*
    @TODO (JPT) 20.09.2021
    Set the number of orbitals to include in the rotations
   */
  template<typename ST>
  inline void SplineR2R<ST>::setOrbitalSetSize( int norbs )
  {
    Identity       = false;
    m_ntot_orbs    = norbs;
    myV.resize(m_ntot_orbs);
    myL.resize(m_ntot_orbs);
    this->checkObject();
  }

  template <typename ST>
  void SplineR2R<ST>::checkObject() const
  {
    // If no rotations, ensure the number of electrons = number of orbitals
    if (Identity)
      {
	if ( m_ntot_orbs != m_nocc_orbs )
	  throw std::runtime_error("SplineR2R::checkObject() Norbs and Nelec must be equal if Identity = true!");
      }
    else
      {
	std::cerr << "SplineR2R::checkObject() case Identity = false is not implemented!\n";
      }
  }
  // End function

  /* 
     The implementation here follow that described in [1,2]
     in each "rotated" orbital, \psi, is composed of a linear combination
     of SPO's \phi:

     \psi_{i} = \psi_i + \sum_{j} \alpha_{i,j}\phi_{j}

     where \alpha is an adjustable parameter, and i runs over the core 
     (occupied) SPOs, and j runs over all virtual (unoccupied) SPOs.
     See [2] for an explanation on why core-core and virtual-virtual rotations
     are not included. This is also what was implemented in the Berkeley branch
     orbopt scheme. Cf. MultiSlateDeterminantFast.cpp lines 154 - 170. 
     
     Note that this implementation requires that NSPOs > Nelec. In other words,
     if you want to run a orbopt wfn, you must include some virtual orbitals!

     Some results from this implementation (Berkeley branch) were published 
     in [3].

     [1] Filippi & Fahy, JCP 112, (2000)
     [2] Toulouse & Umrigar, JCP 126, (2007)
     [3] Townsend et al., PRB 102, (2020)
  */
  template <typename ST>
  void SplineR2R<ST>::applyRotation(const ValueMatrix_t& rot_mat, bool use_stored_copy)
{
  std::cerr << "Entered SplineR2R::applyRotation()...\n";

  // Get the underlying spline coefs via pointer.
  // SplineInst is a MultiBspline. See src/spline2/MultiBspline.hpp
  auto spline_ptr = SplineInst->getSplinePtr();
  if ( spline_ptr == nullptr )
    {
      APP_ABORT("ERROR in SplineR2R::applyRotation()! spline_ptr is nullptr.");
    }

  /* @DEBUG (JPT) 12.10.2021 
     ~~ Spline stuff ~~
     num_splines = Total number of SPOs per spin channel (determinant)
     coefs_size  = Size of the coefs = N_x_coefs * N_y_coefs * N_z_coefs * N_orbs
     spl_coefs   = Raw pointer to a hunk o' memory
     
     Because coefs is a raw pointer, we have to take care to index into it
     with some care. Supposedly, the value of the nth SPO at a point xi, yi, zi
     is given by:

     coefs[offset] = x_stride*xi + y_stride*yi + z_stride*zi + n

     NB: The layout includes padding of the x,y,z, and orbital indices. 
     In the case of the spatial indices, the padding depends on the 
     boundary conditions. For 3d PBC's, the padding seems to be 3 
     (Cf. MultiBspline.hpp). I.e. for a 100x100x100 mesh, the "size" of 
     the x-,y-, and z-dimensions of the coefs is actually 103x103x103. 
     The orbital index is padded in order to achieve 16byte? alignment.

     You can think of spl_coefs as pointing to a matrix of size
     (Nx*Ny*Nz) x Nsplines, with the spline index adjacent in memory.     
  
     const auto x_stride = spline_ptr->x_stride;
     const auto y_stride = spline_ptr->y_stride;
     const auto z_stride = spline_ptr->z_stride;
     
     const auto nx_nopad = spline_ptr->x_grid.num;
     const auto ny_nopad = spline_ptr->y_grid.num;
     const auto nz_nopad = spline_ptr->z_grid.num;
     
     // Padding in x,y,z directions
     // Cf. multi_bspline_create.c lines 290-292.
     const auto nz_pad = y_stride / z_stride;
     const auto ny_pad = x_stride / y_stride;
     const auto nx_pad = coefs_tot_size / ny_pad / nz_pad / z_stride;
  */

  // Figure out the size of the matrix of spline coefficients to rotate
  const auto spl_coefs      = spline_ptr->coefs;
  const auto num_splines    = spline_ptr->num_splines;
  const auto coefs_tot_size = spline_ptr->coefs_size;
  const auto OrbitalSetSize = coefs_tot_size / num_splines;
  const auto BasisSetSize   = num_splines;

  // NOw there is a problem: The spline index above is >= the true
  // number of orbitals because of padding. So rot_mat is too small!
  // Easy fix: Make a new matrix, fill it with rot_mat, and zero otherwise.
  
  
  std::cerr << "Applying the rotation...\n";
  std::cerr << "  Size of rot_mat  = " << rot_mat.size1() << " x " << rot_mat.size2() << "\n";
  std::cerr << "  Size of spl_coefs= " << OrbitalSetSize  << " x " << BasisSetSize << "\n";
  
  /* 
     Need to test:
     1.) Are orbitals adjacent in memory? I.e. is the padding at the end?
     2.) Can we squeeze by with a narrowing conversion on OrbitalSetSize?
  */
  
  //BLAS::gemm('N', 'T', BasisSetSize, OrbitalSetSize, OrbitalSetSize, RealType(1.0), spl_coefs, BasisSetSize, rot_mat.data(), OrbitalSetSize, RealType(0.0), spl_coefs, BasisSetSize);

  std::cerr << "Exited SplineR2R::applyRotation()...\n";
}



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
inline void SplineR2R<ST>::assign_v(int bc_sign, const vContainer_type& myV, ValueVector_t& psi, int first, int last)
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
void SplineR2R<ST>::evaluateValue(const ParticleSet& P, const int iat, ValueVector_t& psi)
{
  const PointType& r = P.activeR(iat);
  PointType ru;
  int bc_sign = convertPos(r, ru);

#pragma omp parallel
  {
    int first, last;
    FairDivideAligned(myV.size(), getAlignment<ST>(), omp_get_num_threads(), omp_get_thread_num(), first, last);

    spline2::evaluate3d(SplineInst->getSplinePtr(), ru, myV, first, last);
    assign_v(bc_sign, myV, psi, first, last);
  }
}

template<typename ST>
void SplineR2R<ST>::evaluateDetRatios(const VirtualParticleSet& VP,
                                      ValueVector_t& psi,
                                      const ValueVector_t& psiinv,
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
    FairDivideAligned(myV.size(), getAlignment<ST>(), omp_get_num_threads(), tid, first, last);
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
                                      ValueVector_t& psi,
                                      GradVector_t& dpsi,
                                      ValueVector_t& d2psi,
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
inline void SplineR2R<ST>::assign_vgl_from_l(int bc_sign, ValueVector_t& psi, GradVector_t& dpsi, ValueVector_t& d2psi)
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
                                ValueVector_t& psi,
                                GradVector_t& dpsi,
                                ValueVector_t& d2psi)
{
  const PointType& r = P.activeR(iat);
  PointType ru;
  int bc_sign = convertPos(r, ru);

#pragma omp parallel
  {
    int first, last;
    FairDivideAligned(myV.size(), getAlignment<ST>(), omp_get_num_threads(), omp_get_thread_num(), first, last);

    spline2::evaluate3d_vgh(SplineInst->getSplinePtr(), ru, myV, myG, myH, first, last);
    assign_vgl(bc_sign, psi, dpsi, d2psi, first, last);
  }
}

template<typename ST>
void SplineR2R<ST>::assign_vgh(int bc_sign,
                               ValueVector_t& psi,
                               GradVector_t& dpsi,
                               HessVector_t& grad_grad_psi,
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
                                ValueVector_t& psi,
                                GradVector_t& dpsi,
                                HessVector_t& grad_grad_psi)
{
  const PointType& r = P.activeR(iat);
  PointType ru;
  int bc_sign = convertPos(r, ru);

#pragma omp parallel
  {
    int first, last;
    FairDivideAligned(myV.size(), getAlignment<ST>(), omp_get_num_threads(), omp_get_thread_num(), first, last);

    spline2::evaluate3d_vgh(SplineInst->getSplinePtr(), ru, myV, myG, myH, first, last);
    assign_vgh(bc_sign, psi, dpsi, grad_grad_psi, first, last);
  }
}

template<typename ST>
void SplineR2R<ST>::assign_vghgh(int bc_sign,
                                 ValueVector_t& psi,
                                 GradVector_t& dpsi,
                                 HessVector_t& grad_grad_psi,
                                 GGGVector_t& grad_grad_grad_psi,
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
                                  ValueVector_t& psi,
                                  GradVector_t& dpsi,
                                  HessVector_t& grad_grad_psi,
                                  GGGVector_t& grad_grad_grad_psi)
{
  const PointType& r = P.activeR(iat);
  PointType ru;
  int bc_sign = convertPos(r, ru);

#pragma omp parallel
  {
    int first, last;
    FairDivideAligned(myV.size(), getAlignment<ST>(), omp_get_num_threads(), omp_get_thread_num(), first, last);

    spline2::evaluate3d_vghgh(SplineInst->getSplinePtr(), ru, myV, myG, myH, mygH, first, last);
    assign_vghgh(bc_sign, psi, dpsi, grad_grad_psi, grad_grad_grad_psi, first, last);
  }
}

template class SplineR2R<float>;
template class SplineR2R<double>;

} // namespace qmcplusplus
