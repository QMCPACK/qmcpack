//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by:
//
// File created by: Jeongnim Kim, jeongnim.kim@intel.com, Intel Corp.
//////////////////////////////////////////////////////////////////////////////////////

/** @file SoaAtomicBasisSet.h
 */
#ifndef QMCPLUSPLUS_SOA_SPHERICALORBITAL_BASISSET_H
#define QMCPLUSPLUS_SOA_SPHERICALORBITAL_BASISSET_H

#include "CPU/math.hpp"
#include "OptimizableObject.h"
#include <ResourceCollection.h>

namespace qmcplusplus
{
/* A basis set for a center type 
   *
   * @tparam ROT : radial function type, e.g.,NGFunctor<T>
   * @tparam SH : spherical or carteisan Harmonics for (l,m) expansion
   *
   * \f$ \phi_{n,l,m}({\bf r})=R_{n,l}(r) Y_{l,m}(\theta) \f$
   */
template<typename ROT, typename SH>
class SoaAtomicBasisSet
{
public:
  using RadialOrbital_t  = ROT;
  using RealType         = typename ROT::RealType;
  using GridType         = typename ROT::GridType;
  using ValueType        = typename QMCTraits::ValueType;
  using OffloadArray4D   = Array<RealType, 4, OffloadPinnedAllocator<RealType>>;
  using OffloadArray3D   = Array<RealType, 3, OffloadPinnedAllocator<RealType>>;
  using OffloadArray2D   = Array<RealType, 2, OffloadPinnedAllocator<RealType>>;
  using OffloadVector    = Vector<ValueType, OffloadPinnedAllocator<ValueType>>;
  using OffloadIntVector = Vector<int, OffloadPinnedAllocator<int>>;

  ///the constructor
  explicit SoaAtomicBasisSet(int lmax, bool addsignforM = false)
      : Ylm(lmax, addsignforM),
        periodic_image_phase_factors_ptr_(std::make_shared<OffloadVector>()),
        periodic_image_displacements_ptr_(std::make_shared<OffloadArray2D>()),
        periodic_image_phase_factors_(*periodic_image_phase_factors_ptr_),
        periodic_image_displacements_(*periodic_image_displacements_ptr_),
        NL_ptr_(std::make_shared<OffloadIntVector>()),
        LM_ptr_(std::make_shared<OffloadIntVector>()),
        NL(*NL_ptr_),
        LM(*LM_ptr_),
        ylm_timer_(createGlobalTimer("SoaAtomicBasisSet::Ylm", timer_level_fine)),
        rnl_timer_(createGlobalTimer("SoaAtomicBasisSet::Rnl", timer_level_fine)),
        pbc_timer_(createGlobalTimer("SoaAtomicBasisSet::pbc_images", timer_level_fine)),
        nelec_pbc_timer_(createGlobalTimer("SoaAtomicBasisSet::nelec_pbc_images", timer_level_fine)),
        phase_timer_(createGlobalTimer("SoaAtomicBasisSet::phase", timer_level_fine)),
        psi_timer_(createGlobalTimer("SoaAtomicBasisSet::psi", timer_level_fine))
  {}

  void checkInVariables(opt_variables_type& active)
  {
    //for(size_t nl=0; nl<Rnl.size(); nl++)
    //  Rnl[nl]->checkInVariables(active);
  }

  void checkOutVariables(const opt_variables_type& active)
  {
    //for(size_t nl=0; nl<Rnl.size(); nl++)
    //  Rnl[nl]->checkOutVariables(active);
  }

  void resetParameters(const opt_variables_type& active)
  {
    //for(size_t nl=0; nl<Rnl.size(); nl++)
    //  Rnl[nl]->resetParameters(active);
  }

  /** return the number of basis functions
      */
  inline int getBasisSetSize() const
  {
    //=NL.size();
    return BasisSetSize;
  }

  /** Set the number of periodic image for the evaluation of the orbitals and the phase factor.
   * In the case of Non-PBC, PBCImages=(1,1,1), SuperTwist(0,0,0) and the PhaseFactor=1.
   */
  void setPBCParams(const TinyVector<int, 3>& pbc_images,
                    const TinyVector<double, 3> supertwist,
                    const OffloadVector& PeriodicImagePhaseFactors,
                    const OffloadArray2D& PeriodicImageDisplacements)
  {
    PBCImages                     = pbc_images;
    periodic_image_phase_factors_ = PeriodicImagePhaseFactors;
    periodic_image_displacements_ = PeriodicImageDisplacements;
    SuperTwist                    = supertwist;

    periodic_image_phase_factors_.updateTo();
    periodic_image_displacements_.updateTo();
  }


  /** implement a BasisSetBase virtual function
   *
   * Set Rmax and BasisSetSize
   * @todo Should be able to overwrite Rmax to be much smaller than the maximum grid
   */
  inline void finalize()
  {
    BasisSetSize = LM.size();
    NL.updateTo();
    LM.updateTo();
    tempS.resize(std::max(Ylm.size(), RnlID.size()));
  }

  /** Set Rmax */
  template<typename T>
  inline void setRmax(T rmax)
  {
    Rmax = (rmax > 0) ? rmax : MultiRnl.rmax();
  }

  ///set the current offset
  inline void setCenter(int c, int offset) {}

  /// Sets a boolean vector for S-type orbitals.  Used for cusp correction.
  void queryOrbitalsForSType(std::vector<bool>& s_orbitals) const
  {
    for (int i = 0; i < BasisSetSize; i++)
    {
      s_orbitals[i] = (RnlID[NL[i]][1] == 0);
    }
  }

  /** evaluate VGL
   */
  template<typename LAT, typename T, typename PosType, typename VGL>
  inline void evaluateVGL(const LAT& lattice, const T r, const PosType& dr, const size_t offset, VGL& vgl, PosType Tv)
  {
    int TransX, TransY, TransZ;

    PosType dr_new;
    T r_new;
    // T psi_new, dpsi_x_new, dpsi_y_new, dpsi_z_new,d2psi_new;

#if not defined(QMC_COMPLEX)
    const ValueType correctphase = 1;
#else

    RealType phasearg = SuperTwist[0] * Tv[0] + SuperTwist[1] * Tv[1] + SuperTwist[2] * Tv[2];
    RealType s, c;
    qmcplusplus::sincos(-phasearg, &s, &c);
    const ValueType correctphase(c, s);
#endif

    constexpr T cone(1);
    constexpr T ctwo(2);


    //one can assert the alignment
    RealType* restrict phi   = tempS.data(0);
    RealType* restrict dphi  = tempS.data(1);
    RealType* restrict d2phi = tempS.data(2);

    //V,Gx,Gy,Gz,L
    auto* restrict psi      = vgl.data(0) + offset;
    const T* restrict ylm_v = Ylm[0]; //value
    auto* restrict dpsi_x   = vgl.data(1) + offset;
    const T* restrict ylm_x = Ylm[1]; //gradX
    auto* restrict dpsi_y   = vgl.data(2) + offset;
    const T* restrict ylm_y = Ylm[2]; //gradY
    auto* restrict dpsi_z   = vgl.data(3) + offset;
    const T* restrict ylm_z = Ylm[3]; //gradZ
    auto* restrict d2psi    = vgl.data(4) + offset;
    const T* restrict ylm_l = Ylm[4]; //lap

    for (size_t ib = 0; ib < BasisSetSize; ++ib)
    {
      psi[ib]    = 0;
      dpsi_x[ib] = 0;
      dpsi_y[ib] = 0;
      dpsi_z[ib] = 0;
      d2psi[ib]  = 0;
    }
    //Phase_idx (iter) needs to be initialized at -1 as it has to be incremented first to comply with the if statement (r_new >=Rmax)
    int iter = -1;
    for (int i = 0; i <= PBCImages[0]; i++) //loop Translation over X
    {
      //Allows to increment cells from 0,1,-1,2,-2,3,-3 etc...
      TransX = ((i % 2) * 2 - 1) * ((i + 1) / 2);
      for (int j = 0; j <= PBCImages[1]; j++) //loop Translation over Y
      {
        //Allows to increment cells from 0,1,-1,2,-2,3,-3 etc...
        TransY = ((j % 2) * 2 - 1) * ((j + 1) / 2);
        for (int k = 0; k <= PBCImages[2]; k++) //loop Translation over Z
        {
          //Allows to increment cells from 0,1,-1,2,-2,3,-3 etc...
          TransZ = ((k % 2) * 2 - 1) * ((k + 1) / 2);

          dr_new[0] = dr[0] + (TransX * lattice.R(0, 0) + TransY * lattice.R(1, 0) + TransZ * lattice.R(2, 0));
          dr_new[1] = dr[1] + (TransX * lattice.R(0, 1) + TransY * lattice.R(1, 1) + TransZ * lattice.R(2, 1));
          dr_new[2] = dr[2] + (TransX * lattice.R(0, 2) + TransY * lattice.R(1, 2) + TransZ * lattice.R(2, 2));

          r_new = std::sqrt(dot(dr_new, dr_new));

          iter++;
          if (r_new >= Rmax)
            continue;

          //SIGN Change!!
          const T x = -dr_new[0], y = -dr_new[1], z = -dr_new[2];
          Ylm.evaluateVGL(x, y, z);

          MultiRnl.evaluate(r_new, phi, dphi, d2phi);

          const T rinv = cone / r_new;

          ///Phase for PBC containing the phase for the nearest image displacement and the correction due to the Distance table.
          const ValueType Phase = periodic_image_phase_factors_[iter] * correctphase;

          for (size_t ib = 0; ib < BasisSetSize; ++ib)
          {
            const int nl(NL[ib]);
            const int lm(LM[ib]);
            const T drnloverr = rinv * dphi[nl];
            const T ang       = ylm_v[lm];
            const T gr_x      = drnloverr * x;
            const T gr_y      = drnloverr * y;
            const T gr_z      = drnloverr * z;
            const T ang_x     = ylm_x[lm];
            const T ang_y     = ylm_y[lm];
            const T ang_z     = ylm_z[lm];
            const T vr        = phi[nl];

            psi[ib] += ang * vr * Phase;
            dpsi_x[ib] += (ang * gr_x + vr * ang_x) * Phase;
            dpsi_y[ib] += (ang * gr_y + vr * ang_y) * Phase;
            dpsi_z[ib] += (ang * gr_z + vr * ang_z) * Phase;
            d2psi[ib] += (ang * (ctwo * drnloverr + d2phi[nl]) + ctwo * (gr_x * ang_x + gr_y * ang_y + gr_z * ang_z) +
                          vr * ylm_l[lm]) *
                Phase;
          }
        }
      }
    }
  }

  template<typename LAT, typename T, typename PosType, typename VGH>
  inline void evaluateVGH(const LAT& lattice, const T r, const PosType& dr, const size_t offset, VGH& vgh)
  {
    int TransX, TransY, TransZ;

    PosType dr_new;
    T r_new;

    constexpr T cone(1);

    //one can assert the alignment
    RealType* restrict phi   = tempS.data(0);
    RealType* restrict dphi  = tempS.data(1);
    RealType* restrict d2phi = tempS.data(2);

    //V,Gx,Gy,Gz,L
    auto* restrict psi      = vgh.data(0) + offset;
    const T* restrict ylm_v = Ylm[0]; //value
    auto* restrict dpsi_x   = vgh.data(1) + offset;
    const T* restrict ylm_x = Ylm[1]; //gradX
    auto* restrict dpsi_y   = vgh.data(2) + offset;
    const T* restrict ylm_y = Ylm[2]; //gradY
    auto* restrict dpsi_z   = vgh.data(3) + offset;
    const T* restrict ylm_z = Ylm[3]; //gradZ

    auto* restrict dhpsi_xx  = vgh.data(4) + offset;
    const T* restrict ylm_xx = Ylm[4];
    auto* restrict dhpsi_xy  = vgh.data(5) + offset;
    const T* restrict ylm_xy = Ylm[5];
    auto* restrict dhpsi_xz  = vgh.data(6) + offset;
    const T* restrict ylm_xz = Ylm[6];
    auto* restrict dhpsi_yy  = vgh.data(7) + offset;
    const T* restrict ylm_yy = Ylm[7];
    auto* restrict dhpsi_yz  = vgh.data(8) + offset;
    const T* restrict ylm_yz = Ylm[8];
    auto* restrict dhpsi_zz  = vgh.data(9) + offset;
    const T* restrict ylm_zz = Ylm[9];

    for (size_t ib = 0; ib < BasisSetSize; ++ib)
    {
      psi[ib]      = 0;
      dpsi_x[ib]   = 0;
      dpsi_y[ib]   = 0;
      dpsi_z[ib]   = 0;
      dhpsi_xx[ib] = 0;
      dhpsi_xy[ib] = 0;
      dhpsi_xz[ib] = 0;
      dhpsi_yy[ib] = 0;
      dhpsi_yz[ib] = 0;
      dhpsi_zz[ib] = 0;
      //      d2psi[ib]  = 0;
    }

    for (int i = 0; i <= PBCImages[0]; i++) //loop Translation over X
    {
      //Allows to increment cells from 0,1,-1,2,-2,3,-3 etc...
      TransX = ((i % 2) * 2 - 1) * ((i + 1) / 2);
      for (int j = 0; j <= PBCImages[1]; j++) //loop Translation over Y
      {
        //Allows to increment cells from 0,1,-1,2,-2,3,-3 etc...
        TransY = ((j % 2) * 2 - 1) * ((j + 1) / 2);
        for (int k = 0; k <= PBCImages[2]; k++) //loop Translation over Z
        {
          //Allows to increment cells from 0,1,-1,2,-2,3,-3 etc...
          TransZ    = ((k % 2) * 2 - 1) * ((k + 1) / 2);
          dr_new[0] = dr[0] + TransX * lattice.R(0, 0) + TransY * lattice.R(1, 0) + TransZ * lattice.R(2, 0);
          dr_new[1] = dr[1] + TransX * lattice.R(0, 1) + TransY * lattice.R(1, 1) + TransZ * lattice.R(2, 1);
          dr_new[2] = dr[2] + TransX * lattice.R(0, 2) + TransY * lattice.R(1, 2) + TransZ * lattice.R(2, 2);
          r_new     = std::sqrt(dot(dr_new, dr_new));

          //const size_t ib_max=NL.size();
          if (r_new >= Rmax)
            continue;

          //SIGN Change!!
          const T x = -dr_new[0], y = -dr_new[1], z = -dr_new[2];
          Ylm.evaluateVGH(x, y, z);

          MultiRnl.evaluate(r_new, phi, dphi, d2phi);

          const T rinv = cone / r_new;


          for (size_t ib = 0; ib < BasisSetSize; ++ib)
          {
            const int nl(NL[ib]);
            const int lm(LM[ib]);
            const T drnloverr = rinv * dphi[nl];
            const T ang       = ylm_v[lm];
            const T gr_x      = drnloverr * x;
            const T gr_y      = drnloverr * y;
            const T gr_z      = drnloverr * z;

            //The non-strictly diagonal term in \partial_i \partial_j R_{nl} is
            // \frac{x_i x_j}{r^2}\left(\frac{\partial^2 R_{nl}}{\partial r^2} - \frac{1}{r}\frac{\partial R_{nl}}{\partial r})
            // To save recomputation, I evaluate everything except the x_i*x_j term once, and store it in
            // gr2_tmp.  The full term is obtained by x_i*x_j*gr2_tmp.
            const T gr2_tmp = rinv * rinv * (d2phi[nl] - drnloverr);
            const T gr_xx   = x * x * gr2_tmp + drnloverr;
            const T gr_xy   = x * y * gr2_tmp;
            const T gr_xz   = x * z * gr2_tmp;
            const T gr_yy   = y * y * gr2_tmp + drnloverr;
            const T gr_yz   = y * z * gr2_tmp;
            const T gr_zz   = z * z * gr2_tmp + drnloverr;

            const T ang_x  = ylm_x[lm];
            const T ang_y  = ylm_y[lm];
            const T ang_z  = ylm_z[lm];
            const T ang_xx = ylm_xx[lm];
            const T ang_xy = ylm_xy[lm];
            const T ang_xz = ylm_xz[lm];
            const T ang_yy = ylm_yy[lm];
            const T ang_yz = ylm_yz[lm];
            const T ang_zz = ylm_zz[lm];

            const T vr = phi[nl];

            psi[ib] += ang * vr;
            dpsi_x[ib] += ang * gr_x + vr * ang_x;
            dpsi_y[ib] += ang * gr_y + vr * ang_y;
            dpsi_z[ib] += ang * gr_z + vr * ang_z;


            // \partial_i \partial_j (R*Y) = Y \partial_i \partial_j R + R \partial_i \partial_j Y
            //                             + (\partial_i R) (\partial_j Y) + (\partial_j R)(\partial_i Y)
            dhpsi_xx[ib] += gr_xx * ang + ang_xx * vr + 2.0 * gr_x * ang_x;
            dhpsi_xy[ib] += gr_xy * ang + ang_xy * vr + gr_x * ang_y + gr_y * ang_x;
            dhpsi_xz[ib] += gr_xz * ang + ang_xz * vr + gr_x * ang_z + gr_z * ang_x;
            dhpsi_yy[ib] += gr_yy * ang + ang_yy * vr + 2.0 * gr_y * ang_y;
            dhpsi_yz[ib] += gr_yz * ang + ang_yz * vr + gr_y * ang_z + gr_z * ang_y;
            dhpsi_zz[ib] += gr_zz * ang + ang_zz * vr + 2.0 * gr_z * ang_z;
          }
        }
      }
    }
  }

  template<typename LAT, typename T, typename PosType, typename VGHGH>
  inline void evaluateVGHGH(const LAT& lattice, const T r, const PosType& dr, const size_t offset, VGHGH& vghgh)
  {
    int TransX, TransY, TransZ;

    PosType dr_new;
    T r_new;

    constexpr T cone(1);

    //one can assert the alignment
    RealType* restrict phi   = tempS.data(0);
    RealType* restrict dphi  = tempS.data(1);
    RealType* restrict d2phi = tempS.data(2);
    RealType* restrict d3phi = tempS.data(3);

    //V,Gx,Gy,Gz,L
    auto* restrict psi      = vghgh.data(0) + offset;
    const T* restrict ylm_v = Ylm[0]; //value
    auto* restrict dpsi_x   = vghgh.data(1) + offset;
    const T* restrict ylm_x = Ylm[1]; //gradX
    auto* restrict dpsi_y   = vghgh.data(2) + offset;
    const T* restrict ylm_y = Ylm[2]; //gradY
    auto* restrict dpsi_z   = vghgh.data(3) + offset;
    const T* restrict ylm_z = Ylm[3]; //gradZ

    auto* restrict dhpsi_xx  = vghgh.data(4) + offset;
    const T* restrict ylm_xx = Ylm[4];
    auto* restrict dhpsi_xy  = vghgh.data(5) + offset;
    const T* restrict ylm_xy = Ylm[5];
    auto* restrict dhpsi_xz  = vghgh.data(6) + offset;
    const T* restrict ylm_xz = Ylm[6];
    auto* restrict dhpsi_yy  = vghgh.data(7) + offset;
    const T* restrict ylm_yy = Ylm[7];
    auto* restrict dhpsi_yz  = vghgh.data(8) + offset;
    const T* restrict ylm_yz = Ylm[8];
    auto* restrict dhpsi_zz  = vghgh.data(9) + offset;
    const T* restrict ylm_zz = Ylm[9];

    auto* restrict dghpsi_xxx = vghgh.data(10) + offset;
    const T* restrict ylm_xxx = Ylm[10];
    auto* restrict dghpsi_xxy = vghgh.data(11) + offset;
    const T* restrict ylm_xxy = Ylm[11];
    auto* restrict dghpsi_xxz = vghgh.data(12) + offset;
    const T* restrict ylm_xxz = Ylm[12];
    auto* restrict dghpsi_xyy = vghgh.data(13) + offset;
    const T* restrict ylm_xyy = Ylm[13];
    auto* restrict dghpsi_xyz = vghgh.data(14) + offset;
    const T* restrict ylm_xyz = Ylm[14];
    auto* restrict dghpsi_xzz = vghgh.data(15) + offset;
    const T* restrict ylm_xzz = Ylm[15];
    auto* restrict dghpsi_yyy = vghgh.data(16) + offset;
    const T* restrict ylm_yyy = Ylm[16];
    auto* restrict dghpsi_yyz = vghgh.data(17) + offset;
    const T* restrict ylm_yyz = Ylm[17];
    auto* restrict dghpsi_yzz = vghgh.data(18) + offset;
    const T* restrict ylm_yzz = Ylm[18];
    auto* restrict dghpsi_zzz = vghgh.data(19) + offset;
    const T* restrict ylm_zzz = Ylm[19];

    for (size_t ib = 0; ib < BasisSetSize; ++ib)
    {
      psi[ib] = 0;

      dpsi_x[ib] = 0;
      dpsi_y[ib] = 0;
      dpsi_z[ib] = 0;

      dhpsi_xx[ib] = 0;
      dhpsi_xy[ib] = 0;
      dhpsi_xz[ib] = 0;
      dhpsi_yy[ib] = 0;
      dhpsi_yz[ib] = 0;
      dhpsi_zz[ib] = 0;

      dghpsi_xxx[ib] = 0;
      dghpsi_xxy[ib] = 0;
      dghpsi_xxz[ib] = 0;
      dghpsi_xyy[ib] = 0;
      dghpsi_xyz[ib] = 0;
      dghpsi_xzz[ib] = 0;
      dghpsi_yyy[ib] = 0;
      dghpsi_yyz[ib] = 0;
      dghpsi_yzz[ib] = 0;
      dghpsi_zzz[ib] = 0;
    }

    for (int i = 0; i <= PBCImages[0]; i++) //loop Translation over X
    {
      //Allows to increment cells from 0,1,-1,2,-2,3,-3 etc...
      TransX = ((i % 2) * 2 - 1) * ((i + 1) / 2);
      for (int j = 0; j <= PBCImages[1]; j++) //loop Translation over Y
      {
        //Allows to increment cells from 0,1,-1,2,-2,3,-3 etc...
        TransY = ((j % 2) * 2 - 1) * ((j + 1) / 2);
        for (int k = 0; k <= PBCImages[2]; k++) //loop Translation over Z
        {
          //Allows to increment cells from 0,1,-1,2,-2,3,-3 etc...
          TransZ    = ((k % 2) * 2 - 1) * ((k + 1) / 2);
          dr_new[0] = dr[0] + TransX * lattice.R(0, 0) + TransY * lattice.R(1, 0) + TransZ * lattice.R(2, 0);
          dr_new[1] = dr[1] + TransX * lattice.R(0, 1) + TransY * lattice.R(1, 1) + TransZ * lattice.R(2, 1);
          dr_new[2] = dr[2] + TransX * lattice.R(0, 2) + TransY * lattice.R(1, 2) + TransZ * lattice.R(2, 2);
          r_new     = std::sqrt(dot(dr_new, dr_new));

          //const size_t ib_max=NL.size();
          if (r_new >= Rmax)
            continue;

          //SIGN Change!!
          const T x = -dr_new[0], y = -dr_new[1], z = -dr_new[2];
          Ylm.evaluateVGHGH(x, y, z);

          MultiRnl.evaluate(r_new, phi, dphi, d2phi, d3phi);

          const T rinv = cone / r_new;
          const T xu = x * rinv, yu = y * rinv, zu = z * rinv;
          for (size_t ib = 0; ib < BasisSetSize; ++ib)
          {
            const int nl(NL[ib]);
            const int lm(LM[ib]);
            const T drnloverr = rinv * dphi[nl];
            const T ang       = ylm_v[lm];
            const T gr_x      = drnloverr * x;
            const T gr_y      = drnloverr * y;
            const T gr_z      = drnloverr * z;

            //The non-strictly diagonal term in \partial_i \partial_j R_{nl} is
            // \frac{x_i x_j}{r^2}\left(\frac{\partial^2 R_{nl}}{\partial r^2} - \frac{1}{r}\frac{\partial R_{nl}}{\partial r})
            // To save recomputation, I evaluate everything except the x_i*x_j term once, and store it in
            // gr2_tmp.  The full term is obtained by x_i*x_j*gr2_tmp.  This is p(r) in the notes.
            const T gr2_tmp = rinv * (d2phi[nl] - drnloverr);

            const T gr_xx = x * xu * gr2_tmp + drnloverr;
            const T gr_xy = x * yu * gr2_tmp;
            const T gr_xz = x * zu * gr2_tmp;
            const T gr_yy = y * yu * gr2_tmp + drnloverr;
            const T gr_yz = y * zu * gr2_tmp;
            const T gr_zz = z * zu * gr2_tmp + drnloverr;

            //This is q(r) in the notes.
            const T gr3_tmp = d3phi[nl] - 3.0 * gr2_tmp;

            const T gr_xxx = xu * xu * xu * gr3_tmp + gr2_tmp * (3. * xu);
            const T gr_xxy = xu * xu * yu * gr3_tmp + gr2_tmp * yu;
            const T gr_xxz = xu * xu * zu * gr3_tmp + gr2_tmp * zu;
            const T gr_xyy = xu * yu * yu * gr3_tmp + gr2_tmp * xu;
            const T gr_xyz = xu * yu * zu * gr3_tmp;
            const T gr_xzz = xu * zu * zu * gr3_tmp + gr2_tmp * xu;
            const T gr_yyy = yu * yu * yu * gr3_tmp + gr2_tmp * (3. * yu);
            const T gr_yyz = yu * yu * zu * gr3_tmp + gr2_tmp * zu;
            const T gr_yzz = yu * zu * zu * gr3_tmp + gr2_tmp * yu;
            const T gr_zzz = zu * zu * zu * gr3_tmp + gr2_tmp * (3. * zu);


            //Angular derivatives up to third
            const T ang_x = ylm_x[lm];
            const T ang_y = ylm_y[lm];
            const T ang_z = ylm_z[lm];

            const T ang_xx = ylm_xx[lm];
            const T ang_xy = ylm_xy[lm];
            const T ang_xz = ylm_xz[lm];
            const T ang_yy = ylm_yy[lm];
            const T ang_yz = ylm_yz[lm];
            const T ang_zz = ylm_zz[lm];

            const T ang_xxx = ylm_xxx[lm];
            const T ang_xxy = ylm_xxy[lm];
            const T ang_xxz = ylm_xxz[lm];
            const T ang_xyy = ylm_xyy[lm];
            const T ang_xyz = ylm_xyz[lm];
            const T ang_xzz = ylm_xzz[lm];
            const T ang_yyy = ylm_yyy[lm];
            const T ang_yyz = ylm_yyz[lm];
            const T ang_yzz = ylm_yzz[lm];
            const T ang_zzz = ylm_zzz[lm];

            const T vr = phi[nl];

            psi[ib] += ang * vr;
            dpsi_x[ib] += ang * gr_x + vr * ang_x;
            dpsi_y[ib] += ang * gr_y + vr * ang_y;
            dpsi_z[ib] += ang * gr_z + vr * ang_z;


            // \partial_i \partial_j (R*Y) = Y \partial_i \partial_j R + R \partial_i \partial_j Y
            //                             + (\partial_i R) (\partial_j Y) + (\partial_j R)(\partial_i Y)
            dhpsi_xx[ib] += gr_xx * ang + ang_xx * vr + 2.0 * gr_x * ang_x;
            dhpsi_xy[ib] += gr_xy * ang + ang_xy * vr + gr_x * ang_y + gr_y * ang_x;
            dhpsi_xz[ib] += gr_xz * ang + ang_xz * vr + gr_x * ang_z + gr_z * ang_x;
            dhpsi_yy[ib] += gr_yy * ang + ang_yy * vr + 2.0 * gr_y * ang_y;
            dhpsi_yz[ib] += gr_yz * ang + ang_yz * vr + gr_y * ang_z + gr_z * ang_y;
            dhpsi_zz[ib] += gr_zz * ang + ang_zz * vr + 2.0 * gr_z * ang_z;

            dghpsi_xxx[ib] += gr_xxx * ang + vr * ang_xxx + 3.0 * gr_xx * ang_x + 3.0 * gr_x * ang_xx;
            dghpsi_xxy[ib] +=
                gr_xxy * ang + vr * ang_xxy + gr_xx * ang_y + ang_xx * gr_y + 2.0 * gr_xy * ang_x + 2.0 * ang_xy * gr_x;
            dghpsi_xxz[ib] +=
                gr_xxz * ang + vr * ang_xxz + gr_xx * ang_z + ang_xx * gr_z + 2.0 * gr_xz * ang_x + 2.0 * ang_xz * gr_x;
            dghpsi_xyy[ib] +=
                gr_xyy * ang + vr * ang_xyy + gr_yy * ang_x + ang_yy * gr_x + 2.0 * gr_xy * ang_y + 2.0 * ang_xy * gr_y;
            dghpsi_xyz[ib] += gr_xyz * ang + vr * ang_xyz + gr_xy * ang_z + ang_xy * gr_z + gr_yz * ang_x +
                ang_yz * gr_x + gr_xz * ang_y + ang_xz * gr_y;
            dghpsi_xzz[ib] +=
                gr_xzz * ang + vr * ang_xzz + gr_zz * ang_x + ang_zz * gr_x + 2.0 * gr_xz * ang_z + 2.0 * ang_xz * gr_z;
            dghpsi_yyy[ib] += gr_yyy * ang + vr * ang_yyy + 3.0 * gr_yy * ang_y + 3.0 * gr_y * ang_yy;
            dghpsi_yyz[ib] +=
                gr_yyz * ang + vr * ang_yyz + gr_yy * ang_z + ang_yy * gr_z + 2.0 * gr_yz * ang_y + 2.0 * ang_yz * gr_y;
            dghpsi_yzz[ib] +=
                gr_yzz * ang + vr * ang_yzz + gr_zz * ang_y + ang_zz * gr_y + 2.0 * gr_yz * ang_z + 2.0 * ang_yz * gr_z;
            dghpsi_zzz[ib] += gr_zzz * ang + vr * ang_zzz + 3.0 * gr_zz * ang_z + 3.0 * gr_z * ang_zz;
          }
        }
      }
    }
  }

  /** evaluate V
  */
  template<typename LAT, typename T, typename PosType, typename VT>
  inline void evaluateV(const LAT& lattice, const T r, const PosType& dr, VT* restrict psi, PosType Tv)
  {
    int TransX, TransY, TransZ;

    PosType dr_new;
    T r_new;

#if not defined(QMC_COMPLEX)
    const ValueType correctphase = 1.0;
#else

    RealType phasearg = SuperTwist[0] * Tv[0] + SuperTwist[1] * Tv[1] + SuperTwist[2] * Tv[2];
    RealType s, c;
    qmcplusplus::sincos(-phasearg, &s, &c);
    const ValueType correctphase(c, s);

#endif

    RealType* restrict ylm_v = tempS.data(0);
    RealType* restrict phi_r = tempS.data(1);

    for (size_t ib = 0; ib < BasisSetSize; ++ib)
      psi[ib] = 0;
    //Phase_idx (iter) needs to be initialized at -1 as it has to be incremented first to comply with the if statement (r_new >=Rmax)
    int iter = -1;
    for (int i = 0; i <= PBCImages[0]; i++) //loop Translation over X
    {
      //Allows to increment cells from 0,1,-1,2,-2,3,-3 etc...
      TransX = ((i % 2) * 2 - 1) * ((i + 1) / 2);
      for (int j = 0; j <= PBCImages[1]; j++) //loop Translation over Y
      {
        //Allows to increment cells from 0,1,-1,2,-2,3,-3 etc...
        TransY = ((j % 2) * 2 - 1) * ((j + 1) / 2);
        for (int k = 0; k <= PBCImages[2]; k++) //loop Translation over Z
        {
          //Allows to increment cells from 0,1,-1,2,-2,3,-3 etc...
          TransZ = ((k % 2) * 2 - 1) * ((k + 1) / 2);

          dr_new[0] = dr[0] + (TransX * lattice.R(0, 0) + TransY * lattice.R(1, 0) + TransZ * lattice.R(2, 0));
          dr_new[1] = dr[1] + (TransX * lattice.R(0, 1) + TransY * lattice.R(1, 1) + TransZ * lattice.R(2, 1));
          dr_new[2] = dr[2] + (TransX * lattice.R(0, 2) + TransY * lattice.R(1, 2) + TransZ * lattice.R(2, 2));

          r_new = std::sqrt(dot(dr_new, dr_new));
          iter++;
          if (r_new >= Rmax)
            continue;

          Ylm.evaluateV(-dr_new[0], -dr_new[1], -dr_new[2], ylm_v);
          MultiRnl.evaluate(r_new, phi_r);
          ///Phase for PBC containing the phase for the nearest image displacement and the correction due to the Distance table.
          const ValueType Phase = periodic_image_phase_factors_[iter] * correctphase;
          for (size_t ib = 0; ib < BasisSetSize; ++ib)
            psi[ib] += ylm_v[LM[ib]] * phi_r[NL[ib]] * Phase;
        }
      }
    }
  }

  /**
   * @brief evaluate VGL for multiple electrons
   * 
   * This function should only assign to elements of psi in the range [[0:nElec],[BasisOffset:BasisOffset+BasisSetSize]].
   * These elements are assumed to be zero when passed to this function.
   * This function only uses only one center (center_idx) from displ_list
   * 
   * @param [in] atom_bs_list multi-walker list of SoaAtomicBasisSet [nWalkers]
   * @param [in] lattice crystal lattice
   * @param [in,out] psi_vgl wavefunction vgl for all electrons [5, nElec, nBasTot]
   * @param [in] displ_list displacement from each electron to each center [NumCenters, nElec, 3] (flattened)
   * @param [in] Tv_list translation vectors for computing overall phase factor [NumCenters, nElec, 3] (flattened)
   * @param [in] nElec number of electrons
   * @param [in] nBasTot total number of basis functions represented in psi_vgl
   * @param [in] center_idx current center index (for indexing into displ_list)
   * @param [in] BasisOffset index of first basis function of this center (for indexing into psi_vgl)
   * @param [in] NumCenters total number of centers in system (for indexing into displ_list)
   *  
  */

  template<typename LAT, typename VT>
  inline void mw_evaluateVGL(const RefVectorWithLeader<SoaAtomicBasisSet>& atom_bs_list,
                             const LAT& lattice,
                             Array<VT, 3, OffloadPinnedAllocator<VT>>& psi_vgl,
                             const Vector<RealType, OffloadPinnedAllocator<RealType>>& displ_list,
                             const Vector<RealType, OffloadPinnedAllocator<RealType>>& Tv_list,
                             const size_t nElec,
                             const size_t nBasTot,
                             const size_t center_idx,
                             const size_t BasisOffset,
                             const size_t NumCenters)
  {
    assert(this == &atom_bs_list.getLeader());
    auto& atom_bs_leader = atom_bs_list.template getCastedLeader<SoaAtomicBasisSet<ROT, SH>>();

    int Nx         = PBCImages[0] + 1;
    int Ny         = PBCImages[1] + 1;
    int Nz         = PBCImages[2] + 1;
    const int Nxyz = Nx * Ny * Nz;

    assert(psi_vgl.size(0) == 5);
    assert(psi_vgl.size(1) == nElec);
    assert(psi_vgl.size(2) == nBasTot);


    auto& ylm_vgl = atom_bs_leader.mw_mem_handle_.getResource().ylm_vgl;
    auto& rnl_vgl = atom_bs_leader.mw_mem_handle_.getResource().rnl_vgl;
    auto& dr      = atom_bs_leader.mw_mem_handle_.getResource().dr;
    auto& r       = atom_bs_leader.mw_mem_handle_.getResource().r;

    size_t nRnl = RnlID.size();
    size_t nYlm = Ylm.size();

    ylm_vgl.resize(5, nElec, Nxyz, nYlm);
    rnl_vgl.resize(3, nElec, Nxyz, nRnl);
    dr.resize(nElec, Nxyz, 3);
    r.resize(nElec, Nxyz);


    // TODO: move these outside?
    auto& correctphase = atom_bs_leader.mw_mem_handle_.getResource().correctphase;
    correctphase.resize(nElec);

    auto* dr_ptr = dr.data();
    auto* r_ptr  = r.data();

    auto* correctphase_ptr = correctphase.data();

    auto* Tv_list_ptr    = Tv_list.data();
    auto* displ_list_ptr = displ_list.data();

    constexpr RealType cone(1);
    constexpr RealType ctwo(2);

    //V,Gx,Gy,Gz,L
    auto* restrict psi_ptr    = psi_vgl.data_at(0, 0, 0);
    auto* restrict dpsi_x_ptr = psi_vgl.data_at(1, 0, 0);
    auto* restrict dpsi_y_ptr = psi_vgl.data_at(2, 0, 0);
    auto* restrict dpsi_z_ptr = psi_vgl.data_at(3, 0, 0);
    auto* restrict d2psi_ptr  = psi_vgl.data_at(4, 0, 0);

    {
      ScopedTimer local_timer(phase_timer_);
#if not defined(QMC_COMPLEX)

      PRAGMA_OFFLOAD("omp target teams distribute parallel for map(to:correctphase_ptr[:nElec]) ")
      for (size_t i_e = 0; i_e < nElec; i_e++)
        correctphase_ptr[i_e] = 1.0;

#else
      auto* SuperTwist_ptr = SuperTwist.data();

      PRAGMA_OFFLOAD("omp target teams distribute parallel for map(to:SuperTwist_ptr[:SuperTwist.size()], \
		      Tv_list_ptr[3*nElec*center_idx:3*nElec], correctphase_ptr[:nElec]) ")
      for (size_t i_e = 0; i_e < nElec; i_e++)
      {
        //RealType phasearg = dot(3, SuperTwist.data(), 1, Tv_list.data() + 3 * i_e, 1);
        RealType phasearg = 0;
        for (size_t i_dim = 0; i_dim < 3; i_dim++)
          phasearg += SuperTwist[i_dim] * Tv_list_ptr[i_dim + 3 * (i_e + center_idx * nElec)];
        RealType s, c;
        qmcplusplus::sincos(-phasearg, &s, &c);
        correctphase_ptr[i_e] = ValueType(c, s);
      }
#endif
    }

    {
      ScopedTimer local_timer(nelec_pbc_timer_);
      auto* periodic_image_displacements_ptr = periodic_image_displacements_.data();
      PRAGMA_OFFLOAD("omp target teams distribute parallel for collapse(2) \
                      map(to:periodic_image_displacements_ptr[:3*Nxyz]) \
                      map(to: dr_ptr[:3*nElec*Nxyz], r_ptr[:nElec*Nxyz], displ_list_ptr[3*nElec*center_idx:3*nElec]) ")
      for (size_t i_e = 0; i_e < nElec; i_e++)
        for (int i_xyz = 0; i_xyz < Nxyz; i_xyz++)
        {
          RealType tmp_r2 = 0.0;
          for (size_t i_dim = 0; i_dim < 3; i_dim++)
          {
            dr_ptr[i_dim + 3 * (i_xyz + Nxyz * i_e)] = -(displ_list_ptr[i_dim + 3 * (i_e + center_idx * nElec)] +
                                                         periodic_image_displacements_ptr[i_dim + 3 * i_xyz]);
            tmp_r2 += dr_ptr[i_dim + 3 * (i_xyz + Nxyz * i_e)] * dr_ptr[i_dim + 3 * (i_xyz + Nxyz * i_e)];
          }
          r_ptr[i_xyz + Nxyz * i_e] = std::sqrt(tmp_r2);
          //printf("particle %lu image %d, %lf, %lf\n", i_e, i_xyz, tmp_r2, dr_ptr[3 * (i_xyz + Nxyz * i_e)]);
        }
    }

    {
      ScopedTimer local(rnl_timer_);
      MultiRnl.batched_evaluateVGL(r, rnl_vgl, Rmax);
    }

    {
      ScopedTimer local(ylm_timer_);
      Ylm.batched_evaluateVGL(dr, ylm_vgl);
    }

    {
      ScopedTimer local_timer(psi_timer_);
      auto* phase_fac_ptr = periodic_image_phase_factors_.data();
      auto* LM_ptr        = LM.data();
      auto* NL_ptr        = NL.data();
      const int bset_size = BasisSetSize;

      RealType* restrict phi_ptr   = rnl_vgl.data_at(0, 0, 0, 0);
      RealType* restrict dphi_ptr  = rnl_vgl.data_at(1, 0, 0, 0);
      RealType* restrict d2phi_ptr = rnl_vgl.data_at(2, 0, 0, 0);


      const RealType* restrict ylm_v_ptr = ylm_vgl.data_at(0, 0, 0, 0); //value
      const RealType* restrict ylm_x_ptr = ylm_vgl.data_at(1, 0, 0, 0); //gradX
      const RealType* restrict ylm_y_ptr = ylm_vgl.data_at(2, 0, 0, 0); //gradY
      const RealType* restrict ylm_z_ptr = ylm_vgl.data_at(3, 0, 0, 0); //gradZ
      const RealType* restrict ylm_l_ptr = ylm_vgl.data_at(4, 0, 0, 0); //lap
      PRAGMA_OFFLOAD("omp target teams distribute parallel for collapse(2) \
                      map(to:phase_fac_ptr[:Nxyz], LM_ptr[:BasisSetSize], NL_ptr[:BasisSetSize]) \
		      map(to:ylm_v_ptr[:nYlm*nElec*Nxyz], ylm_x_ptr[:nYlm*nElec*Nxyz], ylm_y_ptr[:nYlm*nElec*Nxyz], ylm_z_ptr[:nYlm*nElec*Nxyz], ylm_l_ptr[:nYlm*nElec*Nxyz], \
                      phi_ptr[:nRnl*nElec*Nxyz], dphi_ptr[:nRnl*nElec*Nxyz], d2phi_ptr[:nRnl*nElec*Nxyz], \
                      psi_ptr[:nBasTot*nElec], dpsi_x_ptr[:nBasTot*nElec], dpsi_y_ptr[:nBasTot*nElec], dpsi_z_ptr[:nBasTot*nElec], d2psi_ptr[:nBasTot*nElec], \
                      correctphase_ptr[:nElec], r_ptr[:nElec*Nxyz], dr_ptr[:3*nElec*Nxyz]) ")
      for (int i_e = 0; i_e < nElec; i_e++)
        for (int ib = 0; ib < bset_size; ++ib)
        {
          const int nl(NL_ptr[ib]);
          const int lm(LM_ptr[ib]);
          VT psi    = 0;
          VT dpsi_x = 0;
          VT dpsi_y = 0;
          VT dpsi_z = 0;
          VT d2psi  = 0;

          for (int i_xyz = 0; i_xyz < Nxyz; i_xyz++)
          {
            const ValueType Phase    = phase_fac_ptr[i_xyz] * correctphase_ptr[i_e];
            const RealType rinv      = cone / r_ptr[i_xyz + Nxyz * i_e];
            const RealType x         = dr_ptr[0 + 3 * (i_xyz + Nxyz * i_e)];
            const RealType y         = dr_ptr[1 + 3 * (i_xyz + Nxyz * i_e)];
            const RealType z         = dr_ptr[2 + 3 * (i_xyz + Nxyz * i_e)];
            const RealType drnloverr = rinv * dphi_ptr[nl + nRnl * (i_xyz + Nxyz * i_e)];
            const RealType ang       = ylm_v_ptr[lm + nYlm * (i_xyz + Nxyz * i_e)];
            const RealType gr_x      = drnloverr * x;
            const RealType gr_y      = drnloverr * y;
            const RealType gr_z      = drnloverr * z;
            const RealType ang_x     = ylm_x_ptr[lm + nYlm * (i_xyz + Nxyz * i_e)];
            const RealType ang_y     = ylm_y_ptr[lm + nYlm * (i_xyz + Nxyz * i_e)];
            const RealType ang_z     = ylm_z_ptr[lm + nYlm * (i_xyz + Nxyz * i_e)];
            const RealType vr        = phi_ptr[nl + nRnl * (i_xyz + Nxyz * i_e)];

            psi += ang * vr * Phase;
            dpsi_x += (ang * gr_x + vr * ang_x) * Phase;
            dpsi_y += (ang * gr_y + vr * ang_y) * Phase;
            dpsi_z += (ang * gr_z + vr * ang_z) * Phase;
            d2psi += (ang * (ctwo * drnloverr + d2phi_ptr[nl + nRnl * (i_xyz + Nxyz * i_e)]) +
                      ctwo * (gr_x * ang_x + gr_y * ang_y + gr_z * ang_z) +
                      vr * ylm_l_ptr[lm + nYlm * (i_xyz + Nxyz * i_e)]) *
                Phase;
          }

          psi_ptr[BasisOffset + ib + i_e * nBasTot]    = psi;
          dpsi_x_ptr[BasisOffset + ib + i_e * nBasTot] = dpsi_x;
          dpsi_y_ptr[BasisOffset + ib + i_e * nBasTot] = dpsi_y;
          dpsi_z_ptr[BasisOffset + ib + i_e * nBasTot] = dpsi_z;
          d2psi_ptr[BasisOffset + ib + i_e * nBasTot]  = d2psi;
        }
    }
  }

  /**
   * @brief evaluate for multiple electrons
   * 
   * This function should only assign to elements of psi in the range [[0:nElec],[BasisOffset:BasisOffset+BasisSetSize]].
   * These elements are assumed to be zero when passed to this function.
   * This function only uses only one center (center_idx) from displ_list
   * 
   * @param [in] atom_bs_list multi-walker list of SoaAtomicBasisSet [nWalkers]
   * @param [in] lattice crystal lattice
   * @param [in,out] psi wavefunction values for all electrons [nElec, nBasTot]
   * @param [in] displ_list displacement from each electron to each center [NumCenters, nElec, 3] (flattened)
   * @param [in] Tv_list translation vectors for computing overall phase factor [NumCenters, nElec, 3] (flattened)
   * @param [in] nElec number of electrons
   * @param [in] nBasTot total number of basis functions represented in psi
   * @param [in] center_idx current center index (for indexing into displ_list)
   * @param [in] BasisOffset index of first basis function of this center (for indexing into psi)
   * @param [in] NumCenters total number of centers in system (for indexing into displ_list)
   *  
  */
  template<typename LAT, typename VT>
  inline void mw_evaluateV(const RefVectorWithLeader<SoaAtomicBasisSet>& atom_bs_list,
                           const LAT& lattice,
                           Array<VT, 2, OffloadPinnedAllocator<VT>>& psi,
                           const Vector<RealType, OffloadPinnedAllocator<RealType>>& displ_list,
                           const Vector<RealType, OffloadPinnedAllocator<RealType>>& Tv_list,
                           const size_t nElec,
                           const size_t nBasTot,
                           const size_t center_idx,
                           const size_t BasisOffset,
                           const size_t NumCenters)
  {
    assert(this == &atom_bs_list.getLeader());
    auto& atom_bs_leader = atom_bs_list.template getCastedLeader<SoaAtomicBasisSet<ROT, SH>>();
    //TODO: use QMCTraits::DIM instead of 3?
    //      DIM==3 is baked into so many parts here that it's probably not worth it for now
    const int Nx   = PBCImages[0] + 1;
    const int Ny   = PBCImages[1] + 1;
    const int Nz   = PBCImages[2] + 1;
    const int Nxyz = Nx * Ny * Nz;
    assert(psi.size(0) == nElec);
    assert(psi.size(1) == nBasTot);


    auto& ylm_v = atom_bs_leader.mw_mem_handle_.getResource().ylm_v;
    auto& rnl_v = atom_bs_leader.mw_mem_handle_.getResource().rnl_v;
    auto& dr    = atom_bs_leader.mw_mem_handle_.getResource().dr;
    auto& r     = atom_bs_leader.mw_mem_handle_.getResource().r;

    const size_t nRnl = RnlID.size();
    const size_t nYlm = Ylm.size();

    ylm_v.resize(nElec, Nxyz, nYlm);
    rnl_v.resize(nElec, Nxyz, nRnl);
    dr.resize(nElec, Nxyz, 3);
    r.resize(nElec, Nxyz);

    // TODO: move these outside?
    auto& correctphase = atom_bs_leader.mw_mem_handle_.getResource().correctphase;
    correctphase.resize(nElec);

    auto* dr_ptr = dr.data();
    auto* r_ptr  = r.data();

    auto* correctphase_ptr = correctphase.data();

    auto* Tv_list_ptr    = Tv_list.data();
    auto* displ_list_ptr = displ_list.data();

    // need to map Tensor<T,3> vals to device
    auto* latR_ptr = lattice.R.data();


    {
      ScopedTimer local_timer(phase_timer_);
#if not defined(QMC_COMPLEX)

      PRAGMA_OFFLOAD("omp target teams distribute parallel for map(to:correctphase_ptr[:nElec]) ")
      for (size_t i_e = 0; i_e < nElec; i_e++)
        correctphase_ptr[i_e] = 1.0;

#else
      auto* SuperTwist_ptr = SuperTwist.data();

      PRAGMA_OFFLOAD("omp target teams distribute parallel for map(to:SuperTwist_ptr[:SuperTwist.size()], \
		      Tv_list_ptr[3*nElec*center_idx:3*nElec], correctphase_ptr[:nElec]) ")
      for (size_t i_e = 0; i_e < nElec; i_e++)
      {
        //RealType phasearg = dot(3, SuperTwist.data(), 1, Tv_list.data() + 3 * i_e, 1);
        RealType phasearg = 0;
        for (size_t i_dim = 0; i_dim < 3; i_dim++)
          phasearg += SuperTwist[i_dim] * Tv_list_ptr[i_dim + 3 * (i_e + center_idx * nElec)];
        RealType s, c;
        qmcplusplus::sincos(-phasearg, &s, &c);
        correctphase_ptr[i_e] = ValueType(c, s);
      }
#endif
    }

    {
      ScopedTimer local_timer(nelec_pbc_timer_);
      auto* periodic_image_displacements_ptr = periodic_image_displacements_.data();
      PRAGMA_OFFLOAD("omp target teams distribute parallel for collapse(2) \
                      map(to:periodic_image_displacements_ptr[:3*Nxyz]) \
                      map(to: dr_ptr[:3*nElec*Nxyz], r_ptr[:nElec*Nxyz], displ_list_ptr[3*nElec*center_idx:3*nElec]) ")
      for (size_t i_e = 0; i_e < nElec; i_e++)
        for (int i_xyz = 0; i_xyz < Nxyz; i_xyz++)
        {
          RealType tmp_r2 = 0.0;
          for (size_t i_dim = 0; i_dim < 3; i_dim++)
          {
            dr_ptr[i_dim + 3 * (i_xyz + Nxyz * i_e)] = -(displ_list_ptr[i_dim + 3 * (i_e + center_idx * nElec)] +
                                                         periodic_image_displacements_ptr[i_dim + 3 * i_xyz]);
            tmp_r2 += dr_ptr[i_dim + 3 * (i_xyz + Nxyz * i_e)] * dr_ptr[i_dim + 3 * (i_xyz + Nxyz * i_e)];
          }
          r_ptr[i_xyz + Nxyz * i_e] = std::sqrt(tmp_r2);
        }
    }


    {
      ScopedTimer local(rnl_timer_);
      MultiRnl.batched_evaluate(r, rnl_v, Rmax);
    }

    {
      ScopedTimer local(ylm_timer_);
      Ylm.batched_evaluateV(dr, ylm_v);
    }

    {
      ScopedTimer local_timer(psi_timer_);
      ///Phase for PBC containing the phase for the nearest image displacement and the correction due to the Distance table.
      auto* phase_fac_ptr = periodic_image_phase_factors_.data();
      auto* LM_ptr        = LM.data();
      auto* NL_ptr        = NL.data();
      auto* psi_ptr       = psi.data();
      const int bset_size = BasisSetSize;

      auto* ylm_ptr = ylm_v.data();
      auto* rnl_ptr = rnl_v.data();
      PRAGMA_OFFLOAD("omp target teams distribute parallel for collapse(2) \
                      map(to:phase_fac_ptr[:Nxyz], LM_ptr[:BasisSetSize], NL_ptr[:BasisSetSize]) \
		      map(to:ylm_ptr[:nYlm*nElec*Nxyz], rnl_ptr[:nRnl*nElec*Nxyz], psi_ptr[:nBasTot*nElec], correctphase_ptr[:nElec])")
      for (int i_e = 0; i_e < nElec; i_e++)
        for (int ib = 0; ib < bset_size; ++ib)
        {
          VT psi = 0;
          for (int i_xyz = 0; i_xyz < Nxyz; i_xyz++)
          {
            const ValueType Phase = phase_fac_ptr[i_xyz] * correctphase_ptr[i_e];
            psi += ylm_ptr[(i_xyz + Nxyz * i_e) * nYlm + LM_ptr[ib]] *
                rnl_ptr[(i_xyz + Nxyz * i_e) * nRnl + NL_ptr[ib]] * Phase;
          }
          psi_ptr[BasisOffset + ib + i_e * nBasTot] = psi;
        }
    }
  }

  void createResource(ResourceCollection& collection) const
  {
    collection.addResource(std::make_unique<SoaAtomicBSetMultiWalkerMem>());
  }

  void acquireResource(ResourceCollection& collection,
                       const RefVectorWithLeader<SoaAtomicBasisSet>& atom_basis_list) const
  {
    assert(this == &atom_basis_list.getLeader());
    atom_basis_list.template getCastedLeader<SoaAtomicBasisSet>().mw_mem_handle_ =
        collection.lendResource<SoaAtomicBSetMultiWalkerMem>();
  }

  void releaseResource(ResourceCollection& collection,
                       const RefVectorWithLeader<SoaAtomicBasisSet>& atom_basis_list) const
  {
    assert(this == &atom_basis_list.getLeader());
    collection.takebackResource(atom_basis_list.template getCastedLeader<SoaAtomicBasisSet>().mw_mem_handle_);
  }

private:
  /// multi walker shared memory buffer
  struct SoaAtomicBSetMultiWalkerMem : public Resource
  {
    SoaAtomicBSetMultiWalkerMem() : Resource("SoaAtomicBasisSet") {}

    SoaAtomicBSetMultiWalkerMem(const SoaAtomicBSetMultiWalkerMem&) : SoaAtomicBSetMultiWalkerMem() {}

    std::unique_ptr<Resource> makeClone() const override
    {
      return std::make_unique<SoaAtomicBSetMultiWalkerMem>(*this);
    }

    OffloadArray4D ylm_vgl;     // [5][Nelec][PBC][NYlm]
    OffloadArray4D rnl_vgl;     // [5][Nelec][PBC][NRnl]
    OffloadArray3D ylm_v;       // [Nelec][PBC][NYlm]
    OffloadArray3D rnl_v;       // [Nelec][PBC][NRnl]
    OffloadArray3D dr;          // [Nelec][PBC][xyz] ion->elec displacement for each image
    OffloadArray2D r;           // [Nelec][PBC]      ion->elec distance for each image
    OffloadVector correctphase; // [Nelec]           overall phase
  };

  /// multi walker resource handle
  ResourceHandle<SoaAtomicBSetMultiWalkerMem> mw_mem_handle_;
  ///size of the basis set
  int BasisSetSize;
  ///Number of Cell images for the evaluation of the orbital with PBC. If No PBC, should be 0;
  TinyVector<int, 3> PBCImages;
  ///Coordinates of SuperTwist
  TinyVector<double, 3> SuperTwist;
  ///maximum radius of this center
  RealType Rmax;
  ///spherical harmonics
  SH Ylm;
  ///radial orbitals
  ROT MultiRnl;
  ///container for the quantum-numbers
  std::vector<QuantumNumberType> RnlID;
  ///temporary storage
  VectorSoaContainer<RealType, 4> tempS;
  ///Phase Factor array of images
  std::shared_ptr<OffloadVector> periodic_image_phase_factors_ptr_;
  ///Displacements of images
  std::shared_ptr<OffloadArray2D> periodic_image_displacements_ptr_;
  ///reference to the phase factor array of images
  OffloadVector& periodic_image_phase_factors_;
  ///reference to the displacements of images
  OffloadArray2D& periodic_image_displacements_;
  /**index of the corresponding radial orbital with quantum numbers \f$ (n,l) \f$ */
  const std::shared_ptr<OffloadIntVector> NL_ptr_;
  ///index of the corresponding real Spherical Harmonic with quantum numbers \f$ (l,m) \f$
  const std::shared_ptr<OffloadIntVector> LM_ptr_;
  /// reference to NL_ptr_
  OffloadIntVector& NL;
  /// reference to LM_ptr_
  OffloadIntVector& LM;
  // timers
  NewTimer& ylm_timer_;
  NewTimer& rnl_timer_;
  NewTimer& pbc_timer_;
  NewTimer& nelec_pbc_timer_;
  NewTimer& phase_timer_;
  NewTimer& psi_timer_;

  template<typename COT>
  friend class AOBasisBuilder;
  template<typename COT>
  friend class RadialOrbitalSetBuilder;
};

} // namespace qmcplusplus
#endif
