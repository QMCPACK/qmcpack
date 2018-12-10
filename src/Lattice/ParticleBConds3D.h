//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    



#ifndef QMCPLUSPLUS_PARTICLE_BCONDS_3D_H
#define QMCPLUSPLUS_PARTICLE_BCONDS_3D_H

#include <config.h>
#include <Lattice/CrystalLattice.h>

namespace qmcplusplus
{
/** specialization for a periodic 3D, orthorombic cell
*/
template<class T>
struct DTD_BConds<T,3,PPPO>
{
  T Linv0,L0,Linv1,L1,Linv2,L2,r2max,dummy;

  inline DTD_BConds(const CrystalLattice<T,3>& lat)
    : Linv0(lat.OneOverLength[0]), L0(lat.Length[0])
    ,Linv1(lat.OneOverLength[1]), L1(lat.Length[1])
    ,Linv2(lat.OneOverLength[2]), L2(lat.Length[2])
    ,r2max(lat.CellRadiusSq),dummy(T())
  {
  }

  /** apply BC on a displacement vector
   * @param displ
   */
  inline T apply_bc(TinyVector<T,3>& displ)  const
  {
    T x=displ[0]*Linv0;
    displ[0]=L0*(x-round(x));
    T y=displ[1]*Linv1;
    displ[1]=L1*(y-round(y));
    T z=displ[2]*Linv2;
    displ[2]=L2*(z-round(z));
    return displ[0]*displ[0]+displ[1]*displ[1]+displ[2]*displ[2];
  }

  /** evaluate displacement data for a vector
   */
  inline void apply_bc(std::vector<TinyVector<T,3> >& dr
                       , std::vector<T>& r
                       , std::vector<T>& rinv) const
  {
    const int n=r.size();
    //use rinv as temporary rr
    for(int i=0; i<n; ++i)
      rinv[i]=apply_bc(dr[i]);
    simd::sqrt(&rinv[0],&r[0],n);
    simd::inv(&r[0],&rinv[0],n);
  }

  inline void apply_bc(std::vector<TinyVector<T,3> >& dr
                       , std::vector<T>& r) const
  {
    for(int i=0; i<dr.size(); ++i)
      r[i]=dot(dr[i],dr[i]);
  }

  inline void evaluate_rsquared(TinyVector<T,3>* restrict dr, T* restrict rr, int n)
  {
    for(int i=0; i<n; ++i)
      rr[i]=apply_bc(dr[i]);
  }
};

/** specialization for a periodic 3D general cell with wigner-seitz==simulation cell
 *
 * Skip image cells.
 */
template<class T>
struct DTD_BConds<T,3,PPPS>
{
  T r00,r10,r20,r01,r11,r21,r02,r12,r22;
  T g00,g10,g20,g01,g11,g21,g02,g12,g22;
  DTD_BConds(const CrystalLattice<T,3>& lat)
    : r00(lat.R(0)),r10(lat.R(3)),r20(lat.R(6))
    ,r01(lat.R(1)),r11(lat.R(4)),r21(lat.R(7))
    ,r02(lat.R(2)),r12(lat.R(5)),r22(lat.R(8))
    ,g00(lat.G(0)),g10(lat.G(3)),g20(lat.G(6))
    ,g01(lat.G(1)),g11(lat.G(4)),g21(lat.G(7))
    ,g02(lat.G(2)),g12(lat.G(5)),g22(lat.G(8))
  { }

  /** apply BC to a displacement vector a and return the minimum-image distance
   * @param lat lattice
   * @param a displacement vector
   * @return the minimum-image distance
   */
  inline T apply_bc(TinyVector<T,3>& displ)  const
  {
    //cart2unit
    TinyVector<T,3>
    ar(displ[0]*g00+displ[1]*g10+displ[2]*g20
       ,displ[0]*g01+displ[1]*g11+displ[2]*g21
       ,displ[0]*g02+displ[1]*g12+displ[2]*g22);
    //put them in the box
    ar[0]-=round(ar[0]);
    ar[1]-=round(ar[1]);
    ar[2]-=round(ar[2]);
    //unit2cart
    displ[0]=ar[0]*r00+ar[1]*r10+ar[2]*r20;
    displ[1]=ar[0]*r01+ar[1]*r11+ar[2]*r21;
    displ[2]=ar[0]*r02+ar[1]*r12+ar[2]*r22;
    return displ[0]*displ[0]+displ[1]*displ[1]+displ[2]*displ[2];
  }

  inline void apply_bc(std::vector<TinyVector<T,3> >& dr
                       , std::vector<T>& r
                       , std::vector<T>& rinv) const
  {
    const int n=dr.size();
    for(int i=0; i<n; ++i)
      rinv[i]=apply_bc(dr[i]);
    simd::sqrt(&rinv[0],&r[0],n);
    simd::inv(&r[0],&rinv[0],n);
  }

  inline void apply_bc(std::vector<TinyVector<T,3> >& dr
                       , std::vector<T>& r) const
  {
    for(int i=0; i<dr.size(); ++i)
      r[i]=apply_bc(dr[i]);
  }

  inline void evaluate_rsquared(TinyVector<T,3>* restrict dr, T* restrict rr, int n)
  {
    for(int i=0; i<n; ++i)
      rr[i]=apply_bc(dr[i]);
  }
};

/** specialization for a periodic 3D general cell
 *
 * Wigner-Seitz cell radius > simulation cell radius
 * Need to check image cells
*/
template<class T>
struct DTD_BConds<T,3,PPPG>
{
#ifdef BGQPX
  vector4double g0, g1, g2;
  vector4double rb0, rb1, rb2;
  vector4double corner_x1_vec, corner_y1_vec, corner_z1_vec;
  vector4double corner_x2_vec, corner_y2_vec, corner_z2_vec;
#else
  T g00,g10,g20,g01,g11,g21,g02,g12,g22;
#endif
  TinyVector<TinyVector<T,3>,3> rb;
  std::vector<TinyVector<T,3> > corners;

  DTD_BConds(const CrystalLattice<T,3>& lat)
  {
    rb[0]=lat.a(0);
    rb[1]=lat.a(1);
    rb[2]=lat.a(2);
    find_reduced_basis(rb);
    Tensor<T,3> rbt;
    for(int i=0; i<3; ++i)
      for(int j=0; j<3; ++j)
        rbt(i,j)=rb[i][j];
    Tensor<T,3> g=inverse(rbt);
    T minusone=-1.0;
    corners.resize(8);
    corners[0]=0.0;
    corners[1]=minusone*(rb[0]);
    corners[2]=minusone*(rb[1]);
    corners[3]=minusone*(rb[2]);
    corners[4]=minusone*(rb[0]+rb[1]);
    corners[5]=minusone*(rb[0]+rb[2]);
    corners[6]=minusone*(rb[1]+rb[2]);
    corners[7]=minusone*(rb[0]+rb[1]+rb[2]);

#ifdef BGQPX
    g0[0]=g(0);
    g0[1]=g(1);
    g0[2]=g(2);
    g0[3]=0.0;
    g1[0]=g(3);
    g1[1]=g(4);
    g1[2]=g(5);
    g1[3]=0.0;
    g2[0]=g(6);
    g2[1]=g(7);
    g2[2]=g(8);
    g2[3]=0.0;

    for(int i=0;i<3;i++){
      rb0[i]=rb[0][i];
      rb1[i]=rb[1][i];
      rb2[i]=rb[2][i];
    }

    // corners 0~3
    corner_x1_vec[0]=0.0      ;
    corner_x1_vec[1]=-rb[0][0];
    corner_x1_vec[2]=-rb[1][0];
    corner_x1_vec[3]=-rb[2][0];
    corner_y1_vec[0]=0.0      ;
    corner_y1_vec[1]=-rb[0][1];
    corner_y1_vec[2]=-rb[1][1];
    corner_y1_vec[3]=-rb[2][1];
    corner_z1_vec[0]=0.0      ;
    corner_z1_vec[1]=-rb[0][2];
    corner_z1_vec[2]=-rb[1][2];
    corner_z1_vec[3]=-rb[2][2];
    // corners 4~7
    corner_x2_vec[0]=-(rb[0][0]+rb[1][0]);
    corner_x2_vec[1]=-(rb[0][0]+rb[2][0]);
    corner_x2_vec[2]=-(rb[1][0]+rb[2][0]);
    corner_x2_vec[3]=-(rb[0][0]+rb[1][0]+rb[2][0]);
    corner_y2_vec[0]=-(rb[0][1]+rb[1][1]);
    corner_y2_vec[1]=-(rb[0][1]+rb[2][1]);
    corner_y2_vec[2]=-(rb[1][1]+rb[2][1]);
    corner_y2_vec[3]=-(rb[0][1]+rb[1][1]+rb[2][1]);
    corner_z2_vec[0]=-(rb[0][2]+rb[1][2]);
    corner_z2_vec[1]=-(rb[0][2]+rb[2][2]);
    corner_z2_vec[2]=-(rb[1][2]+rb[2][2]);
    corner_z2_vec[3]=-(rb[0][2]+rb[1][2]+rb[2][2]);
#else
    g00=g(0);
    g10=g(3);
    g20=g(6);
    g01=g(1);
    g11=g(4);
    g21=g(7);
    g02=g(2);
    g12=g(5);
    g22=g(8);
#endif
  }

  /** apply BC to a displacement vector a and return the minimum-image distance
   * @param lat lattice
   * @param a displacement vector
   * @return the minimum-image distance
   */
#ifdef BGQPX
  inline T apply_bc(TinyVector<T,3>& displ)  const
  {
    vector4double disp0_vec = vec_splats(displ[0]);
    vector4double disp1_vec = vec_splats(displ[1]);
    vector4double disp2_vec = vec_splats(displ[2]);
    vector4double ar_vec;

    ar_vec=vec_mul (disp0_vec, g0);
    ar_vec=vec_madd(disp1_vec, g1, ar_vec);
    ar_vec=vec_madd(disp2_vec, g2, ar_vec);
    ar_vec=-vec_floor(ar_vec);

    vector4double ar0_vec = vec_splats(ar_vec[0]);
    vector4double ar1_vec = vec_splats(ar_vec[1]);
    vector4double ar2_vec = vec_splats(ar_vec[2]);

    vector4double displ_vec = {displ[0], displ[1], displ[2], 0.0};

    displ_vec=vec_madd(ar0_vec, rb0, displ_vec);
    displ_vec=vec_madd(ar1_vec, rb1, displ_vec);
    displ_vec=vec_madd(ar2_vec, rb2, displ_vec);

    displ[0]=displ_vec[0];
    displ[1]=displ_vec[1];
    displ[2]=displ_vec[2];

    vector4double r2_03_vec, r2_47_vec;

    vector4double x1_vec=vec_splats(displ_vec[0]);
    vector4double x2_vec=vec_splats(displ_vec[0]);
    vector4double y1_vec=vec_splats(displ_vec[1]);
    vector4double y2_vec=vec_splats(displ_vec[1]);
    vector4double z1_vec=vec_splats(displ_vec[2]);
    vector4double z2_vec=vec_splats(displ_vec[2]);

    x1_vec=vec_add(x1_vec, corner_x1_vec);
    y1_vec=vec_add(y1_vec, corner_y1_vec);
    z1_vec=vec_add(z1_vec, corner_z1_vec);
    r2_03_vec=vec_mul(x1_vec,x1_vec);
    r2_03_vec=vec_madd(y1_vec,y1_vec,r2_03_vec);
    r2_03_vec=vec_madd(z1_vec,z1_vec,r2_03_vec);

    x2_vec=vec_add(x2_vec, corner_x2_vec);
    y2_vec=vec_add(y2_vec, corner_y2_vec);
    z2_vec=vec_add(z2_vec, corner_z2_vec);
    r2_47_vec=vec_mul(x2_vec,x2_vec);
    r2_47_vec=vec_madd(y2_vec,y2_vec,r2_47_vec);
    r2_47_vec=vec_madd(z2_vec,z2_vec,r2_47_vec);

    T rmin2=r2_03_vec[0];
    int imin=0;
    for(int i=1; i<4; ++i)
    {
      if(r2_03_vec[i]<rmin2)
      {
        rmin2=r2_03_vec[i];
        imin=i;
      }
    }
    for(int i=0; i<4; ++i)
    {
      if(r2_47_vec[i]<rmin2)
      {
        rmin2=r2_47_vec[i];
        imin=i+4;
      }
    }
    displ+=corners[imin];
    return rmin2;
  }
#else
  inline T apply_bc(TinyVector<T,3>& displ)  const
  {
    //cart2unit
    TinyVector<T,3>
    ar(displ[0]*g00+displ[1]*g10+displ[2]*g20
       ,displ[0]*g01+displ[1]*g11+displ[2]*g21
       ,displ[0]*g02+displ[1]*g12+displ[2]*g22);
    ar[0]=-std::floor(ar[0]);
    ar[1]=-std::floor(ar[1]);
    ar[2]=-std::floor(ar[2]);
    displ+=ar[0]*rb[0]+ar[1]*rb[1]+ar[2]*rb[2];
    T rmin2=dot(displ,displ);
    int imin=0;
    for(int i=1; i<corners.size(); ++i)
    {
      TinyVector<T,3> tv=displ+corners[i];
      T r2=dot(tv,tv);
      if(r2<rmin2)
      {
        rmin2=r2;
        imin=i;
      }
    }
    if(imin>0)
      displ+=corners[imin];
    return rmin2;
  }
#endif

  inline void apply_bc(std::vector<TinyVector<T,3> >& dr
                       , std::vector<T>& r
                       , std::vector<T>& rinv) const
  {
    const int n=dr.size();
    for(int i=0; i<n; ++i)
      rinv[i]=apply_bc(dr[i]);
    simd::sqrt(&rinv[0],&r[0],n);
    simd::inv(&r[0],&rinv[0],n);
  }

  inline void apply_bc(std::vector<TinyVector<T,3> >& dr
                       , std::vector<T>& r) const
  {
    for(int i=0; i<dr.size(); ++i)
      r[i]=apply_bc(dr[i]);
  }

  inline void evaluate_rsquared(TinyVector<T,3>* restrict dr, T* restrict rr, int n)
  {
    for(int i=0; i<n; ++i)
      rr[i]=apply_bc(dr[i]);
  }
};


/** specialization for a slab, general cell
*/
template<class T>
struct DTD_BConds<T,3,PPNG>
{
  T g00,g10,g01,g11;
  TinyVector<TinyVector<T,3>,3> rb;
  std::vector<TinyVector<T,3> > corners;

  DTD_BConds(const CrystalLattice<T,3>& lat)
  {
    rb[0]=lat.a(0);
    rb[1]=lat.a(1);
    rb[2]=lat.a(2); //rb[2]=0.0;
    g00=lat.G(0);
    g10=lat.G(3);
    g01=lat.G(1);
    g11=lat.G(4);
    T minusone=-1.0;
    corners.resize(4);
    corners[0]=0.0;
    corners[1]=minusone*(rb[0]);
    corners[2]=minusone*(rb[1]);
    corners[3]=minusone*(rb[0]+rb[1]);
  }

  inline T apply_bc(TinyVector<T,3>& displ)  const
  {
    //cart2unit
    TinyVector<T,2>
    ar(displ[0]*g00+displ[1]*g10
       ,displ[0]*g01+displ[1]*g11);
    //put them in the box
    ar[0]-=std::floor(ar[0]);
    ar[1]-=std::floor(ar[1]);
    displ+=ar[0]*rb[0]+ar[1]*rb[1];
    T rmin2=dot(displ,displ);
    int imin=0;
    for(int i=1; i<corners.size(); ++i)
    {
      TinyVector<T,3> tv=displ+corners[i];
      T r2=dot(tv,tv);
      if(r2<rmin2)
      {
        rmin2=r2;
        imin=i;
      }
    }
    if(imin>0)
      displ+=corners[imin];
    return rmin2;
  }

  inline void apply_bc(std::vector<TinyVector<T,3> >& dr
                       , std::vector<T>& r
                       , std::vector<T>& rinv) const
  {
    const int n=dr.size();
    for(int i=0; i<n; ++i)
      rinv[i]=apply_bc(dr[i]);
    simd::sqrt(&rinv[0],&r[0],n);
    simd::inv(&r[0],&rinv[0],n);
  }

  inline void apply_bc(std::vector<TinyVector<T,3> >& dr
                       , std::vector<T>& r) const
  {
    for(int i=0; i<dr.size(); ++i)
      r[i]=apply_bc(dr[i]);
  }

  inline void evaluate_rsquared(TinyVector<T,3>* restrict dr, T* restrict rr, int n)
  {
    for(int i=0; i<n; ++i)
      rr[i]=apply_bc(dr[i]);
  }
};

/** specialization for a slab, orthorombic cell
*/
template<class T>
struct DTD_BConds<T,3,PPNO>
{

  T Linv0,L0,Linv1,L1;

  inline DTD_BConds(const CrystalLattice<T,3>& lat)
    : Linv0(lat.OneOverLength[0]), L0(lat.Length[0])
    ,Linv1(lat.OneOverLength[1]), L1(lat.Length[1])
  {
  }

  /** evaluate |a| and apply boundary conditions on a
   */
  inline T apply_bc(TinyVector<T,3>& displ)  const
  {
    T x=displ[0]*Linv0;
    displ[0]=L0*(x-round(x));
    T y=displ[1]*Linv1;
    displ[1]=L1*(y-round(y));
    return displ[0]*displ[0]+displ[1]*displ[1]+displ[2]*displ[2];
  }

  /** evaluate displacement data for a vector
   */
  inline void apply_bc(std::vector<TinyVector<T,3> >& dr
                       , std::vector<T>& r
                       , std::vector<T>& rinv) const
  {
    const int n=r.size();
    //use rinv as temporary rr
    for(int i=0; i<n; ++i)
      rinv[i]=apply_bc(dr[i]);
    simd::sqrt(&rinv[0],&r[0],n);
    simd::inv(&r[0],&rinv[0],n);
  }

  inline void apply_bc(std::vector<TinyVector<T,3> >& dr
                       , std::vector<T>& r) const
  {
    for(int i=0; i<dr.size(); ++i)
      r[i]=dot(dr[i],dr[i]);
  }

  inline void evaluate_rsquared(TinyVector<T,3>* restrict dr, T* restrict rr, int n)
  {
    for(int i=0; i<n; ++i)
      rr[i]=apply_bc(dr[i]);
  }
};

template<class T>
struct DTD_BConds<T,3,PPNS>
{
  T r00,r10,r01,r11;
  T g00,g10,g01,g11;
  DTD_BConds(const CrystalLattice<T,3>& lat)
    : r00(lat.R(0)),r10(lat.R(3))
    ,r01(lat.R(1)),r11(lat.R(4))
    ,g00(lat.G(0)),g10(lat.G(3))
    ,g01(lat.G(1)),g11(lat.G(4))
  { }

  /** apply BC to a displacement vector a and return the minimum-image distance
   * @param lat lattice
   * @param a displacement vector
   * @return the minimum-image distance
   */
  inline T apply_bc(TinyVector<T,3>& displ)  const
  {
    //cart2unit
    TinyVector<T,2>
    ar(displ[0]*g00+displ[1]*g10
       ,displ[0]*g01+displ[1]*g11);
    //put them in the box
    ar[0]-=round(ar[0]);
    ar[1]-=round(ar[1]);
    //unit2cart
    displ[0]=ar[0]*r00+ar[1]*r10;
    displ[1]=ar[0]*r01+ar[1]*r11;
    return displ[0]*displ[0]+displ[1]*displ[1]+displ[2]*displ[2];
  }

  inline void apply_bc(std::vector<TinyVector<T,3> >& dr
                       , std::vector<T>& r
                       , std::vector<T>& rinv) const
  {
    const int n=dr.size();
    for(int i=0; i<n; ++i)
      rinv[i]=apply_bc(dr[i]);
    simd::sqrt(&rinv[0],&r[0],n);
    simd::inv(&r[0],&rinv[0],n);
  }

  inline void apply_bc(std::vector<TinyVector<T,3> >& dr
                       , std::vector<T>& r) const
  {
    for(int i=0; i<dr.size(); ++i)
      r[i]=apply_bc(dr[i]);
  }

  inline void evaluate_rsquared(TinyVector<T,3>* restrict dr, T* restrict rr, int n)
  {
    for(int i=0; i<n; ++i)
      rr[i]=apply_bc(dr[i]);
  }
};


/** specialization for a wire
*/
template<class T>
struct DTD_BConds<T,3,SUPERCELL_WIRE>
{

  T Linv0,L0;

  inline DTD_BConds(const CrystalLattice<T,3>& lat)
    : Linv0(lat.OneOverLength[0]), L0(lat.Length[0])
  {
  }

  /** evaluate |a| and apply boundary conditions on a
   * @param lat lattice
   * @param a Cartesian vector
   * @return |a|^2
   */
  inline T apply_bc(TinyVector<T,3>& displ)  const
  {
    T x=displ[0]*Linv0;
    displ[0]=L0*(x-round(x));
    return displ[0]*displ[0]+displ[1]*displ[1]+displ[2]*displ[2];
  }

  /** evaluate displacement data for a vector
   */
  inline void apply_bc(std::vector<TinyVector<T,3> >& dr
                       , std::vector<T>& r
                       , std::vector<T>& rinv) const
  {
    const int n=r.size();
    //use rinv as temporary rr
    for(int i=0; i<n; ++i)
      rinv[i]=apply_bc(dr[i]);
    simd::sqrt(&rinv[0],&r[0],n);
    simd::inv(&r[0],&rinv[0],n);
  }

  inline void apply_bc(std::vector<TinyVector<T,3> >& dr
                       , std::vector<T>& r) const
  {
    for(int i=0; i<dr.size(); ++i)
      r[i]=dot(dr[i],dr[i]);
  }

  inline void evaluate_rsquared(TinyVector<T,3>* restrict dr, T* restrict rr, int n)
  {
    for(int i=0; i<n; ++i)
      rr[i]=apply_bc(dr[i]);
  }
};

/** specialization for a periodic 3D general cell
 *
 * Slow method and not used unless one needs to check if faster methods fail
 */
template<class T>
struct DTD_BConds<T,3,PPPX>
{
  T r00,r10,r20,r01,r11,r21,r02,r12,r22;
  T g00,g10,g20,g01,g11,g21,g02,g12,g22;
  T r2max;
  std::vector<TinyVector<T,3> > nextcells;

  DTD_BConds(const CrystalLattice<T,3>& lat)
    : r00(lat.R(0)),r10(lat.R(3)),r20(lat.R(6))
    ,r01(lat.R(1)),r11(lat.R(4)),r21(lat.R(7))
    ,r02(lat.R(2)),r12(lat.R(5)),r22(lat.R(8))
    ,g00(lat.G(0)),g10(lat.G(3)),g20(lat.G(6))
    ,g01(lat.G(1)),g11(lat.G(4)),g21(lat.G(7))
    ,g02(lat.G(2)),g12(lat.G(5)),g22(lat.G(8))
    ,r2max(lat.CellRadiusSq)
  {
    nextcells.resize(26);
    int ic=0;
    for(int i=-1; i<=1; ++i)
      for(int j=-1; j<=1; ++j)
        for(int k=-1; k<=1; ++k)
        {
          if(!(i || j || k ))
            continue;//exclude zero
          nextcells[ic][0]=i*r00+j*r10+k*r20;
          nextcells[ic][1]=i*r01+j*r11+k*r21;
          nextcells[ic][2]=i*r02+j*r12+k*r22;
          ++ic;
        }
  }

  /** evaluate the minimum distance
   * @param lat lattice
   * @param a displacement vector [-0.5,0.5)x[-0.5,0.5)x[-0.5,0.5)
   * @param r2max square of the maximum cutoff
   * @return square of the minimum-image distance
   *
   * Search the ghost cells to match Wigner-Seitz cell
   */
  inline T get_min_distance(TinyVector<T,3>& a) const
  {
    T d2 = a[0]*a[0]+a[1]*a[1]+a[2]*a[2];
    if (d2 < r2max)
      return d2;
    else
    {
      T d2min = d2;
      int ic=-1;
      for(int i=0; i<nextcells.size(); ++i)
      {
        TinyVector<T,3> c(a+nextcells[i]);
        d2=c[0]*c[0]+c[1]*c[1]+c[2]*c[2];
        if(d2<d2min)
        {
          d2min=d2;
          ic=i;
        }
      }
      if(ic>=0)
        a += nextcells[ic];
      return d2min;
    }
  }

  /** apply BC to a displacement vector a and return the minimum-image distance
   * @param lat lattice
   * @param a displacement vector
   * @return the minimum-image distance
   */
  inline T apply_bc(TinyVector<T,3>& displ)  const
  {
    //cart2unit
    TinyVector<T,3>
    ar(displ[0]*g00+displ[1]*g10+displ[2]*g20
       ,displ[0]*g01+displ[1]*g11+displ[2]*g21
       ,displ[0]*g02+displ[1]*g12+displ[2]*g22);
    //put them in the box
    ar[0]-=round(ar[0]);
    ar[1]-=round(ar[1]);
    ar[2]-=round(ar[2]);
    //unit2cart
    displ[0]=ar[0]*r00+ar[1]*r10+ar[2]*r20;
    displ[1]=ar[0]*r01+ar[1]*r11+ar[2]*r21;
    displ[2]=ar[0]*r02+ar[1]*r12+ar[2]*r22;
    //return |displ|^2 after checking the ghost cells
    return get_min_distance(displ);
  }

  /** out = prod (in ,lattice)
   * @param lattice 3x3 tensor to for conversion, either CrystalLattice::R or CrystalLattice::G
   * @param in start address of input vectors,  in[n][3]
   * @param out start address of output vectors,  out[n][3]
   * @param n number of 3d vectors
   */
  inline void convert2Cart(const T* restrict in, T* restrict out, int n) const
  {
    for(int i=0,i3=0; i<n; ++i,i3+=3)
    {
      out[i3  ]=in[i3]*r00+in[i3+1]*r10+in[i3+2]*r20;
      out[i3+1]=in[i3]*r01+in[i3+1]*r11+in[i3+2]*r21;
      out[i3+2]=in[i3]*r02+in[i3+1]*r12+in[i3+2]*r22;
    }
  }

  inline void convert2Unit(const T* restrict in, T* restrict out, int n) const
  {
    for(int i=0,i3=0; i<n; ++i,i3+=3)
    {
      out[i3  ]=in[i3]*g00+in[i3+1]*g10+in[i3+2]*g20;
      out[i3+1]=in[i3]*g01+in[i3+1]*g11+in[i3+2]*g21;
      out[i3+2]=in[i3]*g02+in[i3+1]*g12+in[i3+2]*g22;
    }
  }

  inline void apply_bc(std::vector<TinyVector<T,3> >& dr
                       , std::vector<T>& r
                       , std::vector<T>& rinv) const
  {
    const int n=dr.size();
    for(int i=0; i<n; ++i)
      rinv[i]=apply_bc(dr[i]);
    //using inline function but is not better
    //T drnew[n*3];
    //convert2Unit(&dr[0][0],drnew,n);
    //for(int i=0; i<n*3;++i) drnew[i]-= round(drnew[i]);
    //convert2Cart(drnew,&dr[0][0],n);
    //for(int i=0; i<n; ++i) rinv[i]=get_min_distance(dr[i]);
    simd::sqrt(&rinv[0],&r[0],n);
    simd::inv(&r[0],&rinv[0],n);
  }

  inline void apply_bc(std::vector<TinyVector<T,3> >& dr
                       , std::vector<T>& r) const
  {
    for(int i=0; i<dr.size(); ++i)
      r[i]=apply_bc(dr[i]);
  }

  inline void evaluate_rsquared(TinyVector<T,3>* restrict dr, T* restrict rr, int n)
  {
    for(int i=0; i<n; ++i)
      rr[i]=apply_bc(dr[i]);
  }
};

/** specialization for a slab, general cell
*/
template<class T>
struct DTD_BConds<T,3,PPNX>
{
  T r00,r10,r01,r11;
  T g00,g10,g01,g11;
  T r2max;
  std::vector<TinyVector<T,3> > nextcells;

  DTD_BConds(const CrystalLattice<T,3>& lat)
    : r00(lat.R(0)),r10(lat.R(3))
    ,r01(lat.R(1)),r11(lat.R(4))
    ,g00(lat.G(0)),g10(lat.G(3))
    ,g01(lat.G(1)),g11(lat.G(4))
    ,r2max(lat.CellRadiusSq)
  {
    nextcells.resize(8);
    int ic=0;
    for(int i=-1; i<=1; ++i)
      for(int j=-1; j<=1; ++j)
      {
        if(!(i || j))
          continue;//exclude zero
        nextcells[ic][0]=i*r00+j*r10;
        nextcells[ic][1]=i*r01+j*r11;
        nextcells[ic][2]=0;
        ++ic;
      }
  }

  /** evaluate the minimum distance
   * @param lat lattice
   * @param a displacement vector \f$[-0.5,0.5)\times [-0.5,0.5)\times [-\infty,\infty)\f$
   * @param r2max square of the maximum cutoff
   * @return square of the minimum-image distance
   *
   * Search the ghost cells to match Wigner-Seitz cell
   */
  inline T get_min_distance(TinyVector<T,3>& a) const
  {
    T d2 = a[0]*a[0]+a[1]*a[1]+a[2]*a[2];
    if (d2 < r2max)
      return d2;
    else
    {
      T d2min = d2;
      int ic=-1;
      for(int i=0; i<8; ++i)
      {
        TinyVector<T,3> c(a+nextcells[i]);
        d2=c[0]*c[0]+c[1]*c[1]+c[2]*c[2];
        if(d2<d2min)
        {
          d2min=d2;
          ic=i;
        }
      }
      if(ic>=0)
        a += nextcells[ic];
      return d2min;
    }
  }

  /** apply BC to a displacement vector a and return the minimum-image distance
   * @param lat lattice
   * @param a displacement vector
   * @return the minimum-image distance
   */
  inline T apply_bc(TinyVector<T,3>& displ)  const
  {
    //cart2unit
    TinyVector<T,2> ar(displ[0]*g00+displ[1]*g10, displ[0]*g01+displ[1]*g11);
    //put them in the box
    ar[0]-=round(ar[0]);
    ar[1]-=round(ar[1]);
    //unit2cart
    displ[0]=ar[0]*r00+ar[1]*r10;
    displ[1]=ar[0]*r01+ar[1]*r11;
    //return |displ|^2 after checking the ghost cells
    return get_min_distance(displ);
  }

  inline void apply_bc(std::vector<TinyVector<T,3> >& dr
                       , std::vector<T>& r
                       , std::vector<T>& rinv) const
  {
    const int n=dr.size();
    for(int i=0; i<n; ++i)
      rinv[i]=apply_bc(dr[i]);
    simd::sqrt(&rinv[0],&r[0],n);
    simd::inv(&r[0],&rinv[0],n);
  }

  inline void apply_bc(std::vector<TinyVector<T,3> >& dr
                       , std::vector<T>& r) const
  {
    for(int i=0; i<dr.size(); ++i)
      r[i]=apply_bc(dr[i]);
  }

  inline void evaluate_rsquared(TinyVector<T,3>* restrict dr, T* restrict rr, int n)
  {
    for(int i=0; i<n; ++i)
      rr[i]=apply_bc(dr[i]);
  }
};


}

#endif // OHMMS_PARTICLE_BCONDS_3D_H

