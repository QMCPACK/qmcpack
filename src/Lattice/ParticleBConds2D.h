//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    



#ifndef QMCPLUSPLUS_PARTICLE_BCONDS_2D_H
#define QMCPLUSPLUS_PARTICLE_BCONDS_2D_H

#include <config.h>
#include <Lattice/CrystalLattice.h>

namespace qmcplusplus
{
/** specialization for a periodic 2D cell
*/
template<class T>
struct DTD_BConds<T,2,SUPERCELL_BULK>
{
  T r00,r10,r01,r11;
  T g00,g10,g01,g11;
  T r2max;
  std::vector<TinyVector<T,2> > nextcells;

  DTD_BConds(const CrystalLattice<T,2>& lat)
    : r00(lat.R(0)),r10(lat.R(2))
    ,r01(lat.R(1)),r11(lat.R(3))
    ,g00(lat.G(0)),g10(lat.G(2))
    ,g01(lat.G(1)),g11(lat.G(3))
    ,r2max(lat.CellRadiusSq)
  {
    nextcells.resize(8);
    int ic=0;
    for(int i=-1; i<=1; ++i)
      for(int j=-1; j<=1; ++j)
      {
        if(!(i || j ))
          continue;//exclude zero
        nextcells[ic][0]=i*r00+j*r10;
        nextcells[ic][1]=i*r01+j*r11;
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
  inline T get_min_distance(TinyVector<T,2>& a) const
  {
    T d2 = a[0]*a[0]+a[1]*a[1];
    if (d2 < r2max)
      return d2;
    else
    {
      T d2min = d2;
      int ic=-1;
      for(int i=0; i<8; ++i)
      {
        TinyVector<T,2> c(a+nextcells[i]);
        d2=c[0]*c[0]+c[1]*c[1];
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
  inline T apply_bc(TinyVector<T,2>& displ)  const
  {
    //cart2unit
    TinyVector<T,2> ar(displ[0]*g00+displ[1]*g10 ,displ[0]*g01+displ[1]*g11);
    //put them in the box
    ar[0]-=round(ar[0]);
    ar[1]-=round(ar[1]);
    //unit2cart
    displ[0]=ar[0]*r00+ar[1]*r10;
    displ[1]=ar[0]*r01+ar[1]*r11;
    //return |displ|^2 after checking the ghost cells
    return get_min_distance(displ);
  }

  inline void apply_bc(std::vector<TinyVector<T,2> >& dr
                       , std::vector<T>& r
                       , std::vector<T>& rinv) const
  {
    const int n=dr.size();
    for(int i=0; i<n; ++i)
      rinv[i]=apply_bc(dr[i]);
    simd::sqrt(&rinv[0],&r[0],n);
    simd::inv(&r[0],&rinv[0],n);
  }

  inline void apply_bc(std::vector<TinyVector<T,2> >& dr
                       , std::vector<T>& r) const
  {
    for(int i=0; i<dr.size(); ++i)
      r[i]=apply_bc(dr[i]);
  }
};

/** specialization for a periodic 2D general cell with wigner-seitz==simulation cell
 *
 * Skip image cells.
 */
template<class T>
struct DTD_BConds<T,2,PPPS>
{
  T r00,r10,r01,r11;
  T g00,g10,g01,g11;
  T r2max;
  std::vector<TinyVector<T,2> > nextcells;

  DTD_BConds(const CrystalLattice<T,2>& lat)
    : r00(lat.R(0)),r10(lat.R(2))
    ,r01(lat.R(1)),r11(lat.R(3))
    ,g00(lat.G(0)),g10(lat.G(2))
    ,g01(lat.G(1)),g11(lat.G(3))
    ,r2max(lat.CellRadiusSq)
  { }

  /** apply BC to a displacement vector a and return the minimum-image distance
   * @param lat lattice
   * @param a displacement vector
   * @return the minimum-image distance
   */
  inline T apply_bc(TinyVector<T,2>& displ)  const
  {
    //cart2unit
    TinyVector<T,2> ar(displ[0]*g00+displ[1]*g10 ,displ[0]*g01+displ[1]*g11);
    //put them in the box
    ar[0]-=round(ar[0]);
    ar[1]-=round(ar[1]);
    //unit2cart
    displ[0]=ar[0]*r00+ar[1]*r10;
    displ[1]=ar[0]*r01+ar[1]*r11;
    return displ[0]*displ[0]+displ[1]*displ[1];
  }

  inline void apply_bc(std::vector<TinyVector<T,2> >& dr
                       , std::vector<T>& r
                       , std::vector<T>& rinv) const
  {
    const int n=dr.size();
    for(int i=0; i<n; ++i)
      rinv[i]=apply_bc(dr[i]);
    simd::sqrt(&rinv[0],&r[0],n);
    simd::inv(&r[0],&rinv[0],n);
  }

  inline void apply_bc(std::vector<TinyVector<T,2> >& dr
                       , std::vector<T>& r) const
  {
    for(int i=0; i<dr.size(); ++i)
      r[i]=apply_bc(dr[i]);
  }
};

/** specialization for a periodic 2D orthorombic cell
*/
template<class T>
struct DTD_BConds<T,2,PPPO>
{

  T Linv0,L0,Linv1,L1;

  inline DTD_BConds(const CrystalLattice<T,2>& lat)
    : Linv0(lat.OneOverLength[0]), L0(lat.Length[0])
    ,Linv1(lat.OneOverLength[1]), L1(lat.Length[1])
  {
  }

  /** evaluate |a| and apply boundary conditions on a
   * @param lat lattice
   * @param a Cartesian vector
   * @return |a|^2
   */
  inline T apply_bc(TinyVector<T,2>& displ)  const
  {
    T x=displ[0]*Linv0;
    displ[0]=L0*(x-round(x));
    T y=displ[1]*Linv1;
    displ[1]=L1*(y-round(y));
    return displ[0]*displ[0]+displ[1]*displ[1];
  }

  /** evaluate displacement data for a vector
   */
  inline void apply_bc(std::vector<TinyVector<T,2> >& dr
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

  inline void apply_bc(std::vector<TinyVector<T,2> >& dr
                       , std::vector<T>& r) const
  {
    for(int i=0; i<dr.size(); ++i)
      r[i]=dot(dr[i],dr[i]);
  }
};

/** specialization for a wire in 2D
*/
template<class T>
struct DTD_BConds<T,2,SUPERCELL_WIRE>
{

  T Linv0,L0;

  inline DTD_BConds(const CrystalLattice<T,2>& lat)
    : Linv0(lat.OneOverLength[0]), L0(lat.Length[0])
  { }

  /** evaluate |a| and apply boundary conditions on a
   * @param lat lattice
   * @param a Cartesian vector
   * @return |a|^2
   */
  inline T apply_bc(TinyVector<T,2>& displ)  const
  {
    T x=displ[0]*Linv0;
    displ[0]=L0*(x-round(x));
    return displ[0]*displ[0]+displ[1]*displ[1];
  }

  /** evaluate displacement data for a vector
   */
  inline void apply_bc(std::vector<TinyVector<T,2> >& dr
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

  inline void apply_bc(std::vector<TinyVector<T,2> >& dr
                       , std::vector<T>& r) const
  {
    for(int i=0; i<dr.size(); ++i)
      r[i]=dot(dr[i],dr[i]);
  }
};

}
#endif // OHMMS_PARTICLE_BCONDS_2D_H


