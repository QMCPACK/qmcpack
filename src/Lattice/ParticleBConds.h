//////////////////////////////////////////////////////////////////
// (c) Copyright 2003- by Jeongnim Kim and Kenneth P. Esler
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//   National Center for Supercomputing Applications &
//   Materials Computation Center
//   University of Illinois, Urbana-Champaign
//   Urbana, IL 61801
//   e-mail: jnkim@ncsa.uiuc.edu
//
// Supported by 
//   National Center for Supercomputing Applications, UIUC
//   Materials Computation Center, UIUC
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
#ifndef QMCPLUSPLUS_PARTICLE_BCONDS_H
#define QMCPLUSPLUS_PARTICLE_BCONDS_H

#include <config.h>
#include <Lattice/CrystalLattice.h>

namespace APPNAMESPACE 
{
  template<int N,unsigned D>
    struct PowerOfN
    {
      enum {value=N*PowerOfN<N,D-1>::value};
    };

  template<int N>
    struct PowerOfN<N,0>
    {
      enum {value=1};
    };

  const int TwoPowerD=PowerOfN<2,OHMMS_DIM>::value;

  /** generic Boundary condition handler
   *
   * Implements method for any dimension with open boundary condition.
   * Specialization of DTD_BConds should implement
   * - apply_bc(TinyVector<T,D>& displ): apply BC on displ, Cartesian displacement vector, and returns |displ|^2
   * - apply_bc(dr,r,rinv): apply BC on displacements
   * - apply_bc(dr,r): apply BC without inversion calculations
   */
  template<class T, unsigned D, int SC>
    struct DTD_BConds 
    {

      /** constructor: doing nothing */
      inline DTD_BConds(const CrystalLattice<T,D>& lat){}

      /** apply BC on displ and return |displ|^2 
       * @param displ a displacement vector in the Cartesian coordinate
       * @return \f$|displ|^2\f$
       */
      inline T apply_bc(TinyVector<T,D>& displ) const
      {
        return dot(displ,displ);
      }

      /** apply BC on dr and evaluate r and rinv
       * @param dr vector of displacements, in and out
       * @param r vector of distances
       * @param rinv vector of 1/r
       *
       * The input displacement vectors are not modified with the open boundary conditions.
       */
      inline void apply_bc(std::vector<TinyVector<T,D> >& dr
          , std::vector<T>& r
          , std::vector<T>& rinv) const
      {
        const int n=dr.size();
        for(int i=0; i<n; ++i) rinv[i]=dot(dr[i],dr[i]);
        vec_sqrt(n,&rinv[0],&r[0]);
        vec_inv(n,&r[0],&rinv[0]);
      }

      inline void apply_bc(std::vector<TinyVector<T,3> >& dr
          , std::vector<T>& r) const
      {
        for(int i=0;i<dr.size();++i) r[i]=std::sqrt(dot(dr[i],dr[i]));
      }
    };

  /** specialization for a periodic 3D cell
  */
  template<class T>
    struct DTD_BConds<T,3,SUPERCELL_BULK> 
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
        for(int i=-1;i<=1;++i)
          for(int j=-1;j<=1;++j)
            for(int k=-1;k<=1;++k)
            {
              if(!(i || j || k )) continue;//exclude zero
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
          for(int i=0; i<26; ++i)
          {
            TinyVector<T,3> c(a+nextcells[i]);
            d2=c[0]*c[0]+c[1]*c[1]+c[2]*c[2];
            if(d2<d2min) { d2min=d2;ic=i;}
          }
          if(ic>=0) a += nextcells[ic];
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
        ar[0]-=round(ar[0]); ar[1]-=round(ar[1]); ar[2]-=round(ar[2]);

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
        for(int i=0;i<n;++i) rinv[i]=apply_bc(dr[i]);
        //using inline function but is not better
        //T drnew[n*3];
        //convert2Unit(&dr[0][0],drnew,n);
        //for(int i=0; i<n*3;++i) drnew[i]-= round(drnew[i]);
        //convert2Cart(drnew,&dr[0][0],n);
        //for(int i=0; i<n; ++i) rinv[i]=get_min_distance(dr[i]);
        vec_sqrt(n,&rinv[0],&r[0]);
        vec_inv(n,&r[0],&rinv[0]);
      }

      inline void apply_bc(std::vector<TinyVector<T,3> >& dr
          , std::vector<T>& r) const
      {
        for(int i=0;i<dr.size();++i) r[i]=apply_bc(dr[i]);
      }
    };

  /** specialization for a periodic 3D orthorombic cell
  */
  template<class T>
    struct DTD_BConds<T,3,SUPERCELL_BULK+TwoPowerD> 
    {

      T Linv0,L0,Linv1,L1,Linv2,L2,r2max,dummy;

      inline DTD_BConds(const CrystalLattice<T,3>& lat)
        : Linv0(lat.OneOverLength[0]), L0(lat.Length[0])
          ,Linv1(lat.OneOverLength[1]), L1(lat.Length[1])
          ,Linv2(lat.OneOverLength[2]), L2(lat.Length[2])
          ,r2max(lat.CellRadiusSq),dummy(T())
      {}

      /** evaluate |a| and apply boundary conditions on a
       * @param lat lattice 
       * @param a Cartesian vector
       * @return |a|^2
       */
      inline T apply_bc(TinyVector<T,3>& displ)  const
      {
        T x=displ[0]*Linv0; displ[0]=L0*(x-round(x));
        T y=displ[1]*Linv1; displ[1]=L1*(y-round(y));
        T z=displ[2]*Linv2; displ[2]=L2*(z-round(z));
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
        for(int i=0; i<n; ++i) rinv[i]=apply_bc(dr[i]);
        vec_sqrt(n,&rinv[0],&r[0]);
        vec_inv(n,&r[0],&rinv[0]);
      }

      inline void apply_bc(std::vector<TinyVector<T,3> >& dr
          , std::vector<T>& r) const
      {
        for(int i=0;i<dr.size();++i) r[i]=dot(dr[i],dr[i]);
      }
    };

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
        for(int i=-1;i<=1;++i)
          for(int j=-1;j<=1;++j)
            {
              if(!(i || j )) continue;//exclude zero
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
            if(d2<d2min) { d2min=d2;ic=i;}
          }
          if(ic>=0) a += nextcells[ic];
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
        TinyVector<T,3> 
          ar(displ[0]*g00+displ[1]*g10
            ,displ[0]*g01+displ[1]*g11);

        //put them in the box
        ar[0]-=round(ar[0]); ar[1]-=round(ar[1]);

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
        for(int i=0;i<n;++i) rinv[i]=apply_bc(dr[i]);
        vec_sqrt(n,&rinv[0],&r[0]);
        vec_inv(n,&r[0],&rinv[0]);
      }

      inline void apply_bc(std::vector<TinyVector<T,2> >& dr
          , std::vector<T>& r) const
      {
        for(int i=0;i<dr.size();++i) r[i]=apply_bc(dr[i]);
      }
    };

  /** specialization for a periodic 2D orthorombic cell
  */
  template<class T>
    struct DTD_BConds<T,2,SUPERCELL_BULK+TwoPowerD> 
    {

      T Linv0,L0,Linv1,L1,r2max,dummy;

      inline DTD_BConds(const CrystalLattice<T,2>& lat)
        : Linv0(lat.OneOverLength[0]), L0(lat.Length[0])
          ,Linv1(lat.OneOverLength[1]), L1(lat.Length[1])
          ,r2max(lat.CellRadiusSq),dummy(T())
      {}

      /** evaluate |a| and apply boundary conditions on a
       * @param lat lattice 
       * @param a Cartesian vector
       * @return |a|^2
       */
      inline T apply_bc(TinyVector<T,2>& displ)  const
      {
        T x=displ[0]*Linv0; displ[0]=L0*(x-round(x));
        T y=displ[1]*Linv1; displ[1]=L1*(y-round(y));
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
        for(int i=0; i<n; ++i) rinv[i]=apply_bc(dr[i]);
        vec_sqrt(n,&rinv[0],&r[0]);
        vec_inv(n,&r[0],&rinv[0]);
      }

      inline void apply_bc(std::vector<TinyVector<T,2> >& dr
          , std::vector<T>& r) const
      {
        for(int i=0;i<dr.size();++i) r[i]=dot(dr[i],dr[i]);
      }
    };

}
#endif // OHMMS_PARTICLE_BCONDS_H


/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
