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

#include "Lattice/CrystalLattice.h"
#include <limits>

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

  /** common definition of a distance 
   *
   * @param lat supercell
   * @param a in/out displacement
   * @return dot(a,a)
   */
  template<class T, unsigned D, int SC>
    struct DTD_BConds 
    {
      inline static T apply(const CrystalLattice<T,D>& lat, TinyVector<T,D>& a) {
        return dot(a,a);
      }
    };

  /** specialization for a periodic 3D cell
  */
  template<class T>
    struct DTD_BConds<T,3,SUPERCELL_BULK> {
      inline static T apply(const CrystalLattice<T,3>& lat, TinyVector<T,3>& a) {
        TinyVector<T,3> ar(lat.toUnit(a));
        /*
           if(ar[0]<-0.5) ar[0]+=1.0; 
           else if(ar[0]>=0.5) ar[0]-=1.0;
           if(ar[1]<-0.5) ar[1]+=1.0; 
           else if(ar[1]>=0.5) ar[1]-=1.0;
           if(ar[2]<-0.5) ar[2]+=1.0; 
           else if(ar[2]>=0.5) ar[2]-=1.0;
           */
        //T x=fmod(ar[0],1.0); ar[0]=x-static_cast<int>(x*2.0);
        //T y=fmod(ar[1],1.0); ar[1]=y-static_cast<int>(y*2.0);
        //T z=fmod(ar[2],1.0); ar[2]=z-static_cast<int>(z*2.0);
#if defined(HAVE_STD_ROUND)
        ar[0]=ar[0]-round(ar[0]);
        ar[1]=ar[1]-round(ar[1]);
        ar[2]=ar[2]-round(ar[2]);
#else
        T dmy0,dmy1,dmy2;
        T x=modf(ar[0],&dmy0); ar[0]=x-static_cast<int>(x*2.0);
        T y=modf(ar[1],&dmy1); ar[1]=y-static_cast<int>(y*2.0);
        T z=modf(ar[2],&dmy2); ar[2]=z-static_cast<int>(z*2.0);
#endif
        a=lat.toCart(ar);
#if defined(DISABLE_WS_CELL)
        return a[0]*a[0]+a[1]*a[1]+a[2]*a[2];
#else
        T d2 = a[0]*a[0]+a[1]*a[1]+a[2]*a[2];
        if (d2 < lat.SimulationCellRadius * lat.SimulationCellRadius)
          return d2;
        else 
        {
          T d2min = d2;
          TinyVector<T,3> u = ar;
          for (double i=-1.0; i<1.001; i+=1.0) {
            u[0] = ar[0] + i;
            for (double j=-1.0; j<1.001; j+=1.0) {
              u[1] = ar[1] + j;
              for (double k=-1.0; k<1.001; k+=1.0) {
                u[2] = ar[2] + k;
                TinyVector<T,3> anew = lat.toCart(u);
                d2 = anew[0]*anew[0] + anew[1]*anew[1] + anew[2]*anew[2];
                if (d2 < d2min) {
                  a = anew;
                  d2min = d2;
                }
              }
            }
          }
          return d2min;
        }
#endif
      }
    };

  /** specialization for a periodic 3D orthorombic cell
  */
  template<class T>
    struct DTD_BConds<T,3,SUPERCELL_BULK+TwoPowerD> {
      inline static T apply(const CrystalLattice<T,3>& lat, TinyVector<T,3>& a) 
      {
#if defined(HAVE_STD_ROUND)
        T x=a[0]*lat.OneOverLength[0]; a[0]=lat.Length[0]*(x-round(x));
        T y=a[1]*lat.OneOverLength[1]; a[1]=lat.Length[1]*(y-round(y));
        T z=a[2]*lat.OneOverLength[2]; a[2]=lat.Length[2]*(z-round(z));
#else
        T dmy0,dmy1,dmy2;
        T x=modf(a[0]*lat.OneOverLength[0],&dmy0); a[0]=lat.Length[0]*(x-static_cast<int>(x*2.0));
        T y=modf(a[1]*lat.OneOverLength[1],&dmy1); a[1]=lat.Length[1]*(y-static_cast<int>(y*2.0));
        T z=modf(a[2]*lat.OneOverLength[2],&dmy2); a[2]=lat.Length[2]*(z-static_cast<int>(z*2.0));
#endif
        return a[0]*a[0]+a[1]*a[1]+a[2]*a[2];
      }
    };

  /** specialization for a periodic 2D cell
  */
  template<class T>
    struct DTD_BConds<T,2,SUPERCELL_BULK> {
      inline static T apply(const CrystalLattice<T,2>& lat, TinyVector<T,2>& a) {
        TinyVector<T,2> ar(lat.toUnit(a));
#if defined(HAVE_STD_ROUND)
        ar[0]=ar[0]-round(ar[0]);
        ar[1]=ar[1]-round(ar[1]);
#else
        T dmy0,dmy1;
        T x=modf(ar[0],&dmy0); ar[0]=x-static_cast<int>(x*2.0);
        T y=modf(ar[1],&dmy1); ar[1]=y-static_cast<int>(y*2.0);
#endif
        a=lat.toCart(ar);
        return a[0]*a[0]+a[1]*a[1];
      }
    };

  /** specialization for a periodic 2D orthorombic cell
  */
  template<class T>
    struct DTD_BConds<T,2,SUPERCELL_BULK+TwoPowerD> {
      inline static T apply(const CrystalLattice<T,2>& lat, TinyVector<T,2>& a) 
      {
#if defined(HAVE_STD_ROUND)
        T x=a[0]*lat.OneOverLength[0]; a[0]=lat.Length[0]*(x-round(x));
        T y=a[1]*lat.OneOverLength[1]; a[1]=lat.Length[1]*(y-round(y));
#else
        T dmy0,dmy1;
        T x=modf(a[0]*lat.OneOverLength[0],&dmy0); a[0]=lat.Length[0]*(x-static_cast<int>(x*2.0));
        T y=modf(a[1]*lat.OneOverLength[1],&dmy1); a[1]=lat.Length[1]*(y-static_cast<int>(y*2.0));
#endif
        return a[0]*a[0]+a[1]*a[1];
      }
    };

  ///**@ingroup nnlist
  // * @brief class to apply No Boundary conditions for distance evaluation
  // *
  // *NoBConds stands for No Boundary Conditions and is intended for
  // finite systems with open or vanishing boundary conditions.
  // Use a simple dot product assuming cartesian coordinates
  // */
  //template<class T, unsigned D>
  //  struct NoBConds {
  //    inline static T apply(const CrystalLattice<T,D>& lat, TinyVector<T,D>& a) {
  //      return dot(a,a);
  //    }
  //  };

  ///**@ingroup nnlist
  // * @brief class to apply Periodic Boundary conditions for distance evaluation
  // *
  // *PeriodicBConds stands for Periodic Boundary Conditions and is intended 
  // for periodic systems. The Cartesian distances are evaluated
  // according to the minimum-image convention.
  // */
  //template<class T, unsigned D>
  //  struct PeriodicBConds {
  //    inline static T apply(const CrystalLattice<T,D>& lat, TinyVector<T,D>& a) {
  //      TinyVector<T,D> ar(lat.toUnit(a));
  //      for(int idim=0; idim<D; idim++) {
  //        if(lat.BoxBConds[idim]) {
  //          if(ar[idim]<-0.5) ar[idim]+=1.0; 
  //          else if(ar[idim]>=0.5) ar[idim]-=1.0;
  //        }
  //      }
  //      a=lat.toCart(ar);
  //      return dot(a,a);
  //    }
  //  };

  ///**generic ParticleBConds
  //*/
  //template<class T, unsigned D> struct ParticleBConds {};

  ///**class ParticleBConds<T,3> 
  // *brief specialization for 3-dimensional supercells
  // */
  //template<class T>
  //  class ParticleBConds<T,3> 
  //  {

  //    public:

  //      // typedef for a pointer to boundary condition function
  //      typedef T (*ParticleBCond)(const T, const T, const T);

  //    public:
  //      // constructor: initialize all BC's to periodic ones
  //      ParticleBConds() { }

  //      ///apply -> displacement
  //      inline TinyVector<T,3> apply(const TinyVector<T,3>& x) const {
  //        TinyVector<T,3> dr=x;
  //        if(dr[0]<-0.5) 
  //          dr[0]+=1.0;
  //        else if(dr[0]>=0.5) 
  //          dr[0]-=1.0;
  //        if(dr[1]<-0.5) 
  //          dr[1]+=1.0;
  //        else if(dr[1]>=0.5) 
  //          dr[1]-=1.0;
  //        if(dr[2]<-0.5) 
  //          dr[2]+=1.0;
  //        else if(dr[2]>=0.5) 
  //          dr[2]-=1.0;
  //        return dr;
  //      }

  //      ///wrap -> apply
  //      inline TinyVector<T,3> wrap(const TinyVector<T,3>& rin) const {
  //        register T x(rin[0]),y(rin[1]),z(rin[2]);
  //        const T epsilon = -std::numeric_limits<T>::epsilon();
  //        const T plus_one = 1.0;
  //        if(x<epsilon) x+=plus_one;
  //        else if(x>=plus_one) x-=plus_one;
  //        if(y<epsilon) y +=plus_one;
  //        else if(y>=plus_one) y-= plus_one;
  //        if(z<epsilon) z +=plus_one;
  //        else if(z >=plus_one) z -= plus_one;
  //        return TinyVector<T,3>(x,y,z);
  //      }

  //      ///applyBC
  //      inline void applyBC(const TinyVector<T,3>& rin, TinyVector<T,3>& rout) const {
  //        const T epsilon = -std::numeric_limits<T>::epsilon();
  //        const T plus_one = 1.0;
  //        rout=rin;
  //        if(rout[0]<epsilon)        rout[0] += plus_one;
  //        else if(rout[0]>=plus_one) rout[0] -= plus_one;
  //        if(rout[1]<epsilon)        rout[1] += plus_one;
  //        else if(rout[1]>=plus_one) rout[1] -= plus_one;
  //        if(rout[2]<epsilon)        rout[2] += plus_one;
  //        else if(rout[2]>=plus_one) rout[2] -= plus_one;
  //      }
  //  };

}
#endif // OHMMS_PARTICLE_BCONDS_H


/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
