//////////////////////////////////////////////////////////////////
// (c) Copyright 2011- by Jeongnim Kim
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//   Jeongnim Kim
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
#ifndef QMCPLUSPLUS_FAST_MATRIX_HANLDERS_HPP
#define QMCPLUSPLUS_FAST_MATRIX_HANLDERS_HPP

#define USE_DIRAC_FAST_OPERATORS

namespace qmcplusplus {
  template<typename T, unsigned D>
    struct PartialTrace
    {
      static inline void trace(T* restrict target, const T* restrict a, const T* restrict b,int n, int m)
      {
        for(int i=0,ij=0; i<n; ++i)
        {
          T res=T();
          for(int j=0; j<m; ++j,++ij) res += a[ij]*b[ij];
          target[i]=res;
        }
      }

      static inline void traceG(T* restrict target, const T* restrict a, const T* restrict b,int n, int m)
      {
        T res[D];
        for(int i=0,ij=0; i<n*D;)
        {
          for(int k=0; k<D; ++k)
          {
            res[k]=T();
            for(int j=0; j<m; ++j,++ij) res[k] += a[ij]*b[ij];
          }
          for(int k=0; k<D; ++k) target[i++]=res[k];
        }
      }
    };

  template<typename T>
    struct PartialTrace<T,3>
    {
      static inline void traceS(T* restrict target, const T* restrict a, const T* restrict b,int n, int m)
      {
        for(int i=0,ij=0; i<n; ++i)
        {
          T res=T();
          for(int j=0; j<m; ++j,++ij) res += a[ij]*b[ij];
          target[i]=res;
        }
      }

      static inline void traceG(T* restrict target, const T* restrict a, const T* restrict b,int n, int m)
      {
        T v[4];
        for(int i=0,ij=0; i<n;++i)
        {
          v[0]=0;v[1]=0;v[2]=0;
          //v[0]=T();v[1]=T();v[2]=T();
          for(int j=0; j<m; ++j,++ij)
          {
            v[0]+=a[ij]*b[3*ij  ];
            v[1]+=a[ij]*b[3*ij+1];
            v[2]+=a[ij]*b[3*ij+2];
          }
          target[3*i  ]=v[0];
          target[3*i+1]=v[1];
          target[3*i+2]=v[2];
        }
      }

      static inline void computeGL(T* restrict gtarget, T* restrict ltarget
          , const T* restrict psiinv, const T* restrict dpsi, const T* restrict d2psi
          , int n, int m)
      {
        register T v[4];
        for(int i=0,ij=0; i<n;++i)
        {
          v[0]=0;v[1]=0;v[2]=0;v[3]=0;
          //v[0]=T();v[1]=T();v[2]=T();
          for(int j=0; j<m; ++j,++ij)
          {
            register T t=psiinv[ij];
            v[0]+=t*dpsi[3*ij  ];
            v[1]+=t*dpsi[3*ij+1];
            v[2]+=t*dpsi[3*ij+2];
            v[3]+=t*d2psi[ij];
          }
          gtarget[3*i  ]=v[0];
          gtarget[3*i+1]=v[1];
          gtarget[3*i+2]=v[2];
          ltarget[i]=v[3]-v[0]*v[0]-v[1]*v[1]-v[2]*v[2];
        }
      }
    };
}
#endif
/***************************************************************************
 * $RCSfile$   $Author: jmcminis $
 * $Revision: 5077 $   $Date: 2010-12-09 03:14:51 -0600 (Thu, 09 Dec 2010) $
 * $Id: matrixoperators.h 5077 2010-12-09 09:14:51Z jmcminis $ 
 ***************************************************************************/
