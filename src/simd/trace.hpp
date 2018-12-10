//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    


#ifndef QMCPLUSPLUS_VECTOROPERATORS_3D_HPP
#define QMCPLUSPLUS_VECTOROPERATORS_3D_HPP

#include <OhmmsPETE/TinyVector.h>
#include <OhmmsPETE/OhmmsMatrix.h>

namespace qmcplusplus {
  namespace simd {

    /** evaluate gradients and laplacians in determinant classes
     * @param psi starting address of a matrix(n,m)
     * @param psi_g starting address of a gradient, TinyVector<T,D>, matrix (n,m)
     * @param g \f$g_i = \sum_{j} \psi(i,j)\times \psi_g(i,j)\f$
     * @param psi_l starting address of a laplacian matrix (n,m)
     * @param l \f$l_i = \sum_{j} \psi(i,j)\times \psi_l(i,j)-dot(g_i,g_i)\f$
     * @param n rows
     * @param m columns
     *
     * Typically, the index i denotes the particle index with an offset.
     */
    template<typename T, unsigned D>
      inline void evaluate_gl(const T* restrict psi 
          , const TinyVector<T,D>* restrict psi_g , TinyVector<T,D>* restrict g 
          , const T* restrict psi_l, T* restrict l
          , int n, int m)
      {
        for(int i=0; i<n; ++i) 
        {
          g[i]=dot(psi,psi_g,m);
          l[i]=dot(psi,psi_l,m)-dot(g[i],g[i]);
          psi+=m;
          psi_g+=m;
          psi_l+=m;
        }
      }

    /** evaluate gradients and laplacians in determinant classes
     * @param psi starting address of a matrix(n,m)
     * @param psi_gl starting address of a gradient, TinyVector<T,D+1>, matrix (n,m)
     * @param g \f$g_i = \sum_{j} \psi(i,j)\times \psi_gl(i,j)\f$
     * @param l \f$l_i = \sum_{j} \psi(i,j)\times \psi_l(i,j)-dot(g_i,g_i)\f$
     * @param n rows
     * @param m columns
     *
     * Typically, the index i denotes the particle index with an offset.
     */
    template<typename T, unsigned D>
      inline void evaluate_gl(const T* restrict psi , const TinyVector<T,D+1>* restrict psi_gl
          , TinyVector<T,D>* restrict g , T* restrict l , int n, int m)
      {
        register TinyVector<T,D+1> res;
        for(int i=0; i<n; ++i) 
        {
          res=dot(psi,psi_gl,m);
          T sum=0.0;
          for(int k=0;k<D; k) {g[i][k]=res[k]; sum += res[k]*res[k];}
          l[i]=res[D]-sum;
          psi+=m;
          psi_gl+=m;
        }
      }

    /** evaluate gradients and laplacians in determinant classes specialized for 3D
     */
    template<typename T>
      inline void evaluate_gl(const T* restrict psi , const TinyVector<T,4>* restrict psi_gl
          , TinyVector<T,3>* restrict g , T* restrict l , int n, int m)
      {
        register TinyVector<T,4> res;
        for(int i=0; i<n; ++i) 
        {
          res=dot(psi,psi_gl,m);
          g[i][0]=res[0]; g[i][1]=res[1]; g[i][2]=res[2];
          l[i]=res[3]-res[0]*res[0]-res[1]*res[1]-res[2]*res[2];
          psi+=m;
          psi_gl+=m;
        }
      }

    template<typename T>
      inline void trace(const Matrix<T>& c, const Matrix<T>& x, T* restrict y)
      {
        const int m=c.cols();
        for(int i=0; i<c.rows(); ++i) y[i]=dot(c[i],x[i],m);
      }

    template<typename T, unsigned D>
      inline void trace(const Matrix<T>& c, const Matrix<TinyVector<T,D> >& x, TinyVector<T,D>* restrict y)
      {
        const int m=c.cols();
        for(int i=0; i<c.rows(); ++i) y[i]=dot(c[i],x[i],m);
      }

    template<typename T, unsigned D>
      inline void trace(const Matrix<T>& psiinv 
          , const Matrix<TinyVector<T,D> >& psi_g , TinyVector<T,D>* restrict g 
          , const Matrix<T>& psi_l, T* restrict l)
      {
        trace(psiinv.data(), psi_g.data(), g, psi_l.data(),l,psiinv.rows(),psiinv.cols());
      }
  }
}
#endif
