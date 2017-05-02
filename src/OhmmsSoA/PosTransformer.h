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
// -*- C++ -*-
/** @file VectorOperators.h
 * @brief Support funtions to handle position type data manged by soa
 */
#ifndef QMCPLUSPLUS_SOA_FAST_PARTICLE_OPERATORS_H
#define QMCPLUSPLUS_SOA_FAST_PARTICLE_OPERATORS_H

#include <simd/blas1.hpp>

namespace qmcplusplus
{
  //Need to reorg
#if 0
  /** Dummy template class to be specialized
   *
   * - T1 the datatype to be transformed
   * - D dimension
   * - ORTHO true, if only Diagonal Elements are used
   */
  template<class T1, unsigned D, bool ORTHO> struct PosTransformer { };

  /** Specialized PosTransformer<T,3,true> using only the diagonal elements
  */
  template<class T>
    struct PosTransformer<T,3,true>
    {
      using Array_t=VectorSoaContainer<T,3>;
      using Transformer_t=Tensor<T,3>;

      //index for the tensor
      enum {iXX=0, iXY=1, iXZ=2, iYX=3, iYY=4, iYZ=5, iZX=6, iZY=7, iZZ=8};

      inline static void
        apply(const Array_t& pin, const Transformer_t& X, Array_t& pout, int first, int last)
        {
          const int n=last-first;
          blas::axpy(X[iXX],pin.data(0),pout.data(0),n);
          blas::axpy(X[iYY],pin.data(1),pout.data(1),n);
          blas::axpy(X[iZZ],pin.data(2),pout.data(2),n);
        }

      inline static void
        apply(const Transformer_t& X, const Array_t& pin,  Array_t& pout, int first, int last)
        {
          ::apply(pin,X,pout,first,last);
        }

      inline static void
        apply(Array_t& pinout, const Transformer_t& X,int first, int last)
        {
          const int n=last-first;
          blas::scal(X[iXX],pinout.data(0),n);
          blas::scal(X[iYY],pinout.data(1),n);
          blas::scal(X[iZZ],pinout.data(2),n);
        }

      inline static void
        apply(const Transformer_t& X, Array_t& pinout, int first, int last)
        {
          ::apply(pinout,X,first,last);
        }
    };

  template<class T>
    struct PosTransformer<T,3,false>
    {
      using Array_t=VectorSoaContainer<T,3>;
      using Transformer_t=Tensor<T,3>;

      inline static void
        apply(const Array_t& pin, const Transformer_t& X, Array_t& pout, int first, int last)
        {
          const int n=last-first;
          register T x00=X[0],x01=X[1],x02=X[2],
                   x10=X[3],x11=X[4],x12=X[5],
                   x20=X[6],x21=X[7],x22=X[8];
          const T* restrict x_in=pin.data(0)+first; ASSUME_ALIGNED(x_in);
          const T* restrict y_in=pin.data(1)+first; ASSUME_ALIGNED(y_in);
          const T* restrict z_in=pin.data(2)+first; ASSUME_ALIGNED(z_in);
          T* restrict x_out=pout.data(0)+first; ASSUME_ALIGNED(x_out);
          T* restrict y_out=pout.data(1)+first; ASSUME_ALIGNED(y_out);
          T* restrict z_out=pout.data(2)+first; ASSUME_ALIGNED(z_out);
#pragma ivdep
          for(int i=0; i<n; i++)
          {
            x_out[i]=x_in[i]*x00+y_in[i]*x10+z_in[i]*x20;
            y_out[i]=x_in[i]*x01+y_in[i]*x11+z_in[i]*x21;
            z_out[i]=x_in[i]*x02+y_in[i]*x12+z_in[i]*x22;
          }
        }

      inline static void
        apply(const Transformer_t& X, const Array_t& pin,  Array_t& pout, int first, int last)
        {
          ::apply(pin,X,pout,first,last);
        }

      inline static void
        apply(Array_t& pinout, const Transformer_t& X,int first, int last)
        {
          const int n=last-first;
          register T x00=X[0],x01=X[1],x02=X[2],
                   x10=X[3],x11=X[4],x12=X[5],
                   x20=X[6],x21=X[7],x22=X[8];
          T* restrict x_inout=pinout.data(0)+first; ASSUME_ALIGNED(x_inout);
          T* restrict y_inout=pinout.data(1)+first; ASSUME_ALIGNED(y_inout);
          T* restrict z_inout=pinout.data(2)+first; ASSUME_ALIGNED(z_inout);
#pragma ivdep
          for(int i=0; i<n; i++)
          {
            T x=x_inout[i]*x00+y_inout[i]*x10+z_inout[i]*x20;
            T y=x_inout[i]*x01+y_inout[i]*x11+z_inout[i]*x21;
            T z=x_inout[i]*x02+y_inout[i]*x12+z_inout[i]*x22;
            x_inout[i]=x;
            y_inout[i]=y;
            z_inout[i]=z;
          }
        }

      inline static void
        apply(const Transformer_t& X, Array_t& pinout, int first, int last)
        {
          ::apply(X,pinout,first,last);
        }
    };
#endif

  /** General conversion function from AoS[nrows][ncols] to SoA[ncols][ldb]
   * @param nrows the first dimension
   * @param ncols the second dimension
   * @param iptr input pointer
   * @param lda stride of iptr
   * @param out output pointer
   * @param lda strided of out
   *
   * Modeled after blas/lapack for lda/ldb
   */
  template<typename T1, typename T2>
    void PosAoS2SoA(int nrows, int ncols, const T1* restrict iptr, int lda, T2* restrict out, int ldb)
    { 
      T2* restrict x=out      ;
      T2* restrict y=out+  ldb;
      T2* restrict z=out+2*ldb;
      #pragma omp simd aligned(x,y,z)
      for(int i=0; i<nrows;++i)
      {
        x[i]=iptr[i*ncols  ]; //x[i]=in[i][0];
        y[i]=iptr[i*ncols+1]; //y[i]=in[i][1];
        z[i]=iptr[i*ncols+2]; //z[i]=in[i][2];
      }
    }

  /** General conversion function from SoA[ncols][ldb] to AoS[nrows][ncols] 
   * @param nrows the first dimension
   * @param ncols the second dimension
   * @param iptr input pointer
   * @param lda stride of iptr
   * @param out output pointer
   * @param lda strided of out
   *
   * Modeled after blas/lapack for lda/ldb
   */
  template<typename T1, typename T2>
    void PosSoA2AoS(int nrows, int ncols, const T1* restrict iptr, int lda, T2* restrict out, int ldb)
    { 
      const T1* restrict x=iptr      ;
      const T1* restrict y=iptr+  lda;
      const T1* restrict z=iptr+2*lda;
      #pragma omp simd aligned(x,y,z)
      for(int i=0; i<nrows;++i)
      {
        out[i*ldb  ]=x[i]; //out[i][0]=x[i];
        out[i*ldb+1]=y[i]; //out[i][1]=y[i];
        out[i*ldb+2]=z[i]; //out[i][2]=z[i];
      }
    }

#if 0
//#if defined(HAVE_MKL)
  ///specialization for double AoS2SoA
  template<>
    void PosAoS2SoA(int nrows, int ncols, const double* restrict in, int lda, double* restrict out, int ldb)
    { 
      const double zone={1.0};
      mkl_domatcopy('R','T',nrows,ncols,zone,in,lda,out,ldb);
    }

  ///specialization for float AoS2SoA
  template<>
    void PosAoS2SoA(int nrows, int ncols, const float* restrict in, int lda, float* restrict out, int ldb)
    { 
      const float zone={1.0f};
      mkl_somatcopy('R','T',nrows,ncols,zone,in,lda,out,ldb);
    }

  ///specialization for double SoA2AoS
  template<>
    void PosSoA2AoS(int nrows, int ncols, const double* restrict in, int lda, double* restrict out, int ldb)
    { 
      const double zone={1.0};
      mkl_domatcopy('R','T',nrows,ncols,zone,in,lda,out,ldb);
    }

  ///specialization for float SoA2AoS
  template<>
    void PosSoA2AoS(int nrows, int ncols, const float* restrict in, int lda, float* restrict out, int ldb)
    { 
      const float zone={1.0f};
      mkl_somatcopy('R','T',nrows,ncols,zone,in,lda,out,ldb);
    }
#endif


}
#endif
