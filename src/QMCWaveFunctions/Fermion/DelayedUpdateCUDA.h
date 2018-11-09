//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2018 QMCPACK developers.
//
// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_DELAYED_UPDATE_CUDA_H
#define QMCPLUSPLUS_DELAYED_UPDATE_CUDA_H

#include "Numerics/Blasf.h"
#include <OhmmsPETE/OhmmsVector.h>
#include <OhmmsPETE/OhmmsMatrix.h>
#include "QMCWaveFunctions/Fermion/DiracMatrix.h"
#include <simd/simd.hpp>
#include <cublas_v2.h>

namespace qmcplusplus {

  template<typename T, typename T_hp>
    struct DelayedUpdateCUDA
    {
      Matrix<T> U, V, B, Binv;
      Matrix<T, CUDAAllocator<T>> U_gpu, V_gpu, Binv_gpu, temp_gpu;
      Vector<T> temp, rcopy;
      Matrix<T_hp> Binv_hp;
      DiracMatrix<T_hp> deteng;
      std::vector<int> delay_list;
      int delay_count;

      const T* Ainv_row_ptr;
      T curRatio;

      // CUDA specific variables
      cublasHandle_t handle;
      cudaStream_t hstream;

      DelayedUpdateCUDA(): delay_count(0), Ainv_row_ptr(nullptr)
      {
        cublasCreate(&handle);
        cudaStreamCreate(&hstream);
        cublasSetStream(handle,hstream);
      }

      ~DelayedUpdateCUDA()
      {
        cublasDestroy(handle);
        cudaStreamDestroy(hstream);
      }

      inline void resize(int norb, int delay)
      {
        U_gpu.resize(delay, norb);
        V_gpu.resize(delay, norb);
        temp_gpu.resize(norb, delay);
        Binv_gpu.resize(delay, delay);
        V.resize(delay, norb);
        U.resize(delay, norb);
        B.resize(delay, delay);
        Binv.resize(delay, delay);
#ifdef MIXED_PRECISION
        Binv_hp.resize(delay, delay);
        deteng.reset(Binv_hp, delay);
#else
        deteng.reset(Binv, delay);
#endif
        delay_count = 0;
        delay_list.resize(delay);
      }

      inline void getInvRow(const Matrix<T>& Ainv, int rowchanged)
      {
        if ( delay_count == 0 )
        {
          Ainv_row_ptr = Ainv[rowchanged];
          return;
        }
        CONSTEXPR T cone(1);
        CONSTEXPR T czero(0);
        const T* AinvRow = Ainv[rowchanged];
        const int norb = Ainv.rows();
        const int lda_Binv = Binv.cols();
        T temp[lda_Binv];
        // save AinvRow to new_AinvRow
        simd::copy_n(AinvRow, norb, V[delay_count]);
        // multiply V (NxK) Binv(KxK) U(KxN) AinvRow right to the left
        BLAS::gemv('T', norb, delay_count, cone, U.data(), norb, AinvRow, 1, czero, B[delay_count], 1);
        BLAS::gemv('N', delay_count, delay_count, cone, Binv.data(), lda_Binv, B[delay_count], 1, czero, temp, 1);
        BLAS::gemv('N', norb, delay_count, -cone, V.data(), norb, temp, 1, cone, V[delay_count], 1);
        Ainv_row_ptr = V[delay_count];
      }

      template<typename VVT>
      inline T ratio(const Matrix<T>& Ainv, int rowchanged, const VVT& psiV)
      {
        getInvRow(Ainv, rowchanged);
        return curRatio = simd::dot(Ainv_row_ptr,psiV.data(),Ainv.cols());
      }

      template<typename GT>
      inline GT evalGrad(const Matrix<T>& Ainv, int rowchanged, const GT* dpsiV)
      {
        getInvRow(Ainv, rowchanged);
        return simd::dot(Ainv_row_ptr,dpsiV,Ainv.cols());
      }

      template<typename VVT, typename GGT, typename GT>
      inline T ratioGrad(const Matrix<T>& Ainv, int rowchanged, const VVT& psiV, const GGT& dpsiV, GT& g)
      {
        if( delay_count == 0 )
        {
          g = simd::dot(Ainv[rowchanged],dpsiV.data(),Ainv.cols());
          return curRatio = simd::dot(Ainv[rowchanged],psiV.data(),Ainv.cols());
        }
        else if( Ainv_row_ptr )
        {
          g = simd::dot(Ainv_row_ptr,dpsiV.data(),Ainv.cols());
          return curRatio = simd::dot(Ainv_row_ptr,psiV.data(),Ainv.cols());
        }
        else
        {
          throw std::runtime_error("DelayedUpdateCUDA : this should never happen!");
          return T(0);
        }
      }

      // SM-1 Fahy immediate update
      template<typename VVT>
      inline void updateRow(Matrix<T>& a, int rowchanged, const VVT& psiV)
      {
        // safe mechanism
        Ainv_row_ptr = nullptr;

        const int m = a.rows();
        const int lda = a.cols();
        CONSTEXPR T cone(1);
        CONSTEXPR T czero(0);
        temp.resize(lda);
        rcopy.resize(lda);
        T c_ratio = cone / curRatio;
        BLAS::gemv('T', m, m, c_ratio, a.data(), lda, psiV.data(), 1, czero, temp.data(), 1);
        temp[rowchanged] = cone-c_ratio;
        simd::copy_n(a[rowchanged],m,rcopy.data());
        BLAS::ger(m,m,-cone,rcopy.data(),1,temp.data(),1,a.data(),lda);
      }

      // accept with the update delayed
      template<typename VVT>
      inline void acceptRow(Matrix<T>& Ainv, int rowchanged, const VVT& psiV)
      {
        // safe mechanism
        Ainv_row_ptr = nullptr;

        CONSTEXPR T cone(1);
        CONSTEXPR T czero(0);
        const int norb = Ainv.rows();
        const int lda_Binv = Binv.cols();
        simd::copy_n(Ainv[rowchanged], norb, V[delay_count]);
        simd::copy_n(psiV.data(), norb, U[delay_count]);
        delay_list[delay_count] = rowchanged;
        delay_count++;
        // the new row in B has been computed, now compute the new column
        if(delay_count==1)
        {
          B[0][0] = curRatio;
          Binv[0][0] = T_hp(1.0) / static_cast<T_hp>(B[0][0]);
        }
        else
        {
          BLAS::gemv('T', norb, delay_count, cone, V.data(), norb, psiV.data(), 1, czero, B.data()+delay_count-1, lda_Binv);
#ifdef MIXED_PRECISION
          for(int i=0; i<delay_count; i++)
            for(int j=0; j<delay_count; j++)
              Binv_hp[i][j] = B[i][j];
          deteng.invert(Binv_hp,false,delay_count);
          for(int i=0; i<delay_count; i++)
            for(int j=0; j<delay_count; j++)
              Binv[i][j] = Binv_hp[i][j];
#else
          for(int i=0; i<delay_count; i++)
            for(int j=0; j<delay_count; j++)
              Binv[i][j] = B[i][j];
          deteng.invert(Binv,false,delay_count);
#endif
        }
        if(delay_count==lda_Binv) updateInvMat(Ainv);
      }

      inline void updateInvMat(Matrix<T>& Ainv)
      {
        if(delay_count==0) return;
        // update the inverse matrix
        CONSTEXPR T cone(1);
        CONSTEXPR T czero(0);
        const int norb=Ainv.rows();
        if(delay_count==1)
        {
          // Only use the first norb elements of temp_gpu as a temporal array
          BLAS::gemv('T', norb, norb, cone, Ainv.data(), norb, U[0], 1, czero, temp.data(), 1);
          temp[delay_list[0]] -= cone;
          BLAS::ger(norb,norb,-Binv[0][0],V[0],1,temp_gpu[0],1,Ainv.data(),norb);
        }
        else
        {
          const int lda_Binv=Binv.cols();
          const T cminusone(-1);
          //BLAS::gemm('T', 'N', delay_count, norb, norb, cone, U_gpu.data(), norb, Ainv.data(), norb, czero, temp_gpu.data(), lda_Binv);
          cublasDgemm(handle, CUBLAS_OP_T, CUBLAS_OP_N, delay_count, norb, norb, &cone, U_gpu.data(), norb, Ainv.data(), norb, &czero, temp_gpu.data(), lda_Binv);
          for(int i=0; i<delay_count; i++) temp_gpu(delay_list[i], i) -= cone;
          //BLAS::gemm('N', 'N', norb, delay_count, delay_count, cone, V_gpu.data(), norb, Binv_gpu.data(), lda_Binv, czero, U_gpu.data(), norb);
          cublasDgemm(handle, CUBLAS_OP_N, CUBLAS_OP_N, norb, delay_count, delay_count, &cone, V_gpu.data(), norb, Binv_gpu.data(), lda_Binv, &czero, U_gpu.data(), norb);
          //BLAS::gemm('N', 'N', norb, norb, delay_count, -cone, U_gpu.data(), norb, temp_gpu.data(), lda_Binv, cone, Ainv.data(), norb);
          cublasDgemm(handle, CUBLAS_OP_N, CUBLAS_OP_N, norb, norb, delay_count, &cminusone, U_gpu.data(), norb, temp_gpu.data(), lda_Binv, &cone, Ainv.data(), norb);
        }
        delay_count = 0;
      }
    };
}

#endif // QMCPLUSPLUS_DELAYED_UPDATE_H

