//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2017 QMCPACK developers.
//
// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_DELAYED_UPDATE_H
#define QMCPLUSPLUS_DELAYED_UPDATE_H

#include "Numerics/Blasf.h"
#include <OhmmsPETE/OhmmsVector.h>
#include <OhmmsPETE/OhmmsMatrix.h>
#include <simd/simd.hpp>

namespace qmcplusplus {

  template<typename T>
    struct DelayedUpdate
    {
      Matrix<T> U, V, Binv, tempMat;
      // temporal scratch space used by SM-1
      Vector<T> temp;
      // auxiliary arrays for B
      Vector<T> p;
      std::vector<int> delay_list;
      int delay_count;

      const T* Ainv_row_ptr;
      T curRatio;

      DelayedUpdate(): delay_count(0), Ainv_row_ptr(nullptr) {}

      inline void resize(int norb, int delay)
      {
        if(delay<=0) delay=1;
        V.resize(delay, norb);
        U.resize(delay, norb);
        p.resize(delay);
        temp.resize(norb);
        tempMat.resize(norb, delay);
        Binv.resize(delay, delay);
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
        const T cone(1);
        const T czero(0);
        const T* AinvRow = Ainv[rowchanged];
        const int norb = Ainv.rows();
        const int lda_Binv = Binv.cols();
        // save AinvRow to new_AinvRow
        simd::copy_n(AinvRow, norb, V[delay_count]);
        // multiply V (NxK) Binv(KxK) U(KxN) AinvRow right to the left
        BLAS::gemv('T', norb, delay_count, cone, U.data(), norb, AinvRow, 1, czero, p.data(), 1);
        BLAS::gemv('N', delay_count, delay_count, cone, Binv.data(), lda_Binv, p.data(), 1, czero, Binv[delay_count], 1);
        BLAS::gemv('N', norb, delay_count, -cone, V.data(), norb, Binv[delay_count], 1, cone, V[delay_count], 1);
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
          throw std::runtime_error("DelayedUpdate : this should never happen!\n");
          return T(0);
        }
      }

      // accept with the update delayed
      template<typename VVT>
      inline void acceptRow(Matrix<T>& Ainv, int rowchanged, const VVT& psiV)
      {
        // safe mechanism
        Ainv_row_ptr = nullptr;

        const T cminusone(-1);
        const T czero(0);
        const int norb = Ainv.rows();
        const int lda_Binv = Binv.cols();
        simd::copy_n(Ainv[rowchanged], norb, V[delay_count]);
        simd::copy_n(psiV.data(), norb, U[delay_count]);
        delay_list[delay_count] = rowchanged;
        // the new Binv is [[X Y] [Z x]]
        BLAS::gemv('T', norb, delay_count+1, cminusone, V.data(), norb, psiV.data(), 1, czero, p.data(), 1);
        // x
        T y = -p[delay_count];
        for(int i=0; i<delay_count; i++)
          y += Binv[delay_count][i] * p[i];
        Binv[delay_count][delay_count] = y = T(1) / y;
        // Y
        BLAS::gemv('T', delay_count, delay_count, y, Binv.data(), lda_Binv, p.data(), 1, czero, Binv.data()+delay_count, lda_Binv);
        // X
        BLAS::ger(delay_count, delay_count, cminusone, Binv[delay_count], 1, Binv.data()+delay_count, lda_Binv, Binv.data(), lda_Binv);
        // Z
        for(int i=0; i<delay_count; i++)
          Binv[delay_count][i] *= -y;
        delay_count++;
        if(delay_count==lda_Binv) updateInvMat(Ainv);
      }

      inline void updateInvMat(Matrix<T>& Ainv)
      {
        if(delay_count==0) return;
        // update the inverse matrix
        const T cone(1);
        const T czero(0);
        const int norb=Ainv.rows();
        if(delay_count==1)
        {
          // Only use the first norb elements of tempMat as a temporal array
          BLAS::gemv('T', norb, norb, cone, Ainv.data(), norb, U[0], 1, czero, temp.data(), 1);
          temp[delay_list[0]] -= cone;
          BLAS::ger(norb,norb,-Binv[0][0],V[0],1,temp.data(),1,Ainv.data(),norb);
        }
        else
        {
          const int lda_Binv=Binv.cols();
          BLAS::gemm('T', 'N', delay_count, norb, norb, cone, U.data(), norb, Ainv.data(), norb, czero, tempMat.data(), lda_Binv);
          for(int i=0; i<delay_count; i++) tempMat(delay_list[i], i) -= cone;
          BLAS::gemm('N', 'N', norb, delay_count, delay_count, cone, V.data(), norb, Binv.data(), lda_Binv, czero, U.data(), norb);
          BLAS::gemm('N', 'N', norb, norb, delay_count, -cone, U.data(), norb, tempMat.data(), lda_Binv, cone, Ainv.data(), norb);
        }
        delay_count = 0;
      }
    };
}

#endif // QMCPLUSPLUS_DELAYED_UPDATE_H

