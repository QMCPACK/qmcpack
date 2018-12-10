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
    class DelayedUpdate
    {
      /// orbital values of delayed electrons
      Matrix<T> U;
      /// rows of Ainv corresponding to delayed electrons
      Matrix<T> V;
      /// Matrix inverse of B, at maximum KxK
      Matrix<T> Binv;
      /// scratch space, used during inverse update
      Matrix<T> tempMat;
      /// temporal scratch space used by SM-1
      Vector<T> temp;
      /// new column of B
      Vector<T> p;
      /// list of delayed electrons
      std::vector<int> delay_list;
      /// current number of delays, increase one for each acceptance, reset to 0 after updating Ainv
      int delay_count;

      /// pointer to the row of up-to-date Ainv
      const T* Ainv_row_ptr;
      /// electron id of the up-to-date Ainv_row, checked by ratioGrad
      int Ainv_row_ind;
      /// current determinant ratio
      T curRatio;

    public:
      /// default constructor
      DelayedUpdate(): delay_count(0), Ainv_row_ptr(nullptr), Ainv_row_ind(-1) {}

      ///resize the internal storage
      /** resize the internal storage
       * @param norb number of electrons/orbitals
       * @param delay, maximum delay 0<delay<=norb
       */
      inline void resize(int norb, int delay)
      {
        V.resize(delay, norb);
        U.resize(delay, norb);
        p.resize(delay);
        temp.resize(norb);
        tempMat.resize(norb, delay);
        Binv.resize(delay, delay);
        delay_list.resize(delay);
      }

      /** compute the row of up-to-date Ainv
       * @param Ainv inverse matrix
       * @param rowchanged the row id corresponding to the proposed electron
       */
      inline void getInvRow(const Matrix<T>& Ainv, int rowchanged)
      {
        Ainv_row_ind = rowchanged;
        if ( delay_count == 0 )
        {
          // Ainv is fresh, directly access Ainv
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

      /** compute determinant ratio of new determinant
       * @param Ainv inverse matrix
       * @param rowchanged the row id corresponding to the proposed electron
       * @param psiV new orbital values
       */
      template<typename VVT>
      inline T ratio(const Matrix<T>& Ainv, int rowchanged, const VVT& psiV)
      {
        getInvRow(Ainv, rowchanged);
        return curRatio = simd::dot(Ainv_row_ptr,psiV.data(),Ainv.cols());
      }

      /** compute the old gradient
       * @param Ainv inverse matrix
       * @param rowchanged the row id corresponding to the proposed electron
       * @param dpsiV old orbital derivatives
       */
      template<typename GT>
      inline GT evalGrad(const Matrix<T>& Ainv, int rowchanged, const GT* dpsiV)
      {
        getInvRow(Ainv, rowchanged);
        return simd::dot(Ainv_row_ptr,dpsiV,Ainv.cols());
      }

      /** compute determinant ratio and gradients of new determinant
       * @param Ainv inverse matrix
       * @param rowchanged the row id corresponding to the proposed electron
       * @param psiV new orbital values
       * @param dpsiV new orbital derivatives
       * @param g new gradients
       */
      template<typename VVT, typename GGT, typename GT>
      inline T ratioGrad(const Matrix<T>& Ainv, int rowchanged, const VVT& psiV, const GGT& dpsiV, GT& g)
      {
        // check Ainv_row_ind against rowchanged to ensure getInvRow() called before ratioGrad()
        if(Ainv_row_ind != rowchanged)
          getInvRow(Ainv, rowchanged);
        g = simd::dot(Ainv_row_ptr,dpsiV.data(),Ainv.cols());
        return curRatio = simd::dot(Ainv_row_ptr,psiV.data(),Ainv.cols());
      }

      /** accept a move with the update delayed
       * @param Ainv inverse matrix
       * @param rowchanged the row id corresponding to the proposed electron
       * @param psiV new orbital values
       *
       * Before delay_count reaches the maximum delay, only Binv is updated with a recursive algorithm
       */
      template<typename VVT>
      inline void acceptRow(Matrix<T>& Ainv, int rowchanged, const VVT& psiV)
      {
        // safe mechanism
        Ainv_row_ind = -1;

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
        // update Ainv when maximal delay is reached
        if(delay_count==lda_Binv) updateInvMat(Ainv);
      }

      /** update the full Ainv and reset delay_count
       * @param Ainv inverse matrix
       */
      inline void updateInvMat(Matrix<T>& Ainv)
      {
        if(delay_count==0) return;
        // update the inverse matrix
        const T cone(1);
        const T czero(0);
        const int norb=Ainv.rows();
        if(delay_count==1)
        {
          // this is a special case invoking the Fahy's variant of Sherman-Morrison update.
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
        Ainv_row_ind = -1;
      }
    };
}

#endif // QMCPLUSPLUS_DELAYED_UPDATE_H

