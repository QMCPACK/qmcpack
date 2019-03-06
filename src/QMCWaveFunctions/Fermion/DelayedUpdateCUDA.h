//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2019 QMCPACK developers.
//
// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_DELAYED_UPDATE_CUDA_H
#define QMCPLUSPLUS_DELAYED_UPDATE_CUDA_H

#include <OhmmsPETE/OhmmsVector.h>
#include <OhmmsPETE/OhmmsMatrix.h>
#include "simd/CUDAallocator.hpp"
#include <Numerics/CUDA/cuBLAS.hpp>
#include "QMCWaveFunctions/Fermion/delayed_update_helper.h"
#include <cuda_runtime_api.h>
#include "CUDA/cudaError.h"

namespace qmcplusplus {

  class Range
  {
    // [first, last) rows of Ainv
    int first, last;
    public:
    Range() : first(0), last(0) { };
    void setRange(int first_in, int last_in)
    {
      first = first_in;
      last  = last_in;
    }
    inline void clear() { first = last = 0; };
    inline int getOffset(int index) const
    {
      if(!checkRange(index)) throw std::runtime_error("index not in range \n");
      return index-first;
    }
    inline bool checkRange(int index) const { return (index>=first) && (index<last); };
  };

  template<typename T>
    struct DelayedUpdateCUDA
    {
      // Data staged during for delayed acceptRows
      Matrix<T, CUDAHostAllocator<T>> U, Binv;
      Matrix<T> V;
      //Matrix<T> tempMat; // for debugging only
      Matrix<T, CUDAAllocator<T>, MemorySpace::CUDA> U_gpu, V_gpu, Binv_gpu, temp_gpu, Ainv_gpu;
      // temporal scratch space used by SM-1
      Vector<T> temp, rcopy;
      // auxiliary arrays for B
      Vector<T> p;
      Vector<int, CUDAHostAllocator<int>> delay_list;
      Vector<int, CUDAAllocator<int>, MemorySpace::CUDA> delay_list_gpu;
      int delay_count;

      // the range of prefetched_Ainv_rows
      Range prefetched_range;
      // Ainv prefetch buffer
      Matrix<T, CUDAHostAllocator<T>> Ainv_buffer;

      // CUDA specific variables
      cublasHandle_t handle;
      cudaStream_t hstream;

      DelayedUpdateCUDA(): delay_count(0)
      {
        // currently hard-coded first device
        // need to control CUDA_VISIBLE_DEVICES when using multiple MPI on a single node.
        cudaErrorCheck( cudaSetDevice(0), "cudaSetDevice failed!" );
        cudaErrorCheck( cudaStreamCreate(&hstream), "cudaStreamCreate failed!" );
        cublasErrorCheck( cublasCreate(&handle), "cublasCreate failed!" );
        cublasErrorCheck( cublasSetStream(handle,hstream), "cublasSetStream failed!" );
      }

      ~DelayedUpdateCUDA()
      {
        cublasErrorCheck( cublasDestroy(handle), "cublasDestroy failed!" );
        cudaErrorCheck( cudaStreamDestroy(hstream), "cudaStreamDestroy failed!" );
      }

      inline void resize(int norb, int delay)
      {
        //tempMat.resize(norb, delay);
        V.resize(delay, norb);
        U.resize(delay, norb);
        p.resize(delay);
        Binv.resize(delay, delay);
        // prefect 8% more rows corresponding to roughly 96% acceptance ratio
        Ainv_buffer.resize(std::min(static_cast<int>(delay*1.08),norb), norb);

        temp_gpu.resize(norb, delay);
        delay_list.resize(delay);
        U_gpu.resize(delay, norb);
        V_gpu.resize(delay, norb);
        Binv_gpu.resize(delay, delay);
        Ainv_gpu.resize(norb, norb);
        delay_list_gpu.resize(delay);
      }

      inline void initializeInv(const Matrix<T>& Ainv)
      {
        cudaMemcpyAsync(Ainv_gpu.data(), Ainv.data(), Ainv.size()*sizeof(T), cudaMemcpyHostToDevice, hstream);
        // safe mechanism
        delay_count = 0;
        prefetched_range.clear();
      }

      template<typename VVT>
      inline void getInvRow(const Matrix<T>& Ainv, int rowchanged, VVT& invRow)
      {
        if(!prefetched_range.checkRange(rowchanged))
        {
          int last_row = std::min(rowchanged+Ainv_buffer.rows(), Ainv.rows());
          cudaMemcpyAsync(Ainv_buffer.data(), Ainv_gpu[rowchanged],
                          invRow.size()*(last_row-rowchanged)*sizeof(T),
                          cudaMemcpyDeviceToHost, hstream);
          prefetched_range.setRange(rowchanged, last_row);
          waitStream();
        }
        // save AinvRow to new_AinvRow
        std::copy_n(Ainv_buffer[prefetched_range.getOffset(rowchanged)], invRow.size(), invRow.data());
        if ( delay_count > 0 )
        {
          const T cone(1);
          const T czero(0);
          const int norb = Ainv.rows();
          const int lda_Binv = Binv.cols();
          // multiply V (NxK) Binv(KxK) U(KxN) AinvRow right to the left
          BLAS::gemv('T', norb, delay_count, cone, U.data(), norb, invRow.data(), 1, czero, p.data(), 1);
          BLAS::gemv('N', delay_count, delay_count, cone, Binv.data(), lda_Binv, p.data(), 1, czero, Binv[delay_count], 1);
          BLAS::gemv('N', norb, delay_count, -cone, V.data(), norb, Binv[delay_count], 1, cone, invRow.data(), 1);
        }
      }

      // accept with the update delayed
      template<typename VVT>
      inline void acceptRow(Matrix<T>& Ainv, int rowchanged, const VVT& psiV)
      {
        // update Binv from delay_count to delay_count+1
        const T cminusone(-1);
        const T czero(0);
        const int norb = Ainv.rows();
        const int lda_Binv = Binv.cols();
        std::copy_n(Ainv_buffer[prefetched_range.getOffset(rowchanged)], norb, V[delay_count]);
        std::copy_n(psiV.data(), norb, U[delay_count]);
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
        if(delay_count==lda_Binv) updateInvMat(Ainv,false);
      }

      inline void updateInvMat(Matrix<T>& Ainv, bool transfer_to_host = true)
      {
        // update the inverse matrix
        if( delay_count>0 )
        {
          const T cone(1);
          const T czero(0);
          const int norb=Ainv.rows();
          const int lda_Binv=Binv.cols();
          const T cminusone(-1);
          cudaError_t error(cudaSuccess);
          error = cudaMemcpyAsync(U_gpu.data(), U.data(), norb*delay_count*sizeof(T), cudaMemcpyHostToDevice, hstream);
              if( error!=cudaSuccess ) std::cout <<"debug error 1 code " << error << std::endl;
          //BLAS::gemm('T', 'N', delay_count, norb, norb, cone, U.data(), norb, Ainv.data(), norb, czero, tempMat.data(), lda_Binv);
          //    cudaMemPrefetchAsync(Ainv.data(), Ainv.size()*sizeof(T), 0, hstream);
          cuBLAS::gemm(handle, CUBLAS_OP_T, CUBLAS_OP_N, delay_count, norb, norb, &cone, U_gpu.data(), norb, Ainv_gpu.data(), norb, &czero, temp_gpu.data(), lda_Binv);
          error = cudaMemcpyAsync(delay_list_gpu.data(), delay_list.data(), delay_count*sizeof(int), cudaMemcpyHostToDevice, hstream);
              if( error!=cudaSuccess ) std::cout <<"debug error 2 code " << error << std::endl;
          //for(int i=0; i<delay_count; i++) tempMat(delay_list[i], i) -= cone;
          applyW_stageV_cuda(delay_list_gpu.data(), delay_count, temp_gpu.data(), norb, temp_gpu.cols(), V_gpu.data(), Ainv_gpu.data(), hstream);
          //error = cudaMemcpyAsync(tempMat.data(), temp_gpu.data(), tempMat.size()*sizeof(T), cudaMemcpyDeviceToHost, hstream);
          //    if( error!=cudaSuccess ) std::cout <<"debug error 3 code " << error << std::endl;
          //error = cudaStreamSynchronize(hstream);
          //    if( error!=cudaSuccess ) std::cout <<"debug error 4 code " << error << std::endl;
              //std::cout << "debug tempMat " << tempMat << std::endl;
          error = cudaMemcpyAsync(Binv_gpu.data(), Binv.data(), lda_Binv*delay_count*sizeof(T), cudaMemcpyHostToDevice, hstream);
              if( error!=cudaSuccess ) std::cout <<"debug error 3 code " << error << std::endl;
          //BLAS::gemm('N', 'N', norb, delay_count, delay_count, cone, V.data(), norb, Binv.data(), lda_Binv, czero, U.data(), norb);
          cuBLAS::gemm(handle, CUBLAS_OP_N, CUBLAS_OP_N, norb, delay_count, delay_count, &cone, V_gpu.data(), norb, Binv_gpu.data(), lda_Binv, &czero, U_gpu.data(), norb);
          //error = cudaMemcpyAsync(U.data(), U_gpu.data(), norb*delay_count*sizeof(T), cudaMemcpyDeviceToHost, hstream);
          //    if( error!=cudaSuccess ) std::cout <<"debug error 5 code " << error << std::endl;
          //error = cudaStreamSynchronize(hstream);
          //    if( error!=cudaSuccess ) std::cout <<"debug error 6 code " << error << std::endl;
              //std::cout << "debug U " << U << std::endl;
          //BLAS::gemm('N', 'N', norb, norb, delay_count, -cone, U.data(), norb, tempMat.data(), lda_Binv, cone, Ainv.data(), norb);
          cuBLAS::gemm(handle, CUBLAS_OP_N, CUBLAS_OP_N, norb, norb, delay_count, &cminusone, U_gpu.data(), norb, temp_gpu.data(), lda_Binv, &cone, Ainv_gpu.data(), norb);
          delay_count = 0;
          // Ainv is invalid, reset range
          prefetched_range.clear();
        }

        // transfer Ainv_gpu to Ainv and wait till completion
        if(transfer_to_host)
        {
          cudaMemcpyAsync(Ainv.data(), Ainv_gpu.data(), Ainv.size()*sizeof(T), cudaMemcpyDeviceToHost, hstream);
          // no need to wait because : For transfers from device memory to pageable host memory, the function will return only once the copy has completed.
          //waitStream();
        }
      }

      inline void waitStream()
      {
        cudaError_t error = cudaStreamSynchronize(hstream);
        if( error!=cudaSuccess ) throw std::runtime_error("DelayedUpdateCUDA::waitStream cudaStreamSynchronize failed!\n");
      }
    };
}

#endif // QMCPLUSPLUS_DELAYED_UPDATE_H

