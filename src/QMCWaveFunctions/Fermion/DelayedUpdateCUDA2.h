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

#ifndef QMCPLUSPLUS_DELAYED_UPDATE_CUDA2_H
#define QMCPLUSPLUS_DELAYED_UPDATE_CUDA2_H

#include "Numerics/Blasf.h"
#include <OhmmsPETE/OhmmsVector.h>
#include <OhmmsPETE/OhmmsMatrix.h>
#include <simd/simd.hpp>
#include "simd/CUDAallocator.hpp"
#include <cublas_v2.h>
#include <Numerics/CUDA/cuBLAS.hpp>
#include "QMCWaveFunctions/Fermion/delayed_update_helper.h"

namespace qmcplusplus {

  template<typename T>
    struct DelayedUpdateCUDA
    {
      //Matrix<T> U, V, B, Binv;
      //Matrix<T> tempMat; // for debugging only
      Matrix<T, CUDAAllocator<T>, MemorySpace::CUDA> U_gpu, V_gpu, Binv_gpu, temp_gpu, Ainv_gpu;
      Vector<T, CUDAHostAllocator<T>> U_row;
      // temporal scratch space used by SM-1
      Vector<T> temp, rcopy;
      // auxiliary arrays for B
      //Vector<T> p;
      Vector<T, CUDAAllocator<T>, MemorySpace::CUDA> p_gpu;
      //std::vector<int> delay_list;
      Vector<int, CUDAAllocator<int>, MemorySpace::CUDA> delay_list_gpu;
      int delay_count;

      Vector<T, CUDAHostAllocator<T>> Ainv_row;
      Vector<T, CUDAAllocator<T>, MemorySpace::CUDA> Ainv_row_gpu;
      // electron id of the up-to-date Ainv_row
      int Ainv_row_ind;
      // current ratio
      T curRatio;

      // CUDA specific variables
      cublasHandle_t handle;
      cudaStream_t hstream;
      // size 3 constant on GPU for 0, 1, -1
      Vector<T, CUDAAllocator<T>, MemorySpace::CUDA> constants_gpu;

      DelayedUpdateCUDA(): delay_count(0), Ainv_row_ind(-1)
      {
        cublasCreate(&handle);
        cublasSetPointerMode(handle, CUBLAS_POINTER_MODE_DEVICE);
        cudaStreamCreate(&hstream);
        cublasSetStream(handle, hstream);
        const T constants[3] = {0, 1, -1};
        constants_gpu.resize(3);
        cudaMemcpyAsync(constants_gpu.data(), constants, 3*sizeof(T), cudaMemcpyHostToDevice, hstream);
        cudaStreamSynchronize(hstream);
      }

      ~DelayedUpdateCUDA()
      {
        cublasDestroy(handle);
        cudaStreamDestroy(hstream);
      }

      inline void resize(int norb, int delay)
      {
        delay_count = 0;

        //tempMat.resize(norb, delay);
        //V.resize(delay, norb);
        //U.resize(delay, norb);
        U_row.resize(norb);
        //p.resize(delay);
        //Binv.resize(delay, delay);
        Ainv_row.resize(norb);
        //delay_list.resize(delay);

        temp_gpu.resize(norb, delay);
        U_gpu.resize(delay, norb);
        V_gpu.resize(delay, norb);
        p_gpu.resize(delay);
        Binv_gpu.resize(delay, delay);
        Ainv_gpu.resize(norb, norb);
        Ainv_row_gpu.resize(norb);
        delay_list_gpu.resize(delay);
      }

      inline void transferAinvH2D(const Matrix<T>& Ainv)
      {
        cudaMemcpyAsync(Ainv_gpu.data(), Ainv.data(), Ainv.size()*sizeof(T), cudaMemcpyHostToDevice, hstream);
      }

      inline void getInvRow(const Matrix<T>& Ainv, int rowchanged)
      {
        if ( delay_count == 0 )
        {
          cudaMemcpyAsync(Ainv_row.data(), Ainv_gpu[rowchanged], Ainv_row.size()*sizeof(T), cudaMemcpyDeviceToHost, hstream);
        }
        else
        {
          const int norb = Ainv_gpu.rows();
          const int lda_Binv = Binv_gpu.cols();
          // save AinvRow to new_AinvRow
          //simd::copy_n(AinvRow, norb, V[delay_count]);
          cublasScopy(handle, norb, Ainv_gpu[rowchanged], 1, Ainv_row_gpu.data(), 1);
          // multiply V (NxK) Binv(KxK) U(KxN) AinvRow right to the left
          cublasSgemv(handle, CUBLAS_OP_T, norb, delay_count, constants_gpu.data()+1, U_gpu.data(), norb, Ainv_gpu[rowchanged], 1, constants_gpu.data(), p_gpu.data(), 1);
          cublasSgemv(handle, CUBLAS_OP_N, delay_count, delay_count, constants_gpu.data()+1, Binv_gpu.data(), lda_Binv, p_gpu.data(), 1, constants_gpu.data(), Binv_gpu[delay_count], 1);
          cublasSgemv(handle, CUBLAS_OP_N, norb, delay_count, constants_gpu.data()+2, V_gpu.data(), norb, Binv_gpu[delay_count], 1, constants_gpu.data()+1, Ainv_row_gpu.data(), 1);
          cudaMemcpyAsync(Ainv_row.data(), Ainv_row_gpu.data(), Ainv_row.size()*sizeof(T), cudaMemcpyDeviceToHost, hstream);
        }
        waitStream();
        Ainv_row_ind = rowchanged;
      }

      template<typename VVT>
      inline T ratio(const Matrix<T>& Ainv, int rowchanged, const VVT& psiV)
      {
        getInvRow(Ainv, rowchanged);
        //std::cout << "report ratio " << std::endl;
        return curRatio = simd::dot(Ainv_row.data(),psiV.data(),Ainv_row.size());
      }

      template<typename GT>
      inline GT evalGrad(const Matrix<T>& Ainv, int rowchanged, const GT* dpsiV)
      {
        getInvRow(Ainv, rowchanged);
        //std::cout << "report evalGrad " << std::endl;
        return simd::dot(Ainv_row.data(),dpsiV,Ainv_row.size());
      }

      template<typename VVT, typename GGT, typename GT>
      inline T ratioGrad(const Matrix<T>& Ainv, int rowchanged, const VVT& psiV, const GGT& dpsiV, GT& g)
      {
        if(Ainv_row_ind != rowchanged)
          getInvRow(Ainv, rowchanged);
        g = simd::dot(Ainv_row.data(),dpsiV.data(),Ainv_row.size());
        curRatio = simd::dot(Ainv_row.data(),psiV.data(),Ainv_row.size());
        //std::cout << "report ratioGrad " << curRatio << std::endl;
        return curRatio;
      }

      // accept with the update delayed
      template<typename VVT>
      inline void acceptRow(Matrix<T>& Ainv, int rowchanged, const VVT& psiV)
      {
        //std::cout << "report acceptRow start " << std::endl;
        // safe mechanism
        Ainv_row_ind = -1;

        const int norb = Ainv_gpu.rows();
        const int lda_Binv = Binv_gpu.cols();
        //simd::copy_n(Ainv[rowchanged], norb, V[delay_count]);
        cublasScopy(handle, norb, Ainv_gpu[rowchanged], 1, V_gpu[delay_count], 1);
        //save psiV to U_row just in case psiV is changed on the host before cudaMemcpyAsync is completed
        //if the getInvRow is called before the SPO evaluation, this can be removed.
        simd::copy_n(psiV.data(), norb, U_row.data());
        //delay_list[delay_count] = rowchanged;
        cudaMemcpyAsync(U_gpu[delay_count], U_row.data(), norb*sizeof(T), cudaMemcpyHostToDevice, hstream);
        // the new Binv is [[X Y] [Z x]]
        cublasSgemv(handle, CUBLAS_OP_T, norb, delay_count+1,
                    constants_gpu.data()+2, V_gpu.data(), norb, U_gpu[delay_count], 1, constants_gpu.data(), p_gpu.data(), 1);
        // x
        updateBinv_x_cuda(delay_list_gpu.data(), delay_count, rowchanged, Binv_gpu[delay_count], p_gpu.data(), hstream);
        // Y
        cublasSgemv(handle, CUBLAS_OP_T, delay_count, delay_count,
                    Binv_gpu[delay_count]+delay_count, Binv_gpu.data(), lda_Binv, p_gpu.data(), 1, constants_gpu.data(), Binv_gpu.data()+delay_count, lda_Binv);
        // X
        cublasSger(handle, delay_count, delay_count, constants_gpu.data()+2, Binv_gpu[delay_count], 1, Binv_gpu.data()+delay_count, lda_Binv, Binv_gpu.data(), lda_Binv);
        // Z
        cublasSscal(handle, delay_count, p_gpu.data()+delay_count, Binv_gpu[delay_count], 1);
        delay_count++;
        if(delay_count==lda_Binv) updateInvMat(Ainv,false);
        /*
        cudaError_t error(cudaSuccess);
        error = cudaMemcpyAsync(Binv.data(), Binv_gpu.data(), lda_Binv*lda_Binv*sizeof(T), cudaMemcpyDeviceToHost, hstream);
        waitStream();
        for(int i=0; i<delay_count; i++)
        {
          for(int j=0; j<delay_count; j++)
            std::cout << Binv[i][j] << " ";
          std::cout << std::endl;
        }
        std::cout << "report acceptRow end " << std::endl;
        */
      }

      inline void updateInvMat(Matrix<T>& Ainv, bool transfer_to_host = true)
      {
        // update the inverse matrix
        if( delay_count>0 )
        {
          const int norb=Ainv_gpu.rows();
          const int lda_Binv=Binv_gpu.cols();
          //    cudaError_t error(cudaSuccess);
          //    int deviceid;
          //    cudaGetDevice(&deviceid);
              //std::cout << "delay_count = " << delay_count << " id = " << deviceid << std::endl;
          //error = cudaMemcpyAsync(U_gpu.data(), U.data(), norb*delay_count*sizeof(T), cudaMemcpyHostToDevice, hstream);
          //    if( error!=cudaSuccess ) std::cout <<"debug error 1 code " << error << std::endl;
          //BLAS::gemm('T', 'N', delay_count, norb, norb, cone, U.data(), norb, Ainv.data(), norb, czero, tempMat.data(), lda_Binv);
          //    cudaMemPrefetchAsync(Ainv.data(), Ainv.size()*sizeof(T), 0, hstream);
          cuBLAS::gemm(handle, CUBLAS_OP_T, CUBLAS_OP_N, delay_count, norb, norb, constants_gpu.data()+1, U_gpu.data(), norb, Ainv_gpu.data(), norb, constants_gpu.data(), temp_gpu.data(), lda_Binv);
          //error = cudaMemcpyAsync(delay_list_gpu.data(), delay_list.data(), delay_count*sizeof(int), cudaMemcpyHostToDevice, hstream);
          //    if( error!=cudaSuccess ) std::cout <<"debug error 2 code " << error << std::endl;
          //for(int i=0; i<delay_count; i++) tempMat(delay_list[i], i) -= cone;
          applyW_cuda(delay_list_gpu.data(),delay_count, temp_gpu.data(), temp_gpu.cols(), hstream);
          //error = cudaMemcpyAsync(tempMat.data(), temp_gpu.data(), tempMat.size()*sizeof(T), cudaMemcpyDeviceToHost, hstream);
          //    if( error!=cudaSuccess ) std::cout <<"debug error 3 code " << error << std::endl;
          //error = cudaStreamSynchronize(hstream);
          //    if( error!=cudaSuccess ) std::cout <<"debug error 4 code " << error << std::endl;
              //std::cout << "debug tempMat " << tempMat << std::endl;
          //cudaMemcpyAsync(Binv_gpu.data(), Binv.data(), lda_Binv*delay_count*sizeof(T), cudaMemcpyHostToDevice, hstream);
          //BLAS::gemm('N', 'N', norb, delay_count, delay_count, cone, V.data(), norb, Binv.data(), lda_Binv, czero, U.data(), norb);
          cuBLAS::gemm(handle, CUBLAS_OP_N, CUBLAS_OP_N, norb, delay_count, delay_count, constants_gpu.data()+1, V_gpu.data(), norb, Binv_gpu.data(), lda_Binv, constants_gpu.data(), U_gpu.data(), norb);
          //error = cudaMemcpyAsync(U.data(), U_gpu.data(), norb*delay_count*sizeof(T), cudaMemcpyDeviceToHost, hstream);
          //    if( error!=cudaSuccess ) std::cout <<"debug error 5 code " << error << std::endl;
          //error = cudaStreamSynchronize(hstream);
          //    if( error!=cudaSuccess ) std::cout <<"debug error 6 code " << error << std::endl;
              //std::cout << "debug U " << U << std::endl;
          //BLAS::gemm('N', 'N', norb, norb, delay_count, -cone, U.data(), norb, tempMat.data(), lda_Binv, cone, Ainv.data(), norb);
          cuBLAS::gemm(handle, CUBLAS_OP_N, CUBLAS_OP_N, norb, norb, delay_count, constants_gpu.data()+2, U_gpu.data(), norb, temp_gpu.data(), lda_Binv, constants_gpu.data()+1, Ainv_gpu.data(), norb);
          delay_count = 0;
          Ainv_row_ind = -1;
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

