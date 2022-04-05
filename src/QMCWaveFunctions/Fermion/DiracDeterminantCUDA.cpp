//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Ying Wai Li, yingwaili@ornl.gov, Oak Ridge National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


/**@file DiracDeterminantCUDA.h
 * @brief Declaration of DiracDeterminantCUDA with a S(ingle)P(article)O(rbital)Set
 */

#include "DiracDeterminantCUDA.h"
#include "CUDA_legacy/cuda_inverse.h"
#include "QMCWaveFunctions/detail/CUDA_legacy/determinant_update.h"
#include "QMCWaveFunctions/detail/CUDA_legacy/delayed_update.h"
#include "Numerics/DeterminantOperators.h"
#include <unistd.h>

// #define DEBUG_DELAYED
// #define USE_TRSM

namespace qmcplusplus
{
DiracDeterminantCUDA::DiracDeterminantCUDA(std::unique_ptr<SPOSet>&& spos, int first, int last)
    : DiracDeterminantBase("DiracDeterminantCUDA", std::move(spos), first, last),
      UpdateJobList_d("DiracDeterminant::UpdateJobList_d"),
      srcList_d("DiracDeterminant::srcList_d"),
      destList_d("DiracDeterminant::destList_d"),
      AList_d("DiracDeterminant::AList_d"),
      AinvList_d("DiracDeterminant::AinvList_d"),
      newRowList_d("DiracDeterminant::newRowList_d"),
      LemmaList_d("DiracDeterminant::LemmaList_d"),
      LemmaLUList_d("DiracDeterminant::LemmaLUList_d"),
      LemmaInvList_d("DiracDeterminant::LemmaInvList_d"),
      AinvUList_d("DiracDeterminant::AinvUList_d"),
      AinvDeltaList_d("DiracDeterminant::AinvDeltaList_d"),
      AinvColkList_d("DiracDeterminant::AinvColkList_d"),
      gradLaplList_d("DiracDeterminant::gradLaplList_d"),
      newGradLaplList_d("DiracDeterminant::newGradLaplList_d"),
      AWorkList_d("DiracDeterminant::AWorkList_d"),
      AinvWorkList_d("DiracDeterminant::AinvWorkList_d"),
      GLList_d("DiracDeterminant::GLList_d"),
      PivotArray_d("DiracDeterminant::PivotArray_d"),
      infoArray_d("DiracDeterminant::infoArray_d"),
      ratio_d("DiracDeterminant::ratio_d"),
      gradLapl_d("DiracDeterminant::gradLapl_d"),
      iatList_d("DiracDeterminant::iatList_d"),
      NLrowBuffer_d("DiracDeterminant::NLrowBuffer_d"),
      SplineRowList_d("DiracDeterminant::SplineRowList_d"),
      RatioRowList_d("DiracDeterminant::RatioRowList_d"),
      NLposBuffer_d("DiracDeterminant::NLposBuffer_d"),
      NLAinvList_d("DiracDeterminant::NLAinvList_d"),
      NLnumRatioList_d("DiracDeterminant::NLnumRatioList_d"),
      NLelecList_d("DiracDeterminant::NLelecList_d"),
      NLratioList_d("DiracDeterminant::NLratioList_d")
{
  for (int i = 0; i < 2; ++i)
    NLratios_d[i] = gpu::device_vector<CTS::ValueType>("DiracDeterminant::NLratios_d");
}

/////////////////////////////////////
// Vectorized evaluation functions //
/////////////////////////////////////

void CheckAlign(void* p, std::string var)
{
  if ((unsigned long)p & 63)
    app_error() << "CUDA alignment error for variable " << var << "!\n";
}

void DiracDeterminantCUDA::det_lookahead(MCWalkerConfiguration& W,
                                         std::vector<ValueType>& psi_ratios,
                                         std::vector<GradType>& grad,
                                         std::vector<ValueType>& lapl,
                                         int iat,
                                         int k,
                                         int kd,
                                         int nw)
{
  bool klinear = W.getklinear();
  if ((!klinear) &&
      (k ==
       0)) // this only needs to be calculated once and then the columns can be manually updated upon rejection (update function below does this in the update_onemove kernel)
  {
    // make sure the evaluate function kernels are finished
    cudaEventRecord(gpu::syncEvent, gpu::kernelStream);
    cudaEventSynchronize(gpu::syncEvent);
    // calculate new lemma matrix
    cublas_lemma_mats(gpu::cublasHandle, AinvList_d.data(), newRowList_d.data(), LemmaList_d.data(), AinvUList_d.data(),
                      kd, W.getkstart(), NumPtcls, nw, RowStride);
  }
  auto& walkers = W.WalkerList;
  // copy lemma (since it's only calculated once every k blocks) to lemma_lu (for an updated LU decomposition)
  // and copy the k-th row of the A inverse matrix into Ainvcolk (apart from the name that's the place to put it)
  copy_delayed(LemmaLUList_d.data(), LemmaList_d.data(), AinvColkList_d.data(), AinvDeltaList_d.data(), k, kd,
               RowStride, nw);
  // Calculate k-th row of updated A^-1
  if (k > 0)
  {
#ifdef USE_TRSM
    multi_row_copy(AWorkList_d.data(), AinvUList_d.data(), k + 1, W.getkstart() + k, 1, RowStride, nw);
#else
    int curr_ainvu_row = AinvUOffset + W.getkstart() + k; // points to the current row in AinvU
    for (int iw = 0; iw < nw; iw++)
      AinvWorkList[iw] = &(walkers[iw]->cuda_DataSet.data()[curr_ainvu_row]); // don't want to overwrite AinvU list
    AinvWorkList_d.asyncCopy(AinvWorkList);
#endif // use_trsm
    cublas_smw_update(gpu::cublasHandle, AinvDeltaList_d.data(), AinvColkList_d.data(), AinvWorkList_d.data(),
                      AWorkList_d.data(), // <- AinvWork takes the place of A^-1*dU (in the USE_TRSM case it's unused)
                      LemmaInvList_d.data(),
                      LemmaLUList_d.data(), // <- LemmaInv is not needed for USE_TRSM
                      NULL,
                      infoArray_d.data(), // <- important not to use pivoting here
                      k + 1, kd, 1, NumPtcls, nw, RowStride);
  }
  // calculate and collect ratios, gradients, and laplacians
  calc_gradlapl_and_collect(LemmaLUList_d.data(), AinvColkList_d.data(), &(newGradLaplList_d.data()[k * nw]),
                            ratio_d.data(), k, kd, NumPtcls, RowStride, nw);
  if (klinear)
  {
    // Copy back to host
    ratio_host.asyncCopy(ratio_d);
    cudaEventRecord(gpu::ratioSyncDiracEvent, gpu::memoryStream);
  }
  else
  {
    // store output
    ratio_host = ratio_d;
    GradType g;
    for (int iw = 0; iw < nw; iw++)
    {
      psi_ratios[iw + k * nw] *= ratio_host[5 * iw];
      g = GradType(ratio_host[5 * iw + 1], ratio_host[5 * iw + 2], ratio_host[5 * iw + 3]);
      if (k == 0)
      {
#ifdef QMC_COMPLEX
        CTS::ValueType invR = std::conj(ratio_host[5 * iw]) / std::norm(ratio_host[5 * iw]);
#else
        CTS::ValueType invR = 1.0 / ratio_host[5 * iw];
#endif
        g *= invR;
      }
      grad[iw + k * nw] += g;
      lapl[iw + k * nw] += psi_ratios[iw + k * nw];
    }
  }
}

void DiracDeterminantCUDA::update(MCWalkerConfiguration* W,
                                  std::vector<Walker_t*>& walkers,
                                  int iat,
                                  std::vector<bool>* acc,
                                  int k)
{
  int kdelay = 0;
  if (W)
    kdelay = W->getkDelay();
  if (kdelay > 1)
  {
    auto& allwalkers = W->WalkerList;
    int nw           = acc->size();
    int kstart       = W->getkstart();
    int kupdate      = W->getkupdate();

    if (UpdateList.size() < nw)
      UpdateList.resize(nw);
    // Since we need to do some work in case of both accepted and rejected moves
    // the idea is to put them in UpdateList sequentially so a single kernel launch
    // can take care of both scenarios
    //
    // Put accepted moves in the list
    int ws       = walkers.size();
    int accepted = ws;
    for (int iw = 0; iw < ws; iw++)
      UpdateList[iw] = walkers[iw]->cuda_DataSet.data();
    // Put rejected moves in the list
    if (nw - accepted > 0)
    {
      for (int iw = 0; iw < nw; iw++)
        if ((*acc)[iw] == false)
          UpdateList[ws++] = allwalkers[iw]->cuda_DataSet.data();
    }
    UpdateList_d = UpdateList;
    // kernel containing behavior for acceptance and rejection, side benefit: no memory copying needed anymore
    // -> also updates A^-1*dU * lemma^-1 for next round's gradient in "linear" (aka w/ drift) case
    update_onemove(UpdateList_d.data(), newRowOffset + k * RowStride, AOffset + (kstart + k) * RowStride,
                   newGradLaplOffset + 4 * k * RowStride, gradLaplOffset + 4 * (kstart + k) * RowStride, AinvUOffset,
                   LemmaOffset + k * kdelay, LemmaInvOffset, (W->getklinear() && (k < kupdate)) * AWorkOffset, accepted,
                   k, kstart, kdelay, RowStride, ws);
    if (k + 1 == kupdate) // time to update the inverse
    {
#ifndef USE_TRSM
      multi_copy(LemmaLUList_d.data(), LemmaList_d.data(), kupdate * kdelay, nw);
#else // using the triangular solver should be the default
      copy_update_block(LemmaLUList_d.data(), LemmaList_d.data(), AWorkList_d.data(), AinvDeltaList_d.data(), k, kdelay,
                        RowStride, nw);
#endif
      cublas_smw_update(gpu::cublasHandle, AinvDeltaList_d.data(), AinvList_d.data(), AinvUList_d.data(),
                        AWorkList_d.data(), LemmaInvList_d.data(), LemmaLUList_d.data(), PivotArray_d.data(),
                        infoArray_d.data(), kupdate, kdelay, NumPtcls, NumPtcls, nw, RowStride);
    }
  }
  else
  {
    if (UpdateList.size() < walkers.size())
      UpdateList.resize(walkers.size());
    for (int iw = 0; iw < walkers.size(); iw++)
      UpdateList[iw] = walkers[iw]->cuda_DataSet.data();
    UpdateList_d.asyncCopy(UpdateList);
    update_inverse_cuda(UpdateList_d.data(), iat - FirstIndex, AOffset, AinvOffset, newRowOffset, AinvDeltaOffset,
                        AinvColkOffset, NumPtcls, RowStride, walkers.size());
    // Copy temporary gradients and laplacians into matrix
    int gradoff = 4 * (iat - FirstIndex) * RowStride;
    multi_copy(UpdateList_d.data(), gradLaplOffset + gradoff, newGradLaplOffset, 4 * RowStride, walkers.size());
  }
  // multi_copy (destList_d.data(), srcList_d.data(),
  // 		4*RowStride, walkers.size());
  // if (AList.size() < walkers.size())
  //   resizeLists(walkers.size());
  // int gradoff = 4*(iat-FirstIndex)*RowStride;
  // if (UpdateJobList.size() != walkers.size()) {
  //   UpdateJobList.resize(walkers.size());
  //   srcList.resize(walkers.size());
  //   destList.resize(walkers.size());
  //   UpdateJobList_d.resize(walkers.size());
  //   srcList_d.resize(walkers.size());
  //   destList_d.resize(walkers.size());
  // }
  // for (int iw=0; iw<walkers.size(); iw++) {
  //   Walker_t::cuda_Buffer_t &data = walkers[iw]->cuda_DataSet;
  //   updateJob &job = UpdateJobList[iw];
  //   job.A            = &(data.data()[AOffset]);
  //   job.Ainv         = &(data.data()[AinvOffset]);
  //   job.newRow       = &(data.data()[newRowOffset]);
  //   job.AinvDelta    = &(data.data()[AinvDeltaOffset]);
  //   job.AinvColk     = &(data.data()[AinvColkOffset]);
  //   job.gradLapl     = &(data.data()[gradLaplOffset+gradoff]);
  //   job.newGradLapl  = &(data.data()[newGradLaplOffset]);
  //   job.iat          = iat-FirstIndex;
  //   CheckAlign (job.A, "A");
  //   CheckAlign (job.Ainv, "Ainv");
  //   CheckAlign (job.newRow, "newRow");
  //   CheckAlign (job.AinvDelta, "AinvDelta");
  //   CheckAlign (job.AinvColk,  "AinvColk");
  //   CheckAlign (job.gradLapl,  "gradLapl");
  //   CheckAlign (job.newGradLapl, "newGradLapl");
  //   destList[iw]    = &(data.data()[gradLaplOffset+gradoff]);
  //   srcList[iw]     = &(data.data()[newGradLaplOffset]);
  // }
  // // Copy pointers to the GPU
  // UpdateJobList_d = UpdateJobList;
  // srcList_d  = srcList;
  // destList_d = destList;
  // // Call kernel wrapper function
  // CudaRealType dummy;
  // update_inverse_cuda(UpdateJobList_d.data(), dummy,
  // 			NumPtcls, RowStride, walkers.size());
  // // Copy temporary gradients and laplacians into matrix
  // multi_copy (destList_d.data(), srcList_d.data(),
  // 		4*RowStride, walkers.size());


// Check matrix inversion + fixing the CUDA debugging code
#ifdef DEBUG_CUDA
  CTS::ValueType Ainv[NumPtcls][RowStride], A[NumPtcls][RowStride];
  CTS::ValueType new_row[RowStride]; //Ainv_delta[NumPtcls];
  //for (int iw=0; iw<walkers.size(); iw++)
  //{
  int iw                        = 0;
  Walker_t::cuda_Buffer_t& data = walkers[iw]->cuda_DataSet;
  cudaMemcpy(A, &(data.data()[AOffset]), NumPtcls * RowStride * sizeof(CTS::ValueType), cudaMemcpyDeviceToHost);
  cudaMemcpy(Ainv, &(data.data()[AinvOffset]), NumPtcls * RowStride * sizeof(CTS::ValueType), cudaMemcpyDeviceToHost);
  cudaMemcpy(new_row, &(data.data()[newRowOffset]), RowStride * sizeof(CTS::ValueType), cudaMemcpyDeviceToHost);
  // for (int i=0; i<NumPtcls; i++)
  //  	cerr << "new_row(" << i << ") = " << new_row[i]
  // 	     << "  old_row = " << A[iat-FirstIndex][i] << std::endl;
  // float Ainv_k[NumPtcls];
  // for (int i=0; i<NumPtcls; i++) {
  // 	Ainv_delta[i] = 0.0f;
  // 	Ainv_k[i] = Ainv[i][iat-FirstIndex];
  // 	for (int j=0; j<NumPtcls; j++)
  // 	  Ainv_delta[i] += Ainv[j][i] * (new_row[j] - A[(iat-FirstIndex)][j]);
  // }
  // double prefact = 1.0/(1.0+Ainv_delta[iat-FirstIndex]);
  // for (int i=0; i<NumPtcls; i++)
  // 	for (int j=0; j<NumPtcls; j++)
  // 	  Ainv[i][j] += prefact * Ainv_delta[j]*Ainv_k[i];
  // for (int j=0; j<NumPtcls; j++)
  // 	A[iat-FirstIndex][j] = new_row[j];
  for (int i = 0; i < NumPtcls; i++)
    for (int j = 0; j < NumPtcls; j++)
    {
      CTS::ValueType val = 0.0;
      for (int k = 0; k < NumPtcls; k++)
        val += Ainv[i][k] * A[k][j];
      if (i == j && (std::abs(std::real(val) - 1.0) > 1.0e-2))
        std::cerr << "Error in inverse at (i,j) = (" << i << ", " << j << ")  val = " << val << "  walker = " << iw
                  << " of " << walkers.size() << std::endl;
      else if ((i != j) && (std::abs(val) > 1.0e-2))
        std::cerr << "Error in inverse at (i,j) = (" << i << ", " << j << ")  val = " << val << "  walker = " << iw
                  << " of " << walkers.size() << std::endl;
    }
    //}
#endif
}


void DiracDeterminantCUDA::update(const std::vector<Walker_t*>& walkers, const std::vector<int>& iat_list)
{
  int N = walkers.size();
  if (UpdateList.size() < N)
    UpdateList.resize(N);
  if (iatList.size() < N)
    iatList.resize(N);
  if (srcList.size() < N)
    srcList.resize(N);
  if (destList.size() < N)
    destList.resize(N);
  for (int iw = 0; iw < walkers.size(); iw++)
  {
    UpdateList[iw] = walkers[iw]->cuda_DataSet.data();
    iatList[iw]    = iat_list[iw] - FirstIndex;
    srcList[iw]    = walkers[iw]->cuda_DataSet.data() + newGradLaplOffset;
    destList[iw]   = walkers[iw]->cuda_DataSet.data() + gradLaplOffset + 4 * iatList[iw] * RowStride;
  }
  UpdateList_d = UpdateList;
  iatList_d    = iatList;
  srcList_d    = srcList;
  destList_d   = destList;
  update_inverse_cuda(UpdateList_d.data(), iatList_d.data(), AOffset, AinvOffset, newRowOffset, AinvDeltaOffset,
                      AinvColkOffset, NumPtcls, RowStride, walkers.size());
  // Copy temporary gradients and laplacians into matrix
  multi_copy(destList_d.data(), srcList_d.data(), 4 * RowStride, walkers.size());
  // if (AList.size() < walkers.size())
  //   resizeLists(walkers.size());
  // if (UpdateJobList.size() != walkers.size()) {
  //   UpdateJobList.resize(walkers.size());
  //   srcList.resize(walkers.size());
  //   destList.resize(walkers.size());
  //   UpdateJobList_d.resize(walkers.size());
  //   srcList_d.resize(walkers.size());
  //   destList_d.resize(walkers.size());
  // }
  // for (int iw=0; iw<walkers.size(); iw++) {
  //   int gradoff = 4*(iat_list[iw]-FirstIndex)*RowStride;
  //   Walker_t::cuda_Buffer_t &data = walkers[iw]->cuda_DataSet;
  //   updateJob &job = UpdateJobList[iw];
  //   job.A            = &(data.data()[AOffset]);
  //   job.Ainv         = &(data.data()[AinvOffset]);
  //   job.newRow       = &(data.data()[newRowOffset]);
  //   job.AinvDelta    = &(data.data()[AinvDeltaOffset]);
  //   job.AinvColk     = &(data.data()[AinvColkOffset]);
  //   job.gradLapl     = &(data.data()[gradLaplOffset+gradoff]);
  //   job.newGradLapl  = &(data.data()[newGradLaplOffset]);
  //   job.iat          = iat_list[iw] - FirstIndex;
  //   CheckAlign (job.A, "A");
  //   CheckAlign (job.Ainv, "Ainv");
  //   CheckAlign (job.newRow, "newRow");
  //   CheckAlign (job.AinvDelta, "AinvDelta");
  //   CheckAlign (job.AinvColk,  "AinvColk");
  //   CheckAlign (job.gradLapl,  "gradLapl");
  //   CheckAlign (job.newGradLapl, "newGradLapl");
  //   destList[iw]    = &(data.data()[gradLaplOffset+gradoff]);
  //   srcList[iw]     = &(data.data()[newGradLaplOffset]);
  // }
  // // Copy pointers to the GPU
  // UpdateJobList_d = UpdateJobList;
  // srcList_d  = srcList;
  // destList_d = destList;
  // // Call kernel wrapper function
  // CudaRealType dummy;
  // update_inverse_cuda(UpdateJobList_d.data(), dummy,
  // 			NumPtcls, RowStride, walkers.size());
  // // Copy temporary gradients and laplacians into matrix
  // multi_copy (destList_d.data(), srcList_d.data(),
  // 		4*RowStride, walkers.size());
#ifdef DEBUG_CUDA
  CTS::ValueType Ainv[NumPtcls][RowStride], A[NumPtcls][RowStride];
  CTS::ValueType new_row[RowStride]; //Ainv_delta[NumPtcls];
  //for (int iw=0; iw<walkers.size(); iw++)
  //{
  int iw                        = 0;
  Walker_t::cuda_Buffer_t& data = walkers[iw]->cuda_DataSet;
  cudaMemcpy(A, &(data.data()[AOffset]), NumPtcls * RowStride * sizeof(CTS::ValueType), cudaMemcpyDeviceToHost);
  cudaMemcpy(Ainv, &(data.data()[AinvOffset]), NumPtcls * RowStride * sizeof(CTS::ValueType), cudaMemcpyDeviceToHost);
  cudaMemcpy(new_row, &(data.data()[newRowOffset]), RowStride * sizeof(CTS::ValueType), cudaMemcpyDeviceToHost);
  // for (int i=0; i<NumPtcls; i++)
  //  	cerr << "new_row(" << i << ") = " << new_row[i]
  // 	     << "  old_row = " << A[iat-FirstIndex][i] << std::endl;
  // float Ainv_k[NumPtcls];
  // for (int i=0; i<NumPtcls; i++) {
  // 	Ainv_delta[i] = 0.0f;
  // 	Ainv_k[i] = Ainv[i][iat-FirstIndex];
  // 	for (int j=0; j<NumPtcls; j++)
  // 	  Ainv_delta[i] += Ainv[j][i] * (new_row[j] - A[(iat-FirstIndex)][j]);
  // }
  // double prefact = 1.0/(1.0+Ainv_delta[iat-FirstIndex]);
  // for (int i=0; i<NumPtcls; i++)
  // 	for (int j=0; j<NumPtcls; j++)
  // 	  Ainv[i][j] += prefact * Ainv_delta[j]*Ainv_k[i];
  // for (int j=0; j<NumPtcls; j++)
  // 	A[iat-FirstIndex][j] = new_row[j];
  for (int i = 0; i < NumPtcls; i++)
    for (int j = 0; j < NumPtcls; j++)
    {
      CTS::ValueType val = 0.0;
      for (int k = 0; k < NumPtcls; k++)
        val += Ainv[i][k] * A[k][j];
      if (i == j && (std::abs(std::real(val) - 1.0) > 1.0e-2))
        std::cerr << "Error in inverse at (i,j) = (" << i << ", " << j << ")  val = " << val << "  walker = " << iw
                  << " of " << walkers.size() << std::endl;
      else if ((i != j) && (std::abs(val) > 1.0e-2))
        std::cerr << "Error in inverse at (i,j) = (" << i << ", " << j << ")  val = " << val << "  walker = " << iw
                  << " of " << walkers.size() << std::endl;
    }
    //}
#endif
}


void DiracDeterminantCUDA::recompute(MCWalkerConfiguration& W, bool firstTime)
{
  std::vector<Walker_t*> walkers;
  walkers.reserve(W.WalkerList.size());
  for (auto& walker_uptr : W.WalkerList)
    walkers.push_back(walker_uptr.get());
  // HACK HACK HACK
  //     app_log() << "Before recompute:\n";
  //     gpu::host_vector<CUDA_PRECISION> host_data;
  //     for (int iw=0; iw<walkers.size(); iw++) {
  //        Walker_t::cuda_Buffer_t& data = walkers[iw]->cuda_DataSet;
  //        host_data = data;
  //        double sum2 = 0.0;
  //        for (int i=0; i<NumPtcls; i++)
  // 	 for (int j=0; j<NumPtcls; j++) {
  // 	   double val = 0.0;
  // 	   for (int k=0; k<NumPtcls; k++)
  // 	     val += host_data[AinvOffset + i*RowStride+k] *
  // 	       host_data[AOffset+k*RowStride+j];
  // 	   val -= (i==j) ? 1.0 : 0.0;
  // 	   sum2 += val*val;
  // 	 }
  //        app_log() << "iw = " << iw << "  RMS deviation = "
  // 		 << std::sqrt(sum2/(double)(NumPtcls*NumPtcls)) << std::endl;
  //     }
  if (AList.size() < walkers.size())
    resizeLists(walkers.size());
  // Only reevalute the orbitals if this is the first time
  if (firstTime)
  {
    //       int iwsave = walkers.size()-1;
    //       gpu::host_vector<float> old_data, new_data;
    //       old_data = walkers[iwsave]->cuda_DataSet;
    // Recompute A matrices;
    std::vector<PosType> R(walkers.size());
    for (int iat = FirstIndex; iat < LastIndex; iat++)
    {
      int off = (iat - FirstIndex) * RowStride;
      for (int iw = 0; iw < walkers.size(); iw++)
      {
        Walker_t::cuda_Buffer_t& data = walkers[iw]->cuda_DataSet;
        newRowList[iw]                = &(data.data()[AOffset + off]);
        newGradLaplList[iw]           = &(data.data()[gradLaplOffset + 4 * off]);
        R[iw]                         = walkers[iw]->R[iat];
      }
      newRowList_d      = newRowList;
      newGradLaplList_d = newGradLaplList;
      Phi->evaluate(walkers, R, newRowList_d, newGradLaplList_d, RowStride);
    }
    //       new_data = walkers[iwsave]->cuda_DataSet;
    //       for (int i=0; i<NumOrbitals; i++)
    // 	for (int j=0; j<NumOrbitals; j++) {
    // 	  int off = i*RowStride+j + AOffset;
    // 	  double oldA = old_data[off];
    // 	  double newA = new_data[off];
    // 	  if (std::abs(oldA-newA) > 1.0e-9) {
    // 	    char buff[200];
    // 	    snprintf (buff, 200, "(%3d, %3d)  old=%10.6e  new=%10.6e\n",
    // 		      i,j, oldA, newA);
    // 	    app_log() << buff;
    // 	  }
    // 	}
  }
  for (int iw = 0; iw < walkers.size(); iw++)
  {
    Walker_t::cuda_Buffer_t& data = walkers[iw]->cuda_DataSet;
    AList[iw]                     = &(data.data()[AOffset]);
    AinvList[iw]                  = &(data.data()[AinvOffset]);
    AWorkList[iw]                 = &(data.data()[AWorkOffset]);
    AinvWorkList[iw]              = &(data.data()[AinvWorkOffset]);
  }
  AList_d        = AList;
  AinvList_d     = AinvList;
  AWorkList_d    = AWorkList;
  AinvWorkList_d = AinvWorkList;

  // Copy A into Ainv
  //multi_copy (AinvList_d.data(), AList_d.data(),
  //            NumPtcls*RowStride, walkers.size());

  // Invert
  bool useDoublePrecision = true;
  cublas_inverse(gpu::cublasHandle, AList_d.data(), AinvList_d.data(), AWorkList_d.data(), AinvWorkList_d.data(),
                 PivotArray_d.data(), infoArray_d.data(), NumPtcls, RowStride, walkers.size(), useDoublePrecision);

  // checking inversion status
  infoArray_host = infoArray_d;
  bool failed    = false;
  for (int iw = 0; iw < walkers.size(); iw++)
    if (infoArray_host[iw] != 0 || infoArray_host[iw + walkers.size()] != 0)
    {
      failed = true;
      fprintf(stderr, "cublas_inverse failed on walker %d, getrf error %d, getri error %d.\n", iw, infoArray_host[iw],
              infoArray_host[iw + walkers.size()]);
      char name[1000];
      gethostname(name, 1000);
      fprintf(stderr, "Offending hostname = %s\n", name);
      int dev;
      cudaGetDevice(&dev);
      fprintf(stderr, "Offending device = %d\n", dev);
    }
  if (failed)
    abort();

#ifdef DEBUG_CUDA
  CTS::ValueType Ainv[NumPtcls][RowStride], A[NumPtcls][RowStride];
  //for (int iw=0; iw<walkers.size(); iw++)
  //{
  int iw                        = 0;
  Walker_t::cuda_Buffer_t& data = walkers[iw]->cuda_DataSet;
  cudaMemcpy(A, &(data.data()[AOffset]), NumPtcls * RowStride * sizeof(CTS::ValueType), cudaMemcpyDeviceToHost);
  cudaMemcpy(Ainv, &(data.data()[AinvOffset]), NumPtcls * RowStride * sizeof(CTS::ValueType), cudaMemcpyDeviceToHost);

  FILE *f1, *f2;
  f1 = fopen("A.dat", "a");
  f2 = fopen("Ainv.dat", "a");
  for (int i = 0; i < NumPtcls; i++)
    for (int j = 0; j < NumPtcls; j++)
    {
#ifdef QMC_COMPLEX
      fprintf(f1, "%5d %5d %20.15e %20.15e\n", i, j, A[i][j].real(), A[i][j].imag());
      fprintf(f2, "%5d %5d %20.15e %20.15e\n", i, j, Ainv[i][j].real(), Ainv[i][j].imag());
#else
      fprintf(f1, "%5d %5d %20.15e\n", i, j, A[i][j]);
      fprintf(f2, "%5d %5d %20.15e\n", i, j, Ainv[i][j]);
#endif
    }
  //}
  fclose(f1);
  fclose(f2);
#endif
  // HACK HACK HACK
  //     app_log() << "After recompute:\n";
  //     double A[NumPtcls*NumPtcls], work[NumPtcls*NumPtcls];
  //     int piv[NumPtcls];
  //     for (int iw=0; iw<walkers.size(); iw++) {
  //        Walker_t::cuda_Buffer_t& data = walkers[iw]->cuda_DataSet;
  //        host_data = data;
  //        double sum2 = 0.0;
  //        for (int i=0; i<NumPtcls; i++)
  // 	 for (int j=0; j<NumPtcls; j++) {
  // 	   A[i*NumPtcls+j] = host_data[AinvOffset + i*RowStride +j];
  // 	   double val = 0.0;
  // 	   for (int k=0; k<NumPtcls; k++)
  // 	     val += host_data[AinvOffset + i*RowStride+k] *
  // 	       host_data[AOffset+k*RowStride+j];
  // 	   val -= (i==j) ? 1.0 : 0.0;
  // 	   sum2 += val*val;
  // 	 }
  //        double phase;
  //        double error = std::sqrt(sum2/(double)(NumPtcls*NumPtcls));
  //        double logdet =InvertWithLog(A,NumPtcls,NumOrbitals,work,piv,phase);
  //        app_log() << "iw = " << iw << "  RMS deviation = " << error
  // 		 << "  LogDet = " << logdet << std::endl;
  //        if (error > 0.05) {
  // 	 FILE *CUDA = fopen ("Ainv_cuda", "w");
  // 	 FILE *CPU  = fopen ("Ainv_CPU","w");
  // 	 for (int i=0; i<NumPtcls; i++) {
  // 	   for (int j=0; j<NumPtcls; j++) {
  // 	     fprintf (CUDA, "%16.12e ", host_data[AinvOffset + i*RowStride+j]);
  // 	     fprintf (CPU,  "%16.12e ", A[i*NumPtcls+j]);
  // 	   }
  // 	   fprintf (CUDA, "\n");
  // 	   fprintf (CPU,  "\n");
  // 	 }
  // 	 fclose(CUDA);
  // 	 fclose(CPU);
  // 	 abort();
  //        }
  //     }
  // cuda_inverse_many(AinvList_d.data(), workList_d.data(),
  // 		      NumOrbitals, walkers.size());
}

void DiracDeterminantCUDA::addLog(MCWalkerConfiguration& W, std::vector<RealType>& logPsi)
{
  std::vector<Walker_t*> walkers;
  walkers.reserve(W.WalkerList.size());
  for (auto& walker_uptr : W.WalkerList)
    walkers.push_back(walker_uptr.get());
  if (AList.size() < walkers.size())
    resizeLists(walkers.size());
  std::vector<PosType> R(walkers.size());
  // Fill in the A matrix row by row
  for (int iat = FirstIndex; iat < LastIndex; iat++)
  {
    int off = (iat - FirstIndex) * RowStride;
    for (int iw = 0; iw < walkers.size(); iw++)
    {
      Walker_t::cuda_Buffer_t& data = walkers[iw]->cuda_DataSet;
      newRowList[iw]                = &(data.data()[AOffset + off]);
      gradLaplList[iw]              = &(data.data()[gradLaplOffset + 4 * off]);
      R[iw]                         = walkers[iw]->R[iat];
    }
    newRowList_d   = newRowList;
    gradLaplList_d = gradLaplList;
    Phi->evaluate(walkers, R, newRowList_d, gradLaplList_d, RowStride);
  }
  // Now, compute determinant
  for (int iw = 0; iw < walkers.size(); iw++)
  {
    Walker_t::cuda_Buffer_t& data = walkers[iw]->cuda_DataSet;
    gpu::host_vector<CTS::ValueType> host_data;
    Vector<CTS::ValueType> A(NumPtcls * NumOrbitals);
    host_data = data;
    for (int i = 0; i < NumPtcls; i++)
      for (int j = 0; j < NumOrbitals; j++)
        A[i * NumOrbitals + j] = host_data[AOffset + i * RowStride + j];
    logPsi[iw] += std::log(std::abs(Invert(A.data(), NumPtcls, NumOrbitals)));
    int N       = NumPtcls;
    bool passed = true;
    for (int i = 0; i < NumPtcls; i++)
    {
      for (int j = 0; j < NumOrbitals; j++)
        host_data[AinvOffset + i * RowStride + j] = A[i * NumOrbitals + j];
      for (int j = NumOrbitals; j < RowStride; j++)
        host_data[AinvOffset + i * RowStride + j] = 0.0;
    }
    data = host_data;
    // for (int i=0; i<N; i++)
    // 	for (int j=0; j<N; j++) {
    // 	  double val = 0.0;
    // 	  for (int k=0; k<N; k++) {
    // 	    double aval = host_data[AOffset+i*N+k];
    // 	    double ainv = host_data[AinvOffset+k*N+j];
    // 	    val += aval * ainv;
    // 	  }
    // 	  if (i == j) {
    // 	    if (std::abs(val - 1.0) > 1.0e-2) {
    // 	      app_error() << "Error in inverse, (i,j) = " << i << ", " << j << ".\n";
    // 	      passed = false;
    // 	    }
    // 	  }
    // 	  else
    // 	    if (std::abs(val) > 1.0e-2) {
    // 	      app_error() << "Error in inverse, (i,j) = " << i << ", " << j << ".\n";
    // 	      passed = false;
    // 	    }
    // 	}
    // if (!passed)
    // 	app_log() << (passed ? "Passed " : "Failed " ) << "inverse test.\n";
  }
}

void DiracDeterminantCUDA::calcGradient(MCWalkerConfiguration& W, int iat, int k, std::vector<GradType>& grad)
{
  auto& walkers = W.WalkerList;
  int kd        = W.getkDelay();
  int kstk      = 0;
  if (kd)
    kstk = iat - FirstIndex;
  int nw = walkers.size();
  if (newRowList.size() < nw * W.getkblocksize())
    resizeLists(nw, W.getkblocksize());
  for (int iw = 0; iw < nw; iw++)
  {
    Walker_t::cuda_Buffer_t& data = walkers[iw]->cuda_DataSet;
    GLList[iw]                    = &(data.data()[gradLaplOffset + 4 * kstk * RowStride]);
    for (int k2 = 0; k2 < W.getkblocksize(); k2++)
      AinvList[iw + k2 * nw] = &(data.data()[AinvOffset]);
    if (kd && (k == 0))
      AinvColkList[iw] = &(data.data()[AinvColkOffset]);
  }
  if (k == 0)
    AinvList_d.asyncCopy(AinvList);
  GLList_d.asyncCopy(GLList);
  if (kd)
  {
    if (k == 0)
    {
      AinvColkList_d = AinvColkList;
      W.setklinear();
      // just copy the k-th row of the (just updated) A inverse matrix into Ainvcolk (apart from the name that's the place to put it)
      multi_row_copy(AinvColkList_d.data(), AinvList_d.data(), RowStride, kstk, 1, RowStride, nw);
    }
    else // k>0, calculate k-th row of updated A^-1
    {
      // just copy the k-th row of the (just updated) A inverse matrix into Ainvcolk (apart from the name that's the place to put it)
      multi_row_copy(AinvColkList_d.data(), AinvList_d.data(), RowStride, kstk, 1, RowStride, nw);
      cublas_ainv_row(gpu::cublasHandle, AinvDeltaList_d.data(), AWorkList_d.data(), AinvColkList_d.data(), k, NumPtcls,
                      nw, RowStride);
    }
    // calculate and collect gradients only
    calc_gradient_delayed(AinvColkList_d.data(), GLList_d.data(), ratio_d.data(), NumPtcls, RowStride, nw);
  }
  else
  {
    calc_gradient(AinvList_d.data(), GLList_d.data(), ratio_d.data(), NumOrbitals, RowStride, iat - FirstIndex, nw);
  }
  gpu::streamsSynchronize();
  ratio_host.asyncCopy(ratio_d);
  cudaEventRecord(gpu::gradientSyncDiracEvent, gpu::memoryStream);
}

void DiracDeterminantCUDA::addGradient(MCWalkerConfiguration& W, int iat, std::vector<GradType>& grad)
{
  auto& walkers = W.WalkerList;
  cudaEventSynchronize(gpu::gradientSyncDiracEvent);
  for (int iw = 0; iw < walkers.size(); iw++)
  {
#ifdef DEBUG_DELAYED
    fprintf(stderr, "walker #%i: grad = (", iw);
#endif
    for (int dim = 0; dim < OHMMS_DIM; dim++)
    {
#ifdef DEBUG_DELAYED
      if (dim > 0)
        fprintf(stderr, ", ");
      fprintf(stderr, "%f", ratio_host[3 * iw + dim]);
#endif
      grad[iw][dim] += ratio_host[3 * iw + dim];
    }
#ifdef DEBUG_DELAYED
    fprintf(stderr, ")\n");
#endif
  }
#ifdef CUDA_DEBUG3
  if (NumOrbitals == 31)
  {
    gpu::host_vector<CTS::RealType> host_data;
    std::vector<CTS::RealType> cpu_ratios(walkers.size(), 0.0f);
    for (int iw = 0; iw < walkers.size(); iw++)
    {
      host_data = walkers[iw]->cuda_DataSet;
      for (int iorb = 0; iorb < NumOrbitals; iorb++)
      {
        cpu_ratios[iw] += host_data[AinvOffset + RowStride * iorb + iat - FirstIndex] *
            host_data[gradLaplOffset + 4 * RowStride * (iat - FirstIndex) + iorb + RowStride];
      }
      fprintf(stderr, "CPU grad = %10.6e   GPU grad = %10.6e\n", cpu_ratios[iw], grad[iw][1]);
    }
  }
#endif
}

void DiracDeterminantCUDA::ratio(MCWalkerConfiguration& W, int iat, std::vector<ValueType>& psi_ratios)
{
  std::vector<Walker_t*> walkers;
  walkers.reserve(W.WalkerList.size());
  for (auto& walker_uptr : W.WalkerList)
    walkers.push_back(walker_uptr.get());
  if (AList.size() < walkers.size())
    resizeLists(walkers.size());
  // First evaluate orbitals
  for (int iw = 0; iw < walkers.size(); iw++)
  {
    Walker_t::cuda_Buffer_t& data = walkers[iw]->cuda_DataSet;
    AinvList[iw]                  = &(data.data()[AinvOffset]);
    newRowList[iw]                = &(data.data()[newRowOffset]);
  }
  newRowList_d = newRowList;
  AinvList_d   = AinvList;
  Phi->evaluate(walkers, W.Rnew, newRowList_d);
  // Now evaluate ratios
  determinant_ratios_cuda(AinvList_d.data(), newRowList_d.data(), ratio_d.data(), NumPtcls, RowStride, iat - FirstIndex,
                          walkers.size());
  // Copy back to host
  ratio_host = ratio_d;
  for (int iw = 0; iw < psi_ratios.size(); iw++)
    psi_ratios[iw] *= ratio_host[iw];
}

void DiracDeterminantCUDA::ratio(MCWalkerConfiguration& W,
                                 int iat,
                                 std::vector<ValueType>& psi_ratios,
                                 std::vector<GradType>& grad)
{}

// The gradient is (\nabla psi(Rnew))/psi(Rnew)
// The laplaican is (\nabla^2 psi(Rnew))/psi(Rew)
void DiracDeterminantCUDA::ratio(MCWalkerConfiguration& W,
                                 int iat,
                                 std::vector<ValueType>& psi_ratios,
                                 std::vector<GradType>& grad,
                                 std::vector<ValueType>& lapl)
{
  const int N      = W.Rnew.size();
  const int kd     = W.getkDelay();
  const int kstart = W.getkstart();
  const int nw     = W.WalkerList.size();
  std::vector<Walker_t*> walkers;
  walkers.reserve(nw);
  for (auto& walker_uptr : W.WalkerList)
    walkers.push_back(walker_uptr.get());
  if (AinvList.size() < N)
    resizeLists(nw, W.getkblocksize());
  //    if (iat-FirstIndex == 0) {
  // First evaluate orbitals
  for (int iw = 0; iw < N; iw++)
  {
    int k = iw / nw;
    Walker_t::cuda_Buffer_t& data =
        walkers[iw % nw]
            ->cuda_DataSet; // in the k-delay scheme there are now N=nw*kblocksize elements, iw mod nw gives the current walker of each k-block
    AinvList[iw]        = &(data.data()[AinvOffset]);
    newRowList[iw]      = &(data.data()[newRowOffset + k * RowStride]);
    newGradLaplList[iw] = &(data.data()[newGradLaplOffset + 4 * k * RowStride]);
    if (kd && (k == 0))
    {
      AList[iw]         = &(data.data()[AOffset + kstart * RowStride]);
      AinvColkList[iw]  = &(data.data()[AinvColkOffset]);
      AinvDeltaList[iw] = &(
          data.data()
              [AinvOffset +
               kstart]); // for delayed update need to only multiply with k x N sub matrix of A^-1, address needs to start at current k block
      AinvUList[iw]    = &(data.data()[AinvUOffset]);
      AWorkList[iw]    = &(data.data()[AWorkOffset]);
      LemmaList[iw]    = &(data.data()[LemmaOffset]);
      LemmaLUList[iw]  = &(data.data()[LemmaLUOffset]);
      LemmaInvList[iw] = &(data.data()[LemmaInvOffset]);
    }
    if (((unsigned long)AinvList[iw % nw] % 64) || ((unsigned long)newRowList[iw % nw] % 64) ||
        ((unsigned long)newGradLaplList[iw % nw] % 64))
      app_log() << "**** CUDA misalignment!!!! ****\n";
  }
  newRowList_d.asyncCopy(newRowList);
  newGradLaplList_d.asyncCopy(newGradLaplList);
  AinvList_d.asyncCopy(AinvList);
  //    }
  if (kd) // delayed update case
  {
#ifdef DEBUG_DELAYED
    fprintf(stderr, "*** Delayed Update (k start: %i) ***\n", kstart);
#endif
    AList_d.asyncCopy(AList);
    AinvColkList_d.asyncCopy(AinvColkList);
    AinvDeltaList_d.asyncCopy(AinvDeltaList);
    AinvUList_d.asyncCopy(AinvUList);
    AWorkList_d.asyncCopy(AWorkList);
    LemmaList_d.asyncCopy(LemmaList);
    LemmaLUList_d.asyncCopy(LemmaLUList);
    LemmaInvList_d.asyncCopy(LemmaInvList);
    Phi->evaluate(walkers, W.Rnew, newRowList_d, newGradLaplList_d,
                  RowStride); // works with k-delayed positions as well
    // further evaluation (determinant look-ahead) happens in TrialWaveFunction class, copy memory needed later
  }
  else
  {
#ifdef DEBUG_DELAYED
    fprintf(stderr, "*** Old Code Path ***\n");
#endif
    Phi->evaluate(walkers, W.Rnew, newRowList_d, newGradLaplList_d, RowStride);
#ifdef CUDA_DEBUG2
    Vector<ValueType> testPhi(NumOrbitals), testLapl(NumOrbitals);
    Vector<GradType> testGrad(NumOrbitals);
    ParticleSet P;
    P.R.resize(NumPtcls);
    gpu::host_vector<CTS::ValueType> host_vec;
    for (int iw = 0; iw < walkers.size(); iw++)
    {
      host_vec              = walkers[iw]->cuda_DataSet;
      P.R[iat - FirstIndex] = W.Rnew[iw];
      //      Phi->evaluate(P, iat-FirstIndex, testPhi);
      Phi->evaluate(P, iat - FirstIndex, testPhi, testGrad, testLapl);
      for (int iorb = 0; iorb < NumOrbitals; iorb++)
      {
        //if (std::abs(host_vec[newRowOffset+iorb]-testPhi[iorb]) > 1.0e-6)
        //   fprintf (stderr, "CUDA = %1.8e    CPU = %1.8e\n",
        // 	   host_vec[newRowOffset+iorb], testPhi[iorb]);
        fprintf(stderr, "CUDA = %1.8e    CPU = %1.8e\n", host_vec[newGradLaplOffset + 2 * NumOrbitals + iorb],
                testGrad[iorb][2]);
      }
    }
#endif
    determinant_ratios_grad_lapl_cuda(AinvList_d.data(), newRowList_d.data(), newGradLaplList_d.data(), ratio_d.data(),
                                      NumPtcls, RowStride, iat - FirstIndex, walkers.size());
    // Copy back to host
    ratio_host = ratio_d;
#ifdef CUDA_DEBUG
    // Now, check against CPU
    gpu::host_vector<CTS::ValueType> host_data;
    std::vector<CTS::ValueType> cpu_ratios(walkers.size(), 0.0f);
    for (int iw = 0; iw < walkers.size(); iw++)
    {
      host_data = walkers[iw]->cuda_DataSet;
      for (int iorb = 0; iorb < NumOrbitals; iorb++)
      {
        cpu_ratios[iw] += host_data[AinvOffset + RowStride * iorb + iat - FirstIndex] * host_data[newRowOffset + iorb];
      }
      fprintf(stderr, "CPU ratio = %10.6e   GPU ratio = %10.6e\n", cpu_ratios[iw], ratio_host[5 * iw + 0]);
    }
#endif
    // Calculate ratio, gradient and laplacian
    for (int iw = 0; iw < walkers.size(); iw++)
    {
      psi_ratios[iw] *= ratio_host[5 * iw + 0];
      GradType g(ratio_host[5 * iw + 1], ratio_host[5 * iw + 2], ratio_host[5 * iw + 3]);
      grad[iw] += g;
      lapl[iw] += ratio_host[5 * iw + 0];
      // grad[iw] += g / ratio_host[5*iw+0];
      // lapl[iw] += (ratio_host[5*iw+4] - dot(g,g)) / ratio_host[5*iw+0];
    }
#ifdef CUDA_DEBUG
    if (NumOrbitals == 31)
    {
      gpu::host_vector<CTS::ValueType> host_data;
      std::vector<CTS::ValueType> cpu_ratios(walkers.size(), 0.0f);
      for (int iw = 0; iw < walkers.size(); iw++)
      {
        host_data = walkers[iw]->cuda_DataSet;
        for (int iorb = 0; iorb < NumOrbitals; iorb++)
        {
          cpu_ratios[iw] += host_data[AinvOffset + RowStride * iorb + iat - FirstIndex] *
              host_data[newGradLaplOffset + iorb] / ratio_host[5 * iw + 0];
        }
        fprintf(stderr, "ratio CPU grad = %10.6e   GPU grad = %10.6e\n", cpu_ratios[iw], grad[iw][0]);
      }
    }
#endif
  }
}

void DiracDeterminantCUDA::calcRatio(MCWalkerConfiguration& W,
                                     int iat,
                                     std::vector<ValueType>& psi_ratios,
                                     std::vector<GradType>& grad,
                                     std::vector<ValueType>& lapl)
{
  std::vector<Walker_t*> walkers;
  walkers.reserve(W.WalkerList.size());
  for (auto& walker_uptr : W.WalkerList)
    walkers.push_back(walker_uptr.get());
  int nw     = walkers.size();
  int kd     = W.getkDelay();
  int kstart = W.getkstart();
  int k      = W.getkcurr() - (kd > 1);
  if (k < 0)
    k += W.getkupdate();
#ifdef DEBUG_DELAYED
  fprintf(stderr, "k: %i, kcurr: %i, kstart: %i, nw: %i, Rnew size: %lu\n", k, W.getkcurr(), kstart, nw, W.Rnew.size());
#endif
  if (newRowList.size() < nw * W.getkblocksize())
    resizeLists(nw, W.getkblocksize());
  for (int iw = 0; iw < nw; iw++)
  {
    Walker_t::cuda_Buffer_t& data = walkers[iw]->cuda_DataSet;
    for (int k2 = 0; k2 < W.getkblocksize(); k2++) // this also works with the original codepath (kblocksize is 1 then)
    {
      newRowList[iw + k2 * nw]      = &(data.data()[newRowOffset + k2 * RowStride]);
      newGradLaplList[iw + k2 * nw] = &(data.data()[newGradLaplOffset + 4 * k2 * RowStride]);
    }
    if (kd && (k == 0))
    {
      //      AList[iw]         =  &(data.data()[AOffset+kstart*RowStride]);
      AinvDeltaList[iw] = &(
          data.data()
              [AinvOffset +
               kstart]); // for delayed update need to only multiply with k x N sub matrix of A^-1, address needs to start at current k block
      AinvUList[iw]    = &(data.data()[AinvUOffset]);
      AWorkList[iw]    = &(data.data()[AWorkOffset]);
      LemmaList[iw]    = &(data.data()[LemmaOffset]);
      LemmaLUList[iw]  = &(data.data()[LemmaLUOffset]);
      LemmaInvList[iw] = &(data.data()[LemmaInvOffset]);
    }
    if (((unsigned long)AinvList[iw] % 64) || ((unsigned long)newRowList[iw] % 64) ||
        ((unsigned long)newGradLaplList[iw] % 64))
      app_log() << "**** CUDA misalignment!!!! ***\n";
  }
  if (k == 0)
  {
    newRowList_d.asyncCopy(newRowList);
    newGradLaplList_d.asyncCopy(newGradLaplList);
  }
  if (kd && (k == 0))
  {
    AinvDeltaList_d.asyncCopy(AinvDeltaList);
    AinvUList_d.asyncCopy(AinvUList);
    AWorkList_d.asyncCopy(AWorkList);
    LemmaList_d.asyncCopy(LemmaList);
    LemmaLUList_d.asyncCopy(LemmaLUList);
    LemmaInvList_d.asyncCopy(LemmaInvList);
  }
  Phi->evaluate(walkers, W.Rnew, newRowList_d, newGradLaplList_d, RowStride, k, W.getklinear());
  if (kd && W.getklinear())
    calc_lemma_column(AinvList_d.data(), newRowList_d.data(), LemmaList_d.data(), AinvUList_d.data(), k, kd,
                      W.getkstart(), NumPtcls, RowStride, nw);
#ifdef CUDA_DEBUG2
  Vector<ValueType> testPhi(NumOrbitals), testLapl(NumOrbitals);
  Vector<GradType> testGrad(NumOrbitals);
  ParticleSet P;
  P.R.resize(NumPtcls);
  gpu::host_vector<CTS::ValueType> host_vec;
  for (int iw = 0; iw < walkers.size(); iw++)
  {
    host_vec              = walkers[iw]->cuda_DataSet;
    P.R[iat - FirstIndex] = W.Rnew[iw];
    Phi->evaluate(P, iat - FirstIndex, testPhi, testGrad, testLapl);
    for (int iorb = 0; iorb < NumOrbitals; iorb++)
    {
      fprintf(stderr, "CUDA = %1.8e    CPU = %1.8e\n", host_vec[newGradLaplOffset + 2 * NumOrbitals + iorb],
              testGrad[iorb][2]);
    }
  }
#endif
  if (kd == 0)
  {
    determinant_ratios_grad_lapl_cuda(AinvList_d.data(), newRowList_d.data(), newGradLaplList_d.data(), ratio_d.data(),
                                      NumPtcls, RowStride, iat - FirstIndex, walkers.size());
    gpu::streamsSynchronize();
    // Copy back to host
    ratio_host.asyncCopy(ratio_d);
    cudaEventRecord(gpu::ratioSyncDiracEvent, gpu::memoryStream);
  }
#ifdef CUDA_DEBUG
  // Now, check against CPU
  gpu::host_vector<CTS::ValueType> host_data;
  std::vector<CTS::ValueType> cpu_ratios(walkers.size(), 0.0f);
  for (int iw = 0; iw < walkers.size(); iw++)
  {
    host_data = walkers[iw]->cuda_DataSet;
    for (int iorb = 0; iorb < NumOrbitals; iorb++)
    {
      cpu_ratios[iw] += host_data[AinvOffset + RowStride * iorb + iat - FirstIndex] * host_data[newRowOffset + iorb];
    }
    fprintf(stderr, "CPU ratio = %10.6e   GPU ratio = %10.6e\n", cpu_ratios[iw], ratio_host[5 * iw + 0]);
  }
#endif
}

void DiracDeterminantCUDA::addRatio(MCWalkerConfiguration& W,
                                    int iat,
                                    int k,
                                    std::vector<ValueType>& psi_ratios,
                                    std::vector<GradType>& grad,
                                    std::vector<ValueType>& lapl)
{
  auto& walkers = W.WalkerList;
  cudaEventSynchronize(gpu::ratioSyncDiracEvent);
  // Calculate ratio, gradient and laplacian
  for (int iw = 0; iw < walkers.size(); iw++)
  {
    psi_ratios[iw] *= ratio_host[5 * iw + 0];
    GradType g(ratio_host[5 * iw + 1], ratio_host[5 * iw + 2], ratio_host[5 * iw + 3]);
    if (W.getkDelay() && (k == 0))
    {
#ifdef QMC_COMPLEX
      CTS::ValueType invR = std::conj(ratio_host[5 * iw]) / std::norm(ratio_host[5 * iw]);
#else
      CTS::ValueType invR = 1.0 / ratio_host[5 * iw];
#endif
      g *= invR;
    }
    grad[iw] += g;
    lapl[iw] += ratio_host[5 * iw + 0];
#ifdef DEBUG_DELAYED
    fprintf(stderr, "-> walker %i: ratio = %f ; grad = (%f,%f,%f) ; lapl = %f\n", iw, psi_ratios[iw], grad[iw][0],
            grad[iw][1], grad[iw][2], lapl[iw]);
#endif
  }
#ifdef CUDA_DEBUG
  if (NumOrbitals == 31)
  {
    gpu::host_vector<CTS::ValueType> host_data;
    std::vector<CTS::ValueType> cpu_ratios(walkers.size(), 0.0f);
    for (int iw = 0; iw < walkers.size(); iw++)
    {
      host_data = walkers[iw]->cuda_DataSet;
      for (int iorb = 0; iorb < NumOrbitals; iorb++)
      {
        cpu_ratios[iw] += host_data[AinvOffset + RowStride * iorb + iat - FirstIndex] *
            host_data[newGradLaplOffset + iorb] / ratio_host[5 * iw + 0];
      }
      fprintf(stderr, "ratio CPU grad = %10.6e   GPU grad = %10.6e\n", cpu_ratios[iw], grad[iw][0]);
    }
  }
#endif
}


// The gradient is (\nabla psi(Rnew))/psi(Rnew)
// The laplaican is (\nabla^2 psi(Rnew))/psi(Rew)
void DiracDeterminantCUDA::ratio(std::vector<Walker_t*>& walkers,
                                 std::vector<int>& iat_list,
                                 std::vector<PosType>& rNew,
                                 std::vector<ValueType>& psi_ratios,
                                 std::vector<GradType>& grad,
                                 std::vector<ValueType>& lapl)
{
  if (AList.size() < walkers.size())
    resizeLists(walkers.size());
  for (int iw = 0; iw < walkers.size(); iw++)
  {
    Walker_t::cuda_Buffer_t& data = walkers[iw]->cuda_DataSet;
    AinvList[iw]                  = &(data.data()[AinvOffset]);
    newRowList[iw]                = &(data.data()[newRowOffset]);
    newGradLaplList[iw]           = &(data.data()[newGradLaplOffset]);
    if (((unsigned long)AinvList[iw] % 64) || ((unsigned long)newRowList[iw] % 64) ||
        ((unsigned long)newGradLaplList[iw] % 64))
      app_log() << "**** CUDA misalignment!!!! ***\n";
    iatList[iw] = iat_list[iw] - FirstIndex;
  }
  newRowList_d      = newRowList;
  newGradLaplList_d = newGradLaplList;
  AinvList_d        = AinvList;
  iatList_d         = iatList;
  Phi->evaluate(walkers, rNew, newRowList_d, newGradLaplList_d, RowStride);
  determinant_ratios_grad_lapl_cuda(AinvList_d.data(), newRowList_d.data(), newGradLaplList_d.data(), ratio_d.data(),
                                    NumPtcls, RowStride, iatList_d.data(), walkers.size());
  // Copy back to host
  ratio_host = ratio_d;
  // Calculate ratio, gradient and laplacian
  for (int iw = 0; iw < walkers.size(); iw++)
  {
    psi_ratios[iw] *= ratio_host[5 * iw + 0];
    GradType g(ratio_host[5 * iw + 1], ratio_host[5 * iw + 2], ratio_host[5 * iw + 3]);
    grad[iw] += g;
    lapl[iw] += ratio_host[5 * iw + 0];
  }
}

void DiracDeterminantCUDA::gradLapl(MCWalkerConfiguration& W, GradMatrix& grads, ValueMatrix& lapl)
{
  auto& walkers = W.WalkerList;
  int nw        = walkers.size();
#ifdef DEBUG_DELAYED
  fprintf(stderr, "grad/lapl, nw = %i\n", nw);
#endif
  if (AList.size() < nw)
    resizeLists(nw);
  // First evaluate orbitals
  for (int iw = 0; iw < nw; iw++)
  {
    Walker_t::cuda_Buffer_t& data = walkers[iw]->cuda_DataSet;
    AinvList[iw]                  = &(data.data()[AinvOffset]);
    gradLaplList[iw]              = &(data.data()[gradLaplOffset]);
    newGradLaplList[iw]           = &(gradLapl_d.data()[4 * RowStride * iw]);
  }
  AinvList_d        = AinvList;
  gradLaplList_d    = gradLaplList;
  newGradLaplList_d = newGradLaplList;
  calc_grad_lapl(AinvList_d.data(), gradLaplList_d.data(), newGradLaplList_d.data(), NumOrbitals, RowStride, nw);
  // Now copy data into the output matrices
  gradLapl_host = gradLapl_d;
  for (int iw = 0; iw < nw; iw++)
  {
    for (int iat = 0; iat < NumPtcls; iat++)
    {
      GradType g(gradLapl_host[4 * (iw * RowStride + iat) + 0], gradLapl_host[4 * (iw * RowStride + iat) + 1],
                 gradLapl_host[4 * (iw * RowStride + iat) + 2]);
      grads(iw, iat + FirstIndex) += g;
#ifdef QMC_COMPLEX
      //YingWai's trial fix; need to revise the dot product
      // should be converted into "#if defined QMC_COMPLEX" later
      lapl(iw, iat + FirstIndex) += gradLapl_host[4 * (iw * RowStride + iat) + 3] - (CTS::ValueType)dot(g, g);
      ValueType lapl_temp = lapl(iw, iat + FirstIndex);
      if (std::isnan(std::real(lapl_temp)) || std::isnan(std::imag(lapl_temp)))
#else
      lapl(iw, iat + FirstIndex) += gradLapl_host[4 * (iw * RowStride + iat) + 3] - dot(g, g);
      if (std::isnan(lapl(iw, iat + FirstIndex)))
#endif
      {
        fprintf(stderr, "Offending walker = %d\n", iw);
        char name[1000];
        gethostname(name, 1000);
        fprintf(stderr, "Offending hostname = %s\n", name);
        int dev;
        cudaGetDevice(&dev);
        fprintf(stderr, "Offending device = %d\n", dev);
        gpu::host_vector<CTS::ValueType> host_data;
        host_data  = walkers[iw]->cuda_DataSet;
        FILE* Amat = fopen("Amat.dat", "w");
        FILE* Ainv = fopen("Ainv.dat", "w");
        FILE* Lmat = fopen("Alapl.dat", "w");
        FILE* Gmat = fopen("Agrad.dat", "w");
        for (int i = 0; i < NumPtcls; i++)
        {
          for (int j = 0; j < NumPtcls; j++)
          {
#ifdef QMC_COMPLEX
            fprintf(Amat, "%14.8e+%14.8ei ", host_data[AOffset + i * RowStride + j].real(),
                    host_data[AOffset + i * RowStride + j].imag());
            fprintf(Ainv, "%14.8e+%14.8ei ", host_data[AinvOffset + i * RowStride + j].real(),
                    host_data[AinvOffset + i * RowStride + j].imag());
            fprintf(Lmat, "%14.8e+%14.8ei ", host_data[gradLaplOffset + (4 * i + 3) * RowStride + j].real(),
                    host_data[gradLaplOffset + (4 * i + 3) * RowStride + j].imag());
            for (int k = 0; k < 3; k++)
              fprintf(Gmat, "%14.8e+%14.8ei ", host_data[gradLaplOffset + (4 * i + k) * RowStride + j].real(),
                      host_data[gradLaplOffset + (4 * i + k) * RowStride + j].imag());
#else
            fprintf(Amat, "%14.8e ", host_data[AOffset + i * RowStride + j]);
            fprintf(Ainv, "%14.8e ", host_data[AinvOffset + i * RowStride + j]);
            fprintf(Lmat, "%14.8e ", host_data[gradLaplOffset + (4 * i + 3) * RowStride + j]);
            for (int k = 0; k < 3; k++)
              fprintf(Gmat, "%14.8e ", host_data[gradLaplOffset + (4 * i + k) * RowStride + j]);
#endif
          }
          fprintf(Amat, "\n");
          fprintf(Ainv, "\n");
          fprintf(Lmat, "\n");
          fprintf(Gmat, "\n");
        }
        fclose(Amat);
        fclose(Ainv);
        fclose(Lmat);
        fclose(Gmat);
        abort();
        // 	  if (std::isnan(host_data[AinvOffset+i*RowStride+j]))
        // 	    std::cerr << "NAN in inverse at (" << i << "," << j << ")\n";
        std::cerr << "NAN in walker " << iw << ", iat " << iat + FirstIndex
                  << "  grad = " << grads(iw, iat + FirstIndex)
                  << "  lapl = " << gradLapl_host[4 * (iw * RowStride + iat) + 3] << std::endl;
        fprintf(stderr, "grad-lapl row:\n");
        fprintf(stderr, "r = %1.8f %1.8f %1.8f\n", walkers[iw]->R[iat + FirstIndex][0],
                walkers[iw]->R[iat + FirstIndex][1], walkers[iw]->R[iat + FirstIndex][2]);
        for (int orb = 0; orb < NumPtcls; orb++)
#ifdef QMC_COMPLEX
          fprintf(stderr, "%1.10e+%1.10ei %1.10e+%1.10ei %1.10e+%1.10ei %1.10e+%1.10ei \n",
                  host_data[gradLaplOffset + (4 * iat + 0) * RowStride + orb].real(),
                  host_data[gradLaplOffset + (4 * iat + 0) * RowStride + orb].imag(),
                  host_data[gradLaplOffset + (4 * iat + 1) * RowStride + orb].real(),
                  host_data[gradLaplOffset + (4 * iat + 1) * RowStride + orb].imag(),
                  host_data[gradLaplOffset + (4 * iat + 2) * RowStride + orb].real(),
                  host_data[gradLaplOffset + (4 * iat + 2) * RowStride + orb].imag(),
                  host_data[gradLaplOffset + (4 * iat + 3) * RowStride + orb].real(),
                  host_data[gradLaplOffset + (4 * iat + 3) * RowStride + orb].imag());
#else
          fprintf(stderr, "%1.10e %1.10e %1.10e %1.10e \n", host_data[gradLaplOffset + (4 * iat + 0) * RowStride + orb],
                  host_data[gradLaplOffset + (4 * iat + 1) * RowStride + orb],
                  host_data[gradLaplOffset + (4 * iat + 2) * RowStride + orb],
                  host_data[gradLaplOffset + (4 * iat + 3) * RowStride + orb]);
#endif
      }
    }
  }
#ifdef CUDA_DEBUG
  // Now do it on the CPU
  gpu::host_vector<CTS::ValueType> host_data;
  GradMatrix cpu_grads(grads.rows(), grads.cols());
  ValueMatrix cpu_lapl(grads.rows(), grads.cols());
  for (int iw = 0; iw < walkers.size(); iw++)
  {
    fprintf(stderr, "walker #%i:\n", iw);
    host_data = walkers[iw]->cuda_DataSet;
    for (int iat = 0; iat < NumPtcls; iat++)
    {
      cpu_grads(iw, iat + FirstIndex) = GradType();
      cpu_lapl(iw, iat + FirstIndex)  = ValueType();
      for (int iorb = 0; iorb < NumOrbitals; iorb++)
      {
        cpu_lapl(iw, iat + FirstIndex) += host_data[AinvOffset + NumPtcls * iorb + iat] *
            host_data[gradLaplOffset + (4 * iat + 3) * NumOrbitals + iorb];
      }
      fprintf(stderr, "CPU lapl = %10.6e   GPU lapl = %10.6e\n", cpu_lapl(iw, iat + FirstIndex),
              lapl(iw, iat + FirstIndex));
    }
  }
#endif
}

void DiracDeterminantCUDA::NLratios_CPU(MCWalkerConfiguration& W,
                                        std::vector<NLjob>& jobList,
                                        std::vector<PosType>& quadPoints,
                                        std::vector<ValueType>& psi_ratios)
{
  // Phi->evaluate needs to be replaced
  APP_ABORT("DiracDeterminantCUDA::NLratios_CPU is currently disabled.\n");
  auto& walkers = W.WalkerList;
  std::vector<ValueMatrix> Ainv_host;
  int nw = walkers.size();
  Ainv_host.resize(nw);
  int mat_size = NumOrbitals * NumOrbitals * sizeof(CTS::ValueType);
  for (int iw = 0; iw < nw; iw++)
  {
    Ainv_host[iw].resize(NumOrbitals, NumOrbitals);
    ValueType* dest     = &(Ainv_host[iw](0, 0));
    CTS::ValueType* src = &(walkers[iw]->cuda_DataSet.data()[AinvOffset]);
    cudaMemcpy(dest, src, mat_size, cudaMemcpyDeviceToHost);
  }
  std::vector<RealType> phi(NumOrbitals);
  int index = 0;
  for (int ijob = 0; ijob < jobList.size(); ijob++)
  {
    NLjob& job  = jobList[ijob];
    int numQuad = job.numQuadPoints;
    int elec    = job.elec;
    int iw      = job.walker;
    // Check if this electron belongs to this determinant
    if (elec < FirstIndex || elec >= LastIndex)
      index += numQuad;
    else
    {
      for (int iq = 0; iq < numQuad; iq++)
      {
        //The following line should be replaced with Phi->evaluateDetRatios
        //Phi->evaluate(W, quadPoints[index], phi);
        ValueType ratio = 0.0;
        for (int i = 0; i < NumOrbitals; i++)
          ratio += Ainv_host[iw](i, elec - FirstIndex) * phi[i];
        psi_ratios[index] *= ratio;
        index++;
      }
    }
  }
}

void DiracDeterminantCUDA::NLratios(MCWalkerConfiguration& W,
                                    std::vector<NLjob>& jobList,
                                    std::vector<PosType>& quadPoints,
                                    std::vector<ValueType>& psi_ratios)
{
  // DEBUG!
  // std::vector<ValueType> cpu_ratios;
  // cpu_ratios = psi_ratios;
  // NLratios_CPU (W, jobList, quadPoints, cpu_ratios);
  auto& walkers = W.WalkerList;
  int posIndex = 0, numJobs = 0;
  std::vector<PosType> posBuffer[2];
  int rowIndex = 0;
  std::vector<ValueType*> ratio_pointers[2];
  bool hasResults = false;
  for (int i = 0; i < 2; ++i)
  {
    posBuffer[i].clear();
    ratio_pointers[i].clear();
    NLAinvList_host[i].clear();
    NLnumRatioList_host[i].clear();
    NLelecList_host[i].clear();
    NLratioList_host[i].clear();
    RatioRowList_host[i].clear();
  }
  int ijob    = 0;
  int counter = 0;
  //First batch of data (PT)
  for (; ijob < jobList.size(); ++ijob)
  {
    NLjob& job  = jobList[ijob];
    int numQuad = job.numQuadPoints;
    int elec    = job.elec;
    // Check if this electron belongs to this determinant
    if (elec < FirstIndex || elec >= LastIndex)
    {
      posIndex += numQuad;
      continue;
    }
    // Check to see if the buffer is full
    if (rowIndex + numQuad > NLrowBufferRows)
    {
      break;
    }
    int iw = job.walker;
    NLAinvList_host[counter].push_back(&(walkers[iw]->cuda_DataSet.data()[AinvOffset]));
    NLnumRatioList_host[counter].push_back(numQuad);
    NLelecList_host[counter].push_back(job.elec - FirstIndex);
    NLratioList_host[counter].push_back(&(NLratios_d[counter].data()[rowIndex]));
    RatioRowList_host[counter].push_back(&(NLrowBuffer_d.data()[rowIndex * RowStride]));
    for (int iq = 0; iq < numQuad; iq++)
    {
      posBuffer[counter].push_back(quadPoints[posIndex]);
      ratio_pointers[counter].push_back(&(psi_ratios[posIndex]));
      posIndex++;
    }
    rowIndex += numQuad;
    numJobs++;
  }
  //Next batches (PT)
  while (ijob < jobList.size())
  {
    NLjob& job  = jobList[ijob];
    int numQuad = job.numQuadPoints;
    int elec    = job.elec;
    if (rowIndex + numQuad > NLrowBufferRows)
    {
      // Compute orbital rows
      Phi->evaluate(posBuffer[counter], SplineRowList_d);
      // Compute ratios
      NLAinvList_d.asyncCopy(NLAinvList_host[counter]);
      NLnumRatioList_d.asyncCopy(NLnumRatioList_host[counter]);
      NLelecList_d.asyncCopy(NLelecList_host[counter]);
      NLratioList_d.asyncCopy(NLratioList_host[counter]);
      RatioRowList_d.asyncCopy(RatioRowList_host[counter]);
      calc_many_ratios(NLAinvList_d.data(), RatioRowList_d.data(), NLratioList_d.data(), NLnumRatioList_d.data(),
                       NumOrbitals, RowStride, NLelecList_d.data(), numJobs);
      // Write ratios out output vector
      rowIndex = 0;
      numJobs  = 0;
      counter  = (counter + 1) % 2;
    }
    if (hasResults)
    {
      cudaStreamSynchronize(gpu::memoryStream);
      for (int i = 0; i < ratio_pointers[counter].size(); i++)
        *(ratio_pointers[counter][i]) *= NLratios_host[i];
      hasResults = false;
      ratio_pointers[counter].clear();
    }
    for (; ijob < jobList.size(); ijob++)
    {
      NLjob& job  = jobList[ijob];
      int numQuad = job.numQuadPoints;
      int elec    = job.elec;
      // Check if this electron belongs to this determinant
      if (elec < FirstIndex || elec >= LastIndex)
      {
        posIndex += numQuad;
        continue;
      }
      // Check to see if the buffer is full
      if (rowIndex + numQuad > NLrowBufferRows)
      {
        break;
      }
      int iw = job.walker;
      NLAinvList_host[counter].push_back(&(walkers[iw]->cuda_DataSet.data()[AinvOffset]));
      NLnumRatioList_host[counter].push_back(numQuad);
      NLelecList_host[counter].push_back(job.elec - FirstIndex);
      NLratioList_host[counter].push_back(&(NLratios_d[counter].data()[rowIndex]));
      RatioRowList_host[counter].push_back(&(NLrowBuffer_d.data()[rowIndex * RowStride]));
      for (int iq = 0; iq < numQuad; iq++)
      {
        posBuffer[counter].push_back(quadPoints[posIndex]);
        ratio_pointers[counter].push_back(&(psi_ratios[posIndex]));
        posIndex++;
      }
      rowIndex += numQuad;
      numJobs++;
    }
    NLratios_host.asyncCopy(NLratios_d[(counter + 1) % 2]);
    hasResults = true;
    // Reset counters
    posBuffer[(counter + 1) % 2].clear();
    //ratio_pointers[(counter+1)%2].clear();
    NLAinvList_host[(counter + 1) % 2].clear();
    NLnumRatioList_host[(counter + 1) % 2].clear();
    NLelecList_host[(counter + 1) % 2].clear();
    NLratioList_host[(counter + 1) % 2].clear();
    RatioRowList_host[(counter + 1) % 2].clear();
  }
  if (hasResults)
  {
    cudaStreamSynchronize(gpu::memoryStream);
    for (int i = 0; i < ratio_pointers[(counter + 1) % 2].size(); i++)
      *(ratio_pointers[(counter + 1) % 2][i]) *= NLratios_host[i];
  }
  if (posBuffer[counter].size())
  {
    // Compute whatever remains in the buffer
    // Compute orbital rows
    Phi->evaluate(posBuffer[counter], SplineRowList_d);
    // Compute ratios
    NLAinvList_d.asyncCopy(NLAinvList_host[counter]);
    NLnumRatioList_d.asyncCopy(NLnumRatioList_host[counter]);
    NLelecList_d.asyncCopy(NLelecList_host[counter]);
    NLratioList_d.asyncCopy(NLratioList_host[counter]);
    RatioRowList_d.asyncCopy(RatioRowList_host[counter]);
    calc_many_ratios(NLAinvList_d.data(), RatioRowList_d.data(), NLratioList_d.data(), NLnumRatioList_d.data(),
                     NumOrbitals, RowStride, NLelecList_d.data(), numJobs);
    // Write ratios out output vector
    NLratios_host = NLratios_d[counter];
    for (int i = 0; i < ratio_pointers[counter].size(); i++)
      *(ratio_pointers[counter][i]) *= NLratios_host[i];
  }
  // DEBUG DEBUG DEBUG
  // for (int i=0; i<psi_ratios.size(); i++) {
  //   double diff = psi_ratios[i] - cpu_ratios[i];
  //   if (std::abs(diff) > 1.0e-8)
  // 	fprintf (stderr, "i=%d  GPU=%1.12f  CPU=%1.12f  FirstIndex=%d\n",
  // 		 i, psi_ratios[i], cpu_ratios[i], FirstIndex);
  // }
}
} // namespace qmcplusplus
