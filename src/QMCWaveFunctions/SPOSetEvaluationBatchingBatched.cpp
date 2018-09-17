//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2018 QMCPACK developers.
//
// File developed by: Peter Doak, Oak Ridge National Laboratory
//                      refactored from SPOSet.cpp
//
// File created by: Peter Doak, Oak Ridge National Laboratory
//////////////////////////////////////////////////////////////////////////////////////

#include "QMCWaveFunctions/SPOSetEvaluationBatched.h"
#include "QMCWaveFunctions/SPOSet.h"
#include "QMCWaveFunctions/SPOSetTypeAliases.h"

namespace qmcplusplus
{
/** batched walker implementation */

void
SPOSetEvaluation<Batching::BATCHED>::evaluate (std::vector<SSTA::Walker_t*> &walkers, int iat,
					       gpu::device_vector<QMCT::CudaValueType*> &phi)
{
  app_error() << "Need specialization of vectorized evaluate in SPOSet.\n";
  app_error() << "Required CUDA functionality not implemented. Contact developers.\n";
  abort();
}

void
SPOSetEvaluation<Batching::BATCHED>::evaluate (std::vector<SSTA::Walker_t*> &walkers, std::vector<QMCT::PosType> &new_pos,
					       gpu::device_vector<QMCT::CudaValueType*> &phi)
{
  app_error() << "Need specialization of vectorized evaluate in SPOSet.\n";
  app_error() << "Required CUDA functionality not implemented. Contact developers.\n";
  abort();
}

void
SPOSetEvaluation<Batching::BATCHED>::evaluate (std::vector<SSTA::Walker_t*> &walkers,
					       std::vector<QMCT::PosType> &new_pos,
					       gpu::device_vector<QMCT::CudaValueType*> &phi,
					       gpu::device_vector<QMCT::CudaValueType*> &grad_lapl_list,
                           int row_stride)
{
  app_error() << "Need specialization of vectorized eval_grad_lapl in SPOSet.\n";
  app_error() << "Required CUDA functionality not implemented. Contact developers.\n";
  abort();
}

void
SPOSetEvaluation<Batching::BATCHED>::evaluate (std::vector<QMCT::PosType> &pos, gpu::device_vector<QMCT::CudaRealType*> &phi)
{
  app_error() << "Need specialization of vectorized evaluate "
              << "in SPOSet.\n";
  app_error() << "Required CUDA functionality not implemented. Contact developers.\n";
  abort();
}

void
SPOSetEvaluation<Batching::BATCHED>::evaluate (std::vector<QMCT::PosType> &pos, gpu::device_vector<QMCT::CudaComplexType*> &phi)
{
  app_error() << "Need specialization of vectorized evaluate "
              << "in SPOSet.\n";
  app_error() << "Required CUDA functionality not implemented. Contact developers.\n";
  abort();
}


}
