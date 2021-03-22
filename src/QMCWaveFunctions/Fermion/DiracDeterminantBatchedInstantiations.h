//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2021 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Lab
//
// File created by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Lab
//////////////////////////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_DIRAC_DETERMINANT_BATCHED_INSTANTIATIONS_HPP
#define QMCPLUSPLUS_DIRAC_DETERMINANT_BATCHED_INSTANTIATIONS_HPP

#include "DiracDeterminantBatched.h"

namespace qmcplusplus
{
#if defined(ENABLE_CUDA) && defined(ENABLE_OFFLOAD)

extern template void DiracDeterminantDetails::mw_recomputeDispatch<
    MatrixDelayedUpdateCUDA<QMCTraits::ValueType, QMCTraits::QTFull::ValueType>>(
    const RefVectorWithLeader<WaveFunctionComponent>& wfc_list,
    const RefVectorWithLeader<ParticleSet>& p_list,
    const std::vector<bool>& recompute);

extern template struct DiracDeterminantBatchedMultiWalkerResource<
    MatrixDelayedUpdateCUDA<QMCTraits::ValueType, QMCTraits::QTFull::ValueType>>;
extern template class DiracDeterminantBatched<
    MatrixDelayedUpdateCUDA<QMCTraits::ValueType, QMCTraits::QTFull::ValueType>>;

#endif

#endif


extern template struct DiracDeterminantBatchedMultiWalkerResource<
    MatrixUpdateOMPTarget<QMCTraits::ValueType, QMCTraits::QTFull::ValueType>>;
extern template class DiracDeterminantBatched<
    MatrixUpdateOMPTarget<QMCTraits::ValueType, QMCTraits::QTFull::ValueType>>;

#endif
}
