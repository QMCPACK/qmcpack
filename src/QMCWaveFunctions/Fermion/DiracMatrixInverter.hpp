//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2024 QMCPACK developers.
//
// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_DIRAC_MATRIX_INVERTER_H
#define QMCPLUSPLUS_DIRAC_MATRIX_INVERTER_H

#include "type_traits/template_types.hpp"
#include "type_traits/complex_help.hpp"
#include "OhmmsPETE/OhmmsVector.h"
#include "OhmmsPETE/OhmmsMatrix.h"
#include "Queue.hpp"
#include "OMPTarget/OffloadAlignedAllocators.hpp"
#include "ResourceCollection.h"

namespace qmcplusplus
{

template<typename VALUE_INVERTER, typename VALUE>
class DiracMatrixInverter : public Resource
{
  using InverterPrecReal = RealAlias<VALUE_INVERTER>;
  using LogValue         = std::complex<InverterPrecReal>;

  template<typename T>
  using OffloadPinnedMatrix = Matrix<T, OffloadPinnedAllocator<T>>;
  template<typename T>
  using OffloadPinnedVector = Vector<T, OffloadPinnedAllocator<T>>;

public:
  DiracMatrixInverter(const std::string& name) : Resource(name) {}

  virtual void mw_invert_transpose(compute::QueueBase& queue,
                                   const RefVector<const OffloadPinnedMatrix<VALUE>>& a_mats,
                                   const RefVector<OffloadPinnedMatrix<VALUE>>& inv_a_mats,
                                   OffloadPinnedVector<LogValue>& log_values) = 0;
};

} // namespace qmcplusplus

#endif // QMCPLUSPLUS_DIRAC_MATRIX_INVERTER_H
