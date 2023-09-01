//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source
// License. See LICENSE file in top directory for details.
//
// Copyright (c) 2020 QMCPACK developers.
//
// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////

#include "Particle/DynamicCoordinatesT.h"

#include "Particle/RealSpacePositionsT.h"
#include "Particle/RealSpacePositionsTOMPTarget.h"

namespace qmcplusplus
{
/** create DynamicCoordinates based on kind
 */
template <typename T>
std::unique_ptr<DynamicCoordinatesT<T>>
createDynamicCoordinatesT(const DynamicCoordinateKind kind)
{
    if (kind == DynamicCoordinateKind::DC_POS)
        return std::make_unique<RealSpacePositionsT<T>>();
    else if (kind == DynamicCoordinateKind::DC_POS_OFFLOAD)
        return std::make_unique<RealSpacePositionsTOMPTarget<T>>();
    // dummy return
    return std::unique_ptr<RealSpacePositionsT<T>>();
}

template std::unique_ptr<DynamicCoordinatesT<double>>
createDynamicCoordinatesT<double>(const DynamicCoordinateKind kind);
template std::unique_ptr<DynamicCoordinatesT<float>>
createDynamicCoordinatesT<float>(const DynamicCoordinateKind kind);
template std::unique_ptr<DynamicCoordinatesT<std::complex<double>>>
createDynamicCoordinatesT<std::complex<double>>(
    const DynamicCoordinateKind kind);
template std::unique_ptr<DynamicCoordinatesT<std::complex<float>>>
createDynamicCoordinatesT<std::complex<float>>(
    const DynamicCoordinateKind kind);
} // namespace qmcplusplus
