//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source
// License. See LICENSE file in top directory for details.
//
// Copyright (c) 2022 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//////////////////////////////////////////////////////////////////////////////////////

#include "MCCoordsT.hpp"

namespace qmcplusplus
{
template <typename T>
void
MCCoordsT<T, CoordsType::POS>::getSubset(const std::size_t offset,
    const std::size_t size, MCCoordsT<T, CoordsType::POS>& out) const
{
    std::copy_n(positions.begin() + offset, size, out.positions.begin());
}

template <typename T>
MCCoordsT<T, CoordsType::POS>&
MCCoordsT<T, CoordsType::POS>::operator+=(
    const MCCoordsT<T, CoordsType::POS>& rhs)
{
    assert(positions.size() == rhs.positions.size());
    std::transform(positions.begin(), positions.end(), rhs.positions.begin(),
        positions.begin(),
        [](const PosType& x, const PosType& y) { return x + y; });
    return *this;
}

template <typename T>
void
MCCoordsT<T, CoordsType::POS_SPIN>::getSubset(const std::size_t offset,
    const std::size_t size, MCCoordsT<T, CoordsType::POS_SPIN>& out) const
{
    std::copy_n(positions.begin() + offset, size, out.positions.begin());
    std::copy_n(spins.begin() + offset, size, out.spins.begin());
}

template <typename T>
MCCoordsT<T, CoordsType::POS_SPIN>&
MCCoordsT<T, CoordsType::POS_SPIN>::operator+=(
    const MCCoordsT<T, CoordsType::POS_SPIN>& rhs)
{
    assert(positions.size() == rhs.positions.size());
    std::transform(positions.begin(), positions.end(), rhs.positions.begin(),
        positions.begin(),
        [](const PosType& x, const PosType& y) { return x + y; });
    std::transform(spins.begin(), spins.end(), rhs.spins.begin(), spins.begin(),
        [](const FullPrecRealType& x, const FullPrecRealType& y) {
            return x + y;
        });
    return *this;
}

template struct MCCoordsT<double, CoordsType::POS>;
template struct MCCoordsT<double, CoordsType::POS_SPIN>;
template struct MCCoordsT<float, CoordsType::POS>;
template struct MCCoordsT<float, CoordsType::POS_SPIN>;
template struct MCCoordsT<std::complex<double>, CoordsType::POS>;
template struct MCCoordsT<std::complex<double>, CoordsType::POS_SPIN>;
template struct MCCoordsT<std::complex<float>, CoordsType::POS>;
template struct MCCoordsT<std::complex<float>, CoordsType::POS_SPIN>;
} // namespace qmcplusplus
