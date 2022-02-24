//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2022 QMCPACK developers.
//
// File developed by: Cody A. Melton, cmelton@sandia.gov, Sandia National Laboratories
//
//
// File created by: Cody A. Melton, cmelton@sandia.gov, Sandia National Laboratories
//////////////////////////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_DRIVEROPERATIONS_HPP
#define QMCPLUSPLUS_DRIVEROPERATIONS_HPP

#include "Particle/MCCoords.hpp"
#include "TauParams.hpp"
#include <algorithm>

namespace qmcplusplus
{
/// multiplies MCCoords by sqrt(tau)
template<typename RT, CoordsType CT>
MCCoords<CT> scaleBySqrtTau(const TauParams<RT, CT>& taus, const MCCoords<CT>& coords)
{
  MCCoords<CT> out(coords.positions.size());
  std::transform(coords.positions.begin(), coords.positions.end(), out.positions.begin(),
                 [st = taus.sqrttau](const QMCTraits::PosType& pos) { return st * pos; });
  if constexpr (CT == CoordsType::POS_SPIN)
    std::transform(coords.spins.begin(), coords.spins.end(), out.spins.begin(),
                   [st = taus.spin_sqrttau](const QMCTraits::FullPrecRealType& spin) { return st * spin; });
  return out;
}

/** calulates forward GreenFunction from displacements stored in MCCoords
 * [param,out] log_gf
 */
template<CoordsType CT>
void updateForwardLogGreensFunction(const MCCoords<CT>& coords, std::vector<QMCTraits::RealType>& log_gf)
{
  assert(coords.positions.size() == log_gf.size());
  std::transform(coords.positions.begin(), coords.positions.end(), log_gf.begin(), [](const QMCTraits::PosType& pos) {
    constexpr QMCTraits::RealType mhalf(-0.5);
    return mhalf * dot(pos, pos);
  });
  if constexpr (CT == CoordsType::POS_SPIN)
    std::transform(coords.spins.begin(), coords.spins.end(), log_gf.begin(), log_gf.begin(),
                   [](const QMCTraits::FullPrecRealType& spin, const QMCTraits::RealType& loggf) {
                     constexpr QMCTraits::RealType mhalf(-0.5);
                     return loggf + mhalf * spin * spin;
                   });
}

/** calculates reverse Green Function from displacements stored in MCCoords
 * [param, out] log_gb
 */
template<typename RT, CoordsType CT>
void updateReverseLogGreensFunction(const MCCoords<CT>& coords,
                                    const TauParams<RT, CT>& taus,
                                    std::vector<QMCTraits::RealType>& log_gb)
{
  assert(coords.positions.size() == log_gb.size());
  std::transform(coords.positions.begin(), coords.positions.end(), log_gb.begin(),
                 [halfovertau = taus.oneover2tau](const QMCTraits::PosType& pos) {
                   return -halfovertau * dot(pos, pos);
                 });
  if constexpr (CT == CoordsType::POS_SPIN)
    std::transform(coords.spins.begin(), coords.spins.end(), log_gb.begin(), log_gb.begin(),
                   [halfovertau = taus.spin_oneover2tau](const QMCTraits::FullPrecRealType& spin,
                                                         const QMCTraits::RealType& loggb) {
                     return loggb - halfovertau * spin * spin;
                   });
}

#endif
} // namespace qmcplusplus
