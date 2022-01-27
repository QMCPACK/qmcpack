//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2022 developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//////////////////////////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_MCCOORDS_HPP
#define QMCPLUSPLUS_MCCOORDS_HPP

#include "Configuration.h"
#include "type_traits/complex_help.hpp"

#include <vector>

namespace qmcplusplus
{

enum class CoordsTypes
{
  RS,
  RSSPINS
};

template<CoordsTypes MCT = CoordsTypes::RS>
struct MCCoords
{
  static constexpr CoordsTypes mct = MCT;
  // This cleans up some other code.
  void resize(const std::size_t size);
  std::vector<QMCTraits::PosType> positions;
};

template<>
struct MCCoords<CoordsTypes::RSSPINS>
{
  static constexpr CoordsTypes mct = CoordsTypes::RSSPINS;
  // This cleans up some other code.
  void resize(const std::size_t size)
  {
    positions.resize(size);
    spins.resize(size);
  }
  std::vector<QMCTraits::PosType> positions;
  std::vector<QMCTraits::FullPrecRealType> spins;
};

template<CoordsTypes CT = CoordsTypes::RS>
struct MCCIt
{
  std::vector<QMCTraits::PosType>::iterator it_positions;
};

template<>
struct MCCIt<CoordsTypes::RSSPINS>
{
  std::vector<QMCTraits::PosType>::iterator it_positions;
  std::vector<std::complex<QMCTraits::RealType>>::iterator it_spins;
};

/** Object to encapsulate appropriate tau derived values
 *  for a particular MCCoords specialization
 */
template<typename Real, CoordsTypes CT = CoordsTypes::RS>
struct Taus
{
  Real tauovermass;
  Real oneover2tau;
  Real sqrttau;
  Taus(Real tau, Real grp_inv_mass)
  {
    tauovermass = tau * grp_inv_mass;
    oneover2tau = 0.5 / (tauovermass);
    sqrttau     = std::sqrt(tauovermass);
  }
};

template<typename Real>
struct Taus<Real, CoordsTypes::RSSPINS> : public Taus<Real, CoordsTypes::RS>
{
  using Base = Taus<Real, CoordsTypes::RS>;
  Real spin_tauovermass;
  Real spin_oneover2tau;
  Real spin_sqrttau;
  Taus(Real tau, Real grp_inv_mass, Real spin_mass) : Base(tau, grp_inv_mass)
  {
    spin_tauovermass = Base::tauovermass / spin_mass;
    spin_oneover2tau = 0.5 / (spin_tauovermass);
    spin_sqrttau     = std::sqrt(spin_tauovermass);
  }
};

/** This is just one possible factory for Taus based on MCCoordsTypes
 */
template<CoordsTypes CT, typename... ARGS>
auto makeTaus(MCCoords<CT>& mc_coords, const ARGS&... args)
{
  using Real = double;
  if constexpr (CT == CoordsTypes::RS)
  {
    return Taus<Real, CT>(args...);
  }
  else if constexpr (CT == CoordsTypes::RSSPINS)
  {
    return Taus<Real, CT>(args...);
  }
}

extern template struct MCCoords<CoordsTypes::RS>;
extern template struct MCCoords<CoordsTypes::RSSPINS>;
} // namespace qmcplusplus

#endif
