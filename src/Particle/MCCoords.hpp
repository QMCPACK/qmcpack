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

enum class CoordsType
{
  POS,
  POS_SPIN
};

template<CoordsType MCT = CoordsType::POS>
struct MCCoords
{
  // This cleans up some other code.
  void resize(const std::size_t size);
  std::vector<QMCTraits::PosType> positions;
};

template<>
struct MCCoords<CoordsType::POS_SPIN>
{
  // This cleans up some other code.
  void resize(const std::size_t size)
  {
    positions.resize(size);
    spins.resize(size);
  }
  std::vector<QMCTraits::PosType> positions;
  std::vector<QMCTraits::FullPrecRealType> spins;
};

/** Object to encapsulate appropriate tau derived values
 *  for a particular MCCoords specialization
 */
template<typename Real, CoordsType CT = CoordsType::POS>
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
struct Taus<Real, CoordsType::POS_SPIN> : public Taus<Real, CoordsType::POS>
{
  using Base = Taus<Real, CoordsType::POS>;
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

/** Factory function for Taus based on MCCoordsTypes
 *  Note as in previous code value of tau derived values is not full precision.
 */
template<CoordsType CT, typename... ARGS>
auto makeTaus(MCCoords<CT>& mc_coords, const ARGS&... args)
{
  return Taus<QMCTraits::RealType, CT>(args...);
}

extern template struct MCCoords<CoordsType::POS>;
extern template struct MCCoords<CoordsType::POS_SPIN>;
} // namespace qmcplusplus

#endif
