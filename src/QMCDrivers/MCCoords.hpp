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
#include "ParticleBase/RandomSeqGenerator.h"

#include <vector>

namespace qmcplusplus
{

enum class MCCoordsTypes
{
  RS,
  RSSPINS
};

template<MCCoordsTypes MCT = MCCoordsTypes::RS>
struct MCCoords
{
  static constexpr MCCoordsTypes mct = MCT;
  // This cleans up some other code.
  void resize(const std::size_t size);
  std::vector<QMCTraits::PosType> rs;
};

template<>
struct MCCoords<MCCoordsTypes::RSSPINS>
{
  static constexpr MCCoordsTypes mct = MCCoordsTypes::RSSPINS;
  // This cleans up some other code.
  void resize(const std::size_t size)
  {
    rs.resize(size);
    spins.resize(size);
  }
  std::vector<QMCTraits::PosType> rs;
  std::vector<QMCTraits::FullPrecRealType> spins;
};

template<MCCoordsTypes CT = MCCoordsTypes::RS>
struct MCCIt
{
  std::vector<QMCTraits::PosType>::iterator irs;
};

template<>
struct MCCIt<MCCoordsTypes::RSSPINS>
{
  std::vector<QMCTraits::PosType>::iterator irs;
  std::vector<std::complex<QMCTraits::RealType>>::iterator spins;
};

/** Object to encapsulate appropriate tau derived values
 *  for a particular MCCoords specialization
 */
template<typename Real, MCCoordsTypes CT = MCCoordsTypes::RS>
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
struct Taus<Real, MCCoordsTypes::RSSPINS> : public Taus<Real, MCCoordsTypes::RS>
{
  using Base = Taus<Real, MCCoordsTypes::RS>;
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
template<class MCCOORDS, typename... ARGS>
auto makeTaus(MCCOORDS& mc_coords, const ARGS&... args)
{
  using Real = double;
  constexpr MCCoordsTypes mct{std::decay_t<decltype(mc_coords)>::mct};
  if constexpr (mct == MCCoordsTypes::RS)
  {
    return Taus<Real, mct>(args...);
  }
  else if constexpr (mct == MCCoordsTypes::RSSPINS)
  {
    return Taus<Real, mct>(args...);
  }
}

template<class MCC, class RNG>
void generateDeltas(RNG rng, MCC& deltas)
{
  makeGaussRandomWithEngine(deltas.rs, rng);
  if constexpr (std::is_same<MCC, MCCoords<MCCoordsTypes::RSSPINS>>::value)
  {
    makeGaussRandomWithEngine(deltas.spins, rng);
  }
}

extern template struct MCCoords<MCCoordsTypes::RS>;
extern template struct MCCoords<MCCoordsTypes::RSSPINS>;
} // namespace qmcplusplus

#endif
