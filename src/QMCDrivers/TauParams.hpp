//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2022 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//                    Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////


#ifndef QMCPLUSPLUS_TAU_PARAMS_HPP
#define QMCPLUSPLUS_TAU_PARAMS_HPP

#include <MCCoords.hpp>

namespace qmcplusplus
{
/** Object to encapsulate appropriate tau derived parameters
 *  for a particular CoordsType specialization
 */
template<typename RT, CoordsType CT>
struct TauParams;

template<typename RT>
struct TauParams<RT, CoordsType::POS>
{
  RT tauovermass;
  RT oneover2tau;
  RT sqrttau;

  TauParams(RT tau, RT grp_inv_mass, RT spin_mass)
  {
    tauovermass = tau * grp_inv_mass;
    oneover2tau = 0.5 / (tauovermass);
    sqrttau     = std::sqrt(tauovermass);
  }
};

template<typename RT>
struct TauParams<RT, CoordsType::POS_SPIN>
{
  RT tauovermass;
  RT oneover2tau;
  RT sqrttau;

  RT spin_tauovermass;
  RT spin_oneover2tau;
  RT spin_sqrttau;

  TauParams(RT tau, RT grp_inv_mass, RT spin_mass)
  {
    tauovermass = tau * grp_inv_mass;
    oneover2tau = 0.5 / (tauovermass);
    sqrttau     = std::sqrt(tauovermass);

    spin_tauovermass = tauovermass / spin_mass;
    spin_oneover2tau = 0.5 / (spin_tauovermass);
    spin_sqrttau     = std::sqrt(spin_tauovermass);
  }
};

template<CoordsType CT, typename RT>
TauParams<RT, CT> operator*(const TauParams<RT, CT>& taus, RT val)
{
  TauParams<RT, CT> out = taus;
  out.tauovermass *= val;
  out.oneover2tau *= val;
  out.sqrttau *= std::sqrt(val);
  if constexpr (CT == CoordsType::POS_SPIN)
  {
    out.spin_tauovermass *= val;
    out.spin_oneover2tau *= val;
    out.spin_sqrttau *= std::sqrt(val);
  }
  return out;
}

template<CoordsType CT, typename RT>
TauParams<RT, CT> operator*(RT val, const TauParams<RT, CT>& taus)
{
  return taus * val;
}

} // namespace qmcplusplus

#endif
