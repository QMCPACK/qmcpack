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
} // namespace qmcplusplus

#endif
