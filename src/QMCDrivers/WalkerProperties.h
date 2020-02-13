//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2019 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Lab
//
// File created by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Lab
//////////////////////////////////////////////////////////////////////////////////////
#ifndef QMCPLUSPLUS_WALKERPROPERTIES_H
#define QMCPLUSPLUS_WALKERPROPERTIES_H

#include <cstdint>
#include "config.h"

namespace qmcplusplus
{
struct WalkerProperties
{
  /** an enum denoting index of physical properties
 *
 * \todo: this enum and the handling of "Properties" through the class hierarchy
 * is hot garbage. Replace, for now we're just making it safe.
 *
 * LOCALPOTENTIAL being defined last of the Walker "Properties" is special there
 * is likely to be some code lurking that depends on it. Keep until this mess is removed.

 * Believe it or not,you need to modify ParticleSet::initPropertyList to match the walker enum if you add to it.
 * _Don't_, actually work on fixing this.
 */
  enum Indexes : int16_t
  {
    LOGPSI = 0,      /*!< log(std::abs(psi)) instead of square of the many-body wavefunction \f$|\Psi|^2\f$ */
    SIGN,            /*!< value of the many-body wavefunction \f$\Psi(\{R\})\f$ */
    UMBRELLAWEIGHT,  /*!< sum of wavefunction ratios for multiple H and Psi */
    R2ACCEPTED,      /*!< r^2 for accepted moves */
    R2PROPOSED,      /*!< r^2 for proposed moves */
    DRIFTSCALE,      /*!< scaling value for the drift */
    ALTERNATEENERGY, /*!< alternatelocal energy, the sum of all the components */
    LOCALENERGY,     /*!< local energy, the sum of all the components */
    LOCALPOTENTIAL,  /*!< local potential energy = local energy - kinetic energy */
    NUMPROPERTIES,
    MAXPROPERTIES = WALKER_MAX_PROPERTIES /*!< the number of properties */
  };
};
} // namespace qmcplusplus

#endif /* QMCPLUSPLUS_WALKERPROPERTIES_H */
