//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2021 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//////////////////////////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_ESTIMATOR_TESTING_H
#define QMCPLUSPLUS_ESTIMATOR_TESTING_H

#include "ParticleSet.h"
#include "OperatorEstBase.h"

namespace qmcplusplus
{

class SpeciesSet;

namespace testing
{
using POLT    = PtclOnLatticeTraits;
using Lattice = POLT::ParticleLayout;

enum class SpeciesCases
{
  GOOD,
  NO_MEMBERSIZE
};

Lattice makeTestLattice();
SpeciesSet makeSpeciesSet(const SpeciesCases species_case);

/** break encapsulation of data_ by  OperatorEstBase
 *  only for testing!
 */
class OEBAccessor
{
public:
  // break naming rule to make std::vector which we assume is the type of OperatorEstBase::Data
  using value_type = OperatorEstBase::Data::value_type;
  OEBAccessor(OperatorEstBase& oeb);
  value_type& operator[](size_t pos);

private:
  OperatorEstBase& oeb_;
};

} // namespace testing
} // namespace qmcplusplus
#endif
