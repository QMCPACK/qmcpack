/////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2025 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Lab
//////////////////////////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_SPACEGRIDTEST_H
#define QMCPLUSPLUS_SPACEGRIDTEST_H

#include "ValidReferencePointsInput.h"
#include "ValidSpaceGridInput.h"
#include <MinimalParticlePool.h>
#include "OperatorEstBase.h"
#include <NESpaceGrid.h>

namespace qmcplusplus
{
class ReferencePointsInput;
class ParticleSetPool;
class SpaceGridInput;
class NEReferencePoints;
class ParticleSetPool;
namespace testing
{

template<ValidSpaceGridInput::valid VALID>
class SpaceGridEnv
{
public:
  using Input = ValidSpaceGridInput;
  SpaceGridEnv(Communicate* comm);

  SpaceGridEnv(const SpaceGridEnv& env);

  void setDefaultStartingPos();

  /// Canned test data so rng differences don't cause fails.
  std::vector<ParticleSet::ParticlePos> deterministic_rs_;

  UPtr<ReferencePointsInput> rpi_;
  UPtr<SpaceGridInput> sgi_;
  UPtr<NEReferencePoints> ref_points_;
  RefVector<ParticleSet> ref_psets_;
  ParticleSetPool particle_pool_;
  ParticleSet pset_elec_;
  ParticleSet pset_ions_;
};


template<typename REAL>
class NESpaceGridTests
{
public:
  static const auto& getData(const NESpaceGrid<REAL>& nesg) { return nesg.data_; }
  static int getBufferStart(const NESpaceGrid<REAL>& nesg) { return nesg.buffer_start_; }
  static int getBufferEnd(const NESpaceGrid<REAL>& nesg) { return nesg.buffer_end_; }
  static auto* getOdu(const NESpaceGrid<REAL>& nesg) { return nesg.odu_; }
  static int getNDomains(const NESpaceGrid<REAL>& nesg) { return nesg.ndomains_; }
  static int getAxisGridSizes(const NESpaceGrid<REAL>& nesg) { return nesg.ndomains_; }
};

template<typename REAL, ValidSpaceGridInput::valid VALID>
class SpaceGridTest
{
public:
  using MCPWalker = typename OperatorEstBase::MCPWalker;
  SpaceGridTest(const SpaceGridEnv<VALID>& env, int num_walkers, bool generate_test_data = false);
  RefVectorWithLeader<ParticleSet> getPSetList();

private:
  std::vector<ParticleSet> psets_;
  const SpaceGridEnv<VALID>& sgenv_;
  
};

extern template class SpaceGridTest<float, ValidSpaceGridInput::valid::DEFAULT>;
extern template class SpaceGridTest<double, ValidSpaceGridInput::valid::DEFAULT>;
extern template class SpaceGridTest<float, ValidSpaceGridInput::valid::CYLINDRICAL>;
extern template class SpaceGridTest<double, ValidSpaceGridInput::valid::CYLINDRICAL>;

} // namespace testing
} // namespace qmcplusplus

#endif
