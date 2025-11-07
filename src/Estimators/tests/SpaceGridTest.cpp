/////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2025 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Lab
//////////////////////////////////////////////////////////////////////////////////////

#include "SpaceGridTest.hpp"
#include "catch.hpp"
#include "ValidSpaceGridInput.h"
#include <ParticleSetPool.h>
#include "GenerateRandomParticleSets.h"
#include "MockGoldWalkerElements.h"

namespace qmcplusplus
{
namespace testing
{

template<ValidSpaceGridInput::valid VALID>
void SpaceGridEnv<VALID>::setDefaultStartingPos()
{
  deterministic_rs_ = {
      {
          {0.5784564681, -0.7039106488, -1.577452572},
          {0.1348299253, -0.4989691164, -0.04496135769},
          {0.6176573272, 1.249270236, 1.424806025},
          {0.7662747402, 1.008570699, 2.579472429},
          {-0.9434967865, -1.191380189, 1.512721148},
          {0.1492291146, -1.135160203, -1.0874863},
          {4.041292876, 0.2656179458, 4.717080159},
          {0.8551090328, 1.974297541, 2.262264067},
      },
      {
          {-1.908520753, -0.4984436363, -0.5514718161},
          {2.712796299, 0.6255089594, -1.318439055},
          {-0.4716097595, 1.366919373, 1.498546443},
          {4.382191721, 1.209536342, 3.090714299},
          {0.7769128459, 1.093712452, 0.9271239813},
          {-0.7492559163, -0.2648530733, -0.5142864663},
          {3.021315343, 2.232213776, 3.297755805},
          {3.409625477, 1.769879933, 1.223916846},
      },
      {
          {0.4319681346, -1.406019016, 0.5802752151},
          {-1.000397723, 1.340041298, 0.3961784416},
          {0.8525314777, 3.244903308, 0.9471815766},
          {2.138718055, 0.6618596072, 3.518989693},
          {0.1339941753, -1.422134011, -1.574110295},
          {-1.10359286, -0.6955904109, 0.4248506388},
          {2.931111605, 1.601985236, 0.6234291488},
          {1.661872231, 1.923736727, 2.501715031},
      },
      {
          {0.4720635372, 0.2900674782, 1.034044291},
          {-0.595477182, -0.3307560041, -0.3155522407},
          {1.302329658, -0.07223011972, 2.670965585},
          {2.053418512, 2.439064166, 1.918027246},
          {-0.04573295481, -0.8756146251, -0.5210703527},
          {-1.775062015, -0.3894774667, -0.1562153414},
          {1.061934015, -0.2976955493, 1.890430724},
          {1.133206854, 1.987782224, 0.7621546904},
      },
  };
}

template<>
void SpaceGridEnv<ValidSpaceGridInput::valid::CYLINDRICAL>::setDefaultStartingPos()
{}

template<ValidSpaceGridInput::valid VALID>
SpaceGridEnv<VALID>::SpaceGridEnv(Communicate* comm)
    : particle_pool_(MinimalParticlePool::make_diamondC_1x1x1(comm)),
      pset_elec_(*(particle_pool_.getParticleSet("e"))),
      pset_ions_(*(particle_pool_.getParticleSet("ion")))
{
  // Setup particleset
  // particle positions must be inside the unit cell
  setDefaultStartingPos();

  Libxml2Document doc;
  bool okay = doc.parseFromString(Input::xml[VALID]);
  CHECK(okay);

  xmlNodePtr node = doc.getRoot();
  sgi_            = std::make_unique<SpaceGridInput>(node);

  using RPInput = ValidReferencePointsInputs;
  Libxml2Document doc2;
  bool okay2       = doc.parseFromString(RPInput::xml[RPInput::CELL]);
  xmlNodePtr node2 = doc.getRoot();
  rpi_             = std::make_unique<ReferencePointsInput>(node2);
  ref_psets_.push_back(pset_ions_);
  ref_points_ = std::make_unique<NEReferencePoints>(*rpi_, pset_elec_, ref_psets_);
}

template<>
SpaceGridEnv<ValidSpaceGridInput::valid::CYLINDRICAL>::SpaceGridEnv(Communicate* comm)
    : particle_pool_(MinimalParticlePool::make_H2(comm)),
      pset_elec_(*(particle_pool_.getParticleSet("e"))),
      pset_ions_(*(particle_pool_.getParticleSet("ion")))
{
  setDefaultStartingPos();

  Libxml2Document doc;
  bool okay       = doc.parseFromString(Input::xml[ValidSpaceGridInput::valid::CYLINDRICAL]);
  xmlNodePtr node = doc.getRoot();
  sgi_            = std::make_unique<SpaceGridInput>(node);

  using RPInput = ValidReferencePointsInputs;
  Libxml2Document doc2;
  bool okay2 = doc.parseFromString(RPInput::xml[RPInput::CELL]);

  xmlNodePtr node2 = doc.getRoot();
  rpi_             = std::make_unique<ReferencePointsInput>(node2);
  ref_psets_.push_back(pset_ions_);
  ref_points_ = std::make_unique<NEReferencePoints>(*rpi_, pset_elec_, ref_psets_);
}

template<typename REAL, ValidSpaceGridInput::valid VALID>
SpaceGridTest<REAL, VALID>::SpaceGridTest(const SpaceGridEnv<VALID>& env, int num_walkers, bool generate_test_data)
    : sgenv_(env)
{
  psets_ = testing::generateRandomParticleSets(sgenv_.pset_elec_, sgenv_.pset_ions_, sgenv_.deterministic_rs_,
                                               num_walkers, generate_test_data);
}

template<typename REAL, ValidSpaceGridInput::valid VALID>
RefVectorWithLeader<ParticleSet> SpaceGridTest<REAL, VALID>::getPSetList()
{
  return {psets_[0], makeRefVector<ParticleSet>(psets_)};
}


template class SpaceGridEnv<ValidSpaceGridInput::valid::DEFAULT>;
template class SpaceGridEnv<ValidSpaceGridInput::valid::ORIGIN>;
template class SpaceGridEnv<ValidSpaceGridInput::valid::CYLINDRICAL>;
template class SpaceGridEnv<ValidSpaceGridInput::valid::SPHERICAL>;
template class SpaceGridEnv<ValidSpaceGridInput::valid::WEIRD_CARTESIAN>;

template class SpaceGridTest<float, ValidSpaceGridInput::valid::DEFAULT>;
template class SpaceGridTest<double, ValidSpaceGridInput::valid::DEFAULT>;
template class SpaceGridTest<float, ValidSpaceGridInput::valid::CYLINDRICAL>;
template class SpaceGridTest<double, ValidSpaceGridInput::valid::CYLINDRICAL>;
template class SpaceGridTest<float, ValidSpaceGridInput::valid::ORIGIN>;
template class SpaceGridTest<double, ValidSpaceGridInput::valid::ORIGIN>;

} // namespace testing
} // namespace qmcplusplus
