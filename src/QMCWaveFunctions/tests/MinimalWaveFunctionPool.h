//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2019 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//////////////////////////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_MINIMALWAVEFUNCTIONPOOL_H
#define QMCPLUSPLUS_MINIMALWAVEFUNCTIONPOOL_H

#include "catch.hpp"
#include "Message/Communicate.h"
#include "OhmmsData/Libxml2Doc.h"
#include "QMCWaveFunctions/WaveFunctionPool.h"

namespace qmcplusplus
{
class MinimalWaveFunctionPool
{
  static constexpr const char* const wf_input = R"(
<wavefunction target='e'>
  <sposet_collection type="bspline" source="ion" href="diamondC_1x1x1.pwscf.h5" tilematrix="1 0 0 0 1 0 0 0 1" twistnum="0" meshfactor="0.8" twist="0 0 0" precision="double">
    <sposet name="spo_for_dets" size="4" spindataset="0"/>
  </sposet_collection>
  <sposet_collection type="bspline" source="ion" href="diamondC_1x1x1.pwscf.h5" tilematrix="1 0 0 0 1 0 0 0 1" twistnum="0" gpu="no" meshfactor="0.8" twist="0 0 0" precision="double">
    <sposet name="spo_ud" size="4" spindataset="0"/>
    <sposet name="spo_dm" index_min="4" index_max="8" spindataset="0"/>
  </sposet_collection>
  <determinantset>
    <slaterdeterminant>
      <determinant sposet='spo_for_dets'/>
      <determinant sposet='spo_for_dets'/>
    </slaterdeterminant>
  </determinantset>
</wavefunction>
  )";

public:
  static WaveFunctionPool make_diamondC_1x1x1(Communicate* comm, ParticleSetPool& particle_pool)
  {
    WaveFunctionPool wp(particle_pool, comm);

    Libxml2Document doc;
    bool okay = doc.parseFromString(wf_input);
    REQUIRE(okay);

    xmlNodePtr root = doc.getRoot();

    wp.put(root);

    return wp;
  }
};

} // namespace qmcplusplus
#endif
