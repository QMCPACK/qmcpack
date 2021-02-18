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

#include "Message/Communicate.h"
#include "OhmmsData/Libxml2Doc.h"
#include "QMCWaveFunctions/WaveFunctionPool.h"

namespace qmcplusplus
{
class MinimalWaveFunctionPool
{
  const char* wf_input = R"(
<wavefunction target='e'>
  <determinantset type='einspline' href='pwscf.pwscf.h5' tilematrix='1 0 0 0 1 0 0 0 1' twistnum='0' source='ion' meshfactor='1.0' precision='float'>
    <slaterdeterminant>
      <determinant id='updet' size='4'>
        <occupation mode='ground' spindataset='0'/>
      </determinant>
      <determinant id='downdet' size='4'>
        <occupation mode='ground' spindataset='0'/>
      </determinant>
    </slaterdeterminant>
  </determinantset>
</wavefunction>
  )";

public:
  MinimalWaveFunctionPool() : comm_(nullptr) {}

  WaveFunctionPool operator()(Communicate* comm, ParticleSetPool& particle_pool)
  {
    comm_ = comm;
    WaveFunctionPool wp(particle_pool, comm_);

    Libxml2Document* doc = new Libxml2Document;
    bool okay            = doc->parseFromString(wf_input);
    REQUIRE(okay);

    xmlNodePtr root = doc->getRoot();

    wp.put(root);

    delete doc;
    return wp;
  }

private:
  Communicate* comm_;
};

} // namespace qmcplusplus
#endif
