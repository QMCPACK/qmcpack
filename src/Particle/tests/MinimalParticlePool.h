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

#ifndef QMCPLUSPLUS_MINIMALPARTICLEPOOL_H
#define QMCPLUSPLUS_MINIMALPARTICLEPOOL_H

#include "Message/Communicate.h"
#include "OhmmsData/Libxml2Doc.h"
#include "Particle/ParticleSetPool.h"

namespace qmcplusplus
{
/** This should be the minimal ParticleSetPool for integration tests
 *
 */
class MinimalParticlePool
{
  // See ParticleIO/tests/test_xml_io.cpp for particle parsing
  const char* particles_xml = R"(
<tmp>
 <simulationcell>
      <parameter name='lattice' units='bohr'>
          3.37316115        3.37316115        0.00000000
          0.00000000        3.37316115        3.37316115
          3.37316115        0.00000000        3.37316115
      </parameter>
      <parameter name='bconds'>
         p p p
      </parameter>
        <parameter name='LR_dim_cutoff'>15 </parameter>
</simulationcell>
<particleset name="ion" size="2">
  <group name="C">
    <parameter name="charge">4</parameter>
  </group>
  <attrib name="position" datatype="posArray">
    0.00000000  0.00000000  0.00000000
    1.68658058  1.68658058  1.68658058
  </attrib>
</particleset>
 <particleset name="e" random="yes" >
  <group name="u" size="4">
    <parameter name="charge">-1</parameter>
  </group>
  <group name="d" size="4">
    <parameter name="charge">-1</parameter>
  </group>
</particleset>
</tmp>
)";

public:
  ParticleSetPool operator()(Communicate* c)
  {
    Libxml2Document doc;

    doc.parseFromString(particles_xml);

    xmlNodePtr root     = doc.getRoot();
    xmlNodePtr sim_cell = xmlFirstElementChild(root);

    ParticleSetPool pp(c);

    // Need to set up simulation cell lattice before reading particle sets
    pp.readSimulationCellXML(sim_cell);

    xmlNodePtr part_ion = xmlNextElementSibling(sim_cell);
    pp.put(part_ion);
    xmlNodePtr part_elec = xmlNextElementSibling(part_ion);
    pp.put(part_elec);
    pp.randomize();

    return pp;
  }
};

} // namespace qmcplusplus

#endif
