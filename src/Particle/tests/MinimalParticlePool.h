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
  static constexpr const char* const particles_xml = R"(
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

  static constexpr const char* const particles_xml_spinor = R"(
<tmp>
   <simulationcell>
      <parameter name="lattice" units="bohr">
               5.10509515       -3.23993545        0.00000000
               5.10509515        3.23993545       -0.00000000
              -6.49690625        0.00000000        7.08268015
      </parameter>
      <parameter name="bconds">
         p p p
      </parameter>
      <parameter name="LR_dim_cutoff"       >    15                 </parameter>
   </simulationcell>
   <particleset name="ion0">
      <group name="O" size="2" mass="29164.3928678">
         <parameter name="charge"              >    6                     </parameter>
         <parameter name="valence"             >    6                     </parameter>
         <parameter name="atomicnumber"        >    8                     </parameter>
         <parameter name="mass"                >    29164.3928678            </parameter>
         <attrib name="position" datatype="posArray" condition="0">
                 -0.00000000       -0.00000000        1.08659253
                  0.00000000        0.00000000       -1.08659253
         </attrib>
      </group>
   </particleset>
   <particleset name="e" spinor="yes" random="yes">
      <group name="u" size="12" mass="1.0">
         <parameter name="charge"              >    -1                    </parameter>
         <parameter name="mass"                >    1.0                   </parameter>
      </group>
   </particleset>
</tmp>
)";

  static constexpr const char* const particles_xml_NiO_a4 = R"(
<tmp>
  <simulationcell>
    <parameter name="lattice">
      3.94055   3.94055   7.8811
      7.8811   3.94055   3.94055
      3.94055   7.8811   3.94055
    </parameter>
    <parameter name="bconds">p p p </parameter>
    <parameter name="LR_dim_cutoff">15</parameter>
  </simulationcell>
  <particleset name="i" size="4">
    <group name="O">
      <parameter name="charge">6.000000</parameter>
      <parameter name="valence">6.000000</parameter>
      <parameter name="atomicnumber">8.000000</parameter>
    </group>
    <group name="Ni">
      <parameter name="charge">18</parameter>
      <parameter name="valence">18</parameter>
      <parameter name="atomicnumber">28</parameter>
    </group>
    <attrib name="position" datatype="posArray" condition="1">
    0.25   0.25   0.25
    0.75   0.75   0.75
    0   0   0
    0.5   0.5   0.5
    </attrib>
    <attrib name="ionid" datatype="stringArray">
      O  O  Ni  Ni
    </attrib>
  </particleset>
  <particleset name="e" random="yes" randomsrc="i">
    <group name="u" size="24">
      <parameter name="charge">-1</parameter>
    </group>
    <group name="d" size="24">
      <parameter name="charge">-1</parameter>
    </group>
  </particleset>
</tmp>
)";

  static constexpr const char* const particles_xml_H2 = R"(
<tmp>
  <particleset name="ion" size="2">
    <group name="H">
      <parameter name="charge">1</parameter>
      <parameter name="valence">1</parameter>
      <parameter name="atomicnumber">1</parameter>
    </group>
    <attrib name="position" datatype="posArray">
  0.0000000000e+00  0.0000000000e+00 -9.4486299390e-01
  0.0000000000e+00  0.0000000000e+00  9.4486299390e-01
</attrib>
    <attrib name="ionid" datatype="stringArray">
 H H 
</attrib>
  </particleset>
  <particleset name="e">
    <group name="u" size="1">
      <parameter name="charge">-1</parameter>
      <attrib name="position" datatype="posArray">
  5.2261545687e-01  4.6234802699e-01 -1.2355134293e+00
</attrib>
    </group>
    <group name="d" size="1">
      <parameter name="charge">-1</parameter>
      <attrib name="position" datatype="posArray">
  6.4451072699e-01  1.3453681375e-01  1.3161796485e+00
</attrib>
    </group>
  </particleset>
</tmp>
)";

public:
  static void parseParticleSetXML(const char* xml_string, ParticleSetPool& pp)
  {
    Libxml2Document doc;

    doc.parseFromString(xml_string);

    xmlNodePtr root     = doc.getRoot();
    xmlNodePtr sim_cell = xmlFirstElementChild(root);

    // Need to set up simulation cell lattice before reading particle sets
    pp.readSimulationCellXML(sim_cell);

    xmlNodePtr part_ion = xmlNextElementSibling(sim_cell);
    pp.put(part_ion);
    xmlNodePtr part_elec = xmlNextElementSibling(part_ion);
    pp.put(part_elec);
    pp.randomize();
  }

  static ParticleSetPool make_diamondC_1x1x1(Communicate* c)
  {
    ParticleSetPool pp(c);
    parseParticleSetXML(particles_xml, pp);
    return pp;
  }

  static ParticleSetPool make_O2_spinor(Communicate* c)
  {
    ParticleSetPool pp(c);
    parseParticleSetXML(particles_xml_spinor, pp);
    return pp;
  }

  static ParticleSetPool make_NiO_a4(Communicate* c)
  {
    ParticleSetPool pp(c);
    parseParticleSetXML(particles_xml_NiO_a4, pp);
    return pp;
  }

  static ParticleSetPool make_H2(Communicate* c)
  {
    Libxml2Document doc;

    doc.parseFromString(particles_xml_H2);

    xmlNodePtr root     = doc.getRoot();
    xmlNodePtr part_ion = xmlFirstElementChild(root);

    ParticleSetPool pp(c);
    // Don't set up simulation cell lattice before reading particle sets
    // convention is to leave out when doing open boundary cond.

    pp.put(part_ion);
    xmlNodePtr part_elec = xmlNextElementSibling(part_ion);
    pp.put(part_elec);
    pp.randomize();

    return pp;
  }
};

} // namespace qmcplusplus

#endif
