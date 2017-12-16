//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2017 Jeongnim Kim and QMCPACK developers.
//
// File developed by:  Mark Dewing, mdewing@anl.gov Argonne National Laboratory
//
// File created by: Mark Dewing, mdewing@anl.gov Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////


#include "catch.hpp"


#include "Message/Communicate.h"
#include "OhmmsData/Libxml2Doc.h"
#include "QMCApp/ParticleSetPool.h"
#include "QMCApp/WaveFunctionPool.h"


#include <stdio.h>
#include <string>
#include <sstream>



namespace qmcplusplus
{

void setupParticleSetPool(ParticleSetPool &pp)
{

// See ParticleIO/tests/test_xml_io.cpp for particle parsing
const char *particles = \
" \
<tmp> \
 <simulationcell> \
      <parameter name='lattice' units='bohr'> \
          3.37316115        3.37316115        0.00000000 \
          0.00000000        3.37316115        3.37316115 \
          3.37316115        0.00000000        3.37316115 \
      </parameter> \
      <parameter name='bconds'> \
         p p p \
      </parameter> \
        <parameter name='LR_dim_cutoff'>15 </parameter> \
</simulationcell> \
<particleset name=\"ion\" size=\"2\"> \
  <group name=\"C\"> \
    <parameter name=\"charge\">4</parameter> \
  </group> \
  <attrib name=\"position\" datatype=\"posArray\"> \
    0.00000000  0.00000000  0.00000000 \
    1.68658058  1.68658058  1.68658058 \
  </attrib> \
</particleset> \
 <particleset name=\"e\" random=\"yes\" > \
  <group name=\"u\" size=\"4\"> \
    <parameter name=\"charge\">-1</parameter> \
  </group> \
  <group name=\"d\" size=\"4\"> \
    <parameter name=\"charge\">-1</parameter> \
  </group> \
</particleset> \
</tmp> \
";
  Libxml2Document doc;

  bool okay = doc.parseFromString(particles);
  REQUIRE(okay);

  xmlNodePtr root = doc.getRoot();
  xmlNodePtr sim_cell= xmlFirstElementChild(root);
  // Need to set up simulation cell lattice before reading particle sets
  pp.putLattice(sim_cell);

  xmlNodePtr part_ion = xmlNextElementSibling(sim_cell);
  pp.put(part_ion);
  xmlNodePtr part_elec = xmlNextElementSibling(part_ion);
  pp.put(part_elec);
  pp.randomize();
}


TEST_CASE("WaveFunctionPool", "[qmcapp]")
{

  Communicate *c;
  OHMMS::Controller->initialize(0, NULL);
  c = OHMMS::Controller;

  WaveFunctionPool wp(c);

  REQUIRE(wp.empty() == true);


  ParticleSetPool pp(c);
  setupParticleSetPool(pp);
  wp.setParticleSetPool(&pp);

  const char *wf_input = \
  "<wavefunction target='e'>\
     <determinantset type='einspline' href='pwscf.pwscf.h5' tilematrix='1 0 0 0 1 0 0 0 1' twistnum='0' source='ion' meshfactor='1.0' precision='float'> \
         <slaterdeterminant> \
            <determinant id='updet' size='4'> \
              <occupation mode='ground' spindataset='0'/>\
             </determinant>\
              <determinant id='downdet' size='4'>\
                <occupation mode='ground' spindataset='0'/>\
             </determinant>\
         </slaterdeterminant>\
     </determinantset> \
   </wavefunction> \
  ";

  Libxml2Document *doc = new Libxml2Document;
  bool okay = doc->parseFromString(wf_input);
  REQUIRE(okay);

  xmlNodePtr root = doc->getRoot();

  wp.put(root);

  TrialWaveFunction *psi = wp.getWaveFunction("psi0");
  REQUIRE(psi != NULL);
}
}
