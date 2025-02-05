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

  static constexpr const char* const wf_input_spinor = R"(
<wavefunction name="psi0" target="e">
   <sposet_collection name="A" type="einspline" href="o2_45deg_spins.pwscf.h5" tilematrix="1 0 0 0 1 0 0 0 1" twistnum="0" source="ion0" size="12"> 
      <sposet name="spo_ud" size="12"/> 
      <sposet name="spo_dm" size="8" index_min="12" index_max="20"/> 
   </sposet_collection> 
   <determinantset>
      <slaterdeterminant>
         <determinant sposet="spo_ud"/>
      </slaterdeterminant>
   </determinantset>
</wavefunction>
  )";

  static constexpr const char* const wf_input_spinor_J12 = R"(
<wavefunction name="psi0" target="e">
   <sposet_collection name="A" type="einspline" href="o2_45deg_spins.pwscf.h5" tilematrix="1 0 0 0 1 0 0 0 1" twistnum="0" source="ion0" size="12"> 
      <sposet name="spo_ud" size="12"/> 
      <sposet name="spo_dm" size="8" index_min="12" index_max="20"/> 
   </sposet_collection> 
   <determinantset>
      <slaterdeterminant>
         <determinant sposet="spo_ud"/>
      </slaterdeterminant>
   </determinantset>
   <jastrow type="One-Body" name="J1" function="bspline" source="ion0" print="yes">
      <correlation elementType="O" size="9" rcut="2.336894584512495" cusp="0.0">
         <coefficients id="eO" type="Array">                  
-0.51632 -0.1591167977 -0.172367432 -0.1238310413 -0.09792672786 
-0.91785 -0.05476753103 -0.03482448615 -0.01864350288
         </coefficients>
      </correlation>
   </jastrow>
   <jastrow type="Two-Body" name="J2" function="bspline" print="yes">
      <correlation speciesA="u" speciesB="u" size="9" rcut="2.336894584512495" cusp="-0.5">
         <coefficients id="uu" type="Array">                  
0.7554 0.5342428628 0.3861610501 0.2724177345 0.186010153 0.1213795099 
0.04796 0.04068638111 0.01968948012
         </coefficients>
      </correlation>
   </jastrow>
</wavefunction>
  )";

public:
  static WaveFunctionPool make_diamondC_1x1x1(const RuntimeOptions& runtime_options,
                                              Communicate* comm,
                                              ParticleSetPool& particle_pool)
  {
    WaveFunctionPool wp(runtime_options, particle_pool, comm);

    Libxml2Document doc;
    bool okay = doc.parseFromString(wf_input);
    REQUIRE(okay);

    xmlNodePtr root = doc.getRoot();

    wp.put(root);

    return wp;
  }

  static WaveFunctionPool make_O2_spinor(const RuntimeOptions& runtime_options,
                                         Communicate* comm,
                                         ParticleSetPool& particle_pool)
  {
    WaveFunctionPool wp(runtime_options, particle_pool, comm);

    Libxml2Document doc;
    bool okay = doc.parseFromString(wf_input_spinor);
    REQUIRE(okay);

    xmlNodePtr root = doc.getRoot();

    wp.put(root);

    return wp;
  }

  static WaveFunctionPool make_O2_spinor_J12(const RuntimeOptions& runtime_options,
                                             Communicate* comm,
                                             ParticleSetPool& particle_pool)
  {
    WaveFunctionPool wp(runtime_options, particle_pool, comm);

    Libxml2Document doc;
    bool okay = doc.parseFromString(wf_input_spinor_J12);
    REQUIRE(okay);

    xmlNodePtr root = doc.getRoot();

    wp.put(root);

    return wp;
  }
};

} // namespace qmcplusplus
#endif
