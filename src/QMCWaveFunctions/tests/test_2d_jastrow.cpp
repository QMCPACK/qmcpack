//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2022 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Yubo "Paul" Yang, yubo.paul.yang@gmail.com, CCQ @ Flatiron
//
// File created by: Yubo "Paul" Yang, yubo.paul.yang@gmail.com, CCQ @ Flatiron
//////////////////////////////////////////////////////////////////////////////////////


#include "catch.hpp"

#include "OhmmsData/Libxml2Doc.h"
#include "ParticleIO/XMLParticleIO.h" // XMLParticleParser
#include "QMCWaveFunctions/Jastrow/RadialJastrowBuilder.h"
#include "QMCWaveFunctions/Jastrow/TwoBodyJastrow.h"

using std::real;
using std::imag;

namespace qmcplusplus
{
using RealType     = WaveFunctionComponent::RealType;
using PsiValueType = WaveFunctionComponent::PsiValueType;
using ValueType    = QMCTraits::ValueType; // for real and mixed precision
using PosType      = QMCTraits::PosType;

TEST_CASE("Jastrow 2D", "[wavefunction]")
{
  Communicate* c = OHMMS::Controller;
  Libxml2Document doc;
  xmlNodePtr root, node;
  bool okay;

  // Step 1: create Jastrow
  //   TwoBodyJastrow<BsplineFunctor<RealType>> j2;
  //   RadialJastrowBuilder jastrow(Communicator, ParticleSet);
  //   ParticleSet( SimulationCell( CrystalLattice ) )
  CrystalLattice<OHMMS_PRECISION, OHMMS_DIM> lattice;
  lattice.BoxBConds = true;
  lattice.R.diagonal(70.89815403622065); 
  lattice.ndim = 2;
  lattice.reset();
  const SimulationCell cell(lattice);
  ParticleSet elec(cell);
  //   read ParticleSet from xml text
  const char* particle_text = R"(<tmp>
    <particleset name="e" random="no">
      <group name="u" size="2">
        <parameter name="charge"> -1 </parameter>
        <attrib name="position" datatype="posArray" condition="0">
54.66209032978565 2.018663420381362 0.0
23.566489035927912 3.443259712257945 0.0
</attrib>
      </group>
      <group name="d" size="2">
        <parameter name="charge"> -1 </parameter>
        <attrib name="position" datatype="posArray" condition="0">
67.26206197993817 30.29582561496142 0.0
37.00657847142635 1.4508035033146867 0.0
</attrib>
      </group>
    </particleset>
  </tmp>)";
  okay = doc.parseFromString(particle_text);
  REQUIRE(okay);
  root = doc.getRoot();
  node = xmlFirstElementChild(root);
  XMLParticleParser parse_electrons(elec);
  parse_electrons.readXML(node);
  int itab; elec.addTable(elec);
  elec.update(); // update distance tables
  const int nelec = elec.getTotalNum();
  //    read Jastrow component from xml text
  const char* jastrow_text = R"(<tmp>
   <jastrow name="J2" type="Two-Body" function="Bspline">
        <correlation speciesA="u" speciesB="u" size="8">
          <coefficients id="uu" type="Array" optimize="yes">4.868951397 3.154235815 1.719776072 0.9676536301 0.6044866223 0.3368526364 0.1566214572 0.06031539785</coefficients>
        </correlation>
        <correlation speciesA="u" speciesB="d" size="8">
          <coefficients id="ud" type="Array" optimize="yes">6.991319036 3.93760887 2.077967513 1.115208829 0.6946729632 0.3826149045 0.1705411558 0.06155742938</coefficients>
        </correlation>
      </jastrow>
</tmp>)";
  okay = doc.parseFromString(jastrow_text);
  REQUIRE(okay);
  root = doc.getRoot();
  node = xmlFirstElementChild(root);
  RadialJastrowBuilder jastrow(c, elec);
  using J2Type = TwoBodyJastrow<BsplineFunctor<RealType>>;
  auto j2_uptr = jastrow.buildComponent(node);
  J2Type* j2   = dynamic_cast<J2Type*>(j2_uptr.get());
  REQUIRE(j2);

  // Step 2: test Jastrow vglh
  double logpsi = real(j2->evaluateLog(elec, elec.G, elec.L));
  CHECK(real(logpsi) == Approx(-1.51908));
  CHECK(real(elec.G[0][0]) == Approx(0.08490220791));
  CHECK(real(elec.G[0][1]) == Approx(-0.006958062051));
  CHECK(real(elec.G[1][0]) == Approx(-0.1325957095));
  CHECK(real(elec.G[1][1]) == Approx(0.01872445072));
  CHECK(real(elec.G[2][0]) == Approx(0.004059085343));
  CHECK(real(elec.G[2][1]) == Approx(0.009109497845));
  CHECK(real(elec.G[3][0]) == Approx(0.04363441629));
  CHECK(real(elec.G[3][1]) == Approx(-0.02087588652));
  for (int i=0;i<nelec;i++)
  {
    for (int l=0;l<OHMMS_DIM;l++)
      CHECK(imag(elec.G[i][l]) == Approx(0));
    CHECK(real(elec.G[i][2]) == Approx(0)); // ndim=2
  }
  const std::vector<RealType> lap_values = {-0.00916449, -0.0166369, -0.00351783, -0.0153977};
  for (int m=0;m<nelec;m++)
    CHECK(real(elec.L[m]) == Approx(lap_values[m]));

  WaveFunctionComponent::HessVector grad_grad_psi;
  grad_grad_psi.resize(nelec);
  grad_grad_psi = 0.0;

  CHECK_THROWS(j2->evaluateHessian(elec, grad_grad_psi));
  //std::vector<double> hess_values = {
  //  -0.0108098, -0.00172466, 0,
  //  -0.00172466, 0.00164531, 0,
  //  0, 0, 0.00513802,
  //  -0.0254307, 0.00476379, 0,
  //  0.00476379, 0.00879376, 0,
  //  0, 0, 0.00948111,
  //  -0.000367339, -0.00154737, 0,
  //  -0.00154737, -0.00315049, 0,
  //  0, 0, 0.00032215,
  //  -0.0284186, 0.00421815, 0,
  //  0.00421815, 0.0130209, 0,
  //  0, 0, 0.0137115
  //}; // why do the zz components have non-zero values?

  //int m = 0;
  //for (int n = 0; n < nelec; n++)
  //  for (int i = 0; i < OHMMS_DIM; i++)
  //    for (int j = 0; j < OHMMS_DIM; j++, m++)
  //      CHECK(real(grad_grad_psi[n](i, j)) == Approx(hess_values[m]));

  // Step 3: test Jastrow ratio for pbyp move
  PosType newpos(0.3, 0.2, 0.5);

  elec.makeVirtualMoves(newpos);
  std::vector<ValueType> ratios(nelec);
  j2->evaluateRatiosAlltoOne(elec, ratios);
  std::vector<double> ratio_values = {1.46023, 1.46559, 0.444258, 1.90226};
  for (int i=0;i<ratios.size();i++)
    CHECK(real(ratios[i]) == Approx(ratio_values[i]));

  for(int i=0;i<nelec;i++)
  {
    elec.makeMove(i, newpos-elec.R[i]);
    PsiValueType rat1 = j2->ratio(elec, i);
    CHECK(real(rat1) == Approx(ratio_values[i]));
    elec.rejectMove(i);
  }

  // Step 4: test Jastrow parameter derivatives
  UniqueOptObjRefs opt_obj_refs;
  j2->extractOptimizableObjectRefs(opt_obj_refs);
  REQUIRE(opt_obj_refs.size() == 2);

  opt_variables_type optvars;
  Vector<ValueType> dlogpsi;
  Vector<ValueType> dhpsioverpsi;

  for (OptimizableObject& obj : opt_obj_refs)
    obj.checkInVariablesExclusive(optvars);
  optvars.resetIndex();
  const int NumOptimizables(optvars.size());
  j2->checkOutVariables(optvars);
  dlogpsi.resize(NumOptimizables);
  dhpsioverpsi.resize(NumOptimizables);
  j2->evaluateDerivatives(elec, optvars, dlogpsi, dhpsioverpsi);
  app_log() << std::endl << "reporting dlogpsi and dhpsioverpsi" << std::scientific << std::endl;
  for (int iparam = 0; iparam < NumOptimizables; iparam++)
    app_log() << "param=" << iparam << " : " << dlogpsi[iparam] << "  " << dhpsioverpsi[iparam] << std::endl;
  app_log() << std::endl;
  const std::vector<RealType> dlogpsi_values = {0, 0, 0, 0, 0, 0, -1.521258e-04, -2.194165e-01, 0, 0, -2.779981e-02, -5.327999e-01, -9.356640e-01, -4.847466e-01, -1.945081e-02, -2.453297e-01};
  const std::vector<RealType> dhpsi_values = {0, 0, 0, 0, 0, 0, 5.953288e-03, 8.618836e-03, 0, 0, 2.572195e-02, -5.048126e-02, -3.139861e-02, 2.197638e-02, 4.319522e-02, 3.512834e-02};
  for (int iparam = 0; iparam < NumOptimizables; iparam++)
  {
    CHECK(real(dlogpsi[iparam]) == Approx(dlogpsi_values[iparam]));
    CHECK(real(dhpsioverpsi[iparam]) == Approx(dhpsi_values[iparam]));
  }
}
}
