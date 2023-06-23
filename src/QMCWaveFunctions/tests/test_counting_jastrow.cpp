//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Brett Van Der Goetz, bvdg@berkeley.edu, University of California at Berkeley
//
// File created by: Brett Van Der Goetz, bvdg@berkeley.edu, University of California at Berkeley
//////////////////////////////////////////////////////////////////////////////////////

#include "catch.hpp"

#include "Configuration.h"
#include "OhmmsData/Libxml2Doc.h"
#include "Particle/ParticleSet.h"
#include "VariableSet.h"

#include "QMCWaveFunctions/Jastrow/CountingGaussian.h"
#include "QMCWaveFunctions/Jastrow/CountingGaussianRegion.h"
#include "QMCWaveFunctions/Jastrow/CountingJastrow.h"
#include "QMCWaveFunctions/Jastrow/CountingJastrowBuilder.h"

#include <stdio.h>

namespace qmcplusplus
{

// CountingGaussian unit tests
TEST_CASE("Gaussian Functor","[wavefunction]")
{
  using RealType = QMCTraits::RealType;
  using PosType = QMCTraits::PosType;
  using TensorType = QMCTraits::TensorType;

  CountingGaussian gf_abc("gf_abc");
  CountingGaussian gf_adk("gf_adk");

  // test parse/put for both input modes
  const char* gaussian_xml_abc = R"(<function id="g0">
          <var name="A" opt="true">-1.5 0.353553390593273 -0.353553390593273 -2.25 -0.75 -2.25</var>
          <var name="B" opt="true">-1.5 0.353553390593273 -0.353553390593273</var>
          <var name="C" opt="true">-1.5</var>
        </function>)";
  const char* gaussian_xml_adk = R"(<function id="g1">
          <var name="A" opt="true">-1.5 0.353553390593273 -0.353553390593273 -2.25 -0.75 -2.25</var>
          <var name="D" opt="true">1.0 0.0 0.0</var>
          <var name="K" opt="true">0.0</var>
        </function>)";
  Libxml2Document doc_abc, doc_adk;
  bool parse_abc = doc_abc.parseFromString(gaussian_xml_abc);
  bool parse_adk = doc_adk.parseFromString(gaussian_xml_adk);
  xmlNodePtr root_abc = doc_abc.getRoot();
  xmlNodePtr root_adk = doc_adk.getRoot();
  bool put_abc = gf_abc.put(root_abc);
  bool put_adk = gf_adk.put(root_adk);
  REQUIRE( (parse_abc && parse_adk && put_abc && put_adk) == true);

  // test points 
  PosType r1(1,0,0);
  PosType r2(1.707106781186547, 0.5, -0.5);
  PosType r3(0.4714617144631338, 0.1499413068889379, 0.2932213074999387);

  // return variables
  RealType fval, lval, llap, flap;
  PosType fgrad, lgrad;
  opt_variables_type opt_vars;

  std::vector<RealType> dfval;
  std::vector<PosType> dfgrad;
  std::vector<RealType> dflap;
  dfval.resize(10);
  dfgrad.resize(10);
  dflap.resize(10);

  // value tests for ADK input
  gf_adk.evaluate(r1, fval, fgrad, flap);
  gf_adk.evaluateLog(r1, lval, lgrad, llap);
  CHECK( fval     == Approx(1) );
  CHECK( fgrad[0] == Approx(0) );
  CHECK( fgrad[1] == Approx(0) );
  CHECK( fgrad[2] == Approx(0) );
  CHECK( flap     == Approx(-12) );
  CHECK( lval     == Approx(0) );
  CHECK( lgrad[0] == Approx(0) );
  CHECK( lgrad[1] == Approx(0) );
  CHECK( lgrad[2] == Approx(0) );
  CHECK( llap     == Approx(-12) );
  // value tests for ABC input
  gf_abc.evaluate(r1, fval, fgrad, flap);
  gf_abc.evaluateLog(r1, lval, lgrad, llap);
  CHECK( fval     == Approx(1) );
  CHECK( fgrad[0] == Approx(0) );
  CHECK( fgrad[1] == Approx(0) );
  CHECK( fgrad[2] == Approx(0) );
  CHECK( flap     == Approx(-12) );
  CHECK( lval     == Approx(0) );
  CHECK( lgrad[0] == Approx(0) );
  CHECK( lgrad[1] == Approx(0) );
  CHECK( lgrad[2] == Approx(0) );
  CHECK( llap     == Approx(-12) );
  // evaluateDerivatives
  gf_abc.evaluateDerivatives(r3, dfval, dfgrad, dflap);
  CHECK( dfval[0] == Approx(0.113120472934));
  CHECK( dfval[1] == Approx(0.071952529875));
  CHECK( dfval[2] == Approx(0.140708490047));
  CHECK( dfval[3] == Approx(0.011441709933));
  CHECK( dfval[4] == Approx(0.044750218821));
  CHECK( dfval[5] == Approx(0.043756180154));
  CHECK( dfval[6] == Approx(-0.47987130010));
  CHECK( dfval[7] == Approx(-0.15261584911));
  CHECK( dfval[8] == Approx(-0.29845157248));
  CHECK( dfval[9] == Approx(0.508918630484));
  // evaluateLogDerivatives
  gf_abc.evaluateLogDerivatives(r3, dfval, dfgrad, dflap);
  CHECK( dfval[0] == Approx(0.222276148205));
  CHECK( dfval[1] == Approx(0.141383171229));
  CHECK( dfval[2] == Approx(0.276485240702));
  CHECK( dfval[3] == Approx(0.022482395511));
  CHECK( dfval[4] == Approx(0.087931972108));
  CHECK( dfval[5] == Approx(0.085978735172));
  CHECK( dfval[6] == Approx(-0.94292342892));
  CHECK( dfval[7] == Approx(-0.29988261377));
  CHECK( dfval[8] == Approx(-0.58644261500));
  CHECK( dfval[9] == Approx(1));
}

TEST_CASE("CountingJastrow","[wavefunction]")
{
  using PosType = QMCTraits::PosType;
  using GradType = QMCTraits::GradType;
  using RealType = QMCTraits::RealType;
  using ValueType = QMCTraits::ValueType;
  using VariableSet = optimize::VariableSet;
  using LogValueType = std::complex<QMCTraits::QTFull::RealType>;

  Communicate* c = OHMMS::Controller;

  // initialize particle sets
  const SimulationCell simulation_cell;
  ParticleSet elec(simulation_cell);
  std::vector<int> egroup(1);
  int num_els = 4;
  egroup[0] = num_els;
  elec.setName("e");
  elec.create(egroup);
  PosType Re[] = { PosType(2.4601162537, 6.7476360528,-1.9073129953),
                   PosType(2.2585811248, 2.1282254384,0.051545776028),
                   PosType(0.84796873937,5.1735597110,0.84642416761),
                   PosType(3.1597337850, 5.1079432473,1.0545953717)};
  for(int i = 0; i < num_els; ++i)
    for(int k = 0; k < 3; ++k)
      elec.R[i][k] = Re[i][k];

  ParticleSet ion0(simulation_cell);
  std::vector<int> igroup(1);
  int num_ion = 4;
  igroup[0] = num_ion;
  ion0.setName("ion0");
  ion0.create(igroup);
  PosType Ri[] = {PosType(2.61363352510301,  5.01928226905281,  0.0),
                  PosType(3.74709851167814,  3.70007145722224,  0.0),
                  PosType(6.11011670934565,  1.66504047681825,  0.0),
                  PosType(8.10584803421091,  7.78608266172472,  0.0)};
  
  for(int i = 0; i < num_ion; ++i)
    for(int k = 0; k < 3; ++k)
      ion0.R[i][k] = Ri[i][k];

  const char* cj_normgauss_xml = R"(<jastrow name="ncjf_normgauss" type="Counting">
      <var name="F" opt="true">
        4.4903e-01 5.3502e-01 5.2550e-01 6.8081e-01
                   5.1408e-01 4.8658e-01 6.2182e-01
                              2.7189e-01 9.4951e-01
                                         0.0000e+00
      </var>
      <region type="normalized_gaussian" reference_id="g0" opt="true" >
        <function id="g0">
          <var name="A" opt="False">-1.0 -0.0 -0.0 -1.0 -0.0 -1.0</var>
          <var name="B" opt="False">-2.6136335251 -5.01928226905 0.0</var>
          <var name="C" opt="False">-32.0242747</var>
        </function>
        <function id="g1">
          <var name="A" opt="true">-1.0 -0.0 -0.0 -1.0 -0.0 -1.0</var>
          <var name="B" opt="true">-3.74709851168 -3.70007145722 0.0</var>
          <var name="C" opt="true">-27.7312760448</var>
        </function>
        <function id="g2">
          <var name="A" opt="true">-1.0 -0.0 -0.0 -1.0 -0.0 -1.0</var>
          <var name="B" opt="true">-6.11011670935 -1.66504047682 0.0</var>
          <var name="C" opt="true">-40.1058859913</var>
        </function>
        <function id="g3">
          <var name="A" opt="true">-1.0 -0.0 -0.0 -1.0 -0.0 -1.0</var>
          <var name="B" opt="true">-8.10584803421 -7.78608266172 0.0</var>
          <var name="C" opt="true">-126.327855569</var>
        </function>
      </region>
    </jastrow>)";
  // test put for normalized_gaussian
  Libxml2Document doc;
  bool parse_cj = doc.parseFromString(cj_normgauss_xml);
  REQUIRE( parse_cj );

  xmlNodePtr cj_root = doc.getRoot();
  CountingJastrowBuilder cjb(c, elec);

  auto cj_uptr                                = cjb.buildComponent(cj_root);
  CountingJastrow<CountingGaussianRegion>* cj = dynamic_cast<CountingJastrow<CountingGaussianRegion>*>(cj_uptr.get());

  // reference for evaluateLog, evalGrad
  RealType Jval_exact = 7.8100074447e+00;
  PosType Jgrad_exact[] = {PosType(3.6845037054e-04, -4.2882992861e-04, 0),
                            PosType(2.2032083234e-02, -2.5647245917e-02, 0),
                            PosType(6.0625112202e-04, -7.0560012380e-04, 0),
                            PosType(1.0622511249e-01, -1.2363268199e-01, 0)};
  RealType Jlap_exact[] = {1.9649428566e-03, -1.1385706794e-01, 3.2312809658e-03, 4.0668060285e-01};


  // test evaluateLog for cj
  LogValueType logval = cj->evaluateLog(elec, elec.G, elec.L);
  for(int i = 0; i < num_els; ++i)
  {
    for(int k = 0; k < 3; ++k)
      CHECK( Jgrad_exact[i][k] == Approx( std::real(elec.G[i][k])) );
    CHECK( Jlap_exact[i] == Approx( std::real(elec.L[i])) );
  }
  CHECK( ComplexApprox(Jval_exact) == logval );
  
  // test automatic/minimal voronoi generator
  const char* cj_voronoi_xml = R"(<jastrow name="ncjf_voronoi" type="Counting" source="ion0">
      <var name="F" opt="true">
        4.4903e-01 5.3502e-01 5.2550e-01 6.8081e-01
                   5.1408e-01 4.8658e-01 6.2182e-01
                              2.7189e-01 9.4951e-01
                                         0.0000e+00
      </var>
      <region type="voronoi" opt="true">
        <var name="alpha">1.0</var>
      </region>
    </jastrow>)";
  // test put
  Libxml2Document doc2;
  bool parse_cjv = doc2.parseFromString(cj_voronoi_xml);
  REQUIRE( parse_cjv );

  xmlNodePtr cjv_root = doc2.getRoot();
  CountingJastrowBuilder cjvb(c, elec, ion0);

  // test evaluateLog for cjv
  auto cjv_uptr                                = cjvb.buildComponent(cjv_root);
  CountingJastrow<CountingGaussianRegion>* cjv = dynamic_cast<CountingJastrow<CountingGaussianRegion>*>(cjv_uptr.get());

  for(int i = 0; i < num_els; ++i)
  {
    for(int k = 0; k < 3; ++k)
      elec.G[i][k] = 0;
    elec.L[i] = 0;
  }

  logval = cjv->evaluateLog(elec, elec.G, elec.L);
  for(int i = 0; i < num_els; ++i)
  {
    for(int k = 0; k < 3; ++k)
      CHECK( Jgrad_exact[i][k] == Approx( std::real(elec.G[i][k])) );
    CHECK( Jlap_exact[i] == Approx( std::real(elec.L[i])) );
  }
  CHECK( ComplexApprox(Jval_exact) == logval );
  
  // test evalGrad
  for(int iat = 0; iat < num_els; ++iat)
  {
    GradType Jgrad_iat = cj->evalGrad(elec, iat);
    for(int k = 0; k < 3; ++k)
      CHECK( Jgrad_exact[iat][k] == Approx( std::real(Jgrad_iat[k]) ));
  }

  // reference for ratio, ratioGrad, acceptMove
  PosType dr[] = {PosType(0.0984629815,  0.0144420719, 0.1334309321),
                  PosType(-0.1026409581, 0.2289767772, 0.490138058592),
                  PosType(0.03517477469, 0.2693941041,  0.16922043039),
                  PosType(0.3851116893,  -0.1387760973,  0.1980082182)};

  RealType ratioval_exact[] = { 1.00003304765, 0.987624289443, 0.999873210738, 1.09014860194};

  PosType Jgrad_t_exact[] = {PosType(4.4329270315e-04, -5.1593699287e-04, 0),
                              PosType(4.8722465115e-02, -5.6707785196e-02, 0),
                              PosType(3.2691265772e-04, -3.8048525335e-04, 0),
                              PosType(2.0373800011e-01, -2.3712542045e-01, 0)};

  // test ratio, ratioGrad, acceptMove
  for(int iat = 0; iat < num_els; ++iat)
  {
    elec.makeMoveAndCheck(iat,dr[iat]);
  
    RealType ratioval = std::real( cj->ratio(elec, iat) );
    CHECK( ratioval_exact[iat] == Approx(std::real(ratioval)) );

    GradType grad_iat(0, 0, 0);
    RealType gradratioval = std::real( cj->ratioGrad(elec, iat, grad_iat) );
    CHECK( ratioval_exact[iat] == Approx(gradratioval));
    for(int k = 0; k < 3; ++k)
      CHECK(Jgrad_t_exact[iat][k] == Approx(std::real(grad_iat[k])));

    cj->acceptMove(elec, iat);
  }

#ifndef QMC_COMPLEX
  // setup and reference for evaluateDerivatives
  PosType R2[] = { PosType( 4.3280064837, 2.4657709845,  6.3466520181e-01),
                   PosType( 9.7075155012, 7.2984775093, -8.1975111678e-01),
                   PosType( 5.7514912378, 5.2001615327,  6.6673589235e-01),
                   PosType( 8.3805220665, 7.9424368608, -3.5106422506e-02)};
  PosType G2[] = {PosType( 3.3480105792e-01,  2.1316369526e-01, -4.1812914940e-01),
                   PosType(-7.0561066397e-01,  1.7863210008e-01,  3.5677994771e-01),
                   PosType(-9.2302398033e-01, -5.0740272660e-01, -2.6078603626e-01),
                   PosType(-8.3764545416e-01, -4.9181684009e-01,  1.0675382607e-01)};
  for(int i = 0; i < num_els; ++i)
    for(int k = 0; k < 3; ++k)
    {
      elec.R[i][k] = R2[i][k];
      elec.G[i][k] = G2[i][k];
    }

  int num_derivs = 39;
  RealType dlogpsi_exact[] = { 
     7.0982172306e-04, 9.8329357367e-02, 6.6879065207e-03, 1.0670293004e-01, 3.4053136887e+00,
     4.6322726464e-01, 7.3906096412e+00, 1.5753284303e-02, 5.0267496641e-01,-1.3874695168e+00,
    -2.6249136239e+00,-4.2223002567e+00,-3.0798330637e+00, 3.7905326800e+00, 8.4038996349e+00, 
     2.2816901707e+00, 4.1911712810e+00,-9.3658177215e+00,-1.2434457046e+01, 1.6771424507e+00, 
     2.3712452266e+00, 3.6980955070e+00, 2.5407601111e+00, 1.8976924460e-01,-1.0446470315e+00,
    -1.2992491105e+00,-8.5624882767e-01, 9.8686287993e+00, 1.1847431541e+01,-2.5193792283e-01,
    -3.0763224769e-01, 1.2429858182e-01, 1.3295440602e-02, 6.4178676394e-02, 1.2758462324e-01, 
     7.5131761426e-02, 1.1629004831e-01, 3.9639061816e-01, 6.7088705514e-01};
  RealType dhpsioverpsi_exact[] = {
    -1.6695881381e-02,-4.8902571790e-01,-1.2725397012e-01,-6.6714806635e-01, 6.9379259933e+00,
    -4.8393983437e+00, 7.4947765640e+00,-8.8373306290e-01,-6.8244030879e+00, 7.9150085031e-01,
    -1.4313255643e+00, 3.7225112718e+01, 1.7787191536e+01,-1.6672327906e+01,-4.1705496948e+01,
    -9.9674671566e+00,-2.0150790757e+01, 1.1226368249e+02, 1.2744525474e+02,-1.5801247401e+01,
    -1.3618595564e+01,-2.8161585388e+01,-1.4057266918e+01, 1.7626748997e+00, 7.8913447811e+00,
     9.2144952390e+00, 4.6068416473e+00,-9.3975889104e+01,-8.8298321426e+01, 1.5097063606e+01,
     1.8605794463e+01,-7.3647009565e+00,-5.9114663448e-01,-3.9243955679e+00,-7.8630886487e+00,
    -4.4437106408e+00,-7.0313362338e+00,-2.3986142270e+01,-4.0724297500e+01};
  Vector<ValueType> dlogpsi;
  Vector<ValueType> dhpsioverpsi;
  dlogpsi.resize(num_derivs);
  dhpsioverpsi.resize(num_derivs);
  std::fill(dlogpsi.begin(), dlogpsi.end(), 0);
  std::fill(dhpsioverpsi.begin(), dhpsioverpsi.end(), 0);

  // prepare variable set
  VariableSet optVars;
  optVars.clear();
  cj->checkInVariablesExclusive(optVars);
  optVars.resetIndex();
  cj->checkInVariablesExclusive(optVars);
  cj->checkOutVariables(optVars);
  optVars.print(std::cout);

  // test evaluateDerivatives
  cj->evaluateDerivatives(elec, optVars, dlogpsi, dhpsioverpsi);
  for(int p = 0; p < num_derivs; ++p)
  {
    CHECK ( dlogpsi_exact[p] == Approx(std::real(dlogpsi[p])).epsilon(1e-3) );
    CHECK ( dhpsioverpsi_exact[p] == Approx(std::real(dhpsioverpsi[p])).epsilon(1e-3) );
  }




  // test makeClone
  auto cj2_uptr                                = cj->makeClone(elec);
  CountingJastrow<CountingGaussianRegion>* cj2 = dynamic_cast<CountingJastrow<CountingGaussianRegion>*>(cj2_uptr.get());
  std::fill(dlogpsi.begin(), dlogpsi.end(), 0);
  std::fill(dhpsioverpsi.begin(), dhpsioverpsi.end(), 0);

  // prepare variable set
  VariableSet optVars2;
  optVars2.clear();
  cj2->checkInVariablesExclusive(optVars2);
  optVars2.resetIndex();
  cj2->checkInVariablesExclusive(optVars2);
  cj2->checkOutVariables(optVars2);
  optVars2.print(std::cout);

  cj2->evaluateDerivatives(elec, optVars2, dlogpsi, dhpsioverpsi);
  for(int p = 0; p < num_derivs; ++p)
  {
    CHECK ( dlogpsi_exact[p] == Approx(std::real(dlogpsi[p])).epsilon(1e-3) );
    CHECK ( dhpsioverpsi_exact[p] == Approx(std::real(dhpsioverpsi[p])).epsilon(1e-3) );
  }

  // test resetParameters, recompute
  for(int p = 0; p < num_derivs; ++p)
    optVars[p] = 0;
  cj->resetParametersExclusive(optVars);
  cj->recompute(elec);
  REQUIRE( cj->get_log_value() == LogValueType(0) );
#endif

}

} //namespace qmcplusplus
