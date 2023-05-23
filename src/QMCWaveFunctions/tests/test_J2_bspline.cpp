//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Mark Dewing, markdewing@gmail.com, University of Illinois at Urbana-Champaign
//
// File created by: Mark Dewing, markdewing@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


#include "catch.hpp"

#include "OhmmsData/Libxml2Doc.h"
#include "OhmmsPETE/OhmmsMatrix.h"
#include "Particle/ParticleSet.h"
#include "QMCWaveFunctions/WaveFunctionComponent.h"
#include "QMCWaveFunctions/Jastrow/BsplineFunctor.h"
#include "QMCWaveFunctions/Jastrow/RadialJastrowBuilder.h"
#include "ParticleBase/ParticleAttribOps.h"
#include "QMCWaveFunctions/Jastrow/TwoBodyJastrow.h"

#include <cstdio>
#include <string>


// Uncomment to print information and values from the underlying functor
//#define PRINT_SPLINE_DATA

using std::string;

namespace qmcplusplus
{
using RealType     = WaveFunctionComponent::RealType;
using PsiValueType = WaveFunctionComponent::PsiValueType;

TEST_CASE("BSpline builder Jastrow J2", "[wavefunction]")
{
  Communicate* c = OHMMS::Controller;

  const SimulationCell simulation_cell;
  ParticleSet ions_(simulation_cell);
  ParticleSet elec_(simulation_cell);

  ions_.setName("ion");
  ions_.create({1});
  ions_.R[0] = {2.0, 0.0, 0.0};
  elec_.setName("elec");
  elec_.create({1, 1});
  elec_.R[0]                   = {1.00, 0.0, 0.0};
  elec_.R[1]                   = {0.0, 0.0, 0.0};
  SpeciesSet& tspecies         = elec_.getSpeciesSet();
  int upIdx                    = tspecies.addSpecies("u");
  int downIdx                  = tspecies.addSpecies("d");
  int chargeIdx                = tspecies.addAttribute("charge");
  tspecies(chargeIdx, upIdx)   = -1;
  tspecies(chargeIdx, downIdx) = -1;
  elec_.resetGroups();

  const char* particles = R"(<tmp>
<jastrow name="J2" type="Two-Body" function="Bspline" print="yes" gpu="no">
   <correlation rcut="10" size="10" speciesA="u" speciesB="d">
      <coefficients id="ud" type="Array"> 0.02904699284 -0.1004179 -0.1752703883 -0.2232576505 -0.2728029201 -0.3253286875 -0.3624525145 -0.3958223107 -0.4268582166 -0.4394531176</coefficients>
    </correlation>
</jastrow>
</tmp>
)";
  Libxml2Document doc;
  bool okay = doc.parseFromString(particles);
  REQUIRE(okay);

  xmlNodePtr root = doc.getRoot();

  xmlNodePtr jas1 = xmlFirstElementChild(root);

  RadialJastrowBuilder jastrow(c, elec_);

  using J2Type = TwoBodyJastrow<BsplineFunctor<RealType>>;
  auto j2_uptr = jastrow.buildComponent(jas1);
  J2Type* j2   = dynamic_cast<J2Type*>(j2_uptr.get());
  REQUIRE(j2);

  // update all distance tables
  elec_.update();

  double logpsi_real = std::real(j2->evaluateLog(elec_, elec_.G, elec_.L));
  CHECK(logpsi_real == Approx(0.1012632641)); // note: number not validated

  double KE = -0.5 * (Dot(elec_.G, elec_.G) + Sum(elec_.L));
  CHECK(KE == Approx(-0.1616624771)); // note: number not validated

  UniqueOptObjRefs opt_obj_refs;
  j2->extractOptimizableObjectRefs(opt_obj_refs);
  REQUIRE(opt_obj_refs.size() == 1);

  opt_variables_type optvars;
  Vector<WaveFunctionComponent::ValueType> dlogpsi;
  Vector<WaveFunctionComponent::ValueType> dhpsioverpsi;

  for (OptimizableObject& obj : opt_obj_refs)
    obj.checkInVariablesExclusive(optvars);
  optvars.resetIndex();
  const int NumOptimizables(optvars.size());
  j2->checkOutVariables(optvars);
  dlogpsi.resize(NumOptimizables);
  dhpsioverpsi.resize(NumOptimizables);
  j2->evaluateDerivatives(elec_, optvars, dlogpsi, dhpsioverpsi);

  app_log() << std::endl << "reporting dlogpsi and dhpsioverpsi" << std::scientific << std::endl;
  for (int iparam = 0; iparam < NumOptimizables; iparam++)
    app_log() << "param=" << iparam << " : " << dlogpsi[iparam] << "  " << dhpsioverpsi[iparam] << std::endl;
  app_log() << std::endl;

  CHECK(std::real(dlogpsi[2]) == Approx(-0.2211666667));
  CHECK(std::real(dhpsioverpsi[3]) == Approx(0.1331717179));


  // now test evaluateHessian
  WaveFunctionComponent::HessVector grad_grad_psi;
  grad_grad_psi.resize(elec_.getTotalNum());
  grad_grad_psi = 0.0;

  app_log() << "eval hess" << std::endl;
  j2->evaluateHessian(elec_, grad_grad_psi);
  std::vector<double> hess_values = {
      -0.0627236, 0, 0, 0, 0.10652, 0, 0, 0, 0.10652, -0.0627236, 0, 0, 0, 0.10652, 0, 0, 0, 0.10652,
  };

  int m = 0;
  for (int n = 0; n < elec_.getTotalNum(); n++)
    for (int i = 0; i < OHMMS_DIM; i++)
      for (int j = 0; j < OHMMS_DIM; j++, m++)
      {
        CHECK(std::real(grad_grad_psi[n](i, j)) == Approx(hess_values[m]));
      }


  struct JValues
  {
    double r;
    double u;
    double du;
    double ddu;
  };

  // Cut and paste from output of gen_bspline_jastrow.py
  const int N     = 20;
  JValues Vals[N] = {{0.00, 0.1374071801, -0.5, 0.7866949593},
                     {0.60, -0.04952403966, -0.1706645865, 0.3110897524},
                     {1.20, -0.121361995, -0.09471371432, 0.055337302},
                     {1.80, -0.1695590431, -0.06815900213, 0.0331784053},
                     {2.40, -0.2058414025, -0.05505192964, 0.01049597156},
                     {3.00, -0.2382237097, -0.05422744821, -0.002401552969},
                     {3.60, -0.2712606182, -0.05600918024, -0.003537553803},
                     {4.20, -0.3047843679, -0.05428535477, 0.0101841028},
                     {4.80, -0.3347515004, -0.04506573714, 0.01469003611},
                     {5.40, -0.3597048574, -0.03904232165, 0.005388015505},
                     {6.00, -0.3823503292, -0.03657502025, 0.003511355265},
                     {6.60, -0.4036800017, -0.03415678101, 0.007891305516},
                     {7.20, -0.4219818468, -0.02556305518, 0.02075444724},
                     {7.80, -0.4192355508, 0.06799438701, 0.3266190181},
                     {8.40, -0.3019238309, 0.32586994, 0.2880861726},
                     {9.00, -0.09726352421, 0.2851358014, -0.4238666348},
                     {9.60, -0.006239062395, 0.04679296796, -0.2339648398},
                     {10.20, 0, 0, 0},
                     {10.80, 0, 0, 0},
                     {11.40, 0, 0, 0}};


  BsplineFunctor<RealType>* bf = j2->getPairFunctions()[0];

  for (int i = 0; i < N; i++)
  {
    RealType dv  = 0.0;
    RealType ddv = 0.0;
    RealType val = bf->evaluate(Vals[i].r, dv, ddv);
    CHECK(Vals[i].u == Approx(val));
    CHECK(Vals[i].du == Approx(dv));
    CHECK(Vals[i].ddu == Approx(ddv));
  }

#ifdef PRINT_SPLINE_DATA
  // write out values of the Bspline functor
  //BsplineFunctor<double> *bf = j2->F[0];
  printf("NumParams = %d\n", bf->NumParams);
  printf("CuspValue = %g\n", bf->CuspValue);
  printf("DeltaR = %g\n", bf->DeltaR);
  printf("SplineCoeffs size = %d\n", bf->SplineCoefs.size());
  for (int j = 0; j < bf->SplineCoefs.size(); j++)
  {
    printf("%d %g\n", j, bf->SplineCoefs[j]);
  }
  printf("\n");

  for (int i = 0; i < 20; i++)
  {
    double r      = 0.6 * i;
    elec_.R[0][0] = r;
    elec_.update();
    double logpsi_real = std::real(j2->evaluateLog(elec_, elec_.G, elec_.L));
    //double alt_val = bf->evaluate(r);
    double dv      = 0.0;
    double ddv     = 0.0;
    double alt_val = bf->evaluate(r, dv, ddv);
    printf("%g %g %g %g %g\n", r, logpsi_real, alt_val, dv, ddv);
  }
#endif

  using ValueType = QMCTraits::ValueType;
  using PosType   = QMCTraits::PosType;

  // set virtutal particle position
  PosType newpos(0.3, 0.2, 0.5);

  elec_.makeVirtualMoves(newpos);
  std::vector<ValueType> ratios(elec_.getTotalNum());
  j2->evaluateRatiosAlltoOne(elec_, ratios);

  CHECK(std::real(ratios[0]) == Approx(0.9522052017));
  CHECK(std::real(ratios[1]) == Approx(0.9871985577));

  elec_.makeMove(0, newpos - elec_.R[0]);
  PsiValueType ratio_0 = j2->ratio(elec_, 0);
  elec_.rejectMove(0);

  CHECK(std::real(ratio_0) == Approx(0.9522052017));

  VirtualParticleSet VP(elec_, 2);
  std::vector<PosType> newpos2(2);
  std::vector<ValueType> ratios2(2);
  newpos2[0] = newpos - elec_.R[1];
  newpos2[1] = PosType(0.2, 0.5, 0.3) - elec_.R[1];
  VP.makeMoves(elec_, 1, newpos2);
  j2->evaluateRatios(VP, ratios2);

  CHECK(std::real(ratios2[0]) == Approx(0.9871985577));
  CHECK(std::real(ratios2[1]) == Approx(0.9989268241));

  //test acceptMove
  elec_.makeMove(1, newpos - elec_.R[1]);
  PsiValueType ratio_1 = j2->ratio(elec_, 1);
  j2->acceptMove(elec_, 1);
  elec_.acceptMove(1);

  CHECK(std::real(ratio_1) == Approx(0.9871985577));
  CHECK(std::real(j2->get_log_value()) == Approx(0.0883791773));
}
} // namespace qmcplusplus
