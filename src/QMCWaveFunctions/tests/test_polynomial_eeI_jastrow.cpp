//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////
    
    

#include "catch.hpp"

#include "OhmmsData/Libxml2Doc.h"
#include "OhmmsPETE/OhmmsMatrix.h"
#include "Lattice/ParticleBConds.h"
#include "Particle/ParticleSet.h"
#include "Particle/DistanceTableData.h"
#include "Particle/DistanceTable.h"
#include "Particle/SymmetricDistanceTableData.h"
#include "QMCWaveFunctions/OrbitalBase.h"
#include "QMCWaveFunctions/TrialWaveFunction.h"
#include "QMCWaveFunctions/Jastrow/PolynomialFunctor3D.h"
#include "QMCWaveFunctions/Jastrow/eeI_JastrowOrbital.h"
#include "QMCWaveFunctions/Jastrow/JeeIOrbitalSoA.h"
#include "QMCWaveFunctions/Jastrow/eeI_JastrowBuilder.h"
#include "ParticleBase/ParticleAttribOps.h"


#include <stdio.h>
#include <string>

using std::string;

namespace qmcplusplus
{

TEST_CASE("PolynomialFunctor3D functor zero", "[wavefunction]")
{

  PolynomialFunctor3D functor;

  double r = 1.2;
  double u = functor.evaluate(r,r,r);
  REQUIRE(u == 0.0);
}

TEST_CASE("PolynomialFunctor3D Jastrow", "[wavefunction]")
{

  Communicate *c;
  OHMMS::Controller->initialize(0, NULL);
  c = OHMMS::Controller;

  ParticleSet ions_;
  ParticleSet elec_;

  ions_.setName("ion");
  ions_.create(2);
  ions_.R[0][0] = 2.0;
  ions_.R[0][1] = 0.0;
  ions_.R[0][2] = 0.0;
  ions_.R[1][0] = -2.0;
  ions_.R[1][1] = 0.0;
  ions_.R[1][2] = 0.0;
  SpeciesSet& source_species(ions_.getSpeciesSet());
  source_species.addSpecies("O");
  ions_.RSoA=ions_.R;
  //ions_.resetGroups();

  elec_.setName("elec");
  std::vector<int> ud(2); ud[0]=ud[1]=2;
  elec_.create(ud);
  elec_.R[0][0] = 1.00;
  elec_.R[0][1] = 0.0;
  elec_.R[0][2] = 0.0;
  elec_.R[1][0] = 0.0;
  elec_.R[1][1] = 0.0;
  elec_.R[1][2] = 0.0;
  elec_.R[2][0] = -1.00;
  elec_.R[2][1] = 0.0;
  elec_.R[2][2] = 0.0;
  elec_.R[3][0] = 0.0;
  elec_.R[3][1] = 0.0;
  elec_.R[3][2] = 2.0;

  SpeciesSet& target_species(elec_.getSpeciesSet());
  int upIdx = target_species.addSpecies("u");
  int downIdx = target_species.addSpecies("d");
  int chargeIdx = target_species.addAttribute("charge");
  target_species(chargeIdx, upIdx) = -1;
  target_species(chargeIdx, downIdx) = -1;
  //elec_.resetGroups();

#ifdef ENABLE_SOA
  elec_.addTable(ions_,DT_SOA);
#else
  elec_.addTable(ions_,DT_AOS);
#endif
  elec_.update();

  TrialWaveFunction psi = TrialWaveFunction(c);

const char *particles = \
"<tmp> \
    <jastrow name=\"J3\" type=\"eeI\" function=\"polynomial\" source=\"ion\" print=\"yes\"> \
      <correlation ispecies=\"O\" especies=\"u\" isize=\"3\" esize=\"3\" rcut=\"10\"> \
        <coefficients id=\"uuO\" type=\"Array\" optimize=\"yes\"> 8.227710241e-06 2.480817653e-06 -5.354068112e-06 -1.112644787e-05 -2.208006078e-06 5.213121933e-06 -1.537865869e-05 8.899030233e-06 6.257255156e-06 3.214580988e-06 -7.716743107e-06 -5.275682077e-06 -1.778457637e-06 7.926231121e-06 1.767406868e-06 5.451359059e-08 2.801423724e-06 4.577282736e-06 7.634608083e-06 -9.510673173e-07 -2.344131575e-06 -1.878777219e-06 3.937363358e-07 5.065353773e-07 5.086724869e-07 -1.358768154e-07</coefficients> \
      </correlation> \
      <correlation ispecies=\"O\" especies1=\"u\" especies2=\"d\" isize=\"3\" esize=\"3\" rcut=\"10\"> \
        <coefficients id=\"udO\" type=\"Array\" optimize=\"yes\"> -6.939530224e-06 2.634169299e-05 4.046077477e-05 -8.002682388e-06 -5.396795988e-06 6.697370507e-06 5.433953051e-05 -6.336849668e-06 3.680471431e-05 -2.996059772e-05 1.99365828e-06 -3.222705626e-05 -8.091669063e-06 4.15738535e-06 4.843939112e-06 3.563650208e-07 3.786332474e-05 -1.418336941e-05 2.282691374e-05 1.29239286e-06 -4.93580873e-06 -3.052539228e-06 9.870288001e-08 1.844286407e-06 2.970561871e-07 -4.364303677e-08</coefficients> \
      </correlation> \
    </jastrow> \
</tmp> \
";
  Libxml2Document doc;
  bool okay = doc.parseFromString(particles);
  REQUIRE(okay);

  xmlNodePtr root = doc.getRoot();

  xmlNodePtr jas_eeI = xmlFirstElementChild(root);

  eeI_JastrowBuilder jastrow(elec_, psi, ions_);
  bool build_okay = jastrow.put(jas_eeI);
  REQUIRE(build_okay);

  OrbitalBase *orb = psi.getOrbitals()[0];

#ifdef ENABLE_SOA
  typedef JeeIOrbitalSoA<PolynomialFunctor3D> J3Type;
#else
  typedef eeI_JastrowOrbital<PolynomialFunctor3D> J3Type;
#endif
  J3Type *j3 = dynamic_cast<J3Type *>(orb);
  REQUIRE(j3 != NULL);

  double logpsi = psi.evaluateLog(elec_);
  REQUIRE(logpsi == Approx(-1.193457749)); // note: number not validated

  double KE = -0.5*(Dot(elec_.G,elec_.G)+Sum(elec_.L));
  REQUIRE(KE == Approx(-0.058051245)); // note: number not validated

  typedef QMCTraits::ValueType ValueType;
  typedef QMCTraits::PosType PosType;

  // set virtutal particle position
  PosType newpos(0.3,0.2,0.5);

  elec_.makeVirtualMoves(newpos);
  std::vector<ValueType> ratios(elec_.getTotalNum());
  j3->evaluateRatiosAlltoOne(elec_,ratios);

  REQUIRE(ratios[0] == ComplexApprox(0.8744938582).compare_real_only());
  REQUIRE(ratios[1] == ComplexApprox(1.0357541137).compare_real_only());
  REQUIRE(ratios[2] == ComplexApprox(0.8302245609).compare_real_only());
  REQUIRE(ratios[3] == ComplexApprox(0.7987703724).compare_real_only());

  elec_.makeMove(0, newpos-elec_.R[0]);
  ValueType ratio_0 = j3->ratio(elec_,0);
  elec_.rejectMove(0);

  elec_.makeMove(1, newpos-elec_.R[1]);
  ValueType ratio_1 = j3->ratio(elec_,1);
  elec_.rejectMove(1);

  elec_.makeMove(2, newpos-elec_.R[2]);
  ValueType ratio_2 = j3->ratio(elec_,2);
  elec_.rejectMove(2);

  elec_.makeMove(3, newpos-elec_.R[3]);
  ValueType ratio_3 = j3->ratio(elec_,3);
  elec_.rejectMove(3);

  REQUIRE(ratio_0 == ComplexApprox(0.8744938582).compare_real_only());
  REQUIRE(ratio_1 == ComplexApprox(1.0357541137).compare_real_only());
  REQUIRE(ratio_2 == ComplexApprox(0.8302245609).compare_real_only());
  REQUIRE(ratio_3 == ComplexApprox(0.7987703724).compare_real_only());

  opt_variables_type optvars;
  std::vector<OrbitalBase::RealType> dlogpsi;
  std::vector<OrbitalBase::RealType> dhpsioverpsi;

  psi.checkInVariables(optvars);
  optvars.resetIndex();
  const int NumOptimizables(optvars.size());
  psi.checkOutVariables(optvars);
  dlogpsi.resize(NumOptimizables);
  dhpsioverpsi.resize(NumOptimizables);
  psi.evaluateDerivatives(elec_, optvars, dlogpsi, dhpsioverpsi);

  std::cout << std::endl << "reporting dlogpsi and dhpsioverpsi" << std::scientific << std::endl;
  for(int iparam=0; iparam<NumOptimizables; iparam++)
    std::cout << "param=" << iparam << " : " << dlogpsi[iparam] << "  " << dhpsioverpsi[iparam] << std::endl;

  REQUIRE(dlogpsi[43] == Approx(1.3358726814e+05));
  REQUIRE(dhpsioverpsi[43] == Approx(-2.3246270644e+05));

}
}

