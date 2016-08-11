
#include "catch.hpp"

#include "OhmmsData/Libxml2Doc.h"
#include "OhmmsPETE/OhmmsMatrix.h"
#include "Utilities/OhmmsInfo.h"
#include "Lattice/ParticleBConds.h"
#include "Particle/ParticleSet.h"
#include "Particle/DistanceTableData.h"
#include "Particle/DistanceTable.h"
#include "Particle/SymmetricDistanceTableData.h"
#include "QMCApp/ParticleSetPool.h"
#include "QMCWaveFunctions/OrbitalBase.h"
#include "QMCWaveFunctions/TrialWaveFunction.h"
#include "QMCWaveFunctions/Jastrow/TwoBodyJastrowOrbital.h"
#include "QMCWaveFunctions/Jastrow/BsplineFunctor.h"
#include "QMCWaveFunctions/Jastrow/BsplineJastrowBuilder.h"


#include <stdio.h>
#include <string>

using std::string;

namespace qmcplusplus
{

TEST_CASE("BSpline functor zero", "[wavefunction]")
{

  BsplineFunctor<double> bf;

  double r = 1.2;
  double u = bf.evaluate(r);
  REQUIRE(u == 0.0);
}

TEST_CASE("BSpline functor one", "[wavefunction]")
{

  BsplineFunctor<double> bf;

  bf.resize(1);

  double r = 1.2;
  double u = bf.evaluate(r);
  REQUIRE(u == 0.0);
}

TEST_CASE("BSpline builder Jastrow", "[wavefunction]")
{

  Communicate *c;
  OHMMS::Controller->initialize(0, NULL);
  c = OHMMS::Controller;
  OhmmsInfo("testlogfile");

  ParticleSet ions_;
  ParticleSet elec_;

  ions_.setName("ion");
  ions_.R.resize(1);
  ions_.R[0][0] = 2.0;
  ions_.R[0][1] = 0.0;
  ions_.R[0][2] = 0.0;
  ions_.create(1);

  elec_.setName("elec");
  elec_.R.resize(2);
  elec_.R[0][0] = 1.00;
  elec_.R[0][1] = 0.0;
  elec_.R[0][2] = 0.0;
  elec_.R[1][0] = 0.0;
  elec_.R[1][1] = 0.0;
  elec_.R[1][2] = 0.0;
  elec_.create(2);

  SpeciesSet &tspecies =  elec_.getSpeciesSet();
  int upIdx = tspecies.addSpecies("u");
  int downIdx = tspecies.addSpecies("d");
  int chargeIdx = tspecies.addAttribute("charge");
  tspecies(chargeIdx, upIdx) = -1;
  tspecies(chargeIdx, downIdx) = -1;

  elec_.addTable(ions_);
  elec_.update();


  TrialWaveFunction psi = TrialWaveFunction(c);
  // Need 1 electron and 1 proton, somehow
  //ParticleSet target = ParticleSet();
  ParticleSetPool ptcl = ParticleSetPool(c);

const char *particles = \
"<tmp> \
<jastrow name=\"J2\" type=\"Two-Body\" function=\"Bspline\" print=\"yes\"> \
   <correlation rcut=\"10\" size=\"10\" speciesA=\"u\" speciesB=\"u\"> \
      <coefficients id=\"uu\" type=\"Array\"> 0.02904699284 -0.1004179 -0.1752703883 -0.2232576505 -0.2728029201 -0.3253286875 -0.3624525145 -0.3958223107 -0.4268582166 -0.4394531176</coefficients> \
    </correlation> \
</jastrow> \
</tmp> \
";
  Libxml2Document doc;
  bool okay = doc.parseFromString(particles);
  REQUIRE(okay);

  xmlNodePtr root = doc.getRoot();

  xmlNodePtr jas1 = xmlFirstElementChild(root);

  BsplineJastrowBuilder jastrow(elec_, psi);
  bool build_okay = jastrow.put(jas1);
  REQUIRE(build_okay);

  OrbitalBase *orb = psi.getOrbitals()[0];

  typedef TwoBodyJastrowOrbital<BsplineFunctor<double> > J2Type;
  J2Type *j2 = dynamic_cast<J2Type *>(orb);
  REQUIRE(j2 != NULL);

  double logpsi = psi.evaluateLog(elec_);
  REQUIRE(logpsi == Approx(0.1012632641)); // note: number not validated
  
#if 0
  // write out values of the Bspline functor
  BsplineFunctor<double> *bf = j2->F[0];
  printf("NumParams = %d\n",bf->NumParams);
  printf("CuspValue = %g\n",bf->CuspValue);
  printf("DeltaR = %g\n",bf->DeltaR);
  printf("SplineCoeffs size = %d\n",bf->SplineCoefs.size());
  for (int j = 0; j < bf->SplineCoefs.size(); j++)
  {
    printf("%d %g\n",j,bf->SplineCoefs[j]);
  }
  printf("\n");

  for (int i = 0; i < 50; i++) {
    double r = 0.2*i;
    elec_.R[0][0] = r;
    elec_.update();
    double logpsi = psi.evaluateLog(elec_);
    double alt_val = bf->evaluate(r);
    printf("%g %g %g\n",r,logpsi,alt_val);
  }
#endif


}
}

