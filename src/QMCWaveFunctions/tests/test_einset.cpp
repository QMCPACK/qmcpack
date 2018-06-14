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
#include "Lattice/ParticleBConds.h"
#include "Particle/ParticleSet.h"
#include "Particle/DistanceTableData.h"
#include "Particle/DistanceTable.h"
#include "Particle/SymmetricDistanceTableData.h"
#include "QMCApp/ParticleSetPool.h"
#include "QMCWaveFunctions/OrbitalBase.h"
#include "QMCWaveFunctions/TrialWaveFunction.h"
#include "QMCWaveFunctions/EinsplineSetBuilder.h"


#include <stdio.h>
#include <string>

using std::string;

namespace qmcplusplus
{

TEST_CASE("Einspline SPO from HDF", "[wavefunction]")
{

  Communicate *c;
  OHMMS::Controller->initialize(0, NULL);
  c = OHMMS::Controller;

  ParticleSet ions_;
  ParticleSet elec_;

  ions_.setName("ion");
  ions_.create(2);
  ions_.R[0][0] = 0.0;
  ions_.R[0][1] = 0.0;
  ions_.R[0][2] = 0.0;
  ions_.R[1][0] = 1.68658058;
  ions_.R[1][1] = 1.68658058;
  ions_.R[1][2] = 1.68658058;


  elec_.setName("elec");
  elec_.create(2);
  elec_.R[0][0] = 0.00;
  elec_.R[0][1] = 0.0;
  elec_.R[0][2] = 0.0;
  elec_.R[1][0] = 0.0;
  elec_.R[1][1] = 1.0;
  elec_.R[1][2] = 0.0;

  // monoO
 /*
  elec_.Lattice.R(0,0) = 5.10509515;
  elec_.Lattice.R(0,1) = -3.23993545;
  elec_.Lattice.R(0,2) = 0.0;
  elec_.Lattice.R(1,0) = 5.10509515;
  elec_.Lattice.R(1,1) = 3.23993545;
  elec_.Lattice.R(1,2) = 0.0;
  elec_.Lattice.R(2,0) = -6.49690625;
  elec_.Lattice.R(2,1) = 0.0;
  elec_.Lattice.R(2,2) = 7.08268015;
 */

  // diamondC_1x1x1
  elec_.Lattice.R(0,0) = 3.37316115;
  elec_.Lattice.R(0,1) = 3.37316115;
  elec_.Lattice.R(0,2) = 0.0;
  elec_.Lattice.R(1,0) = 0.0;
  elec_.Lattice.R(1,1) = 3.37316115;
  elec_.Lattice.R(1,2) = 3.37316115;
  elec_.Lattice.R(2,0) = 3.37316115;
  elec_.Lattice.R(2,1) = 0.0;
  elec_.Lattice.R(2,2) = 3.37316115;

  SpeciesSet &tspecies =  elec_.getSpeciesSet();
  int upIdx = tspecies.addSpecies("u");
  int downIdx = tspecies.addSpecies("d");
  int chargeIdx = tspecies.addAttribute("charge");
  tspecies(chargeIdx, upIdx) = -1;
  tspecies(chargeIdx, downIdx) = -1;

#ifdef ENABLE_SOA
  elec_.addTable(ions_,DT_SOA);
#else
  elec_.addTable(ions_,DT_AOS);
#endif
  elec_.resetGroups();
  elec_.update();


  TrialWaveFunction psi = TrialWaveFunction(c);
  // Need 1 electron and 1 proton, somehow
  //ParticleSet target = ParticleSet();
  ParticleSetPool ptcl = ParticleSetPool(c);
  ptcl.addParticleSet(&elec_);
  ptcl.addParticleSet(&ions_);

//diamondC_1x1x1
const char *particles = 
"<tmp> \
<determinantset type=\"einspline\" href=\"pwscf.pwscf.h5\" tilematrix=\"1 0 0 0 1 0 0 0 1\" twistnum=\"0\" source=\"ion\" meshfactor=\"1.0\" precision=\"float\" size=\"4\"/> \
</tmp> \
";

// monoO
//<determinantset type=\"einspline\" href=\"pwscf.pwscf.h5\" tilematrix=\"1 0 0 0 1 0 0 0 1\" twistnum=\"0\" source=\"ion\" meshfactor=\"1.0\" precision=\"float\" size=\"6\"/> \

  Libxml2Document doc;
  bool okay = doc.parseFromString(particles);
  REQUIRE(okay);

  xmlNodePtr root = doc.getRoot();

  xmlNodePtr ein1 = xmlFirstElementChild(root);

  EinsplineSetBuilder einSet(elec_, ptcl.getPool(), ein1);
  SPOSetBase *spo = einSet.createSPOSetFromXML(ein1);
  REQUIRE(spo != NULL);

#if !defined(QMC_CUDA) || defined(QMC_COMPLEX)
  // due to the different ordering of bands skip the tests on CUDA+Real builds
  // checking evaluations, reference values are not independently generated.
  // for vgl
  SPOSetBase::ValueMatrix_t psiM(elec_.R.size(),spo->getOrbitalSetSize());
  SPOSetBase::GradMatrix_t dpsiM(elec_.R.size(),spo->getOrbitalSetSize());
  SPOSetBase::ValueMatrix_t d2psiM(elec_.R.size(),spo->getOrbitalSetSize());
  spo->evaluate_notranspose(elec_, 0, elec_.R.size(), psiM, dpsiM, d2psiM);

  // value
  REQUIRE(psiM[1][0] == ComplexApprox(-0.8886948824).compare_real_only());
  REQUIRE(psiM[1][1] == ComplexApprox(1.4194120169).compare_real_only());
  // grad
  REQUIRE(dpsiM[1][0][0] == ComplexApprox(-0.0000183403).compare_real_only());
  REQUIRE(dpsiM[1][0][1] == ComplexApprox(0.1655139178).compare_real_only());
  REQUIRE(dpsiM[1][0][2] == ComplexApprox(-0.0000193077).compare_real_only());
  REQUIRE(dpsiM[1][1][0] == ComplexApprox(-1.3131694794).compare_real_only());
  REQUIRE(dpsiM[1][1][1] == ComplexApprox(-1.1174004078).compare_real_only());
  REQUIRE(dpsiM[1][1][2] == ComplexApprox(-0.8462534547).compare_real_only());
  // lapl
  REQUIRE(d2psiM[1][0] == ComplexApprox(1.3313053846).compare_real_only());
  REQUIRE(d2psiM[1][1] == ComplexApprox(-4.712583065).compare_real_only());

  // for vgh
  SPOSetBase::ValueVector_t psiV(psiM[1],spo->getOrbitalSetSize());
  SPOSetBase::GradVector_t dpsiV(dpsiM[1],spo->getOrbitalSetSize());
  SPOSetBase::HessVector_t ddpsiV(spo->getOrbitalSetSize());
  spo->evaluate(elec_, 1, psiV, dpsiV, ddpsiV);

  //hess
  REQUIRE(ddpsiV[1](0,0) == ComplexApprox(-2.3160984034).compare_real_only());
  REQUIRE(ddpsiV[1](0,1) == ComplexApprox(1.8089479397).compare_real_only());
  REQUIRE(ddpsiV[1](0,2) == ComplexApprox(0.5608575749).compare_real_only());
  REQUIRE(ddpsiV[1](1,0) == ComplexApprox(1.8089479397).compare_real_only());
  REQUIRE(ddpsiV[1](1,1) == ComplexApprox(-0.0799493616).compare_real_only());
  REQUIRE(ddpsiV[1](1,2) == ComplexApprox(0.5237969314).compare_real_only());
  REQUIRE(ddpsiV[1](2,0) == ComplexApprox(0.5608575749).compare_real_only());
  REQUIRE(ddpsiV[1](2,1) == ComplexApprox(0.5237969314).compare_real_only());
  REQUIRE(ddpsiV[1](2,2) == ComplexApprox(-2.316497764).compare_real_only());
#endif

#if 0
  // Dump values of the orbitals
  int orbSize= spo->getOrbitalSetSize();
  int basisSize= spo->getBasisSetSize();
  printf("orb size = %d basis set size = %d\n",orbSize, basisSize);

  FILE *fspo = fopen("spo.dat", "w");
  for (int ix = 0; ix < 30; ix++) {
    for (int iy = 0; iy < 30; iy++) {
      for (int iz = 0; iz < 30; iz++) {
        double x = 0.1*ix - 1.5;
        double y = 0.1*iy - 1.5;
        double z = 0.1*iz - 1.5;
        elec_.R[0][0] = x;
        elec_.R[0][1] = y;
        elec_.R[0][2] = z;
        elec_.update();
        SPOSetBase::ValueVector_t orbs(orbSize);
        spo->evaluate(elec_, 0, orbs);
        fprintf(fspo, "%g %g %g",x,y,z);
        for (int j = 0; j < orbSize; j++) {
          fprintf(fspo, " %g ",orbs[j]);
        }
        fprintf(fspo, "\n");
      }
    }
  }
  fclose(fspo);
#endif

}

TEST_CASE("EinsplineSetBuilder CheckLattice", "[wavefunction]")
{

  Communicate *c;
  OHMMS::Controller->initialize(0, NULL);
  c = OHMMS::Controller;

  ParticleSet *elec = new ParticleSet;

  elec->setName("elec");
  std::vector<int> agroup(2);
  agroup[0] = 1;
  agroup[1] = 1;
  elec->create(agroup);
  elec->R[0][0] = 0.00;
  elec->R[0][1] = 0.0;
  elec->R[0][2] = 0.0;
  elec->R[1][0] = 0.0;
  elec->R[1][1] = 1.0;
  elec->R[1][2] = 0.0;

  elec->Lattice.R = 0.0;
  elec->Lattice.R(0,0) = 1.0;
  elec->Lattice.R(1,1) = 1.0;
  elec->Lattice.R(2,2) = 1.0;

  EinsplineSetBuilder::PtclPoolType ptcl_map;
  ptcl_map["e"] = elec;

  xmlNodePtr cur = NULL;
  EinsplineSetBuilder esb(*elec, ptcl_map, cur);

  esb.SuperLattice = 0.0;
  esb.SuperLattice(0,0) = 1.0;
  esb.SuperLattice(1,1) = 1.0;
  esb.SuperLattice(2,2) = 1.0;

  REQUIRE(esb.CheckLattice());

  esb.SuperLattice(0,0) = 1.1;
  REQUIRE_FALSE(esb.CheckLattice());
}
}

