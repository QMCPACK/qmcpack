//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2023 QMCPACK developers.
//
// File developed by: Raymond Clay, rclay@sandia.gov, Sandia National Laboratories
//
// File created by: Raymond Clay, rclay@sandia.gov, Sandia National Laboratories
//////////////////////////////////////////////////////////////////////////////////////


#include "catch.hpp"
#include "Configuration.h"
#include "QMCWaveFunctions/OrbitalSetTraits.h"
#include "QMCWaveFunctions/ConstantSPOSet.h"

namespace qmcplusplus
{
//Ray:  Figure out how to template me on value type.
TEST_CASE("ConstantSPOSet", "[wavefunction]")
{
  //For now, do a small square case.
  const int nelec = 2;
  const int norb  = 2;

  using RealType    = QMCTraits::RealType;
  using ValueType   = QMCTraits::ValueType;
  using GradType   =  QMCTraits::GradType;
  using ValueVector = OrbitalSetTraits<ValueType>::ValueVector;
  using GradVector  = OrbitalSetTraits<ValueType>::GradVector;
  using ValueMatrix = OrbitalSetTraits<ValueType>::ValueMatrix;
  using GradMatrix = OrbitalSetTraits<ValueType>::GradMatrix;

  ValueVector row0{ValueType(0.92387953), ValueType(0.92387953)};
  ValueVector row1{ValueType(0.29131988), ValueType(0.81078057)};

  GradVector  grow0{GradType({ -2.22222, -1.11111, 0.33333}), GradType({8.795388,-0.816057,-0.9238793})};
  GradVector  grow1{GradType({  2.22222,  1.11111, -0.33333}), GradType({-8.795388,0.816057,0.9238793})};

  ValueVector lrow0{ValueType(-0.2234545), ValueType(0.72340234)};
  ValueVector lrow1{ValueType(-12.291810), ValueType(6.879057)};

  
  ValueMatrix spomat;
  GradMatrix gradspomat;
  ValueMatrix laplspomat;

  spomat.resize(nelec, norb);
  gradspomat.resize(nelec, norb);
  laplspomat.resize(nelec, norb);

  for (int iorb = 0; iorb < norb; iorb++)
  {
    spomat(0, iorb) = row0[iorb];
    spomat(1, iorb) = row1[iorb];
    
    gradspomat(0,iorb) = grow0[iorb];
    gradspomat(1,iorb) = grow1[iorb];
   
    laplspomat(0,iorb) = lrow0[iorb];
    laplspomat(1,iorb) = lrow1[iorb];
  }


  const SimulationCell simulation_cell;
  ParticleSet elec(simulation_cell);

  elec.create({nelec});

  ValueVector psiV = {0.0, 0.0};
  ValueVector psiL = {0.0, 0.0};
  GradVector psiG;
  psiG.resize(norb);

  //Test of value only constructor.
  auto sposet = std::make_unique<ConstantSPOSet>("constant_spo", nelec, norb);
  sposet->setRefVals(spomat);
  sposet->setRefEGrads(gradspomat);
  sposet->setRefELapls(laplspomat);

  sposet->evaluateValue(elec, 0, psiV);

  CHECK(psiV[0] == row0[0]);
  CHECK(psiV[1] == row0[1]);


  psiV = 0.0;

  sposet->evaluateValue(elec, 1, psiV);
  CHECK(psiV[0] == row1[0]);
  CHECK(psiV[1] == row1[1]);

  psiV = 0.0;

  sposet->evaluateVGL(elec,1,psiV,psiG,psiL);
   
  CHECK(psiV[0] == row1[0]);
  CHECK(psiV[1] == row1[1]);
  
  CHECK(psiG[0][0] == grow1[0][0]); 
  CHECK(psiG[0][1] == grow1[0][1]); 
  CHECK(psiG[0][2] == grow1[0][2]); 

  CHECK(psiG[1][0] == grow1[1][0]); 
  CHECK(psiG[1][1] == grow1[1][1]); 
  CHECK(psiG[1][2] == grow1[1][2]); 

  CHECK(psiL[0] == lrow1[0]);
  CHECK(psiL[1] == lrow1[1]);

  //Test of evaluate_notranspose.
  ValueMatrix phimat,lphimat;
  GradMatrix gphimat;
  phimat.resize(nelec,norb);
  gphimat.resize(nelec,norb);
  lphimat.resize(nelec,norb);

  const int first_index=0; //Only 2 electrons in this case.  
  const int last_index=2; 
  sposet->evaluate_notranspose(elec,first_index,last_index,phimat,gphimat,lphimat);
 
  CHECK(phimat[0][0] == row0[0]);
  CHECK(phimat[0][1] == row0[1]);
  CHECK(phimat[1][0] == row1[0]);
  CHECK(phimat[1][1] == row1[1]);

  CHECK(lphimat[0][0] == lrow0[0]);
  CHECK(lphimat[0][1] == lrow0[1]);
  CHECK(lphimat[1][0] == lrow1[0]);
  CHECK(lphimat[1][1] == lrow1[1]);

  //Test of makeClone()
  auto sposet_vgl2 = sposet->makeClone();
  phimat=0.0;
  gphimat=0.0;
  lphimat=0.0;

  sposet_vgl2->evaluate_notranspose(elec,first_index,last_index,phimat,gphimat,lphimat);
  CHECK(phimat[0][0] == row0[0]);
  CHECK(phimat[0][1] == row0[1]);
  CHECK(phimat[1][0] == row1[0]);
  CHECK(phimat[1][1] == row1[1]);

  CHECK(lphimat[0][0] == lrow0[0]);
  CHECK(lphimat[0][1] == lrow0[1]);
  CHECK(lphimat[1][0] == lrow1[0]);
  CHECK(lphimat[1][1] == lrow1[1]);

  //Lastly, check if name is correct.
  std::string myname=sposet_vgl2->getClassName();
  std::string targetstring("ConstantSPOSet");
  CHECK(myname == targetstring); 

}
} // namespace qmcplusplus
