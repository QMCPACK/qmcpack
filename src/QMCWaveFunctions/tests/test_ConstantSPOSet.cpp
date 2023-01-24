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
  app_log()<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n";
  //For now, do a small square case.
  const int nelec = 2;
  const int norb  = 2;

  using RealType          = QMCTraits::RealType;
  using ValueType         = QMCTraits::ValueType;
  using ValueVector = OrbitalSetTraits<ValueType>::ValueVector;
  using GradVector = OrbitalSetTraits<ValueType>::GradVector;

  ValueVector row0{ValueType(0.92387953), ValueType(0.92387953)};
  ValueVector row1{ValueType(0.29131988), ValueType(0.81078057)};
 
  using ValueMatrix = OrbitalSetTraits<ValueType>::ValueMatrix;
  ValueMatrix spomat;
  spomat.resize(nelec, nelec);
  for (int iorb = 0; iorb < norb; iorb++)
  {
    spomat(0, iorb) = row0[iorb];
    spomat(1, iorb) = row1[iorb];
  }
  auto sposet = std::make_unique<ConstantSPOSet>("constant_spo",spomat);

  const SimulationCell simulation_cell;
  ParticleSet elec(simulation_cell);

  elec.create({3});

  ValueVector psiV = {0.0, 0.0};
  ValueVector psiL = {0.0, 0.0};
  GradVector  psiG;
  psiG.resize(norb);

  sposet->evaluateValue(elec, 0, psiV);
 
  CHECK( psiV[0] == row0[0]);
  CHECK( psiV[1] == row0[1]);

  
  psiV=0.0;

  sposet->evaluateValue(elec, 1, psiV);
  CHECK( psiV[0] == row1[0]);
  CHECK( psiV[1] == row1[1]);
 
  psiV=0.0;

  sposet->evaluateVGL(elec,1,psiV,psiG,psiL);
  CHECK( psiV[0] == row1[0]);
  CHECK( psiV[1] == row1[1]);

  CHECK( psiG[0][0] == 0.0 );
  CHECK( psiG[0][1] == 0.0 );
  CHECK( psiG[0][2] == 0.0 );

  CHECK( psiG[1][0] == 0.0 );
  CHECK( psiG[1][1] == 0.0 );
  CHECK( psiG[1][2] == 0.0 );

  CHECK( psiL[0] == 0.0);
  CHECK( psiL[1] == 0.0);
  
}
} // namespace qmcplusplus
