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
#include "QMCWaveFunctions/WaveFunctionTypes.hpp"
#include "QMCWaveFunctions/tests/ConstantSPOSet.h"
#include "Utilities/for_testing/checkMatrix.hpp"
namespace qmcplusplus
{
//Ray:  Figure out how to template me on value type.
TEST_CASE("ConstantSPOSet", "[wavefunction]")
{
  //For now, do a small square case.
  const int nelec   = 2;
  const int norb    = 2;
  using WF          = WaveFunctionTypes<QMCTraits::ValueType, QMCTraits::FullPrecValueType>;
  using Real        = WF::Real;
  using Value       = WF::Value;
  using Grad        = WF::Grad;
  using ValueVector = Vector<Value>;
  using GradVector  = Vector<Grad>;
  using ValueMatrix = Matrix<Value>;
  using GradMatrix  = Matrix<Grad>;

  ValueVector row0{Value(0.92387953), Value(0.92387953)};
  ValueVector row1{Value(0.29131988), Value(0.81078057)};

  GradVector grow0{Grad({-2.22222, -1.11111, 0.33333}), Grad({8.795388, -0.816057, -0.9238793})};
  GradVector grow1{Grad({2.22222, 1.11111, -0.33333}), Grad({-8.795388, 0.816057, 0.9238793})};

  ValueVector lrow0{Value(-0.2234545), Value(0.72340234)};
  ValueVector lrow1{Value(-12.291810), Value(6.879057)};


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

    gradspomat(0, iorb) = grow0[iorb];
    gradspomat(1, iorb) = grow1[iorb];

    laplspomat(0, iorb) = lrow0[iorb];
    laplspomat(1, iorb) = lrow1[iorb];
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

  sposet->evaluateVGL(elec, 1, psiV, psiG, psiL);

  for (int iorb = 0; iorb < norb; iorb++)
  {
    CHECK(psiV[iorb] == row1[iorb]);
    CHECK(psiL[iorb] == lrow1[iorb]);

    for (int idim = 0; idim < OHMMS_DIM; idim++)
      CHECK(psiG[iorb][idim] == grow1[iorb][idim]);
  }
  //Test of evaluate_notranspose.
  ValueMatrix phimat, lphimat;
  GradMatrix gphimat;
  phimat.resize(nelec, norb);
  gphimat.resize(nelec, norb);
  lphimat.resize(nelec, norb);

  const int first_index = 0; //Only 2 electrons in this case.
  const int last_index  = 2;
  sposet->evaluate_notranspose(elec, first_index, last_index, phimat, gphimat, lphimat);

  checkMatrix(phimat, spomat);
  checkMatrix(lphimat, laplspomat);

  //Test of makeClone()
  auto sposet_vgl2 = sposet->makeClone();
  phimat           = 0.0;
  gphimat          = 0.0;
  lphimat          = 0.0;

  sposet_vgl2->evaluate_notranspose(elec, first_index, last_index, phimat, gphimat, lphimat);

  checkMatrix(phimat, spomat);
  checkMatrix(lphimat, laplspomat);

  //Lastly, check if name is correct.
  std::string myname = sposet_vgl2->getClassName();
  std::string targetstring("ConstantSPOSet");
  CHECK(myname == targetstring);
}
} // namespace qmcplusplus
