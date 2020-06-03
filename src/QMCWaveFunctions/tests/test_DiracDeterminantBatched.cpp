//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2020 QMCPACK developers.
//
// File developed by: Mark Dewing, mdewing@anl.gov, Argonne National Laboratory
//                    Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////


#include "catch.hpp"

#include "OhmmsData/Libxml2Doc.h"
#include "OhmmsPETE/OhmmsMatrix.h"
#include "QMCWaveFunctions/WaveFunctionComponent.h"
#include "QMCWaveFunctions/Fermion/DiracDeterminantBatched.h"
#include "QMCWaveFunctions/tests/FakeSPO.h"

#ifdef QMC_COMPLEX //This is for the spinor test.
#include "QMCWaveFunctions/ElectronGas/ElectronGasComplexOrbitalBuilder.h"
#endif

#include "QMCWaveFunctions/SpinorSet.h"

#include <stdio.h>
#include <string>

using std::string;


namespace qmcplusplus
{
using RealType     = QMCTraits::RealType;
using ValueType    = QMCTraits::ValueType;
using LogValueType = std::complex<QMCTraits::QTFull::RealType>;
using PsiValueType = QMCTraits::QTFull::ValueType;

#if defined(ENABLE_CUDA) && defined(ENABLE_OFFLOAD)
typedef DiracDeterminantBatched<MatrixUpdateCUDA<ValueType, QMCTraits::QTFull::ValueType>> DetType;
#else
typedef DiracDeterminantBatched<> DetType;
#endif

template<typename T1, typename ALLOC1, typename T2, typename ALLOC2>
void check_matrix(Matrix<T1, ALLOC1>& a, Matrix<T2, ALLOC2>& b)
{
  REQUIRE(a.rows() >= b.rows());
  REQUIRE(a.cols() >= b.cols());
  for (int i = 0; i < b.rows(); i++)
    for (int j = 0; j < b.cols(); j++)
    {
      REQUIRE(a(i, j) == ValueApprox(b(i, j)));
    }
}

TEST_CASE("DiracDeterminantBatched_first", "[wavefunction][fermion]")
{
  FakeSPO* spo = new FakeSPO();
  spo->setOrbitalSetSize(3);
  DetType ddb(spo);

  int norb = 3;
  ddb.set(0, norb);

  // occurs in call to registerData
  ddb.dpsiV.resize(norb);
  ddb.d2psiV.resize(norb);


  ParticleSet elec;

  elec.create(3);
  ddb.recompute(elec);

  Matrix<ValueType> b;
  b.resize(3, 3);

  b(0, 0) = 0.6159749342;
  b(0, 1) = -0.2408954682;
  b(0, 2) = -0.1646081192;
  b(1, 0) = 0.07923894288;
  b(1, 1) = 0.1496231042;
  b(1, 2) = -0.1428117337;
  b(2, 0) = -0.2974298429;
  b(2, 1) = -0.04586322768;
  b(2, 2) = 0.3927890292;

  check_matrix(ddb.psiMinv, b);


  ParticleSet::GradType grad;
  PsiValueType det_ratio  = ddb.ratioGrad(elec, 0, grad);
  PsiValueType det_ratio1 = 0.178276269185;
  REQUIRE(det_ratio1 == ValueApprox(det_ratio));

  ddb.acceptMove(elec, 0);

  b(0, 0) = 3.455170657;
  b(0, 1) = -1.35124809;
  b(0, 2) = -0.9233316353;
  b(1, 0) = 0.05476311768;
  b(1, 1) = 0.1591951095;
  b(1, 2) = -0.1362710138;
  b(2, 0) = -2.235099338;
  b(2, 1) = 0.7119205298;
  b(2, 2) = 0.9105960265;

  check_matrix(ddb.psiMinv, b);
}

//#define DUMP_INFO

TEST_CASE("DiracDeterminantBatched_second", "[wavefunction][fermion]")
{
  FakeSPO* spo = new FakeSPO();
  spo->setOrbitalSetSize(4);
  DetType ddb(spo);

  int norb = 4;
  ddb.set(0, norb);

  // occurs in call to registerData
  ddb.dpsiV.resize(norb);
  ddb.d2psiV.resize(norb);


  ParticleSet elec;

  elec.create(4);
  ddb.recompute(elec);

  Matrix<ValueType> orig_a;
  orig_a.resize(4, 4);
  orig_a = spo->a2;

  for (int i = 0; i < 3; i++)
  {
    for (int j = 0; j < norb; j++)
    {
      orig_a(j, i) = spo->v2(i, j);
    }
  }

  //check_matrix(ddb.psiMinv, b);
  DiracMatrix<ValueType> dm;

  Matrix<ValueType> a_update1, scratchT;
  a_update1.resize(4, 4);
  scratchT.resize(4, 4);
  a_update1 = spo->a2;
  for (int j = 0; j < norb; j++)
  {
    a_update1(j, 0) = spo->v2(0, j);
  }

  Matrix<ValueType> a_update2;
  a_update2.resize(4, 4);
  a_update2 = spo->a2;
  for (int j = 0; j < norb; j++)
  {
    a_update2(j, 0) = spo->v2(0, j);
    a_update2(j, 1) = spo->v2(1, j);
  }

  Matrix<ValueType> a_update3;
  a_update3.resize(4, 4);
  a_update3 = spo->a2;
  for (int j = 0; j < norb; j++)
  {
    a_update3(j, 0) = spo->v2(0, j);
    a_update3(j, 1) = spo->v2(1, j);
    a_update3(j, 2) = spo->v2(2, j);
  }

  ParticleSet::GradType grad;
  PsiValueType det_ratio = ddb.ratioGrad(elec, 0, grad);

  simd::transpose(a_update1.data(), a_update1.rows(), a_update1.cols(), scratchT.data(), scratchT.rows(),
                  scratchT.cols());
  LogValueType det_update1;
  dm.invert_transpose(scratchT, a_update1, det_update1);
  PsiValueType det_ratio1 = LogToValue<ValueType>::convert(det_update1 - ddb.LogValue);
#ifdef DUMP_INFO
  std::cout << "det 0 = " << std::exp(ddb.LogValue) << std::endl;
  std::cout << "det 1 = " << std::exp(det_update1) << std::endl;
  std::cout << "det ratio 1 = " << det_ratio1 << std::endl;
#endif
  //double det_ratio1 = 0.178276269185;

  REQUIRE(det_ratio1 == ValueApprox(det_ratio));

  ddb.acceptMove(elec, 0);

  PsiValueType det_ratio2 = ddb.ratioGrad(elec, 1, grad);
  LogValueType det_update2;
  simd::transpose(a_update2.data(), a_update2.rows(), a_update2.cols(), scratchT.data(), scratchT.rows(),
                  scratchT.cols());
  dm.invert_transpose(scratchT, a_update2, det_update2);
  PsiValueType det_ratio2_val = LogToValue<ValueType>::convert(det_update2 - det_update1);
#ifdef DUMP_INFO
  std::cout << "det 1 = " << std::exp(ddb.LogValue) << std::endl;
  std::cout << "det 2 = " << std::exp(det_update2) << std::endl;
  std::cout << "det ratio 2 = " << det_ratio2 << std::endl;
#endif
  //double det_ratio2_val = 0.178276269185;
  REQUIRE(det_ratio2 == ValueApprox(det_ratio2_val));

  ddb.acceptMove(elec, 1);

  PsiValueType det_ratio3 = ddb.ratioGrad(elec, 2, grad);
  LogValueType det_update3;
  simd::transpose(a_update3.data(), a_update3.rows(), a_update3.cols(), scratchT.data(), scratchT.rows(),
                  scratchT.cols());
  dm.invert_transpose(scratchT, a_update3, det_update3);
  PsiValueType det_ratio3_val = LogToValue<ValueType>::convert(det_update3 - det_update2);
#ifdef DUMP_INFO
  std::cout << "det 2 = " << std::exp(ddb.LogValue) << std::endl;
  std::cout << "det 3 = " << std::exp(det_update3) << std::endl;
  std::cout << "det ratio 3 = " << det_ratio3 << std::endl;
#endif
  REQUIRE(det_ratio3 == ValueApprox(det_ratio3_val));
  //check_value(det_ratio3, det_ratio3_val);

  ddb.acceptMove(elec, 2);

  simd::transpose(orig_a.data(), orig_a.rows(), orig_a.cols(), scratchT.data(), scratchT.rows(), scratchT.cols());
  dm.invert_transpose(scratchT, orig_a, det_update3);

#ifdef DUMP_INFO
  std::cout << "original " << std::endl;
  std::cout << orig_a << std::endl;
  std::cout << "block update " << std::endl;
  std::cout << ddb.psiMinv << std::endl;
#endif

  check_matrix(ddb.psiMinv, orig_a);
}

TEST_CASE("DiracDeterminantBatched_delayed_update", "[wavefunction][fermion]")
{
  FakeSPO* spo = new FakeSPO();
  spo->setOrbitalSetSize(4);
  DetType ddc(spo);

  int norb = 4;
  // maximum delay 2
  ddc.set(0, norb, 2);

  // occurs in call to registerData
  ddc.dpsiV.resize(norb);
  ddc.d2psiV.resize(norb);


  ParticleSet elec;

  elec.create(4);
  ddc.recompute(elec);

  Matrix<ValueType> orig_a;
  orig_a.resize(4, 4);
  orig_a = spo->a2;

  for (int i = 0; i < 3; i++)
  {
    for (int j = 0; j < norb; j++)
    {
      orig_a(j, i) = spo->v2(i, j);
    }
  }

  //check_matrix(ddc.psiMinv, b);
  DiracMatrix<ValueType> dm;

  Matrix<ValueType> a_update1, scratchT;
  scratchT.resize(4, 4);
  a_update1.resize(4, 4);
  a_update1 = spo->a2;
  for (int j = 0; j < norb; j++)
  {
    a_update1(j, 0) = spo->v2(0, j);
  }

  Matrix<ValueType> a_update2;
  a_update2.resize(4, 4);
  a_update2 = spo->a2;
  for (int j = 0; j < norb; j++)
  {
    a_update2(j, 0) = spo->v2(0, j);
    a_update2(j, 1) = spo->v2(1, j);
  }

  Matrix<ValueType> a_update3;
  a_update3.resize(4, 4);
  a_update3 = spo->a2;
  for (int j = 0; j < norb; j++)
  {
    a_update3(j, 0) = spo->v2(0, j);
    a_update3(j, 1) = spo->v2(1, j);
    a_update3(j, 2) = spo->v2(2, j);
  }


  ParticleSet::GradType grad;
  PsiValueType det_ratio = ddc.ratioGrad(elec, 0, grad);

  simd::transpose(a_update1.data(), a_update1.rows(), a_update1.cols(), scratchT.data(), scratchT.rows(),
                  scratchT.cols());
  LogValueType det_update1;
  dm.invert_transpose(scratchT, a_update1, det_update1);
  PsiValueType det_ratio1 = LogToValue<ValueType>::convert(det_update1 - ddc.LogValue);
#ifdef DUMP_INFO
  std::cout << "det 0 = " << std::exp(ddc.LogValue) << std::endl;
  std::cout << "det 1 = " << std::exp(det_update1) << std::endl;
  std::cout << "det ratio 1 = " << det_ratio1 << std::endl;
#endif
  //double det_ratio1 = 0.178276269185;

  REQUIRE(det_ratio1 == ValueApprox(det_ratio));

  // update of Ainv in ddc is delayed
  ddc.acceptMove(elec, 0, true);
  // force update Ainv in ddc using SM-1 code path
  ddc.completeUpdates();

  check_matrix(ddc.psiMinv, a_update1);

  grad                    = ddc.evalGrad(elec, 1);
  PsiValueType det_ratio2 = ddc.ratioGrad(elec, 1, grad);
  simd::transpose(a_update2.data(), a_update2.rows(), a_update2.cols(), scratchT.data(), scratchT.rows(),
                  scratchT.cols());
  LogValueType det_update2;
  dm.invert_transpose(scratchT, a_update2, det_update2);
  PsiValueType det_ratio2_val = LogToValue<ValueType>::convert(det_update2 - det_update1);
#ifdef DUMP_INFO
  std::cout << "det 1 = " << std::exp(ddc.LogValue) << std::endl;
  std::cout << "det 2 = " << std::exp(det_update2) << std::endl;
  std::cout << "det ratio 2 = " << det_ratio2 << std::endl;
#endif
  // check ratio computed directly and the one computed by ddc with no delay
  //double det_ratio2_val = 0.178276269185;
  REQUIRE(det_ratio2 == ValueApprox(det_ratio2_val));

  // update of Ainv in ddc is delayed
  ddc.acceptMove(elec, 1, true);

  grad                    = ddc.evalGrad(elec, 2);
  PsiValueType det_ratio3 = ddc.ratioGrad(elec, 2, grad);
  simd::transpose(a_update3.data(), a_update3.rows(), a_update3.cols(), scratchT.data(), scratchT.rows(),
                  scratchT.cols());
  LogValueType det_update3;
  dm.invert_transpose(scratchT, a_update3, det_update3);
  PsiValueType det_ratio3_val = LogToValue<ValueType>::convert(det_update3 - det_update2);
#ifdef DUMP_INFO
  std::cout << "det 2 = " << std::exp(ddc.LogValue) << std::endl;
  std::cout << "det 3 = " << std::exp(det_update3) << std::endl;
  std::cout << "det ratio 3 = " << det_ratio3 << std::endl;
#endif
  // check ratio computed directly and the one computed by ddc with 1 delay
  REQUIRE(det_ratio3 == ValueApprox(det_ratio3_val));
  //check_value(det_ratio3, det_ratio3_val);

  // maximal delay reached and Ainv is updated fully
  ddc.acceptMove(elec, 2, true);
  ddc.completeUpdates();

  // fresh invert orig_a
  simd::transpose(orig_a.data(), orig_a.rows(), orig_a.cols(), scratchT.data(), scratchT.rows(), scratchT.cols());
  dm.invert_transpose(scratchT, orig_a, det_update3);

#ifdef DUMP_INFO
  std::cout << "original " << std::endl;
  std::cout << orig_a << std::endl;
  std::cout << "delayed update " << std::endl;
  std::cout << ddc.psiMinv << std::endl;
#endif

  // compare all the elements of psiMinv in ddc and orig_a
  check_matrix(ddc.psiMinv, orig_a);
}

} // namespace qmcplusplus
