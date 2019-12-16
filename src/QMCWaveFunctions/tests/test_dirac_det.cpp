//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2017 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Mark Dewing, mdewing@anl.gov, Argonne National Laboratory
//
// File created by: Mark Dewing, mdewing@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////


#include "catch.hpp"

#include "OhmmsData/Libxml2Doc.h"
#include "OhmmsPETE/OhmmsMatrix.h"
#include "QMCWaveFunctions/WaveFunctionComponent.h"
#include "QMCWaveFunctions/SPOSet.h"
#include "QMCWaveFunctions/Fermion/DiracDeterminant.h"

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

#ifdef ENABLE_CUDA
typedef DiracDeterminant<DelayedUpdateCUDA<ValueType, QMCTraits::QTFull::ValueType>> DetType;
#else
typedef DiracDeterminant<> DetType;
#endif

template<typename T1, typename T2>
void check_matrix(Matrix<T1>& a, Matrix<T2>& b)
{
  REQUIRE(a.size() == b.size());
  for (int i = 0; i < a.rows(); i++)
  {
    for (int j = 0; j < a.cols(); j++)
    {
      REQUIRE(a(i, j) == ValueApprox(b(i, j)));
    }
  }
}

class FakeSPO : public SPOSet
{
public:
  Matrix<ValueType> a;
  Matrix<ValueType> a2;
  Vector<ValueType> v;
  Matrix<ValueType> v2;

  FakeSPO();
  virtual ~FakeSPO() {}

  virtual void report() {}
  virtual void resetParameters(const opt_variables_type& optVariables) {}
  virtual void resetTargetParticleSet(ParticleSet& P) {}
  virtual void setOrbitalSetSize(int norbs);

  virtual void evaluateValue(const ParticleSet& P, int iat, ValueVector_t& psi);

  virtual void evaluateVGL(const ParticleSet& P, int iat, ValueVector_t& psi, GradVector_t& dpsi, ValueVector_t& d2psi);

  virtual void evaluate_notranspose(const ParticleSet& P,
                                    int first,
                                    int last,
                                    ValueMatrix_t& logdet,
                                    GradMatrix_t& dlogdet,
                                    ValueMatrix_t& d2logdet);
};

FakeSPO::FakeSPO()
{
  className = "FakeSPO";
  a.resize(3, 3);

  a(0, 0) = 2.3;
  a(0, 1) = 4.5;
  a(0, 2) = 2.6;
  a(1, 0) = 0.5;
  a(1, 1) = 8.5;
  a(1, 2) = 3.3;
  a(2, 0) = 1.8;
  a(2, 1) = 4.4;
  a(2, 2) = 4.9;

  v.resize(3);
  v[0] = 1.9;
  v[1] = 2.0;
  v[2] = 3.1;


  a2.resize(4, 4);
  a2(0, 0) = 2.3;
  a2(0, 1) = 4.5;
  a2(0, 2) = 2.6;
  a2(0, 3) = 1.2;
  a2(1, 0) = 0.5;
  a2(1, 1) = 8.5;
  a2(1, 2) = 3.3;
  a2(1, 3) = 0.3;
  a2(2, 0) = 1.8;
  a2(2, 1) = 4.4;
  a2(2, 2) = 4.9;
  a2(2, 3) = 2.8;
  a2(3, 0) = 0.8;
  a2(3, 1) = 4.1;
  a2(3, 2) = 3.2;
  a2(3, 3) = 1.1;

  v2.resize(3, 4);

  v2(0, 0) = 3.2;
  v2(0, 1) = 0.5;
  v2(0, 2) = 5.9;
  v2(0, 3) = 3.7;
  v2(1, 0) = 0.3;
  v2(1, 1) = 1.4;
  v2(1, 2) = 3.9;
  v2(1, 3) = 8.2;
  v2(2, 0) = 3.3;
  v2(2, 1) = 5.4;
  v2(2, 2) = 4.9;
  v2(2, 3) = 2.2;
}

void FakeSPO::setOrbitalSetSize(int norbs) { OrbitalSetSize = norbs; }

void FakeSPO::evaluateValue(const ParticleSet& P, int iat, ValueVector_t& psi)
{
  if (OrbitalSetSize == 3)
  {
    for (int i = 0; i < 3; i++)
    {
      psi[i] = a(iat, i);
    }
  }
  else if (OrbitalSetSize == 4)
  {
    for (int i = 0; i < 4; i++)
    {
      psi[i] = a2(iat, i);
    }
  }
}

void FakeSPO::evaluateVGL(const ParticleSet& P, int iat, ValueVector_t& psi, GradVector_t& dpsi, ValueVector_t& d2psi)
{
  if (OrbitalSetSize == 3)
  {
    for (int i = 0; i < 3; i++)
    {
      psi[i] = v[i];
    }
  }
  else if (OrbitalSetSize == 4)
  {
    for (int i = 0; i < 4; i++)
    {
      psi[i] = v2(iat, i);
    }
  }
}

void FakeSPO::evaluate_notranspose(const ParticleSet& P,
                                   int first,
                                   int last,
                                   ValueMatrix_t& logdet,
                                   GradMatrix_t& dlogdet,
                                   ValueMatrix_t& d2logdet)
{
  if (OrbitalSetSize == 3)
  {
    for (int i = 0; i < 3; i++)
    {
      for (int j = 0; j < 3; j++)
      {
        logdet(j, i) = a(i, j);
      }
    }
  }
  else if (OrbitalSetSize == 4)
  {
    for (int i = 0; i < 4; i++)
    {
      for (int j = 0; j < 4; j++)
      {
        logdet(j, i) = a2(i, j);
      }
    }
  }
}

TEST_CASE("DiracDeterminant_first", "[wavefunction][fermion]")
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

  check_matrix(ddb.psiM, b);


  ParticleSet::GradType grad;
  PsiValueType det_ratio  = ddb.ratioGrad(elec, 0, grad);
  PsiValueType det_ratio1 = 0.178276269185;
  REQUIRE(det_ratio1 == ValueApprox(det_ratio));

  ddb.acceptMove(elec, 0);
  ddb.completeUpdates();

  b(0, 0) = 3.455170657;
  b(0, 1) = -1.35124809;
  b(0, 2) = -0.9233316353;
  b(1, 0) = 0.05476311768;
  b(1, 1) = 0.1591951095;
  b(1, 2) = -0.1362710138;
  b(2, 0) = -2.235099338;
  b(2, 1) = 0.7119205298;
  b(2, 2) = 0.9105960265;

  check_matrix(ddb.psiM, b);
}

//#define DUMP_INFO

TEST_CASE("DiracDeterminant_second", "[wavefunction][fermion]")
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

  //check_matrix(ddb.psiM, b);
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
  ddb.completeUpdates();

  simd::transpose(orig_a.data(), orig_a.rows(), orig_a.cols(), scratchT.data(), scratchT.rows(), scratchT.cols());
  dm.invert_transpose(scratchT, orig_a, det_update3);

#ifdef DUMP_INFO
  std::cout << "original " << std::endl;
  std::cout << orig_a << std::endl;
  std::cout << "block update " << std::endl;
  std::cout << ddb.psiM << std::endl;
#endif

  check_matrix(orig_a, ddb.psiM);
}

TEST_CASE("DiracDeterminant_delayed_update", "[wavefunction][fermion]")
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

  //check_matrix(ddc.psiM, b);
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
  ddc.acceptMove(elec, 0);
  // force update Ainv in ddc using SM-1 code path
  ddc.completeUpdates();

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
  ddc.acceptMove(elec, 1);

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
  ddc.acceptMove(elec, 2);
  ddc.completeUpdates();

  // fresh invert orig_a
  simd::transpose(orig_a.data(), orig_a.rows(), orig_a.cols(), scratchT.data(), scratchT.rows(), scratchT.cols());
  dm.invert_transpose(scratchT, orig_a, det_update3);

#ifdef DUMP_INFO
  std::cout << "original " << std::endl;
  std::cout << orig_a << std::endl;
  std::cout << "delayed update " << std::endl;
  std::cout << ddc.psiM << std::endl;
#endif

  // compare all the elements of psiM in ddc and orig_a
  check_matrix(orig_a, ddc.psiM);
}

#ifdef QMC_COMPLEX
TEST_CASE("DiracDeterminant_spinor_update", "[wavefunction][fermion]")
{
  typedef QMCTraits::ValueType ValueType;
  typedef QMCTraits::PosType PosType;
  typedef QMCTraits::GradType GradType;
  typedef WaveFunctionComponent::LogValueType LogValueType;
  typedef ParticleSet::ParticlePos_t ParticlePos_t;
  typedef ParticleSet::ParticleGradient_t ParticleGradient_t;
  typedef ParticleSet::ParticleLaplacian_t ParticleLaplacian_t;
  //Shamelessly stealing this from test_einset.cpp.  3 particles though.
  ParticleSet ions_;
  ParticleSet elec_;
  ions_.setName("ion");
  ions_.create(2);

  ions_.R[0][0] = 0.00000000;
  ions_.R[0][1] = 0.00000000;
  ions_.R[0][2] = 1.08659253;
  ions_.R[1][0] = 0.00000000;
  ions_.R[1][1] = 0.00000000;
  ions_.R[1][2] = -1.08659253;

  elec_.setName("elec");
  elec_.create(3);
  elec_.R[0][0] = 0.1;
  elec_.R[0][1] = -0.3;
  elec_.R[0][2] = 1.0;
  elec_.R[1][0] = -0.1;
  elec_.R[1][1] = 0.3;
  elec_.R[1][2] = 1.0;
  elec_.R[2][0] = 0.1;
  elec_.R[2][1] = 0.2;
  elec_.R[2][2] = 0.3;

  elec_.spins[0] = 0.0;
  elec_.spins[1] = 0.2;
  elec_.spins[2] = 0.4;

  // O2 test example from pwscf non-collinear calculation.
  elec_.Lattice.R(0, 0) = 5.10509515;
  elec_.Lattice.R(0, 1) = -3.23993545;
  elec_.Lattice.R(0, 2) = 0.00000000;
  elec_.Lattice.R(1, 0) = 5.10509515;
  elec_.Lattice.R(1, 1) = 3.23993545;
  elec_.Lattice.R(1, 2) = 0.00000000;
  elec_.Lattice.R(2, 0) = -6.49690625;
  elec_.Lattice.R(2, 1) = 0.00000000;
  elec_.Lattice.R(2, 2) = 7.08268015;

  SpeciesSet& tspecies       = elec_.getSpeciesSet();
  int upIdx                  = tspecies.addSpecies("u");
  int chargeIdx              = tspecies.addAttribute("charge");
  tspecies(chargeIdx, upIdx) = -1;

#ifdef ENABLE_SOA
  elec_.addTable(ions_, DT_SOA);
#else
  elec_.addTable(ions_, DT_AOS);
#endif
  elec_.resetGroups();
  elec_.update();
  // </steal>


  QMCTraits::IndexType nelec = elec_.R.size();
  QMCTraits::IndexType norb  = 3;
  //Our test case is going to be three electron gas orbitals distinguished by 3 different kpoints.
  //Independent SPO's for the up and down channels.
  //
  std::vector<PosType> kup, kdn;
  std::vector<RealType> k2up, k2dn;


  kup.resize(nelec);
  kup[0] = PosType(0, 0, 0);
  kup[1] = PosType(0.1, 0.2, 0.3);
  kup[2] = PosType(0.4, 0.5, 0.6);

  k2up.resize(nelec);
  //For some goofy reason, EGOSet needs to be initialized with:
  //1.) A k-vector list (fine).
  //2.) A list of -|k|^2.  To save on expensive - sign multiplication apparently.
  k2up[0] = -dot(kup[0], kup[0]);
  k2up[1] = -dot(kup[1], kup[1]);
  k2up[2] = -dot(kup[2], kup[2]);

  kdn.resize(nelec);
  kdn[0] = PosType(0, 0, 0);
  kdn[1] = PosType(-0.1, 0.2, -0.3);
  kdn[2] = PosType(0.4, -0.5, 0.6);

  k2dn.resize(nelec);
  k2dn[0] = -dot(kdn[0], kdn[0]);
  k2dn[1] = -dot(kdn[1], kdn[1]);
  k2dn[2] = -dot(kdn[2], kdn[2]);

  std::shared_ptr<EGOSet> spo_up(new EGOSet(kup, k2up));
  std::shared_ptr<EGOSet> spo_dn(new EGOSet(kdn, k2dn));

  SpinorSet* spinor_set = new SpinorSet();
  spinor_set->set_spos(spo_up, spo_dn);

  DetType dd(spinor_set);
  dd.resize(nelec, norb);
  app_log() << " nelec=" << nelec << " norb=" << norb << std::endl;

  ParticleGradient_t G;
  ParticleLaplacian_t L;
  //This is a vector of ValueType, so we're using it..."
  ParticleLaplacian_t SG;

  G.resize(nelec);
  L.resize(nelec);
  SG.resize(nelec);

  G  = 0.0;
  L  = 0.0;
  SG = 0.0;

  PosType dr(0.1, -0.05, 0.2);
  RealType ds = 0.3;

  app_log() << " BEFORE\n";
  app_log() << " R = " << elec_.R << std::endl;
  app_log() << " s = " << elec_.spins << std::endl;

  //In this section, we're going to test that values and various derivatives come out
  //correctly at the reference configuration.

  LogValueType logref = dd.evaluateLog(elec_, G, L);

  REQUIRE(logref == ComplexApprox(ValueType(-1.1619939279564413, 0.8794794652468605)));
  REQUIRE(G[0][0] == ComplexApprox(ValueType(0.13416635, 0.2468612)));
  REQUIRE(G[0][1] == ComplexApprox(ValueType(-1.1165475, 0.71497753)));
  REQUIRE(G[0][2] == ComplexApprox(ValueType(0.0178403, 0.08212244)));
  REQUIRE(G[1][0] == ComplexApprox(ValueType(1.00240841, 0.12371593)));
  REQUIRE(G[1][1] == ComplexApprox(ValueType(1.62679698, -0.41080777)));
  REQUIRE(G[1][2] == ComplexApprox(ValueType(1.81324632, 0.78589013)));
  REQUIRE(G[2][0] == ComplexApprox(ValueType(-1.10994555, 0.15525902)));
  REQUIRE(G[2][1] == ComplexApprox(ValueType(-0.46335602, -0.50809713)));
  REQUIRE(G[2][2] == ComplexApprox(ValueType(-1.751199, 0.10949589)));
  REQUIRE(L[0] == ComplexApprox(ValueType(-2.06554158, 1.18145239)));
  REQUIRE(L[1] == ComplexApprox(ValueType(-5.06340536, 0.82126749)));
  REQUIRE(L[2] == ComplexApprox(ValueType(-4.82375261, -1.97943258)));

  //This is a workaround for the fact that I haven't implemented
  // evaluateLogWithSpin().  Shouldn't be needed unless we do drifted all-electron moves...
  for (int iat = 0; iat < nelec; iat++)
    dd.evalGradWithSpin(elec_, iat, SG[iat]);

  REQUIRE(SG[0] == ComplexApprox(ValueType(-1.05686704, -2.01802154)));
  REQUIRE(SG[1] == ComplexApprox(ValueType(1.18922259, 2.80414598)));
  REQUIRE(SG[2] == ComplexApprox(ValueType(-0.62617675, -0.51093984)));

  GradType g_singleeval(0.0);
  g_singleeval = dd.evalGrad(elec_, 1);

  REQUIRE(g_singleeval[0] == ComplexApprox(G[1][0]));
  REQUIRE(g_singleeval[1] == ComplexApprox(G[1][1]));
  REQUIRE(g_singleeval[2] == ComplexApprox(G[1][2]));


  //And now we're going to propose a trial spin+particle move and check the ratio and gradients at the
  //new location.
  //
  elec_.makeMoveAndCheckWithSpin(1, dr, ds);

  ValueType ratio_new;
  LogValueType spingrad_new;
  GradType grad_new;

  //This tests ratio only evaluation.  Indirectly a call to evaluate(P,iat)
  ratio_new = dd.ratio(elec_, 1);
  REQUIRE(ratio_new == ComplexApprox(ValueType(1.7472917722050971, 1.1900872950904169)));

  ratio_new = dd.ratioGrad(elec_, 1, grad_new);
  REQUIRE(ratio_new == ComplexApprox(ValueType(1.7472917722050971, 1.1900872950904169)));
  REQUIRE(grad_new[0] == ComplexApprox(ValueType(0.5496675534224996, -0.07968022499097227)));
  REQUIRE(grad_new[1] == ComplexApprox(ValueType(0.4927399293808675, -0.29971549854643653)));
  REQUIRE(grad_new[2] == ComplexApprox(ValueType(1.2792642963632226, 0.12110307514989149)));

  grad_new     = 0;
  spingrad_new = 0;
  ratio_new    = dd.ratioGradWithSpin(elec_, 1, grad_new, spingrad_new);
  REQUIRE(ratio_new == ComplexApprox(ValueType(1.7472917722050971, 1.1900872950904169)));
  REQUIRE(grad_new[0] == ComplexApprox(ValueType(0.5496675534224996, -0.07968022499097227)));
  REQUIRE(grad_new[1] == ComplexApprox(ValueType(0.4927399293808675, -0.29971549854643653)));
  REQUIRE(grad_new[2] == ComplexApprox(ValueType(1.2792642963632226, 0.12110307514989149)));
  REQUIRE(spingrad_new == ComplexApprox(ValueType(1.164708841479661, 0.9576425115390172)));


  //Cool.  Now we test the transition between rejecting a move and accepting a move.
  //Reject the move first.  We want to see if everything stays the same.  evalGrad and evalSpinGrad for ease of use.

  elec_.rejectMove(1);
  //Going to check evalGrad and evalGradWithSpin for simplicity.
  g_singleeval = dd.evalGrad(elec_, 1);
  REQUIRE(g_singleeval[0] == ComplexApprox(G[1][0]));
  REQUIRE(g_singleeval[1] == ComplexApprox(G[1][1]));
  REQUIRE(g_singleeval[2] == ComplexApprox(G[1][2]));

  LogValueType spingrad_old_test;
  g_singleeval = dd.evalGradWithSpin(elec_, 1, spingrad_old_test);

  REQUIRE(spingrad_old_test == ComplexApprox(SG[1]));
  REQUIRE(g_singleeval[0] == ComplexApprox(G[1][0]));
  REQUIRE(g_singleeval[1] == ComplexApprox(G[1][1]));
  REQUIRE(g_singleeval[2] == ComplexApprox(G[1][2]));

  //Now we test what happens if we accept a move...
  elec_.makeMoveAndCheckWithSpin(1, dr, ds);
  elec_.acceptMove(1);

  LogValueType lognew(0.0);
  G      = 0.0; //evalauteLog += onto the G and L arguments.  So we zero them out.
  L      = 0.0;
  SG     = 0.0;
  lognew = dd.evaluateLog(elec_, G, L);

  for (int iat = 0; iat < nelec; iat++)
    dd.evalGradWithSpin(elec_, iat, SG[iat]);
  //logval for the new configuration has been computed with python.
  //The others reference values are computed earlier in this section.  New values equal the previous
  // "new values" associated with the previous trial moves.
  REQUIRE(lognew == ComplexApprox(ValueType(-0.41337396772929913, 1.4774106123071726)));
  REQUIRE(G[1][0] == ComplexApprox(grad_new[0]));
  REQUIRE(G[1][1] == ComplexApprox(grad_new[1]));
  REQUIRE(G[1][2] == ComplexApprox(grad_new[2]));
  REQUIRE(SG[1] == ComplexApprox(spingrad_new));
}
#endif

} // namespace qmcplusplus
