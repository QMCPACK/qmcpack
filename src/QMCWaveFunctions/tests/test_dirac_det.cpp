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
#include "QMCWaveFunctions/ElectronGas/ElectronGasComplexOrbitalBuilder.h"

#include "QMCWaveFunctions/SpinorSet.h"

#include <stdio.h>
#include <string>

using std::string;


namespace qmcplusplus
{

using RealType = QMCTraits::RealType;
using ValueType = QMCTraits::ValueType;
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

  virtual void evaluate(const ParticleSet& P, int iat, ValueVector_t& psi);

  virtual void evaluate(const ParticleSet& P, int iat, ValueVector_t& psi, GradVector_t& dpsi, ValueVector_t& d2psi);

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

void FakeSPO::evaluate(const ParticleSet& P, int iat, ValueVector_t& psi)
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

void FakeSPO::evaluate(const ParticleSet& P, int iat, ValueVector_t& psi, GradVector_t& dpsi, ValueVector_t& d2psi)
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
  ValueType det_ratio = ddb.ratioGrad(elec, 0, grad);
  ValueType det_ratio1 = 0.178276269185;
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
  ValueType det_ratio = ddb.ratioGrad(elec, 0, grad);

  simd::transpose(a_update1.data(), a_update1.rows(), a_update1.cols(), scratchT.data(), scratchT.rows(), scratchT.cols());
  LogValueType det_update1;
  dm.invert_transpose(scratchT, a_update1, det_update1);
  ValueType det_ratio1  = LogToValue<ValueType>::convert(det_update1 - ddb.LogValue);
#ifdef DUMP_INFO
  std::cout << "det 0 = " << std::exp(ddb.LogValue) << std::endl;
  std::cout << "det 1 = " << std::exp(det_update1) << std::endl;
  std::cout << "det ratio 1 = " << det_ratio1 << std::endl;
#endif
  //double det_ratio1 = 0.178276269185;

  REQUIRE(det_ratio1 == ValueApprox(det_ratio));

  ddb.acceptMove(elec, 0);


  ValueType det_ratio2 = ddb.ratioGrad(elec, 1, grad);
  LogValueType det_update2;
  simd::transpose(a_update2.data(), a_update2.rows(), a_update2.cols(), scratchT.data(), scratchT.rows(), scratchT.cols());
  dm.invert_transpose(scratchT, a_update2, det_update2);
  ValueType det_ratio2_val = LogToValue<ValueType>::convert(det_update2 - det_update1);
#ifdef DUMP_INFO
  std::cout << "det 1 = " << std::exp(ddb.LogValue) << std::endl;
  std::cout << "det 2 = " << std::exp(det_update2) << std::endl;
  std::cout << "det ratio 2 = " << det_ratio2 << std::endl;
#endif
  //double det_ratio2_val = 0.178276269185;
  REQUIRE(det_ratio2 == ValueApprox(det_ratio2_val));

  ddb.acceptMove(elec, 1);

  ValueType det_ratio3 = ddb.ratioGrad(elec, 2, grad);
  LogValueType det_update3;
  simd::transpose(a_update3.data(), a_update3.rows(), a_update3.cols(), scratchT.data(), scratchT.rows(), scratchT.cols());
  dm.invert_transpose(scratchT, a_update3, det_update3);
  ValueType det_ratio3_val = LogToValue<ValueType>::convert(det_update3 - det_update2);
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
  ValueType det_ratio = ddc.ratioGrad(elec, 0, grad);

  simd::transpose(a_update1.data(), a_update1.rows(), a_update1.cols(), scratchT.data(), scratchT.rows(), scratchT.cols());
  LogValueType det_update1;
  dm.invert_transpose(scratchT, a_update1, det_update1);
  ValueType det_ratio1  = LogToValue<ValueType>::convert(det_update1 - ddc.LogValue);
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

  grad                 = ddc.evalGrad(elec, 1);
  ValueType det_ratio2 = ddc.ratioGrad(elec, 1, grad);
  simd::transpose(a_update2.data(), a_update2.rows(), a_update2.cols(), scratchT.data(), scratchT.rows(), scratchT.cols());
  LogValueType det_update2;
  dm.invert_transpose(scratchT, a_update2, det_update2);
  ValueType det_ratio2_val = LogToValue<ValueType>::convert(det_update2 - det_update1);
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

  grad                 = ddc.evalGrad(elec, 2);
  ValueType det_ratio3 = ddc.ratioGrad(elec, 2, grad);
  simd::transpose(a_update3.data(), a_update3.rows(), a_update3.cols(), scratchT.data(), scratchT.rows(), scratchT.cols());
  LogValueType det_update3;
  dm.invert_transpose(scratchT, a_update3, det_update3);
  ValueType det_ratio3_val = LogToValue<ValueType>::convert(det_update3 - det_update2);
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
  std::cout << " i!!!!!!!!!!!!!WEEEEE OH MY GOD WE'RE GETTING TESTED!!!!!!!!!!!!!\n";
  REQUIRE( 1==1);
  typedef QMCTraits::PosType PosType;
  typedef ParticleSet::ParticlePos_t ParticlePos_t;
  
  //Shamelessly stealing this from test_einset.cpp.  3 particles though.
  ParticleSet ions_;
  ParticleSet elec_;
  ions_.setName("ion");
  ions_.create(2);

  ions_.R[0][0] = 0.00000000 ; 
  ions_.R[0][1] = 0.00000000 ;    
  ions_.R[0][2] = 1.08659253 ;  
  ions_.R[1][0] = 0.00000000 ;   
  ions_.R[1][1] = 0.00000000 ;  
  ions_.R[1][2] =-1.08659253 ;

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

  elec_.spins[0] = 0.0 ;
  elec_.spins[1] = 0.2 ;
  elec_.spins[2] = 0.4 ;

  // O2 test example from pwscf non-collinear calculation.
  elec_.Lattice.R(0, 0) =  5.10509515 ;
  elec_.Lattice.R(0, 1) = -3.23993545 ;
  elec_.Lattice.R(0, 2) =  0.00000000 ;
  elec_.Lattice.R(1, 0) = 5.10509515 ;
  elec_.Lattice.R(1, 1) = 3.23993545 ;
  elec_.Lattice.R(1, 2) = 0.00000000 ;
  elec_.Lattice.R(2, 0) = -6.49690625 ;
  elec_.Lattice.R(2, 1) =  0.00000000 ;
  elec_.Lattice.R(2, 2) =  7.08268015 ; 

  SpeciesSet& tspecies         = elec_.getSpeciesSet();
  int upIdx                    = tspecies.addSpecies("u");
  int chargeIdx                = tspecies.addAttribute("charge");
  tspecies(chargeIdx, upIdx)   = -1;

#ifdef ENABLE_SOA
  elec_.addTable(ions_, DT_SOA);
#else
  elec_.addTable(ions_, DT_AOS);
#endif
  elec_.resetGroups();
  elec_.update();
  // </steal>

  //Our test case is going to be three electron gas orbitals distinguished by 3 different kpoints.
  //Independent SPO's for the up and down channels.
  //
  std::vector<PosType> kup,kdn;
  std::vector<RealType> k2up,k2dn;  

 
  kup.resize(3);
  kup[0]=PosType(0,0,0);
  kup[1]=PosType(0.1,0.2,0.3);
  kup[2]=PosType(0.4,0.5,0.6);
  
  k2up.resize(3);
  k2up[0]=dot(kup[0],kup[0]);
  k2up[1]=dot(kup[1],kup[1]);
  k2up[2]=dot(kup[2],kup[2]);
 
  kdn.resize(3);
  kdn[0]=PosType(0,0,0);
  kdn[1]=PosType(-0.1,0.2,-0.3);
  kdn[2]=PosType(0.4,-0.5,0.6);

  k2dn.resize(3);
  k2dn[0]=dot(kdn[0],kdn[0]);
  k2dn[1]=dot(kdn[1],kdn[1]);
  k2dn[2]=dot(kdn[2],kdn[2]);

  std::shared_ptr<EGOSet> spo_up(new EGOSet(kup,k2up));
  std::shared_ptr<EGOSet> spo_dn(new EGOSet(kdn,k2dn));

  SpinorSet* spinor_set = new SpinorSet();
  spinor_set->set_spos(spo_up,spo_dn);

  DetType dd(spinor_set);


  //So we need a ratio test.
   //make a trial spin move.  

  //a ratioGrad test.  
  
  //an evaluateLog with gradients and laplacians.  
  


}
#endif

} // namespace qmcplusplus
