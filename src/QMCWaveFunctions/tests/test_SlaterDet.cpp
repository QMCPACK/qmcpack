//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2024 QMCPACK developers.
//
// File developed by: Cody A. Melton, cmelton@sandia.gov, Sandia National Laboratories
//
// File created by: Cody A. Melton, cmelton@sandia.gov, Sandia National Laboratories
//////////////////////////////////////////////////////////////////////////////////////

#include "catch.hpp"

#include "QMCWaveFunctions/Fermion/DiracDeterminantBase.h"
#include "QMCWaveFunctions/Fermion/SlaterDet.h"
#include "QMCWaveFunctions/tests/ConstantSPOSet.h"
#include "Particle/tests/MinimalParticlePool.h"
#include <ResourceCollection.h>
#include <algorithm>

namespace qmcplusplus
{

class DummyDiracDetWithoutMW : public DiracDeterminantBase
{
public:
  DummyDiracDetWithoutMW(const std::string& class_name, std::unique_ptr<SPOSet>&& spos, int first, int last)
      : DiracDeterminantBase(getClassName(), std::move(spos), first, last)
  {}
  std::string getClassName() const override { return "DummyDiracDetWithoutMW"; }
  LogValue evaluateLog(const ParticleSet& P,
                       ParticleSet::ParticleGradient& G,
                       ParticleSet::ParticleLaplacian& L) override
  {
    G = 0.0;
    L = 0.0;
    return 0.0;
  }
  void acceptMove(ParticleSet& P, int iat, bool safe_to_delay = false) override {}
  void restore(int iat) override {}
  PsiValue ratio(ParticleSet& P, int iat) override { return 1.0; }
  GradType evalGrad(ParticleSet& P, int iat) override
  {
    GradType grad;
    grad[0] = 123.;
    grad[1] = 456.;
    grad[2] = 789.;
    return grad;
  }
  GradType evalGradWithSpin(ParticleSet& P, int iat, ComplexType& spingrad) override
  {
    GradType grad;
    grad[0]  = 0.123;
    grad[1]  = 0.456;
    grad[2]  = 0.789;
    spingrad = ComplexType(0.1, 0.2);
    return grad;
  }
  PsiValue ratioGrad(ParticleSet& P, int iat, GradType& grad_iat) override
  {
    grad_iat[0] = 123.;
    grad_iat[1] = 456.;
    grad_iat[2] = 789.;
    return 1;
  }
  PsiValue ratioGradWithSpin(ParticleSet& P, int iat, GradType& grad_iat, ComplexType& spingrad_iat) override
  {
    grad_iat[0]  = 0.123;
    grad_iat[1]  = 0.456;
    grad_iat[2]  = 0.789;
    spingrad_iat = ComplexType(0.1, 0.2);
    return 1;
  }
  void registerData(ParticleSet& P, WFBufferType& buf) override {}
  LogValue updateBuffer(ParticleSet& P, WFBufferType& buf, bool fromscratch = false) override { return 0.0; }
  void copyFromBuffer(ParticleSet& P, WFBufferType& buf) override {}
  void evaluateDerivatives(ParticleSet& P,
                           const opt_variables_type& optvars,
                           Vector<ValueType>& dlogpsi,
                           Vector<ValueType>& dhpsioverpsi) override
  {}
  std::unique_ptr<DiracDeterminantBase> makeCopy(std::unique_ptr<SPOSet>&& spo) const override
  {
    return std::make_unique<DummyDiracDetWithoutMW>(getClassName(), std::move(spo), FirstIndex, LastIndex);
  }
};

class DummyDiracDetWithMW : public DummyDiracDetWithoutMW
{
public:
  DummyDiracDetWithMW(const std::string& class_name, std::unique_ptr<SPOSet>&& spos, int first, int last)
      : DummyDiracDetWithoutMW(getClassName(), std::move(spos), first, last)
  {}

  void mw_evalGrad(const RefVectorWithLeader<WaveFunctionComponent>& wfc_list,
                   const RefVectorWithLeader<ParticleSet>& p_list,
                   const int iat,
                   std::vector<GradType>& grad_now) const override
  {
    for (auto& grad : grad_now)
    {
      grad[0] = 321.;
      grad[1] = 654.;
      grad[2] = 987.;
    }
  }
  void mw_evalGradWithSpin(const RefVectorWithLeader<WaveFunctionComponent>& wfc_list,
                           const RefVectorWithLeader<ParticleSet>& p_list,
                           const int iat,
                           std::vector<GradType>& grad_now,
                           std::vector<ComplexType>& spingrad_now) const override
  {
    for (auto& grad : grad_now)
    {
      grad[0] = 0.321;
      grad[1] = 0.654;
      grad[2] = 0.987;
    }
    for (auto& spingrad : spingrad_now)
      spingrad = ComplexType(0.2, 0.1);
  }
  void mw_ratioGrad(const RefVectorWithLeader<WaveFunctionComponent>& wfc_list,
                    const RefVectorWithLeader<ParticleSet>& p_list,
                    int iat,
                    std::vector<PsiValue>& ratios,
                    std::vector<GradType>& grad_new) const override
  {
    for (auto& grad : grad_new)
    {
      grad[0] = 321.;
      grad[1] = 654.;
      grad[2] = 987.;
    }
  }
  void mw_ratioGradWithSpin(const RefVectorWithLeader<WaveFunctionComponent>& wfc_list,
                            const RefVectorWithLeader<ParticleSet>& p_list,
                            int iat,
                            std::vector<PsiValue>& ratios,
                            std::vector<GradType>& grad_new,
                            std::vector<ComplexType>& spingrad_new) const override
  {
    for (auto& grad : grad_new)
    {
      grad[0] = 0.321;
      grad[1] = 0.654;
      grad[2] = 0.987;
    }
    for (auto& spingrad : spingrad_new)
      spingrad = ComplexType(0.2, 0.1);
  }
};

TEST_CASE("SlaterDet mw_ APIs", "[wavefunction]")
{
  Communicate* comm = OHMMS::Controller;

  auto particle_pool = MinimalParticlePool::make_O2_spinor(comm);
  auto& elec0        = *(particle_pool).getParticleSet("e");
  auto& elec1        = *(particle_pool).getParticleSet("e");
  RefVectorWithLeader<ParticleSet> p_list(elec0, {elec0, elec1});

  std::unique_ptr<ConstantSPOSet> spo_ptr0 = std::make_unique<ConstantSPOSet>("dummySPO", 3, 3);
  std::unique_ptr<ConstantSPOSet> spo_ptr1 = std::make_unique<ConstantSPOSet>("dummySPO", 3, 3);
  //Right now, DiracDeterminantBatched has mw_ WithSpin APIs but DiracDeterminant does not.
  //We want to add a test to make sure Slater determinant chooses the mw_ implementation if it has it.

  // First, do without mw_ APIs
  {
    std::unique_ptr<DiracDeterminantBase> det_ptr0 =
        std::make_unique<DummyDiracDetWithoutMW>("dummy", std::move(spo_ptr0), 0, 12);
    std::unique_ptr<DiracDeterminantBase> det_ptr1 =
        std::make_unique<DummyDiracDetWithoutMW>("dummy", std::move(spo_ptr1), 0, 12);

    std::vector<std::unique_ptr<DiracDeterminantBase>> dirac_dets0;
    dirac_dets0.push_back(std::move(det_ptr0));
    std::vector<std::unique_ptr<DiracDeterminantBase>> dirac_dets1;
    dirac_dets1.push_back(std::move(det_ptr1));

    SlaterDet slaterdet0(elec0, std::move(dirac_dets0));
    SlaterDet slaterdet1(elec0, std::move(dirac_dets1));

    RefVectorWithLeader<WaveFunctionComponent> sd_list(slaterdet0, {slaterdet0, slaterdet1});
    ResourceCollection sd_res("test_sd_res");
    slaterdet0.createResource(sd_res);
    ResourceCollectionTeamLock<WaveFunctionComponent> mw_sd_lock(sd_res, sd_list);

    std::vector<SPOSet::GradType> grads(2);
    std::vector<WaveFunctionComponent::PsiValue> ratios(2);
    std::vector<SPOSet::ComplexType> spingrads(2);
    slaterdet0.mw_evalGrad(sd_list, p_list, 0, grads);
    for (auto grad : grads)
    {
      CHECK(grad[0] == ValueApprox(123.));
      CHECK(grad[1] == ValueApprox(456.));
      CHECK(grad[2] == ValueApprox(789.));
    }

    slaterdet0.mw_evalGradWithSpin(sd_list, p_list, 0, grads, spingrads);
    for (auto grad : grads)
    {
      CHECK(grad[0] == ValueApprox(0.123));
      CHECK(grad[1] == ValueApprox(0.456));
      CHECK(grad[2] == ValueApprox(0.789));
    }
    for (auto sgrad : spingrads)
      CHECK(sgrad == ComplexApprox(SPOSet::ComplexType(0.1, 0.2)));

    slaterdet0.mw_ratioGrad(sd_list, p_list, 0, ratios, grads);
    for (auto grad : grads)
    {
      CHECK(grad[0] == ValueApprox(123.));
      CHECK(grad[1] == ValueApprox(456.));
      CHECK(grad[2] == ValueApprox(789.));
    }

    slaterdet0.mw_ratioGradWithSpin(sd_list, p_list, 0, ratios, grads, spingrads);
    for (auto grad : grads)
    {
      CHECK(grad[0] == ValueApprox(0.123));
      CHECK(grad[1] == ValueApprox(0.456));
      CHECK(grad[2] == ValueApprox(0.789));
    }
    for (auto sgrad : spingrads)
      CHECK(sgrad == ComplexApprox(SPOSet::ComplexType(0.1, 0.2)));
  }
  //Now do with MW
  {
    std::unique_ptr<DiracDeterminantBase> det_ptr0 =
        std::make_unique<DummyDiracDetWithMW>("dummy", std::move(spo_ptr0), 0, 12);
    std::unique_ptr<DiracDeterminantBase> det_ptr1 =
        std::make_unique<DummyDiracDetWithMW>("dummy", std::move(spo_ptr1), 0, 12);

    std::vector<std::unique_ptr<DiracDeterminantBase>> dirac_dets0;
    dirac_dets0.push_back(std::move(det_ptr0));
    std::vector<std::unique_ptr<DiracDeterminantBase>> dirac_dets1;
    dirac_dets1.push_back(std::move(det_ptr1));

    SlaterDet slaterdet0(elec0, std::move(dirac_dets0));
    SlaterDet slaterdet1(elec1, std::move(dirac_dets1));

    RefVectorWithLeader<WaveFunctionComponent> sd_list(slaterdet0, {slaterdet0, slaterdet1});

    ResourceCollection sd_res("test_sd_res");
    slaterdet0.createResource(sd_res);
    ResourceCollectionTeamLock<WaveFunctionComponent> mw_sd_lock(sd_res, sd_list);

    std::vector<SPOSet::GradType> grads(2);
    std::vector<WaveFunctionComponent::PsiValue> ratios(2);
    std::vector<SPOSet::ComplexType> spingrads(2);
    slaterdet0.mw_evalGrad(sd_list, p_list, 0, grads);
    for (auto grad : grads)
    {
      CHECK(grad[0] == ValueApprox(321.));
      CHECK(grad[1] == ValueApprox(654.));
      CHECK(grad[2] == ValueApprox(987.));
    }

    slaterdet0.mw_evalGradWithSpin(sd_list, p_list, 0, grads, spingrads);
    for (auto grad : grads)
    {
      CHECK(grad[0] == ValueApprox(0.321));
      CHECK(grad[1] == ValueApprox(0.654));
      CHECK(grad[2] == ValueApprox(0.987));
    }
    for (auto sgrad : spingrads)
      CHECK(sgrad == ComplexApprox(SPOSet::ComplexType(0.2, 0.1)));

    slaterdet0.mw_ratioGrad(sd_list, p_list, 0, ratios, grads);
    for (auto grad : grads)
    {
      CHECK(grad[0] == ValueApprox(321.));
      CHECK(grad[1] == ValueApprox(654.));
      CHECK(grad[2] == ValueApprox(987.));
    }

    slaterdet0.mw_ratioGradWithSpin(sd_list, p_list, 0, ratios, grads, spingrads);
    for (auto grad : grads)
    {
      CHECK(grad[0] == ValueApprox(0.321));
      CHECK(grad[1] == ValueApprox(0.654));
      CHECK(grad[2] == ValueApprox(0.987));
    }
    for (auto sgrad : spingrads)
      CHECK(sgrad == ComplexApprox(SPOSet::ComplexType(0.2, 0.1)));
  }
}

} // namespace qmcplusplus
