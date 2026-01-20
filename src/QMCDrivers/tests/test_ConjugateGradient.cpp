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
#include "QMCDrivers/WFOpt/ConjugateGradient.h"
#include "QMCDrivers/WFOpt/QMCCostFunctionBase.h"

namespace qmcplusplus
{

namespace testing
{

class LinearSystem : public QMCCostFunctionBase
{
public:
  LinearSystem(MCWalkerConfiguration& w, TrialWaveFunction& psi, QMCHamiltonian& h, Communicate* comm)
      : QMCCostFunctionBase(w, psi, h, comm){};
  void GradCost(std::vector<Return_rt>& PGradient, const std::vector<Return_rt>& PM, Return_rt FiniteDiff = 0) override
  {}

  void resetPsi(bool final_reset) override{};
  Return_rt fillOverlapHamiltonianMatrices(Matrix<Return_rt>& Left, Matrix<Return_rt>& Right) override { return 0; }
  void getConfigurations(const std::string& aroot) override {}
  void checkConfigurations(EngineHandle& handle) override {}
  EffectiveWeight correlatedSampling(bool needGrad) override { return 0; }

  void setAMatrix(std::vector<std::vector<Return_rt>>& Amat)
  {
    int np = Amat.size();
    Amat_.resize(np * np);
    for (int i = 0; i < np; i++)
    {
      for (int j = 0; j < np; j++)
        Amat_[i * np + j] = Amat[i][j];
      std::string name = "tmp" + std::to_string(i);
      opt_vars.insert(name, 0);
    }
  }

  void calcOvlParmVec(const std::vector<Return_rt>& param, std::vector<Return_rt>& ovlParmVec) override
  {
    assert(ovlParmVec.size() == param.size());
    assert(param.size() * param.size() == Amat_.size());

    std::fill(ovlParmVec.begin(), ovlParmVec.end(), 0.0);
    for (int i = 0; i < param.size(); i++)
      for (int j = 0; j < param.size(); j++)
        ovlParmVec[i] += Amat_[i * param.size() + j] * param[j];
  }

private:
  std::vector<Return_rt> Amat_;
};

} // namespace testing

using Real = QMCTraits::RealType;
TEST_CASE("ConjugateGradient", "[drivers]")
{
  //Testing that congjuate gradient algorithm correctly solves a particular Ax = b linear system

  //Set A matrix
  std::vector<Real> row0{0.25592304, 0.59014979, 0.46571581, 0.12010389, 0.27168455};
  std::vector<Real> row1{0.59014979, 0.17202548, 0.95357168, 0.57385606, 0.36200427};
  std::vector<Real> row2{0.46571581, 0.95357168, 0.41047519, 0.21899367, 0.80101095};
  std::vector<Real> row3{0.12010389, 0.57385606, 0.21899367, 0.65777236, 0.66479641};
  std::vector<Real> row4{0.27168455, 0.36200427, 0.80101095, 0.66479641, 0.30727352};
  std::vector<std::vector<Real>> Amat{row0, row1, row2, row3, row4};

  //Ainv found from numpy.linalg.inv
  std::vector<Real> invrow0{26.82697284, 14.05843768, -18.80573026, 16.63247607, -27.24381534};
  std::vector<Real> invrow1{14.05843768, 4.38858593, -8.87096117, 7.60630171, -10.93179847};
  std::vector<Real> invrow2{-18.80573026, -8.87096117, 13.65194417, -12.83325154, 19.25546651};
  std::vector<Real> invrow3{16.63247607, 7.60630171, -12.83325154, 11.60545707, -15.32182682};
  std::vector<Real> invrow4{-27.24381534, -10.93179847, 19.25546651, -15.32182682, 23.17523938};
  std::vector<std::vector<Real>> Ainv{invrow0, invrow1, invrow2, invrow3, invrow4};

  std::vector<Real> x{0.79414581, 0.61498605, 0.12401388, 0.33310696, 0.45705058};

  std::vector<Real> b{0.78811034, 1.04932409, 1.44623499, 0.99840586, 0.89960905};

  // Sanity check
  int np = x.size();
  std::vector<Real> tmp(np, 0);
  for (int i = 0; i < np; i++)
    for (int j = 0; j < np; j++)
      tmp[i] += Amat[i][j] * x[j];
  for (int i = 0; i < np; i++)
    CHECK(tmp[i] == Approx(b[i]));

  std::fill(tmp.begin(), tmp.end(), 0);
  for (int i = 0; i < np; i++)
    for (int j = 0; j < np; j++)
      tmp[i] += Ainv[i][j] * b[j];
  for (int i = 0; i < np; i++)
    CHECK(tmp[i] == Approx(x[i]));


  const SimulationCell simulation_cell;
  MCWalkerConfiguration w(simulation_cell);
  QMCHamiltonian h;
  RuntimeOptions runtime_options;
  TrialWaveFunction psi(runtime_options);
  Communicate* comm = OHMMS::Controller;

  testing::LinearSystem ls(w, psi, h, comm);
  ls.setAMatrix(Amat);

  ConjugateGradient cg;
  std::vector<Real> soln(np, 0);
  int niterations = cg.run(ls, b, soln);
  for (int i = 0; i < np; i++)
    CHECK(soln[i] == Approx(x[i]).epsilon(0.001));
}

} // namespace qmcplusplus
