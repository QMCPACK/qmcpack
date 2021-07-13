//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2020 QMCPACK developers.
//
// File developed by: Mark Dewing, mdewing@anl.gov, Argonne National Laboratory
//
// File created by: Mark Dewing, mdewing@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////

#include "catch.hpp"
#include "QMCDrivers/WFOpt/QMCCostFunctionBatched.h"


namespace qmcplusplus
{
void compute_batch_parameters(int sample_size, int batch_size, int& num_batches, int& final_batch_size);

TEST_CASE("compute_batch_parameters", "[drivers]")
{
  int sample_size = 1;
  int batch_size  = 1;

  int num_batches;
  int final_batch_size;

  compute_batch_parameters(sample_size, batch_size, num_batches, final_batch_size);
  CHECK(num_batches == 1);
  CHECK(final_batch_size == 1);


  sample_size = 11;
  batch_size  = 4;

  compute_batch_parameters(sample_size, batch_size, num_batches, final_batch_size);
  CHECK(num_batches == 3);
  CHECK(final_batch_size == 3);
}


class QMCCostFunctionBatchedTest : public QMCCostFunctionBatched
{
public:
  QMCCostFunctionBatchedTest(MCWalkerConfiguration& w,
                             TrialWaveFunction& psi,
                             QMCHamiltonian& h,
                             SampleStack& samples,
                             int opt_num_crowds,
                             int crowd_size,
                             Communicate* comm)
      : QMCCostFunctionBatched(w, psi, h, samples, opt_num_crowds, crowd_size, comm)
  {}
  // For unit testing
  Matrix<Return_rt>& getRecordsOnNode() { return RecordsOnNode_; }
  Matrix<Return_rt>& getDerivRecords() { return DerivRecords_; }
  Matrix<Return_rt>& getHDerivRecords() { return HDerivRecords_; }

  std::vector<Return_rt>& getSumValue() { return SumValue; }

  Matrix<Return_rt>& getRecords() { return Records; }

  opt_variables_type& getOptVariables() { return OptVariables; }

  void setNumOptimizables(int n) { NumOptimizables = n; }
};

struct LinearMethod
{
  int numSamples;
  int numParam;
  SampleStack samples;
  MCWalkerConfiguration w;
  QMCHamiltonian h;
  TrialWaveFunction psi;
  QMCCostFunctionBatchedTest costFn;

  LinearMethod(int num_opt_crowds, int crowd_size, Communicate* comm)
      : costFn(w, psi, h, samples, num_opt_crowds, crowd_size, comm)
  {}

  void set_samples_and_param(int nsamples, int nparam)
  {
    numSamples = nsamples;
    numParam   = nparam;
    samples.setTotalNum(1);
    samples.setMaxSamples(numSamples);

    costFn.getOptVariables().insert("var1", 1.0);
    costFn.setNumOptimizables(numParam);

    costFn.getRecordsOnNode().resize(numSamples, QMCCostFunctionBase::SUM_INDEX_SIZE);
    costFn.getDerivRecords().resize(numSamples, numParam);
    costFn.getHDerivRecords().resize(numSamples, numParam);
  }
};

TEST_CASE("fillOverlapAndHamiltonianMatrices", "[drivers]")
{
  int num_opt_crowds = 1;
  int crowd_size     = 1;

  using Return_rt = qmcplusplus::QMCTraits::RealType;

  Communicate* comm = OHMMS::Controller;

  LinearMethod lin(num_opt_crowds, crowd_size, comm);

  int numSamples = 1;
  int numParam   = 1;
  lin.set_samples_and_param(numSamples, numParam);

  QMCCostFunctionBatchedTest& costFn = lin.costFn;

  std::vector<Return_rt>& SumValue           = costFn.getSumValue();
  SumValue[QMCCostFunctionBase::SUM_WGT]     = 1.0;
  SumValue[QMCCostFunctionBase::SUM_E_WGT]   = -1.3;
  SumValue[QMCCostFunctionBase::SUM_ESQ_WGT] = 1.69;

  auto& RecordsOnNode                               = costFn.getRecordsOnNode();
  RecordsOnNode(0, QMCCostFunctionBase::REWEIGHT)   = 1.0;
  RecordsOnNode(0, QMCCostFunctionBase::ENERGY_NEW) = -1.4;


  auto& derivRecords = costFn.getDerivRecords();
  derivRecords(0, 0) = 1.1;

  auto& HDerivRecords = costFn.getDerivRecords();
  HDerivRecords(0, 0) = -1.2;

  int N = numParam + 1;
  Matrix<Return_rt> ham(N, N);
  Matrix<Return_rt> ovlp(N, N);
  costFn.fillOverlapHamiltonianMatrices(ham, ovlp);

  CHECK(ovlp(0, 0) == Approx(1.0));
  CHECK(ovlp(1, 0) == Approx(0.0));
  CHECK(ovlp(0, 1) == Approx(0.0));
  // With one sample, value is always zero
  CHECK(ovlp(1, 1) == Approx(0.0));

  CHECK(ham(0, 0) == Approx(-1.3));
  // With one sample, values are always zero
  CHECK(ham(1, 0) == Approx(0.0));
  CHECK(ham(0, 1) == Approx(0.0));
  CHECK(ham(1, 1) == Approx(0.0));
}


} // namespace qmcplusplus
