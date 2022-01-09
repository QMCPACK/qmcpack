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
#include "FillData.h"
// Input data and gold data for fillFromText test
#include "diamond_fill_data.h"


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

namespace testing
{
class LinearMethodTestSupport
{
public:
  int numSamples;
  int numParam;
  SampleStack samples;
  const SimulationCell simulation_cell;
  MCWalkerConfiguration w;
  QMCHamiltonian h;
  TrialWaveFunction psi;
  QMCCostFunctionBatched costFn;

  LinearMethodTestSupport(int num_opt_crowds, int crowd_size, Communicate* comm)
      : w(simulation_cell), costFn(w, psi, h, samples, num_opt_crowds, crowd_size, comm)
  {}

  std::vector<QMCCostFunctionBase::Return_rt>& getSumValue() { return costFn.SumValue; }
  Matrix<QMCCostFunctionBase::Return_rt>& getRecordsOnNode() { return costFn.RecordsOnNode_; }
  Matrix<QMCCostFunctionBase::Return_rt>& getDerivRecords() { return costFn.DerivRecords_; }
  Matrix<QMCCostFunctionBase::Return_rt>& getHDerivRecords() { return costFn.HDerivRecords_; }

  void set_samples_and_param(int nsamples, int nparam)
  {
    numSamples = nsamples;
    numParam   = nparam;

    costFn.rank_local_num_samples_ = nsamples;

    for (int i = 0; i < nparam; i++)
    {
      std::string varname = "var" + std::to_string(i);
      costFn.OptVariables.insert(varname, 1.0);
    }

    costFn.NumOptimizables = numParam;

    getRecordsOnNode().resize(numSamples, QMCCostFunctionBase::SUM_INDEX_SIZE);
    getDerivRecords().resize(numSamples, numParam);
    getHDerivRecords().resize(numSamples, numParam);
  }
};

} // namespace testing

TEST_CASE("fillOverlapAndHamiltonianMatrices", "[drivers]")
{
  int num_opt_crowds = 1;
  int crowd_size     = 1;

  using Return_rt = qmcplusplus::QMCTraits::RealType;

  Communicate* comm = OHMMS::Controller;

  testing::LinearMethodTestSupport lin(num_opt_crowds, crowd_size, comm);

  int numSamples = 1;
  int numParam   = 1;
  lin.set_samples_and_param(numSamples, numParam);

  std::vector<Return_rt>& SumValue           = lin.getSumValue();
  SumValue[QMCCostFunctionBase::SUM_WGT]     = 1.0;
  SumValue[QMCCostFunctionBase::SUM_E_WGT]   = -1.3;
  SumValue[QMCCostFunctionBase::SUM_ESQ_WGT] = 1.69;

  auto& RecordsOnNode                               = lin.getRecordsOnNode();
  RecordsOnNode(0, QMCCostFunctionBase::REWEIGHT)   = 1.0;
  RecordsOnNode(0, QMCCostFunctionBase::ENERGY_NEW) = -1.4;

  auto& derivRecords = lin.getDerivRecords();
  derivRecords(0, 0) = 1.1;

  auto& HDerivRecords = lin.getHDerivRecords();
  HDerivRecords(0, 0) = -1.2;

  int N = numParam + 1;
  Matrix<Return_rt> ham(N, N);
  Matrix<Return_rt> ovlp(N, N);
  lin.costFn.fillOverlapHamiltonianMatrices(ham, ovlp);

  CHECK(ovlp(0, 0) == Approx(1.0));
  CHECK(ovlp(1, 0) == Approx(0.0));
  CHECK(ovlp(0, 1) == Approx(0.0));
  // With one sample, value is always zero
  CHECK(ovlp(1, 1) == Approx(0.0));

  CHECK(ham(0, 0) == Approx(-1.3));
  // With one sample, values are always zero
  CHECK(ham(1, 0) == Approx(0.0));
  CHECK(ham(0, 1) == Approx(-1.2));
  CHECK(ham(1, 1) == Approx(0.0));
}

// Test QMCCostFunctionBatched::fillOverlapHamiltonianMatrices
// Inputs are the number of crowds (threads) and
// the input/gold data (from a file created by convert_hdf_to_cpp.py)
void fill_from_text(int num_opt_crowds, FillData& fd)
{
  // Not used in the function under test
  int crowd_size = 1;

  using Return_rt = qmcplusplus::QMCTraits::RealType;

  Communicate* comm = OHMMS::Controller;

  testing::LinearMethodTestSupport lin(num_opt_crowds, crowd_size, comm);

  int numSamples = fd.numSamples;
  int numParam   = fd.numParam;
  lin.set_samples_and_param(numSamples, numParam);

  std::vector<Return_rt>& SumValue           = lin.getSumValue();
  SumValue[QMCCostFunctionBase::SUM_WGT]     = fd.sum_wgt;
  SumValue[QMCCostFunctionBase::SUM_E_WGT]   = fd.sum_e_wgt;
  SumValue[QMCCostFunctionBase::SUM_ESQ_WGT] = fd.sum_esq_wgt;

  auto& RecordsOnNode = lin.getRecordsOnNode();
  for (int iw = 0; iw < numSamples; iw++)
  {
    RecordsOnNode(iw, QMCCostFunctionBase::REWEIGHT)   = fd.reweight[iw];
    RecordsOnNode(iw, QMCCostFunctionBase::ENERGY_NEW) = fd.energy_new[iw];
  }

  auto& derivRecords = lin.getDerivRecords();
  derivRecords       = fd.derivRecords;

  auto& HDerivRecords = lin.getHDerivRecords();
  HDerivRecords       = fd.HDerivRecords;

  int N = numParam + 1;
  Matrix<Return_rt> ham(N, N);
  Matrix<Return_rt> ovlp(N, N);
  lin.costFn.fillOverlapHamiltonianMatrices(ham, ovlp);

  for (int iw = 0; iw < numParam; iw++)
  {
    for (int iw2 = 0; iw2 < numParam; iw2++)
    {
      //app_log() << "iw = " << iw << " iw2 = " << iw2 << " ovlp = " << ovlp(iw,iw2) << " " << ovlp_gold(iw,iw2);
      //app_log() << " ham = " << ham(iw,iw2) << " " << ham_gold(iw,iw2) << std::endl;
      CHECK(ovlp(iw, iw2) == Approx(fd.ovlp_gold(iw, iw2)));
      CHECK(ham(iw, iw2) == Approx(fd.ham_gold(iw, iw2)));
    }
  }
}


// Test fillOverlapHamiltonianMatrices function using gold data
// This can test the parallelization of that function using different numbers of crowds
TEST_CASE("fillfromText", "[drivers]")
{
  FillData fd;

  // Generated from short-diamondC_1x1x1_pp-opt_sdj-1-16 with 1 thread and 10 samples
  // And the bsplines sizes were reduced from 8 to 4 to give a total of 12 parameters
  get_diamond_fill_data(fd);

  // Test for 1 and 2 crowds (threads)
  for (int num_opt_crowds = 1; num_opt_crowds < 3; num_opt_crowds++)
  {
    fill_from_text(num_opt_crowds, fd);
  }
}


} // namespace qmcplusplus
