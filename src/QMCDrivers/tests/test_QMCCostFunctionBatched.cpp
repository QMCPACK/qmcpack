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

namespace testing
{
class LinearMethodTestSupport
{
public:
  int numSamples;
  int numParam;
  SampleStack samples;
  MCWalkerConfiguration w;
  QMCHamiltonian h;
  TrialWaveFunction psi;
  QMCCostFunctionBatched costFn;

  LinearMethodTestSupport(int num_opt_crowds, int crowd_size, Communicate* comm)
      : costFn(w, psi, h, samples, num_opt_crowds, crowd_size, comm)
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

void fill_from_h5(std::string& prefix);

TEST_CASE("fillFromH5", "[drivers]")
{
  std::vector<std::string> prefixes = {"diamondC_1x1x1"};
  for (auto prefix : prefixes)
  {
    fill_from_h5(prefix);
  }
}

void fill_from_h5(std::string& prefix)
{
  int num_opt_crowds = 1;
  int crowd_size     = 1;
  using Return_rt    = qmcplusplus::QMCTraits::RealType;

  Communicate* comm = OHMMS::Controller;
  testing::LinearMethodTestSupport lin(num_opt_crowds, crowd_size, comm);

  std::string matrix_inputs = prefix + "matrix_inputs.h5";
  hdf_archive h5_inputs;
  h5_inputs.open(matrix_inputs, H5F_ACC_RDONLY);

  std::vector<Return_rt>& SumValue = lin.getSumValue();
  h5_inputs.read(SumValue[QMCCostFunctionBase::SUM_WGT], "weight");
  h5_inputs.read(SumValue[QMCCostFunctionBase::SUM_E_WGT], "e_weight");
  h5_inputs.read(SumValue[QMCCostFunctionBase::SUM_ESQ_WGT], "e_sq_weight");


  std::vector<int> sizes(2);
  h5_inputs.getShape<Return_rt>("deriv_records", sizes);
  int numSamples = sizes[0];
  int numParam   = sizes[1];

  lin.set_samples_and_param(numSamples, numParam);

  std::vector<Return_rt> reweight(numSamples);
  std::vector<Return_rt> energy_new(numSamples);
  h5_inputs.read(reweight, "reweight");
  h5_inputs.read(energy_new, "energy_new");
  //app_log() << "size of reweight = " << reweight.size() << std::endl;

  auto& RecordsOnNode = lin.getRecordsOnNode();
  for (int iw = 0; iw < numSamples; iw++)
  {
    RecordsOnNode(iw, QMCCostFunctionBase::REWEIGHT)   = reweight[iw];
    RecordsOnNode(iw, QMCCostFunctionBase::ENERGY_NEW) = energy_new[iw];
  }

  auto& derivRecords = lin.getDerivRecords();
  h5_inputs.read(derivRecords, "deriv_records");
  auto& HDerivRecords = lin.getHDerivRecords();
  h5_inputs.read(HDerivRecords, "H_deriv_records");


  int N = numParam + 1;
  Matrix<Return_rt> ham(N, N);
  Matrix<Return_rt> ovlp(N, N);
  lin.costFn.fillOverlapHamiltonianMatrices(ham, ovlp);


  std::string linear_matrices = prefix + "_linear_matrices.h5";
  hdf_archive h5_matrices;
  h5_matrices.open(linear_matrices, H5F_ACC_RDONLY);
  Matrix<Return_rt> ovlp_gold(N, N);
  Matrix<Return_rt> ham_gold(N, N);
  h5_matrices.read(ovlp_gold, "overlap");
  h5_matrices.read(ham_gold, "Hamiltonian");

  for (int iw = 0; iw < numParam; iw++)
  {
    for (int iw2 = 0; iw2 < numParam; iw2++)
    {
      //app_log() << "iw = " << iw << " iw2 = " << iw2 << " ovlp = " << ovlp(iw,iw2) << " " << ovlp_gold(iw,iw2);
      //app_log() << " ham = " << ham(iw,iw2) << " " << ham_gold(iw,iw2) << std::endl;
      CHECK(ovlp(iw, iw2) == Approx(ovlp_gold(iw, iw2)));
      CHECK(ham(iw, iw2) == Approx(ham_gold(iw, iw2)));
    }
  }
}


} // namespace qmcplusplus
