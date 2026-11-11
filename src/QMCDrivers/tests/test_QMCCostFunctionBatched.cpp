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
#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/matchers/catch_matchers_exception.hpp>
#include "Utilities/for_testing/Catch2Approx.h"
#include "QMCDrivers/WFOpt/QMCCostFunctionBatched.h"
#include "MinimalParticlePool.h"
#include "MinimalWaveFunctionPool.h"
#include "MinimalHamiltonianPool.h"
#include "QMCWaveFunctions/Jastrow/RadialJastrowBuilder.h"
#include "Particle/MCSample.h"
#include "ProjectData.h"
#include "StdRandom.h"
#include "FillData.h"
// Input data and gold data for fillFromText test
#include "diamond_fill_data.h"
#include "Utilities/RuntimeOptions.h"


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

  batch_size = 0;
  compute_batch_parameters(sample_size, batch_size, num_batches, final_batch_size);
  CHECK(num_batches == 0);
  CHECK(final_batch_size == 0);
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
  RuntimeOptions runtime_options_;
  TrialWaveFunction psi = TrialWaveFunction(runtime_options_);
  QMCCostFunctionBatched costFn;

  LinearMethodTestSupport(const std::vector<int>& walkers_per_crowd, Communicate& comm)
      : w(simulation_cell), costFn(w, psi, h, samples, walkers_per_crowd, comm)
  {}

  std::vector<QMCCostFunctionBase::Return_rt>& getSumValue() { return costFn.SumValue; }
  Matrix<QMCCostFunctionBase::Return_rt>& getRecordsOnNode() { return costFn.RecordsOnNode_; }
  Matrix<QMCCostFunctionBase::Return_t>& getDerivRecords() { return costFn.DerivRecords_; }
  Matrix<QMCCostFunctionBase::Return_rt>& getHDerivRecords() { return costFn.HDerivRecords_; }
  const std::vector<int>& getWalkersPerCrowd() const { return costFn.walkers_per_crowd_; }

  void set_samples_and_param(int nsamples, int nparam)
  {
    numSamples = nsamples;
    numParam   = nparam;

    costFn.rank_local_num_samples_ = nsamples;

    for (int i = 0; i < nparam; i++)
    {
      std::string varname = "var" + std::to_string(i);
      costFn.opt_vars.insert(varname, 1.0);
    }

    getRecordsOnNode().resize(numSamples, QMCCostFunctionBase::SUM_INDEX_SIZE);
    getDerivRecords().resize(numSamples, numParam);
    getHDerivRecords().resize(numSamples, numParam);
  }
};

class OptimizationCrowdTestCost : public QMCCostFunctionBatched
{
public:
  using QMCCostFunctionBatched::correlatedSampling;
  using QMCCostFunctionBatched::QMCCostFunctionBatched;

  const auto& getSums() const { return SumValue; }
  const auto& getRecords() const { return RecordsOnNode_; }
  const auto& getDerivatives() const { return DerivRecords_; }
  const auto& getEnergyDerivatives() const { return HDerivRecords_; }
};

} // namespace testing

TEST_CASE("QMCCostFunctionBatched active crowds", "[drivers][optimization_crowds]")
{
  Communicate& comm(*OHMMS::Controller);
  for (const std::vector<int>& capacities : {std::vector<int>{1, 0}, {0, 1}, {0, 2, 0, 1, 0}, {2, 1}})
  {
    CAPTURE(capacities);
    testing::LinearMethodTestSupport test(capacities, comm);
    std::vector<int> expected;
    for (int capacity : capacities)
      if (capacity > 0)
        expected.push_back(capacity);
    CHECK(test.getWalkersPerCrowd() == expected);
  }

  for (const std::vector<int>& capacities : {std::vector<int>{}, {0}, {0, 0}})
  {
    CAPTURE(capacities);
    CHECK_THROWS_MATCHES(testing::LinearMethodTestSupport(capacities, comm), std::runtime_error,
                         Catch::Matchers::Message("QMCCostFunctionBatched requires at least one active crowd."));
  }
  CHECK_THROWS_MATCHES(testing::LinearMethodTestSupport({1, -1}, comm), std::runtime_error,
                       Catch::Matchers::Message("QMCCostFunctionBatched crowd capacities must not be negative."));
}

TEST_CASE("QMCCostFunctionBatched empty crowds do not discard samples", "[drivers][optimization_crowds]")
{
  Communicate& comm(*OHMMS::Controller);
  const bool use_sr = GENERATE(false, true);
  for (const std::vector<int>& capacities : {std::vector<int>{1, 0}, {0, 1}, {0, 2, 0, 1, 0}, {2, 1}})
  {
    CAPTURE(use_sr, capacities);
    ProjectData project;
    auto particles     = MinimalParticlePool::make_diamondC_1x1x1(comm);
    auto wavefunctions = MinimalWaveFunctionPool::make_diamondC_1x1x1(project.getRuntimeOptions(), comm, particles);
    auto hamiltonians  = MinimalHamiltonianPool::make_hamWithEE(comm, particles, wavefunctions);
    auto& pset         = *particles.getParticleSet("e");
    auto& psi          = wavefunctions.getWaveFunction().value().get();
    auto& ham          = hamiltonians.getHamiltonian().value().get();
    Libxml2Document jastrow_doc;
    REQUIRE(jastrow_doc.parseFromString(R"(
      <jastrow name="J2" type="Two-Body" function="Bspline">
        <correlation rcut="1.5" size="4" speciesA="u" speciesB="d">
          <coefficients id="ud" type="Array">0.1 0.08 0.04 0.01</coefficients>
        </correlation>
      </jastrow>)"));
    RadialJastrowBuilder builder(comm, pset);
    psi.addComponent(builder.buildComponent(jastrow_doc.getRoot()));

    SampleStack samples;
    constexpr int num_samples = 10;
    samples.setMaxSamples(num_samples);
    for (int i = 0; i < num_samples; ++i)
    {
      pset.R[0][0] += 0.01;
      samples.appendSample(MCSample(pset));
    }
    StdRandom<QMCTraits::FullPrecRealType> rng(1234);
    RefVector<RandomBase<QMCTraits::FullPrecRealType>> rngs{rng};
    testing::OptimizationCrowdTestCost reference(pset, psi, ham, samples, {1}, comm);
    testing::OptimizationCrowdTestCost cost(pset, psi, ham, samples, capacities, comm);
    Libxml2Document doc;
    REQUIRE(doc.parseFromString("<qmc/>"));
    auto initialize = [&](auto& target) {
      REQUIRE(target.put(doc.getRoot()));
      REQUIRE(target.getNumParams() > 0);
      target.setRng(rngs);
      target.getConfigurations("");
      NullEngineHandle handle;
      if (use_sr)
        target.checkConfigurationsSR(handle);
      else
        target.checkConfigurations(handle);
      REQUIRE(target.getNumSamples() == num_samples * comm.size());
    };
    initialize(reference);
    initialize(cost);

    auto check_results = [&]() {
      for (int i = 0; i < reference.getSums().size(); ++i)
        CHECK(cost.getSums()[i] == Approx(reference.getSums()[i]));
      for (int sample = 0; sample < num_samples; ++sample)
      {
        CAPTURE(sample);
        for (int field = QMCCostFunctionBase::LOGPSI_FIXED; field <= QMCCostFunctionBase::REWEIGHT; ++field)
          CHECK(cost.getRecords()(sample, field) == Approx(reference.getRecords()(sample, field)));
        for (int param = 0; param < reference.getNumParams(); ++param)
        {
          CHECK(std::abs(cost.getDerivatives()(sample, param) - reference.getDerivatives()(sample, param)) ==
                Approx(0.0).margin(1e-10));
          if (!use_sr)
            CHECK(cost.getEnergyDerivatives()(sample, param) ==
                  Approx(reference.getEnergyDerivatives()(sample, param)));
        }
      }
    };
    check_results();

    if (!use_sr)
      for (bool need_grad : {false, true})
      {
        reference.Params(0) += 0.01;
        cost.Params(0) += 0.01;
        const auto reference_weight = reference.correlatedSampling(need_grad);
        CHECK(cost.correlatedSampling(need_grad) == Approx(reference_weight));
        check_results();
      }
  }
}

TEST_CASE("fillOverlapAndHamiltonianMatrices", "[drivers]")
{
  std::vector<int> walkers_per_crowd{1};

  using Return_rt = qmcplusplus::QMCTraits::RealType;

  Communicate& comm(*OHMMS::Controller);
  if (comm.size() > 1)
    return;

  testing::LinearMethodTestSupport lin(walkers_per_crowd, comm);

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
  std::vector<int> walkers_per_crowd(num_opt_crowds, 1);

  using Return_rt = qmcplusplus::QMCTraits::RealType;

  Communicate& comm(*OHMMS::Controller);
  if (comm.size() > 1)
    return;

  testing::LinearMethodTestSupport lin(walkers_per_crowd, comm);

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
  Communicate& comm(*OHMMS::Controller);
  if (comm.size() > 1)
    return;

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
