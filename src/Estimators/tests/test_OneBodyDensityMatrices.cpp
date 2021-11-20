//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2021 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Lab
//
// File created by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Lab
//////////////////////////////////////////////////////////////////////////////////////


#include "catch.hpp"

#include "OneBodyDensityMatrices.h"
#include "ValidOneBodyDensityMatricesInput.h"
#include "InvalidOneBodyDensityMatricesInput.h"
#include "EstimatorTesting.h"
#include "EstimatorInput.h"
#include "ParticleSet.h"
#include "TrialWaveFunction.h"
#include "OhmmsData/Libxml2Doc.h"
#include "Message/UniformCommunicateError.h"
#include "Particle/tests/MinimalParticlePool.h"
#include "QMCWaveFunctions/tests/MinimalWaveFunctionPool.h"
#include "Utilities/StdRandom.h"
#include "Utilities/StlPrettyPrint.hpp"
#include "Utilities/for_testing/NativeInitializerPrint.hpp"
#include <iostream>

namespace qmcplusplus
{

// set to true to regenrate static testing data
constexpr bool generate_test_data = false;

namespace testing
{
using OBDMI = OneBodyDensityMatricesInput;

template<typename T>
class OneBodyDensityMatricesTests
{
public:
  using Evaluators  = OneBodyDensityMatricesInput::Evaluator;
  using Integrators = OneBodyDensityMatricesInput::Integrator;
  using Sampling    = OneBodyDensityMatrices::Sampling;
  using MCPWalker   = OneBodyDensityMatrices::MCPWalker;
  using Data        = OneBodyDensityMatrices::Data;
  using Real        = Data::value_type;

  void testCopyConstructor(const OneBodyDensityMatrices& obdm)
  {
    OneBodyDensityMatrices obdm2(obdm);
    CHECK(obdm.sampling_ == obdm2.sampling_);
    CHECK(obdm.data_ != obdm2.data_);
  }

  OneBodyDensityMatricesTests() = default;
  void testGenerateSamples(onebodydensitymatrices::Inputs input,
                           OneBodyDensityMatrices& obdm,
                           ParticleSet& pset_target,
                           StdRandom<T>& rng)
  {
    using namespace onebodydensitymatrices;
    switch (input)
    {
    case (valid_obdm_input):
      obdm.generateSamples(1.0, pset_target, rng);
      CHECK(obdm.nmoves_ == 64);
      break;
    case (valid_obdm_input_scale):
      obdm.generateSamples(1.0, pset_target, rng);
      CHECK(obdm.nmoves_ == 0);
      break;
    case (valid_obdm_input_grid):
      obdm.generateSamples(1.0, pset_target, rng);
      CHECK(obdm.nmoves_ == 0);
      CHECK(obdm.samples_ == pow(22, OHMMS_DIM));
      break;
    }
  }

  /** Checking approximate equality for complex valued data as
   *  two reals is not consistent with testing practices elsewhere in the code.
   *  Values that are slightly off may now fall in approximation limit properly.
   *
   *  The smell from OperatorEstBase continuing the tradition of
   *  packing everything into a Real buffer is not lost on me.
   */
  void checkData(Real* ref_in, Real* test_in, size_t size)
  {
    if constexpr (IsComplex_t<OneBodyDensityMatrices::Value>::value)
    {
      size /= 2;
      auto* ref_data  = reinterpret_cast<std::complex<Real>*>(ref_in);
      auto* test_data = reinterpret_cast<std::complex<Real>*>(test_in);
      for (size_t id = 0; id < size; id += 2)
        CHECK(ref_data[id] == ComplexApprox(test_data[id]));
    }
    else
    {
      for (size_t id = 0; id < size; ++id)
        CHECK(ref_in[id] == Approx(test_in[id]));
    }
  }

  /** no change test for accumulate
   */
  void testAccumulate(OneBodyDensityMatrices& obdm,
                      RefVector<MCPWalker>& walkers,
                      RefVector<ParticleSet>& psets,
                      RefVector<TrialWaveFunction>& twfcs,
                      StdRandom<T>& rng)
  {
    obdm.implAccumulate(walkers, psets, twfcs, rng);
    Data data(getAccumulateData());
    auto& returned_data = obdm.data_;
    checkData(data.data(), returned_data.data(), data.size());
  }

  /** no change test for evaluateMatrix.
   */
  void testEvaluateMatrix(OneBodyDensityMatrices& obdm,
                          ParticleSet& pset,
                          TrialWaveFunction& trial_wavefunction,
                          MCPWalker& walker,
                          StdRandom<T>& rng)
  {
    obdm.evaluateMatrix(pset, trial_wavefunction, walker, rng);
    Data data(getEvaluateMatrixData(obdm.input_.get_integrator()));
    auto& returned_data = obdm.data_;
    checkData(returned_data.data(), data.data(), data.size());
  }

  void dumpData(OneBodyDensityMatrices& obdm)
  {
    std::cout << "Here is what is in your OneBodyDensityMatrices:\n" << NativePrint(obdm.data_) << '\n';
  }

private:
  Data getEvaluateMatrixData(OBDMI::Integrator integrator);
  Data getAccumulateData();
};

} // namespace testing


TEST_CASE("OneBodyDensityMatrices::OneBodyDensityMatrices", "[estimators]")
{
  using namespace testing;
  using namespace onebodydensitymatrices;
  Libxml2Document doc;
  bool okay = doc.parseFromString(valid_one_body_density_matrices_input_sections[valid_obdm_input]);
  if (!okay)
    throw std::runtime_error("cannot parse OneBodyDensitMatricesInput section");
  xmlNodePtr node = doc.getRoot();
  OneBodyDensityMatricesInput obdmi(node);
  auto lattice     = testing::makeTestLattice();
  auto species_set = testing::makeSpeciesSet(SpeciesCases::GOOD);

  Communicate* comm;
  comm = OHMMS::Controller;
  outputManager.pause();

  MinimalParticlePool mpp;
  ParticleSetPool particle_pool = mpp(comm);
  MinimalWaveFunctionPool wfp;
  WaveFunctionPool wavefunction_pool = wfp(comm, particle_pool);
  auto& pset_target                  = *(particle_pool.getParticleSet("e"));
  auto& wf_factory                   = *(wavefunction_pool.getWaveFunctionFactory("wavefunction"));

  // Good constructor
  OneBodyDensityMatrices obdm(std::move(obdmi), lattice, species_set, wf_factory, pset_target);
  // make sure there is something in obdm's data
  OEBAccessor oeba(obdm);
  oeba[0] = 1.0;
  testing::OneBodyDensityMatricesTests<double> obdmt;
  obdmt.testCopyConstructor(obdm);

  species_set = testing::makeSpeciesSet(SpeciesCases::NO_MEMBERSIZE);
  CHECK_THROWS_AS(OneBodyDensityMatrices(std::move(obdmi), lattice, species_set, wf_factory, pset_target),
                  UniformCommunicateError);

  outputManager.resume();
}

TEST_CASE("OneBodyDensityMatrices::generateSamples", "[estimators]")
{
  using namespace testing;
  using namespace onebodydensitymatrices;

  using MCPWalker = OperatorEstBase::MCPWalker;

  Communicate* comm;
  comm = OHMMS::Controller;
  outputManager.pause();

  MinimalParticlePool mpp;
  ParticleSetPool particle_pool = mpp(comm);
  MinimalWaveFunctionPool wfp;
  WaveFunctionPool wavefunction_pool = wfp(comm, particle_pool);
  wavefunction_pool.setPrimary(wavefunction_pool.getWaveFunction("psi0"));
  auto& pset_target = *(particle_pool.getParticleSet("e"));
  auto& species_set = pset_target.getSpeciesSet();
  auto& wf_factory  = *(wavefunction_pool.getWaveFunctionFactory("wavefunction"));

  auto samplingCaseRunner = [&pset_target, &species_set, &wf_factory](Inputs test_case) {
    Libxml2Document doc;

    bool okay = doc.parseFromString(valid_one_body_density_matrices_input_sections[test_case]);
    if (!okay)
      throw std::runtime_error("cannot parse OneBodyDensitMatricesInput section");
    xmlNodePtr node = doc.getRoot();
    OneBodyDensityMatricesInput obdmi(node);

    OneBodyDensityMatrices obDenMat(std::move(obdmi), pset_target.Lattice, species_set, wf_factory, pset_target);

    OneBodyDensityMatricesTests<double> obdmt;
    //Get control over which rng is used.
    //we don't want FakeRandom.
    StdRandom<double> rng;
    obdmt.testGenerateSamples(test_case, obDenMat, pset_target, rng);
  };

  samplingCaseRunner(valid_obdm_input);
  samplingCaseRunner(valid_obdm_input_scale);
  samplingCaseRunner(valid_obdm_input_grid);

  outputManager.resume();
}

TEST_CASE("OneBodyDensityMatrices::spawnCrowdClone()", "[estimators]")
{
  using namespace testing;
  using namespace onebodydensitymatrices;

  using MCPWalker = OperatorEstBase::MCPWalker;

  Communicate* comm;
  comm = OHMMS::Controller;
  outputManager.pause();

  MinimalParticlePool mpp;
  ParticleSetPool particle_pool = mpp(comm);
  MinimalWaveFunctionPool wfp;
  WaveFunctionPool wavefunction_pool = wfp(comm, particle_pool);
  wavefunction_pool.setPrimary(wavefunction_pool.getWaveFunction("psi0"));
  auto& pset_target = *(particle_pool.getParticleSet("e"));
  auto& species_set = pset_target.getSpeciesSet();
  auto& wf_factory  = *(wavefunction_pool.getWaveFunctionFactory("wavefunction"));

  Libxml2Document doc;
  bool okay = doc.parseFromString(valid_one_body_density_matrices_input_sections[Inputs::valid_obdm_input]);
  if (!okay)
    throw std::runtime_error("cannot parse OneBodyDensitMatricesInput section");
  xmlNodePtr node = doc.getRoot();
  OneBodyDensityMatricesInput obdmi(node);

  OneBodyDensityMatrices original(std::move(obdmi), pset_target.Lattice, species_set, wf_factory, pset_target);
  auto clone = original.spawnCrowdClone();
  REQUIRE(clone != nullptr);
  REQUIRE(clone.get() != &original);
  REQUIRE(dynamic_cast<decltype(&original)>(clone.get()) != nullptr);
}

TEST_CASE("OneBodyDensityMatrices::accumulate", "[estimators]")
{
  using namespace testing;
  using namespace onebodydensitymatrices;
  using MCPWalker = OperatorEstBase::MCPWalker;

  Communicate* comm;
  comm = OHMMS::Controller;
  outputManager.pause();

  Libxml2Document doc;
  bool okay = doc.parseFromString(valid_one_body_density_matrices_input_sections[valid_obdm_input]);
  if (!okay)
    throw std::runtime_error("cannot parse OneBodyDensitMatricesInput section");
  xmlNodePtr node = doc.getRoot();
  OneBodyDensityMatricesInput obdmi(node);
  MinimalParticlePool mpp;
  ParticleSetPool particle_pool = mpp(comm);
  MinimalWaveFunctionPool wfp;
  WaveFunctionPool wavefunction_pool = wfp(comm, particle_pool);
  auto& wf_factory                   = *(wavefunction_pool.getWaveFunctionFactory("wavefunction"));
  wavefunction_pool.setPrimary(wavefunction_pool.getWaveFunction("psi0"));
  auto& pset_target = *(particle_pool.getParticleSet("e"));
  auto& pset_source = *(particle_pool.getParticleSet("ion"));
  auto& species_set = pset_target.getSpeciesSet();
  OneBodyDensityMatrices obdm(std::move(obdmi), pset_target.Lattice, species_set, wf_factory, pset_target);

  std::vector<MCPWalker> walkers;
  int nwalkers = 3;
  for (int iw = 0; iw < nwalkers; ++iw)
    walkers.emplace_back(8);

  std::vector<ParticleSet::ParticlePos_t> deterministic_rs = {
      {
          {0.4538051387, -1.057858504, -0.2663974283},
          {-0.01056867486, -1.705679936, 0.4891466675},
          {1.327399178, 0.9107993414, 1.318479613},
          {2.957728881, 2.394986509, 0.08985640442},
          {-0.03391919965, -0.7795750272, 1.248553217},
          {0.6829973971, 0.8692599761, -1.304274991},
          {2.169392989, 0.3343588146, 0.697054544},
          {2.115217961, -0.07220208998, 2.301371666},
      },
      {
          {-1.391095985, -1.367703815, 1.270925225},
          {-1.364643736, 1.164920615, 0.1394268939},
          {2.738471999, 0.5792334772, 2.86064708},
          {3.394945867, 0.9643879042, 2.405711727},
          {0.8161897766, 0.361988661, 0.3694065569},
          {-0.8830466952, 0.04930234077, -0.05773275004},
          {0.9373587674, 2.687607336, 2.002084488},
          {1.055408731, 2.758634075, 0.1239970225},
      },
      {
          {-1.362662151, -0.2322816531, -1.300150266},
          {0.9228481455, 0.06853710039, -1.244370413},
          {1.213026597, 0.5986030429, 2.578909128},
          {2.621095431, 2.138690329, 2.577521166},
          {0.2229965821, -1.282158307, -0.4278759218},
          {0.7976909814, -0.6759093043, 0.8356683624},
          {3.671887703, 1.783684111, 0.9190620094},
          {-0.206488308, -0.5605483407, -0.957887542},
      },
  };
  if constexpr (generate_test_data)
    std::cout << "Initialize OneBodyDensityMatrices::accumulate psets with:\n{";
  std::vector<ParticleSet> psets;
  for (int iw = 0; iw < nwalkers; ++iw)
  {
    psets.emplace_back(pset_target);
    if constexpr (generate_test_data)
    {
      psets.back().randomizeFromSource(pset_source);
      std::cout << "{";
      for (auto r : psets.back().R)
        std::cout << NativePrint(r) << ",";
      std::cout << "},\n";
    }
    else
      psets.back().R = deterministic_rs[iw];
  }
  if constexpr (generate_test_data)
    std::cout << "}\n";

  auto& trial_wavefunction = *(wavefunction_pool.getPrimary());
  std::vector<UPtr<TrialWaveFunction>> twfcs(nwalkers);
  for (int iw = 0; iw < nwalkers; ++iw)
    twfcs[iw] = trial_wavefunction.makeClone(psets[iw]);

  StdRandom<double> rng;
  rng.init(0, 1, 101);

  auto updateWalker = [](auto& walker, auto& pset_target, auto& trial_wavefunction) {
    pset_target.update(true);
    pset_target.donePbyP();
    trial_wavefunction.evaluateLog(pset_target);
    pset_target.saveWalker(walker);
  };

  for (int iw = 0; iw < nwalkers; ++iw)
    updateWalker(walkers[iw], psets[iw], *(twfcs[iw]));

  auto ref_walkers(makeRefVector<MCPWalker>(walkers));
  auto ref_psets(makeRefVector<ParticleSet>(psets));
  auto ref_twfcs(convertUPtrToRefVector(twfcs));

  OneBodyDensityMatricesTests<double> obdmt;
  obdmt.testAccumulate(obdm, ref_walkers, ref_psets, ref_twfcs, rng);

  if constexpr (generate_test_data)
    obdmt.dumpData(obdm);

  outputManager.resume();
}

TEST_CASE("OneBodyDensityMatrices::evaluateMatrix", "[estimators]")
{
  using namespace testing;
  using namespace onebodydensitymatrices;
  using MCPWalker = OperatorEstBase::MCPWalker;
  using namespace onebodydensitymatrices;

  Communicate* comm;
  comm = OHMMS::Controller;
  outputManager.pause();

  for (auto valid_integrator : std::vector<int>{valid_obdm_input, valid_obdm_input_scale, valid_obdm_input_grid})
  {
    Libxml2Document doc;
    bool okay = doc.parseFromString(valid_one_body_density_matrices_input_sections[valid_integrator]);
    if (!okay)
      throw std::runtime_error("cannot parse OneBodyDensitMatricesInput section");
    xmlNodePtr node = doc.getRoot();
    OneBodyDensityMatricesInput obdmi(node);

    std::string integrator_str =
        InputSection::reverseLookupInputEnumMap(obdmi.get_integrator(), OBDMI::lookup_input_enum_value);
    std::cout << "Test evaluateMatrix for: " << integrator_str << '\n';

    MinimalParticlePool mpp;
    ParticleSetPool particle_pool = mpp(comm);
    MinimalWaveFunctionPool wfp;
    WaveFunctionPool wavefunction_pool = wfp(comm, particle_pool);
    auto& wf_factory                   = *(wavefunction_pool.getWaveFunctionFactory("wavefunction"));
    wavefunction_pool.setPrimary(wavefunction_pool.getWaveFunction("psi0"));
    auto& pset_target = *(particle_pool.getParticleSet("e"));
    if constexpr (generate_test_data)
    {
      std::cout << "Initialize pset_target.R with the following:\n{";
      for (auto r : pset_target.R)
        std::cout << NativePrint(r) << ",";
      std::cout << "}\n";
    }
    auto& species_set = pset_target.getSpeciesSet();
    OneBodyDensityMatrices obdm(std::move(obdmi), pset_target.Lattice, species_set, wf_factory, pset_target);
    auto& trial_wavefunction = *(wavefunction_pool.getPrimary());

    // We can't reason about the state of the global Random in tests. A User can run only some tests,
    // new tests will get added, other tests modified so global Random is called more times or fewer.
    // Also due to use of FakeRandom in unit tests in other tests of this executable its difficult
    // to know which global Random this test will have have access to. So trying to initialize it to
    // a known state is not maintainable.
    // So we must initialize particle positions to known values.
    pset_target.R =
        ParticleSet::ParticlePos_t{{1.751870349, 4.381521229, 2.865202269}, {3.244515371, 4.382273176, 4.21105285},
                                   {3.000459944, 3.329603408, 4.265030556}, {3.748660329, 3.63420622, 5.393637791},
                                   {3.033228526, 3.391869137, 4.654413566}, {3.114198787, 2.654334594, 5.231075822},
                                   {3.657151589, 4.883870516, 4.201243939}, {2.97317591, 4.245644974, 4.284564732}};

    StdRandom<double> rng;
    rng.init(0, 1, 101);
    MCPWalker walker;
    // Now we have to bring the pset, trial_wavefunction and walker to valid state.
    //pset.loadWalker(walker, false);
    pset_target.update(true);
    pset_target.donePbyP();
    trial_wavefunction.evaluateLog(pset_target);
    pset_target.saveWalker(walker);
    OneBodyDensityMatricesTests<double> obdmt;
    obdmt.testEvaluateMatrix(obdm, pset_target, trial_wavefunction, walker, rng);
    // You can use this to regenerate the test data
    if constexpr (generate_test_data)
      obdmt.dumpData(obdm);
  }
  outputManager.resume();
}

namespace testing
{
// The test result data is defined down here for readability of the test code.
template<typename T>
typename OneBodyDensityMatricesTests<T>::Data OneBodyDensityMatricesTests<T>::getEvaluateMatrixData(
    OBDMI::Integrator integrator)
{
  Data data;
  switch (integrator)
  {
  case OBDMI::Integrator::UNIFORM_GRID: {
    if constexpr (IsComplex_t<OneBodyDensityMatrices::Value>::value)
    {
      if constexpr (std::is_same<OneBodyDensityMatrices::Real, double>::value)
        data = {
            0.8479310253,     1.110223025e-16,
            -0.003246774574,  -0.001925348328,
            -0.01697761665,   -0.0003681976742,
            -0.1742565222,    3.700360712e-10,
            0.1992540403,     4.606063586e-10,
            -0.004738188201,  -0.006972389413,
            0.006377855498,   -0.005179873185,
            0.2403578726,     -2.081668171e-16,
            -0.003246774574,  0.001925348328,
            0.6491139457,     -1.110223025e-16,
            0.0008416059524,  -0.0009537904934,
            0.000469580579,   -0.0003005351381,
            -0.001166491073,  0.0006185243955,
            0.01544061242,    -0.02985155826,
            -0.02589818355,   -0.02743137999,
            -0.0008422855246, 0.0004561738209,
            -0.01697761665,   0.0003681976742,
            0.0008416059524,  0.0009537904934,
            0.6574162459,     0,
            0.00265009784,    -5.325351768e-05,
            -0.005454136903,  0.0001322819456,
            -0.02584983289,   -0.02361723534,
            0.02712753804,    -0.01330769562,
            -0.004022521874,  9.551741183e-05,
            -0.1742565222,    -3.700359602e-10,
            0.000469580579,   0.0003005351381,
            0.00265009784,    5.325351768e-05,
            0.6259294634,     -1.665334537e-16,
            -0.3056315893,    -9.288864122e-11,
            -0.0002500001889, -0.0004682526462,
            0.0003365028092,  -0.0003478829106,
            -0.1082773831,    -8.841181259e-11,
            -0.0878135962,    -4.743097312e-08,
            0.8701071598,     0.02319774265,
            0.2045565786,     0.09867468728,
            -0.197541384,     1.012967707e-08,
            0.07482205604,    -1.661263613e-08,
            0.01827034258,    -0.04928728352,
            -0.02459283283,   -0.03661617618,
            -0.003082208891,  -1.476835605e-08,
            -0.754193754,     0.2240498756,
            -0.3115042983,    0.6811980058,
            -0.2233958458,    -0.4031699305,
            0.3782798955,     -0.08517944449,
            -0.2766538428,    0.07068578771,
            0.0242920127,     0.05942351867,
            0.0195332263,     -0.01991019668,
            -0.2362713493,    0.06708718283,
            1.015184,         0.1664496483,
            0.5053365411,     0.5402530165,
            0.1354505239,     -0.2159682122,
            -0.5091844144,    -0.06328095235,
            0.3723904607,     0.05251341435,
            0.01953323876,    -0.01593043434,
            0.01251077105,    -0.05942350205,
            0.3180335199,     0.04983996418,
            -0.2198624567,    -5.870304243e-15,
            -2.088917198,     0.001959818254,
            0.01728161399,    -0.2368940956,
            0.6639046983,     2.079009037e-08,
            -0.3276186073,    1.264326699e-10,
            -0.05888980528,   0.1049113351,
            0.07926872673,    0.07794001871,
            -0.1242315553,    -1.630640067e-15,
            0.8479310253,     -8.881784197e-16,
            -0.003246774574,  -0.001925348328,
            -0.01697761665,   -0.0003681976742,
            -0.1742565222,    3.700351137e-10,
            0.1992540403,     4.606060255e-10,
            -0.004738188201,  -0.006972389413,
            0.006377855498,   -0.005179873185,
            0.2403578726,     -2.775557562e-17,
            -0.003246774574,  0.001925348328,
            0.6491139457,     7.771561172e-16,
            0.0008416059524,  -0.0009537904934,
            0.000469580579,   -0.0003005351381,
            -0.001166491073,  0.0006185243955,
            0.01544061242,    -0.02985155826,
            -0.02589818355,   -0.02743137999,
            -0.0008422855246, 0.0004561738209,
            -0.01697761665,   0.0003681976742,
            0.0008416059524,  0.0009537904934,
            0.6574162459,     3.330669074e-16,
            0.00265009784,    -5.325351768e-05,
            -0.005454136903,  0.0001322819456,
            -0.02584983289,   -0.02361723534,
            0.02712753804,    -0.01330769562,
            -0.004022521874,  9.551741183e-05,
            -0.1742565222,    -3.700360435e-10,
            0.000469580579,   0.0003005351381,
            0.00265009784,    5.325351768e-05,
            0.6259294634,     0,
            -0.3056315893,    -9.288872449e-11,
            -0.0002500001889, -0.0004682526462,
            0.0003365028092,  -0.0003478829106,
            -0.1082773831,    -8.841196525e-11,
            -2.25611399,      -5.094168354e-08,
            -1.496397952,     0.04197205755,
            0.3701067606,     -0.1696992504,
            -3.028538005,     -5.011290716e-08,
            1.036419558,      1.770635139e-08,
            -0.03321625225,   0.09519291876,
            0.04471084154,    0.07072013214,
            -0.2908638254,    -7.614071845e-09,
            -0.683381127,     0.4755374123,
            0.2033652126,     0.4686325908,
            -0.2001582441,    -0.4364641008,
            -0.3388052931,    1.352856329,
            0.05470563214,    -0.5379968647,
            0.03288396914,    0.02810335112,
            -0.005775626092,  -0.02746264737,
            -0.1460550888,    -0.01037877896,
            0.9198664375,     0.3532831852,
            -0.1662607588,    0.3825385354,
            0.1545051752,     -0.3602434421,
            0.4560494115,     1.005055153,
            -0.07363658276,   -0.3996850954,
            -0.005775647383,  0.01051241664,
            0.03636750456,    -0.02810336883,
            0.1965977315,     -0.007710521785,
            1.021599858,      -7.327471963e-15,
            0.1415421361,     0.023055409,
            0.2033013166,     0.01605157874,
            1.821121796,      -5.46969614e-09,
            -0.6716770369,    1.612600609e-09,
            -0.01362750517,   -0.02990320627,
            0.01834330404,    -0.02221551834,
            0.08718674567,    -1.540434447e-15,
        };
      else if constexpr (std::is_same<OneBodyDensityMatrices::Real, float>::value)
        data = {};
    }
    else if constexpr (std::is_floating_point<OneBodyDensityMatrices::Value>::value)
    {
      if constexpr (std::is_same<OneBodyDensityMatrices::Real, double>::value)
        data = {};
      else if constexpr (std::is_same<OneBodyDensityMatrices::Real, float>::value)
        data = {};
    }
    break;
  }
  case OBDMI::Integrator::UNIFORM: {
    if constexpr (IsComplex_t<OneBodyDensityMatrices::Value>::value)
    {
      if constexpr (std::is_same<OneBodyDensityMatrices::Real, double>::value)
        data = {
            0.8207296586,   1.665334537e-16,  0.07548328902,   -0.01020249039,   -0.08996493371,  0.008560200457,
            -0.1317480093,  1.023118146e-09,  0.2491948091,    -3.494558712e-09, 0.02968127153,   0.06452312383,
            -0.03995252482, 0.04793510042,    0.2173928965,    5.551115123e-17,  0.07548328902,   0.01020249039,
            0.6335421171,   -5.551115123e-17, -0.03366469844,  0.0132172667,     -0.05322681298,  -0.001350193778,
            0.08303390974,  0.009749756098,   -0.005002151922, -0.1452516051,    -0.01165302175,  -0.1021548287,
            -0.02561569231, 0.01308015249,    -0.08996493371,  -0.008560200457,  -0.03366469844,  -0.0132172667,
            0.5184920036,   5.551115123e-17,  0.01190593785,   0.006036202016,   -0.08597273722,  -0.009416483968,
            0.008131097685, -0.07804258925,   -0.04470233214,  -0.05707880237,   -0.1153399797,   0.002904949865,
            -0.1317480093,  -1.02311859e-09,  -0.05322681298,  0.001350193778,   0.01190593785,   -0.006036202016,
            0.7212254051,   1.110223025e-16,  -0.3727674491,   -3.194204901e-10, -0.05125175435,  0.07631675717,
            0.06898751753,  0.056696722,      -0.1492371305,   -2.549892672e-09, -0.01459098526,  -5.016695612e-08,
            0.8472382505,   0.01406996187,    0.124068409,     0.09608124191,    -0.3093600708,   1.526067131e-08,
            0.1802012276,   -2.054179687e-08, 0.02641734172,   -0.2560471361,    -0.03555912291,  -0.1902208237,
            -0.05249722171, -1.32974614e-08,  -0.7416709929,   0.3409016998,     -0.3827033145,   0.720966957,
            -0.1103351807,  -0.389630854,     0.4069930313,    -0.14487758,      -0.3699170719,   0.2237867371,
            0.05006187413,  0.07533128018,    0.1452226817,    0.05318051365,    -0.1919875344,   0.09889887838,
            0.9983277158,   0.253260426,      0.5961189214,    0.5422991188,     -0.02303775641,  -0.1896161183,
            -0.5478337909,  -0.1076314866,    0.4979276119,    0.1662541471,     0.145222696,     -0.09861573644,
            -0.03752754636, -0.07533125964,   0.2584252008,    0.07347329964,    -0.4323912902,   -5.440092821e-15,
            -2.105461521,   0.01572812901,    0.1386897412,    -0.2387703133,    0.9366795805,    1.854620557e-08,
            -0.6773405088,  -1.211485423e-08, -0.06712648164,  0.5390751007,     0.09035570759,   0.4004860643,
            -0.08784656675, -1.262878691e-15, 0.8207296586,    -1.276756478e-15, 0.07548328902,   -0.01020249039,
            -0.08996493371, 0.008560200457,   -0.1317480093,   1.023115856e-09,  0.2491948091,    -3.494558268e-09,
            0.02968127153,  0.06452312383,    -0.03995252482,  0.04793510042,    0.2173928965,    7.077671782e-16,
            0.07548328902,  0.01020249039,    0.6335421171,    -1.110223025e-16, -0.03366469844,  0.0132172667,
            -0.05322681298, -0.001350193778,  0.08303390974,   0.009749756098,   -0.005002151922, -0.1452516051,
            -0.01165302175, -0.1021548287,    -0.02561569231,  0.01308015249,    -0.08996493371,  -0.008560200457,
            -0.03366469844, -0.0132172667,    0.5184920036,    3.330669074e-16,  0.01190593785,   0.006036202016,
            -0.08597273722, -0.009416483968,  0.008131097685,  -0.07804258925,   -0.04470233214,  -0.05707880237,
            -0.1153399797,  0.002904949865,   -0.1317480093,   -1.023118479e-09, -0.05322681298,  0.001350193778,
            0.01190593785,  -0.006036202016,  0.7212254051,    1.665334537e-16,  -0.3727674491,   -3.194207399e-10,
            -0.05125175435, 0.07631675717,    0.06898751753,   0.056696722,      -0.1492371305,   -2.549892783e-09,
            -2.623553571,   -5.482606902e-08, -1.468262591,    0.06912332787,    0.6095249145,    -0.1665085377,
            -3.626402182,   -6.917701145e-08, 1.002242017,     2.240606978e-08,  0.1913697865,    -0.4051608567,
            -0.2575936146,  -0.3009993784,    0.04358624391,   -1.440615186e-09, -0.6550296498,   0.6610480811,
            0.1736058217,   0.4426971401,     -0.1060672518,   -0.4251930577,    -0.4725422331,   1.584683643,
            0.09989718465,  -0.5312309328,    -0.192516785,    -0.247450456,     -0.1587886145,   0.04242197737,
            -0.06872100423, -0.07817777616,   0.8817038607,    0.4911015533,     -0.1296286188,   0.3419515806,
            0.03641562731,  -0.3452085377,    0.6360662621,    1.177282788,      -0.1344667649,   -0.3946586024,
            -0.1587886486,  0.1068249126,     -0.09674523302,  0.2474504336,     0.09250204624,   -0.05807931571,
            1.10814767,     -8.54871729e-15,  0.09283970382,   0.00440416402,    0.038835841,     0.01052846149,
            2.221793506,    -1.71057013e-09,  -0.8158278984,   -3.590288289e-09, -0.1100928074,   0.3270731201,
            0.1481906104,   0.242986956,      -0.1424011075,   1.151856388e-15,
        };
      else if constexpr (std::is_same<OneBodyDensityMatrices::Real, float>::value)
        data = {};
    }
    else if constexpr (std::is_floating_point<OneBodyDensityMatrices::Value>::value)
    {
      if constexpr (std::is_same<OneBodyDensityMatrices::Real, double>::value)
        data = {};
      else if constexpr (std::is_same<OneBodyDensityMatrices::Real, float>::value)
        data = {};
    }
    break;
  }
  case OBDMI::Integrator::DENSITY: {
    if constexpr (IsComplex_t<OneBodyDensityMatrices::Value>::value)
    {
      if constexpr (std::is_same<OneBodyDensityMatrices::Real, double>::value)
        data = {0.9972842135,   2.775557562e-16,  -0.1509463392,  0.004894026847,   0.04315523355,  -0.01711810294,
                0.1232433221,   6.700087429e-10,  0.1927144236,   6.442509581e-10,  -0.094787711,   0.1537809336,
                0.1275891946,   0.114245917,      0.009762182978, 1.769417945e-16,  -0.1509463392,  -0.004894026847,
                1.167677748,    -4.440892099e-16, 0.05516205268,  0.03235550535,    0.1969117701,   -0.008414514051,
                0.01633315462,  -0.007457786918,  -0.02730020562, -0.2330227348,    0.03183169144,  -0.162739637,
                -0.2566088424,  0.005950756757,   0.04315523355,  0.01711810294,    0.05516205268,  -0.03235550535,
                0.8860381802,   -2.775557562e-16, 0.07419862606,  -0.02233081948,   0.06576238506,  -0.001852263199,
                0.01793673063,  -0.01792147225,   -0.07817004956, -0.01922402746,   -0.05247343171, 0.02910077141,
                0.1232433221,   -6.700090205e-10, 0.1969117701,   0.008414514051,   0.07419862606,  0.02233081948,
                0.9160994045,   -1.110223025e-16, 0.1678893864,   1.051832649e-10,  0.01637708678,  0.01636964028,
                -0.02204439798, 0.01216122985,    -0.3464414664,  -3.63824329e-09,  -0.4029298437,  -3.912557406e-08,
                1.539625298,    0.03517084686,    0.3101348509,   0.1746015219,     -0.06421021074, -1.950993521e-08,
                -0.05079505994, 3.741992265e-09,  -0.01038711951, -0.347553722,     0.0139815873,   -0.2582023181,
                -0.2398699887,  7.46367293e-09,   -0.6968783912,  0.04616429667,    -0.4092305246,  1.152793152,
                -0.3844659898,  -0.4696152905,    0.1178922745,   0.1425202428,     -0.1194995868,  0.01710804859,
                0.2877854559,   -0.06386091967,   0.03221321673,  0.1106168689,     0.0162332681,   -0.2252878362,
                0.9380345297,   0.03429608874,    0.6498300211,   0.915771426,      0.2376849138,   -0.2407116018,
                -0.1586891256,  0.1058801743,     0.1608526338,   0.01270981038,    0.03221320771,  -0.07209989828,
                0.268356413,    0.06386091592,    -0.02185083227, -0.1673693325,    0.5665475714,   -1.076916334e-14,
                -3.55533077,    -0.009126973382,  -0.08048105243, -0.4031930198,    0.3123355945,   3.756725633e-08,
                0.1134356285,   -2.7655428e-08,   0.1049166466,   0.7517269135,     -0.1412232565,  0.5584679678,
                0.4721033136,   -2.498001805e-16, 0.9972842135,   -2.775557562e-16, -0.1509463392,  0.004894026847,
                0.04315523355,  -0.01711810294,   0.1232433221,   6.700072788e-10,  0.1927144236,   6.442505557e-10,
                -0.094787711,   0.1537809336,     0.1275891946,   0.114245917,      0.009762182978, -1.07813064e-15,
                -0.1509463392,  -0.004894026847,  1.167677748,    -7.771561172e-16, 0.05516205268,  0.03235550535,
                0.1969117701,   -0.008414514051,  0.01633315462,  -0.007457786918,  -0.02730020562, -0.2330227348,
                0.03183169144,  -0.162739637,     -0.2566088424,  0.005950756757,   0.04315523355,  0.01711810294,
                0.05516205268,  -0.03235550535,   0.8860381802,   3.885780586e-16,  0.07419862606,  -0.02233081948,
                0.06576238506,  -0.001852263199,  0.01793673063,  -0.01792147225,   -0.07817004956, -0.01922402746,
                -0.05247343171, 0.02910077141,    0.1232433221,   -6.70009194e-10,  0.1969117701,   0.008414514051,
                0.07419862606,  0.02233081948,    0.9160994045,   -1.665334537e-16, 0.1678893864,   1.051833065e-10,
                0.01637708678,  0.01636964028,    -0.02204439798, 0.01216122985,    -0.3464414664,  -3.638242235e-09,
                -4.218460121,   -9.610451324e-08, -3.272413151,   -0.03429277204,   -0.3023918958,  -0.3711085646,
                -6.325229493,   -7.875119135e-08, -1.746291197,   -4.946045018e-08, 0.3508551411,   -0.1669920235,
                -0.4722693032,  -0.1240606884,    2.589688623,    4.144042354e-08,  -1.120194689,   1.2106985,
                0.2804650255,   1.13361394,       -0.4366230486,  -0.2974182405,    -0.837001073,   2.480582466,
                -0.3370383963,  0.5834726525,     0.0197252187,   -0.3202170206,    -0.1163293998,  -0.01093766396,
                0.2250211263,   -1.000648999,     1.507840126,    0.8994442544,     -0.3005177755,  0.9142309287,
                0.3109934929,   -0.2786655311,    1.126646723,    1.842858089,      0.4536711259,   0.4334696902,
                -0.1163293559,  0.2040729096,     0.08988792882,  0.3202170701,     -0.302890033,   -0.7433956089,
                2.319844279,    -1.043609643e-14, 0.6702898076,   0.0742522338,     0.6547518612,   0.07601428408,
                3.460919978,    -1.978514064e-08, 0.9746423386,   -2.257782517e-09, -0.1160181893,  0.292467088,
                0.1561665529,   0.2172777448,     -1.250567834,   8.659739592e-15};
      else if constexpr (std::is_same<OneBodyDensityMatrices::Real, float>::value)
        data = {0.997284174,    0,
                -0.1509462148,  0.004894062877,
                0.04315539822,  -0.01711797714,
                0.1232431382,   5.960464478e-08,
                0.1927144527,   1.490116119e-08,
                -0.09478760511, 0.1537808031,
                0.1275892109,   0.1142460853,
                0.009762163274, -7.450580597e-09,
                -0.1509461701,  -0.004894018173,
                1.167678118,    -5.960464478e-08,
                0.05516195297,  0.03235545754,
                0.1969118416,   -0.008414544165,
                0.01633344032,  -0.007457806263,
                -0.02730023116, -0.2330225706,
                0.03183176368,  -0.1627395749,
                -0.256608963,   0.005950763822,
                0.04315534234,  0.01711807586,
                0.0551616475,   -0.0323554799,
                0.8860384226,   0,
                0.07419875264,  -0.0223308336,
                0.06576254964,  -0.001852300018,
                0.01793673821,  -0.01792119071,
                -0.07817010581, -0.0192239508,
                -0.05247352645, 0.02910077758,
                0.1232429594,   -7.450580597e-09,
                0.1969116628,   0.008414536715,
                0.07419854403,  0.02233078144,
                0.9160988331,   -2.980232239e-08,
                0.1678893715,   1.490116119e-08,
                0.01637715101,  0.01636958495,
                -0.02204445377, 0.012161172,
                -0.3464412391,  0,
                -0.4029290378,  -5.960464478e-08,
                1.539624691,    0.03517085314,
                0.3101349175,   0.1746013612,
                -0.06420990825, -5.587935448e-08,
                -0.05079455674, -1.303851604e-08,
                -0.01038721204, -0.347553134,
                0.01398165524,  -0.2582020462,
                -0.2398701012,  1.490116119e-08,
                -0.6968790293,  0.04616469145,
                -0.409229666,   1.152794003,
                -0.3844661713,  -0.4696149528,
                0.1178922132,   0.142519787,
                -0.1194998473,  0.01710827276,
                0.2877854109,   -0.06386129558,
                0.032213144,    0.1106166169,
                0.01623325795,  -0.2252878547,
                0.9380354881,   0.03429636359,
                0.6498287916,   0.9157721996,
                0.2376853228,   -0.24071154,
                -0.1586889923,  0.1058801115,
                0.1608530283,   0.01271001995,
                0.03221330047,  -0.07209946215,
                0.2683564723,   0.06386158615,
                -0.02185085416, -0.1673694402,
                0.5665459037,   0,
                -3.555330276,   -0.009126901627,
                -0.08048132062, -0.4031928182,
                0.3123348355,   8.940696716e-08,
                0.1134345308,   0,
                0.104916811,    0.7517259121,
                -0.1412234902,  0.5584673882,
                0.4721037149,   -2.980232239e-08,
                0.9972836971,   -8.940696716e-08,
                -0.1509464681,  0.004893258214,
                0.04315529019,  -0.01711768284,
                0.123244673,    3.725290298e-07,
                0.1927143633,   1.11758709e-07,
                -0.09478767961, 0.1537810266,
                0.1275890619,   0.1142454594,
                0.009762742557, -3.073364496e-08,
                -0.1509454846,  -0.004894219339,
                1.167678595,    -9.536743164e-07,
                0.05516173691,  0.03235505521,
                0.1969116032,   -0.008414916694,
                0.01633333229,  -0.007457929663,
                -0.02730023861, -0.2330227196,
                0.03183183074,  -0.1627394408,
                -0.2566090226,  0.005951091647,
                0.04315596819,  0.01711825281,
                0.05516267568,  -0.03235335648,
                0.8860384226,   2.682209015e-07,
                0.07419607788,  -0.02232901752,
                0.06576249003,  -0.001851793379,
                0.01793645881,  -0.01792129315,
                -0.07816983759, -0.01922356337,
                -0.05247297883, 0.0291005224,
                0.1232430413,   -8.195638657e-08,
                0.1969119757,   0.008415028453,
                0.07419873774,  0.02233074792,
                0.9160985947,   1.788139343e-07,
                0.1678895056,   -5.215406418e-08,
                0.01637711562,  0.01636960916,
                -0.0220443625,  0.01216138527,
                -0.3464415669,  2.980232239e-08,
                -4.218452454,   -4.768371582e-07,
                -3.272411823,   -0.0342912674,
                -0.3023903668,  -0.3711089492,
                -6.325219154,   -7.152557373e-07,
                -1.746289253,   -1.788139343e-07,
                0.3508545458,   -0.166991502,
                -0.4722686708,  -0.1240597963,
                2.58968401,     5.960464478e-07,
                -1.120192409,   1.210695028,
                0.2804673016,   1.133612633,
                -0.436622709,   -0.29741925,
                -0.8369976878,  2.480578899,
                -0.3370373845,  0.5834715366,
                0.01972543076,  -0.3202166855,
                -0.1163287833,  -0.01093763486,
                0.225019455,    -1.000647306,
                1.507837296,    0.8994423151,
                -0.3005181253,  0.9142314196,
                0.3109933138,   -0.2786653638,
                1.126644135,    1.842858195,
                0.4536704123,   0.4334697425,
                -0.1163290516,  0.2040731758,
                0.08988789469,  0.3202165067,
                -0.302887857,   -0.7433953285,
                2.31983757,     1.192092896e-07,
                0.6702869534,   0.0742533803,
                0.6547510028,   0.07601451874,
                3.460909367,    1.072883606e-06,
                0.9746402502,   4.470348358e-07,
                -0.1160178259,  0.2924669087,
                0.156166032,    0.2172774523,
                -1.250563502,   -4.768371582e-07};
    }
    else if constexpr (std::is_floating_point<OneBodyDensityMatrices::Value>::value)
    {
      if constexpr (std::is_same<OneBodyDensityMatrices::Real, double>::value)
        data = {0.9965771993,    -0.1276230838,  0.03958306806,  0.1387017217,   0.1942437768,    0.053929644,
                0.2344135141,    -0.0072116162,  -0.1276230838,  1.14757642,     0.2606661124,    0.1992496192,
                0.01161410961,   -0.2376481391,  -0.1358804612,  -0.2716422407,  0.03958306806,   0.2606661124,
                0.8895496478,    0.09026675397,  0.07482099268,  0.03203129787,  -0.09998410562,  -0.06962064713,
                0.1387017217,    0.1992496192,   0.09026675397,  0.9362099992,   0.1647085609,    0.04014883082,
                -0.008667251236, -0.3387070854,  -0.3816205747,  1.526601118,    0.450628534,     -0.08325125513,
                -0.06505223916,  -0.3367568853,  -0.2337969074,  -0.2501181474,  -0.759979096,    -1.598167941,
                0.001566609973,  -0.02491515452, -0.1152966847,  0.381176093,    -0.07186867215,  0.2844624377,
                0.9034968623,    -0.1833555236,  0.6301141723,   -0.2633959431,  0.1582965722,    0.09111738873,
                0.1645013359,    0.1367509408,   0.5272612767,   -3.474323999,   -0.4137162493,   0.3501207451,
                0.153163578,     0.8376243065,   0.387078839,    0.5159687433,   0.9965771993,    -0.1276230838,
                0.03958306806,   0.1387017217,   0.1942437768,   0.053929644,    0.2344135141,    -0.0072116162,
                -0.1276230838,   1.14757642,     0.2606661124,   0.1992496192,   0.01161410961,   -0.2376481391,
                -0.1358804612,   -0.2716422407,  0.03958306806,  0.2606661124,   0.8895496478,    0.09026675397,
                0.07482099268,   0.03203129787,  -0.09998410562, -0.06962064713, 0.1387017217,    0.1992496192,
                0.09026675397,   0.9362099992,   0.1647085609,   0.04014883082,  -0.008667251236, -0.3387070854,
                -4.341682703,    -3.281905856,   -0.63616415,    -6.494174955,   -1.698130443,    0.157715294,
                -0.6031292071,   2.641093171,    -2.383983684,   -0.9329968953,  -0.08113582861,  -3.414342806,
                -0.9024677642,   -0.08564081593, -0.4186924916,  1.246196012,    0.5913805452,    -1.098837966,
                0.5427940957,    -0.7226756762,  0.04220981851,  0.2642804489,   0.1699938682,    0.4461506245,
                2.379646766,     0.7448243926,   0.7276662244,   3.55662162,     0.9666690056,    0.2069702368,
                0.3616379717,    -1.254351175};
      else if constexpr (std::is_same<OneBodyDensityMatrices::Real, float>::value)
        data = {0.9965772033,    -0.1276224554,   0.03958324343,  0.138701871,    0.194243744,     0.05392972752,
                0.2344133556,    -0.007211369928, -0.1276228428,  1.147576571,    0.2606660724,    0.1992495656,
                0.01161403582,   -0.2376479805,   -0.1358801872,  -0.2716423869,  0.03958294168,   0.2606661916,
                0.8895497322,    0.09026694298,   0.07482103258,  0.03203130513,  -0.0999841243,   -0.06962074339,
                0.1387016773,    0.1992500126,    0.09026675671,  0.9362098575,   0.1647084951,    0.04014879465,
                -0.008667248301, -0.3387069404,   -0.3816198409,  1.526600122,    0.4506285191,    -0.08325134218,
                -0.06505221874,  -0.336756438,    -0.233796373,   -0.2501182556,  -0.7599802017,   -1.598167896,
                0.001566099701,  -0.0249146726,   -0.1152965948,  0.3811755478,   -0.07186914235,  0.2844621241,
                0.9034972191,    -0.1833569407,   0.6301141381,   -0.2633955181,  0.1582967192,    0.09111790359,
                0.1645013839,    0.1367513388,    0.5272595286,   -3.474322319,   -0.4137164652,   0.3501208723,
                0.1531635821,    0.8376233578,    0.3870776892,   0.5159689784,   0.9965775609,    -0.1276229024,
                0.03958255798,   0.1387042105,    0.1942443401,   0.05392966419,  0.234413594,     -0.007211854216,
                -0.1276231557,   1.147577047,     0.260666281,    0.1992495805,   0.01161361579,   -0.2376479208,
                -0.1358803362,   -0.2716422677,   0.03958233446,  0.2606659532,   0.8895499706,    0.09026726335,
                0.07482092828,   0.03203126043,   -0.09998448938, -0.06961926818, 0.1387016624,    0.1992497295,
                0.09026705474,   0.9362094998,    0.1647085547,   0.04014874622,  -0.008667317219, -0.3387072086,
                -4.341678143,    -3.281904936,    -0.6361619234,  -6.494166851,   -1.698127627,    0.1577153355,
                -0.6031289101,   2.641089678,     -2.383980513,   -0.932995379,   -0.081134215,    -3.414337158,
                -0.9024663568,   -0.08564066887,  -0.4186921716,  1.246193886,    0.5913794041,    -1.098839045,
                0.5427934527,    -0.722673595,    0.04221029207,  0.2642802894,   0.1699934751,    0.4461522698,
                2.379641771,     0.74482131,      0.7276645899,   3.556614637,    0.9666671157,    0.2069700211,
                0.3616372049,    -1.254347205};
    }
    break;
  }
  }
  return data;
}

template<typename T>
typename OneBodyDensityMatricesTests<T>::Data OneBodyDensityMatricesTests<T>::getAccumulateData()
{
  Data data;
  if constexpr (IsComplex_t<OneBodyDensityMatrices::Value>::value)
  {
    if constexpr (std::is_same<OneBodyDensityMatrices::Real, double>::value)
      data = {
          2.908407079,    -3.330669074e-16, -0.6090936877,  -0.03566623534,   -0.3145026132,  -0.06907439849,
          0.210018867,    1.240839002e-08,  0.1558723712,   1.16299807e-08,   -0.2592098517,  0.08285488702,
          0.348909943,    0.06155400836,    -0.2766418068,  8.101158633e-16,  -0.6090936877,  0.03566623534,
          3.12886181,     -1.110223025e-15, 0.1256129134,   0.03323621594,    0.2834652262,   -0.01422030447,
          0.3349894174,   0.01710927639,    -0.03193561362, -0.03601825785,   0.07626525472,  -0.02949074837,
          -0.08891826234, -0.02205348601,   -0.3145026132,  0.06907439849,    0.1256129134,   -0.03323621594,
          2.839556058,    1.665334537e-16,  0.1253937038,   -0.03214643392,   -0.1508684612,  -0.03798953922,
          -0.01576665987, 0.1459023256,     0.01238018024,  0.09705332889,    0.1944662261,   0.0100837862,
          0.210018867,    -1.240839329e-08, 0.2834652262,   0.01422030447,    0.1253937038,   0.03214643392,
          2.803045564,    -3.330669074e-16, 0.3900011161,   1.80151058e-08,   -0.2131650988,  -0.2076116,
          0.2869313158,   -0.1542374131,    -0.3794457049,  2.962675631e-09,  -7.023068491,   -4.718818206e-07,
          -6.460503203,   -1.499620153,     -13.22356032,   -0.732654863,     -13.02276263,   -4.262060947e-07,
          2.415996444,    -5.804901433e-08, -0.956212458,   2.72789813,       1.287111704,    2.026591678,
          3.560345445,    -5.395647001e-08, -7.516621334,   2.204375057,      4.637228779,    -0.2603328162,
          -5.028529083,   5.980200315,      -9.638443217,   7.929036465,      -0.65885009,    -0.1899233672,
          3.549012354,    -1.282437437,     -2.648159502,   1.938622271,      -1.167170184,   -0.3187049961,
          10.11776283,    1.637659951,      -7.51766905,    0.9896308454,     6.696171718,    3.199540378,
          12.9738454,     5.890587508,      0.8868463511,   -0.1410964105,    -2.648159412,   -1.165134593,
          5.146219825,    1.282437495,      1.571071889,    -0.2367703412,    -8.212154617,   -9.168078623e-14,
          -5.965144505,   -1.410723553,     -12.439673,     -0.6764782201,    -10.74364679,   3.074102123e-07,
          2.773701782,    -5.531456067e-08, -0.3160472,     2.479063913,      0.4254160718,   1.8417298,
          3.351835501,    1.625435897e-14,  2.908407079,    -5.551115123e-17, -0.6090936877,  -0.03566623534,
          -0.3145026132,  -0.06907439849,   0.210018867,    1.240839311e-08,  0.1558723712,   1.162997969e-08,
          -0.2592098517,  0.08285488702,    0.348909943,    0.06155400836,    -0.2766418068,  1.179611964e-16,
          -0.6090936877,  0.03566623534,    3.12886181,     -1.110223025e-15, 0.1256129134,   0.03323621594,
          0.2834652262,   -0.01422030447,   0.3349894174,   0.01710927639,    -0.03193561362, -0.03601825785,
          0.07626525472,  -0.02949074837,   -0.08891826234, -0.02205348601,   -0.3145026132,  0.06907439849,
          0.1256129134,   -0.03323621594,   2.839556058,    3.885780586e-16,  0.1253937038,   -0.03214643392,
          -0.1508684612,  -0.03798953922,   -0.01576665987, 0.1459023256,     0.01238018024,  0.09705332889,
          0.1944662261,   0.0100837862,     0.210018867,    -1.240839284e-08, 0.2834652262,   0.01422030447,
          0.1253937038,   0.03214643392,    2.803045564,    -1.110223025e-16, 0.3900011161,   1.801510592e-08,
          -0.2131650988,  -0.2076116,       0.2869313158,   -0.1542374131,    -0.3794457049,  2.96267581e-09,
          0.3941080137,   -4.726492245e-08, 2.052722417,    -0.001451978485,  -0.01280291539, 0.2327894166,
          1.16436472,     -4.341952076e-08, -0.1468177385,  2.13771334e-10,   -0.3629857602,  -0.7030273922,
          0.488597706,    -0.5222884379,    -0.8341649642,  1.331856431e-08,  -1.012439432,   -0.26066536,
          -0.264523576,   1.027857623,      0.06141787593,  -0.8607775579,    -0.8968882375,  -0.8798400742,
          -0.5714652082,  0.308993155,      -0.1661014632,  -0.05876212861,   0.04677985043,  0.04953921946,
          0.3514125982,   0.1961046442,     1.362796128,    -0.1936517929,    0.5503747678,   0.7219805028,
          -0.3214497861,  -0.5547816418,    1.20725823,     -0.6536449884,    0.7692218753,   0.2295551372,
          0.04677984276,  -0.01409753328,   -0.1943162002,  0.0587621284,     -0.4730196602,  0.1456887927,
          -4.17192743,    -3.802513859e-15, 2.327346242,    0.3958794534,     3.490840243,    0.2639331793,
          -0.1189745691,  -9.287628714e-08, 0.5912812019,   -4.039935734e-08, 0.2471126461,   0.9240454459,
          -0.3326264459,  0.686485726,      0.7926337257,   2.248201625e-15,
      };
    else if constexpr (std::is_same<OneBodyDensityMatrices::Real, float>::value)
      data = {
          2.908409595,   -2.384185791e-07, -0.6090942621,  -0.03566826135,
          -0.3145041466, -0.06907411665,   0.2100194395,   3.757886589e-07,
          0.1558737606,  -1.639127731e-07, -0.2592100501,  0.08285357803,
          0.348911345,   0.06155376136,    -0.2766413093,  -1.303851604e-07,
          -0.6090935469, 0.03566679358,    3.128861427,    -1.221895218e-06,
          0.1256155074,  0.0332358852,     0.2834667563,   -0.01422011852,
          0.3349895775,  0.0171098616,     -0.03193645179, -0.03601863235,
          0.07626512647, -0.02949096262,   -0.0889185071,  -0.02205353603,
          -0.3145017922, 0.06907459348,    0.1256124973,   -0.0332384482,
          2.839557648,   -4.172325134e-07, 0.1253938824,   -0.03214593604,
          -0.150870502,  -0.03798932582,   -0.0157662034,  0.1459026039,
          0.01237866282, 0.09705364704,    0.1944674253,   0.01008377969,
          0.2100191712,  -1.536682248e-08, 0.28346771,     0.01422321983,
          0.1253924817,  0.0321469456,     2.803043365,    4.470348358e-07,
          0.3899998665,  1.713633537e-07,  -0.2131657749,  -0.2076116204,
          0.286931932,   -0.154237777,     -0.3794481158,  -1.373700798e-07,
          -7.023084641,  -2.801418304e-06, -6.46050787,    -1.499619961,
          -13.22356701,  -0.7326583862,    -13.02278328,   -2.726912498e-06,
          2.415994406,   2.086162567e-07,  -0.9562116265,  2.727905273,
          1.287110806,   2.026596308,      3.56035161,     6.258487701e-07,
          -7.51664114,   2.204387665,      4.637217999,    -0.2603224516,
          -5.028517723,  5.980193138,      -9.638428688,   7.929022312,
          -0.6588373184, -0.1899331212,    3.549021244,    -1.282437444,
          -2.648169279,  1.938628912,      -1.16714859,    -0.318720758,
          10.11778641,   1.637673616,      -7.51765728,    0.9896341562,
          6.696152687,   3.199539423,      12.97381973,    5.890585899,
          0.886828661,   -0.141102612,     -2.64816761,    -1.165143013,
          5.146233559,   1.282441139,      1.571043015,    -0.2367825508,
          -8.212172508,  -3.101537004e-06, -5.965149879,   -1.410723448,
          -12.43968201,  -0.6764817238,    -10.7436657,    -1.817941666e-06,
          2.773700714,   6.351619959e-07,  -0.316046685,   2.479072571,
          0.4254153073,  1.841736436,      3.351844311,    3.483146429e-07,
          2.908406734,   -5.960464478e-08, -0.6090934277,  -0.03566631675,
          -0.3145028353, -0.0690741986,    0.2100202143,   4.516914487e-08,
          0.1558735073,  2.235174179e-08,  -0.2592100203,  0.08285429329,
          0.3489104509,  0.0615535602,     -0.2766424119,  2.793967724e-09,
          -0.6090936065, 0.03566636145,    3.12886095,     -4.172325134e-07,
          0.1256127954,  0.0332361497,     0.2834657431,   -0.01422024984,
          0.3349898756,  0.01710940897,    -0.03193625063, -0.03601834923,
          0.07626627386, -0.02949080616,   -0.08891823888, -0.02205368131,
          -0.3145031035, 0.06907439232,    0.1256127506,   -0.03323595226,
          2.839556694,   5.960464478e-08,  0.1253929287,   -0.03214644641,
          -0.1508699208, -0.03798954934,   -0.01576604694, 0.1459029913,
          0.01237925887, 0.09705369174,    0.1944676042,   0.01008378342,
          0.2100204527,  3.166496754e-08,  0.2834661603,   0.01422035228,
          0.1253925115,  0.03214619309,    2.803045273,    0,
          0.3900001645,  8.195638657e-08,  -0.2131664902,  -0.2076121569,
          0.2869332433,  -0.1542378664,    -0.3794475496,  1.396983862e-09,
          0.3941099048,  9.685754776e-08,  2.052723169,    -0.001452207565,
          -0.0128043294, 0.2327890694,     1.16436553,     -1.564621925e-07,
          -0.1468167156, -7.078051567e-08, -0.3629867733,  -0.7030285001,
          0.488599062,   -0.5222893953,    -0.834166646,   5.960464478e-08,
          -1.01243937,   -0.2606658936,    -0.264524281,   1.027858019,
          0.06141793728, -0.8607769012,    -0.8968876004,  -0.8798406124,
          -0.5714646578, 0.3089928627,     -0.1661020517,  -0.05876155943,
          0.04677910358, 0.04953817278,    0.3514130116,   0.1961062402,
          1.362796068,   -0.1936522126,    0.5503755212,   0.7219808102,
          -0.3214501143, -0.5547809601,    1.207257509,    -0.6536455154,
          0.769221127,   0.2295548022,     0.04677911103,  -0.01409687102,
          -0.1943161637, 0.05876130611,    -0.4730205536,  0.1456899792,
          -4.171929359,  -1.192092896e-07, 2.32734561,     0.3958798349,
          3.490841866,   0.2639330924,     -0.1189776659,  -9.685754776e-08,
          0.5912772417,  -3.725290298e-08, 0.2471138537,   0.9240474701,
          -0.3326281905, 0.6864871979,     0.7926388383,   -3.725290298e-08,
      };
  }
  else if constexpr (std::is_floating_point<OneBodyDensityMatrices::Value>::value)
  {
    if constexpr (std::is_same<OneBodyDensityMatrices::Real, double>::value)
      data = {
          2.81628362,    -0.01775168275, 0.1840875608,    -0.06937192284, 0.2972079123,   -0.518014912,  0.222584989,
          0.6303155211,  -0.01775168275, 3.177122572,     0.7447037649,   0.1186410285,   -0.1685090599, -0.5863662783,
          -0.317914909,  -0.03609288544, 0.1840875608,    0.7447037649,   3.218855407,    0.1396371692,  -0.7107370541,
          -0.8615134233, 0.2398770742,   0.3233909845,    -0.06937192284, 0.1186410285,   0.1396371692,  2.930167779,
          -0.6871870313, -0.3673295739,  0.1229376989,    -0.3746158801,  -5.749109707,   -10.58726466,  -12.71894083,
          -10.86344686,  3.557737246,    7.2022276,       1.739331109,    -2.555554659,   -8.796558043,  0.3637535164,
          -13.70880564,  -15.58312967,   6.313333104,     14.98263762,    -1.764196583,   -7.651551478,  7.696425026,
          -5.460943109,  3.809347506,    6.310960067,     -3.061742306,   -6.285584869,   2.377243215,   4.908643128,
          -6.951595685,  -11.03818442,   -12.23920128,    -8.113342449,   2.756520499,    8.189056561,   1.701740681,
          -3.844500193,  2.81628362,     -0.01775168275,  0.1840875608,   -0.06937192284, 0.2972079123,  -0.518014912,
          0.222584989,   0.6303155211,   -0.01775168275,  3.177122572,    0.7447037649,   0.1186410285,  -0.1685090599,
          -0.5863662783, -0.317914909,   -0.03609288544,  0.1840875608,   0.7447037649,   3.218855407,   0.1396371692,
          -0.7107370541, -0.8615134233,  0.2398770742,    0.3233909845,   -0.06937192284, 0.1186410285,  0.1396371692,
          2.930167779,   -0.6871870313,  -0.3673295739,   0.1229376989,   -0.3746158801,  -0.3228930173, 1.998020536,
          -0.3083396032, 1.295016536,    -0.01908192398,  -0.4548772415,  -0.01128301401, -0.845260842,  -1.310094894,
          -0.926372115,  0.6224620856,   -0.005416125989, 0.2878646673,   0.3025350725,   -0.2815865688, 0.1839608816,
          1.111915377,   -0.1469707379,  -0.3571450662,   1.688192462,    0.7720891953,   0.126319864,   -0.1756465577,
          -0.267816385,  -2.150403951,   1.820412041,     3.943488427,    0.1215705779,   -1.278372043,  -0.4183742809,
          0.3387590886,  -0.2382786398,
      };
    else if constexpr (std::is_same<OneBodyDensityMatrices::Real, float>::value)
      data = {
          2.816282749,   -0.01775467396, 0.1840839088,    -0.06936730444, 0.2972071171,   -0.5180174112, 0.2225853503,
          0.6303144693,  -0.01775258593, 3.177124739,     0.7447066307,   0.1186397672,   -0.168511793,  -0.5863617659,
          -0.3179149628, -0.03609220684, 0.1840886027,    0.7447082996,   3.218853951,    0.1396351308,  -0.7107358575,
          -0.8615118265, 0.2398764491,   0.3233931363,    -0.06936863065, 0.1186418384,   0.1396371722,  2.930168629,
          -0.6871862411, -0.3673290908,  0.1229372546,    -0.3746160269,  -5.749118805,   -10.58727074,  -12.71894264,
          -10.86344337,  3.557730913,    7.202232838,     1.7393291,      -2.555566549,   -8.796613693,  0.3637294769,
          -13.70880127,  -15.58312702,   6.313296318,     14.98267078,    -1.76416862,    -7.651569843,  7.696452141,
          -5.460932732,  3.809342146,    6.31096077,      -3.061723709,   -6.285599709,   2.377227545,   4.908651829,
          -6.951607704,  -11.0381918,    -12.23921108,    -8.113337517,   2.756513834,    8.189065933,   1.701741219,
          -3.844512224,  2.816282988,    -0.01775406115,  0.1840856224,   -0.06936839223, 0.2972068191,  -0.5180151463,
          0.2225854397,  0.6303145885,   -0.01775429584,  3.177122116,    0.7447044849,   0.1186416447,  -0.168511495,
          -0.5863642097, -0.3179140091,  -0.03609220684,  0.1840862632,   0.7447050214,   3.21885252,    0.1396374106,
          -0.7107351422, -0.8615128398,  0.2398774922,    0.3233922422,   -0.06936822832, 0.1186413094,  0.1396369785,
          2.930169582,   -0.6871861815,  -0.3673312664,   0.1229370609,   -0.3746150136,  -0.3228892386, 1.99802053,
          -0.308336854,  1.295014977,    -0.01908075809,  -0.4548793137,  -0.01128436625, -0.8452599645, -1.310091734,
          -0.9263712168, 0.6224637628,   -0.005418270826, 0.2878654003,   0.3025329113,   -0.2815873027, 0.1839606613,
          1.11191988,    -0.1469697803,  -0.3571412563,   1.688190818,    0.7720907331,   0.126316905,   -0.1756481826,
          -0.2678145766, -2.150410652,   1.820416451,     3.943485975,    0.1215699911,   -1.278371572,  -0.4183692932,
          0.3387611508,  -0.238276422,
      };
  }
  return data;
}

} // namespace testing
} // namespace qmcplusplus
