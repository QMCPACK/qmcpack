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

/** \file
 *  Tests for the OneBodyDensityMatrices.
 *  We can't reason about the state of the global Random in tests. A User can run only some tests,
 *  new tests will get added, other tests modified so global Random is called more times or fewer.
 *  Use of FakeRandom in unit tests in other tests of this executable can result in ambiguity
 *  about which global Random this test will have have access to as well.
 *  But several ParticleSet functions use the global Random so we have to avoid the normal
 *  sequence of particleset state transforms and set particle positions explicitly.
 */

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

/** set to true to generate new random R for particles.
 *  EXPECT all tests to fail in this case!
 */
constexpr bool generate_test_data = false;
// set to true to dump obdm for same R but perhaps different RNG effecting the integration.
constexpr bool dump_obdm = false;

std::vector<ParticleSet> generateRandomParticleSets(ParticleSet& pset_target,
                                                    ParticleSet& pset_source,
                                                    std::vector<ParticleSet::ParticlePos>& deterministic_rs,
                                                    int num_psets)
{
  int nwalkers = num_psets;
  std::vector<ParticleSet> psets(num_psets, pset_target);
  if constexpr (generate_test_data)
  {
    std::cout << "Initialize OneBodyDensityMatrices::accumulate psets with:\n{";
    std::vector<ParticleSet> psets;
    for (int iw = 0; iw < nwalkers; ++iw)
    {
      //psets.emplace_back(pset_target);
      psets.back().randomizeFromSource(pset_source);
      std::cout << "{";
      for (auto r : psets.back().R)
        std::cout << NativePrint(r) << ",";
      std::cout << "},\n";
    }
    std::cout << "}\n";
  }
  else
  {
    for (int iw = 0; iw < nwalkers; ++iw)
      psets.back().R = deterministic_rs[iw];
  }
  return psets;
}

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
#if defined(MIXED_PRECISION)
        CHECK(ref_data[id] == ComplexApprox(test_data[id]).epsilon(1e-4));
#else
        CHECK(ref_data[id] == ComplexApprox(test_data[id]));
#endif
    }
    else
    {
      for (size_t id = 0; id < size; ++id)
#if defined(MIXED_PRECISION)
        CHECK(ref_in[id] == Approx(test_in[id]).epsilon(1e-4));
#else
        CHECK(ref_in[id] == Approx(test_in[id]));
#endif
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
    if constexpr (generate_test_data)
      FAIL_CHECK("Test always fails when generating new test reference data.");
    else
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
    if constexpr (generate_test_data)
      FAIL_CHECK("Test always fails when generating new test reference data.");
    else
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

  auto particle_pool     = MinimalParticlePool::make_diamondC_1x1x1(comm);
  auto wavefunction_pool = MinimalWaveFunctionPool::make_diamondC_1x1x1(comm, particle_pool);
  auto& pset_target      = *(particle_pool.getParticleSet("e"));
  auto& wf_factory       = *(wavefunction_pool.getWaveFunctionFactory("wavefunction"));

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

  auto particle_pool     = MinimalParticlePool::make_diamondC_1x1x1(comm);
  auto wavefunction_pool = MinimalWaveFunctionPool::make_diamondC_1x1x1(comm, particle_pool);
  auto& pset_target      = *(particle_pool.getParticleSet("e"));
  auto& species_set      = pset_target.getSpeciesSet();
  auto& wf_factory       = *(wavefunction_pool.getWaveFunctionFactory("wavefunction"));

  auto samplingCaseRunner = [&pset_target, &species_set, &wf_factory](Inputs test_case) {
    Libxml2Document doc;

    bool okay = doc.parseFromString(valid_one_body_density_matrices_input_sections[test_case]);
    if (!okay)
      throw std::runtime_error("cannot parse OneBodyDensitMatricesInput section");
    xmlNodePtr node = doc.getRoot();
    OneBodyDensityMatricesInput obdmi(node);

    OneBodyDensityMatrices obDenMat(std::move(obdmi), pset_target.getLattice(), species_set, wf_factory, pset_target);

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

  auto particle_pool     = MinimalParticlePool::make_diamondC_1x1x1(comm);
  auto wavefunction_pool = MinimalWaveFunctionPool::make_diamondC_1x1x1(comm, particle_pool);
  auto& pset_target      = *(particle_pool.getParticleSet("e"));
  auto& species_set      = pset_target.getSpeciesSet();
  auto& wf_factory       = *(wavefunction_pool.getWaveFunctionFactory("wavefunction"));

  Libxml2Document doc;
  bool okay = doc.parseFromString(valid_one_body_density_matrices_input_sections[Inputs::valid_obdm_input]);
  if (!okay)
    throw std::runtime_error("cannot parse OneBodyDensitMatricesInput section");
  xmlNodePtr node = doc.getRoot();
  OneBodyDensityMatricesInput obdmi(node);

  OneBodyDensityMatrices original(std::move(obdmi), pset_target.getLattice(), species_set, wf_factory, pset_target);
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
  auto particle_pool     = MinimalParticlePool::make_diamondC_1x1x1(comm);
  auto wavefunction_pool = MinimalWaveFunctionPool::make_diamondC_1x1x1(comm, particle_pool);
  auto& wf_factory       = *(wavefunction_pool.getWaveFunctionFactory("wavefunction"));
  auto& pset_target      = *(particle_pool.getParticleSet("e"));
  auto& pset_source      = *(particle_pool.getParticleSet("ion"));
  auto& species_set      = pset_target.getSpeciesSet();
  OneBodyDensityMatrices obdm(std::move(obdmi), pset_target.getLattice(), species_set, wf_factory, pset_target);

  std::vector<MCPWalker> walkers;
  int nwalkers = 3;
  for (int iw = 0; iw < nwalkers; ++iw)
    walkers.emplace_back(8);

  std::vector<ParticleSet::ParticlePos> deterministic_rs = {{
                                                                {-0.6759092808, 0.835668385, 1.985307097},
                                                                {0.09710352868, -0.76751858, -1.89306891},
                                                                {-0.5605484247, -0.9578875303, 1.476860642},
                                                                {2.585144997, 1.862680197, 3.282609463},
                                                                {-0.1961335093, 1.111888766, -0.578481257},
                                                                {1.794641614, 1.6000278, -0.9474347234},
                                                                {2.157717228, 0.9254754186, 2.263158321},
                                                                {1.883366346, 2.136350632, 3.188981533},
                                                            },
                                                            {
                                                                {-0.2079261839, -0.2796236873, 0.5512072444},
                                                                {-0.2823159397, 0.7537326217, 0.01526880637},
                                                                {3.533515453, 2.433290243, 0.9281452894},
                                                                {2.051767349, 2.312927485, 0.7089259624},
                                                                {-1.043096781, 0.8190526962, -0.1958218962},
                                                                {0.9210210443, 0.7726522088, 0.3962054551},
                                                                {2.043324947, 0.3482068777, 3.39059639},
                                                                {0.9103830457, 2.167978764, 2.341906071},
                                                            },
                                                            {
                                                                {-0.466550231, 0.09173964709, -0.3779250085},
                                                                {-0.4211375415, -2.017466068, -1.691870451},
                                                                {2.090800285, 1.88529861, 2.152359247},
                                                                {2.973145723, 1.718174577, 3.822324753},
                                                                {-0.8552014828, -0.3484517336, -0.2870049179},
                                                                {0.2349359095, -0.5025780797, 0.2305756211},
                                                                {-0.03547382355, 2.279159069, 3.057915211},
                                                                {2.535993099, 1.637133598, 3.689830303},
                                                            }};
  std::vector<ParticleSet> psets = generateRandomParticleSets(pset_target, pset_source, deterministic_rs, nwalkers);

  auto& trial_wavefunction = *(wavefunction_pool.getPrimary());
  std::vector<UPtr<TrialWaveFunction>> twfcs(nwalkers);
  for (int iw = 0; iw < nwalkers; ++iw)
    twfcs[iw] = trial_wavefunction.makeClone(psets[iw]);

  StdRandom<double> rng;
  rng.init(101);

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

  if constexpr (dump_obdm)
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

    auto particle_pool     = MinimalParticlePool::make_diamondC_1x1x1(comm);
    auto wavefunction_pool = MinimalWaveFunctionPool::make_diamondC_1x1x1(comm, particle_pool);
    auto& wf_factory       = *(wavefunction_pool.getWaveFunctionFactory("wavefunction"));
    auto& pset_target      = *(particle_pool.getParticleSet("e"));
    auto& species_set      = pset_target.getSpeciesSet();
    OneBodyDensityMatrices obdm(std::move(obdmi), pset_target.getLattice(), species_set, wf_factory, pset_target);
    auto& trial_wavefunction = *(wavefunction_pool.getPrimary());

    // Because we can't control or consistent know the global random state we must initialize particle positions to known values.
    pset_target.R = ParticleSet::ParticlePos{
        {4.120557308, 2.547962427, 2.11555481},   {2.545657158, 2.021627665, 3.17555666},
        {1.251996636, 1.867651463, 0.7268046737}, {4.749059677, 5.845647812, 3.871560574},
        {5.18129015, 4.168475151, 2.748870373},   {6.24560833, 4.087143421, 4.187825203},
        {3.173382998, 3.651777267, 2.970916748},  {1.576967478, 2.874752045, 3.687536716},
    };

    StdRandom<double> rng;
    rng.init(101);
    MCPWalker walker;
    // Now we have to bring the pset, trial_wavefunction and walker to valid state.
    //pset.loadWalker(walker, false);
    pset_target.update(true);
    pset_target.donePbyP();
    trial_wavefunction.evaluateLog(pset_target);
    pset_target.saveWalker(walker);
    OneBodyDensityMatricesTests<double> obdmt;
    obdmt.testEvaluateMatrix(obdm, pset_target, trial_wavefunction, walker, rng);
    if constexpr (dump_obdm)
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
      data = {
          0.8559739634,     -2.775557562e-16, -0.0009312358165, -0.0002137464399, -0.00188481782,   -0.00010560363,
          -0.1690957318,    4.68545619e-11,   0.2135706414,     4.229001871e-11,  -0.0003234768609, -0.0008471008737,
          0.0004354250735,  -0.000629313951,  0.2632686421,     -4.718447855e-16, -0.0009312358165, 0.0002137464399,
          0.6904253737,     -5.551115123e-17, -2.29895663e-05,  -6.759071163e-05, 8.571154293e-05,  -2.029129991e-05,
          -0.0002755957962, 6.15929991e-05,   0.01997752724,    -0.03959823579,   -0.03441868427,   -0.03703684606,
          -0.000228978625,  5.163016816e-05,  -0.00188481782,   0.00010560363,    -2.29895663e-05,  6.759071163e-05,
          0.6910137193,     1.110223025e-15,  0.0001789239062,  -9.720928471e-06, -0.0005431524679, 3.125021981e-05,
          -0.0370646689,    -0.03445096274,   0.03964615719,    -0.02000771091,   -0.0004552815242, 2.596515681e-05,
          -0.1690957318,    -4.685549171e-11, 8.571154293e-05,  2.029129991e-05,  0.0001789239062,  9.720928472e-06,
          0.6510614284,     4.996003611e-16,  -0.3264829535,    -1.120145643e-11, -1.628936668e-05, -5.104526467e-05,
          2.191572411e-05,  -3.7934521e-05,   -0.1113259098,    -1.13588236e-11,  0.9949076208,     5.254666641e-08,
          0.09716631627,    -0.008698974489,  -0.07670746649,   0.01101924465,    -1.054191923,     -1.009885029e-07,
          0.6429960606,     5.478465426e-08,  0.006654585037,   -0.003422468497,  -0.008957381204,  -0.002542568933,
          0.3883611988,     2.492687778e-08,  -0.2567403319,    -2.274337275,     -0.070057294,     -2.130173974,
          -0.185057423,     -2.774684629,     0.494079399,      4.39071663,       -0.2679771972,    -2.380684367,
          -0.2537442403,    0.1010428549,     -0.2026493778,    -0.04440749192,   -0.1216357873,    -1.078602297,
          0.345585874,      -1.689635575,     0.7414178343,     -1.48274947,      0.742418168,      -1.971367257,
          -0.6650567573,    3.261922039,      0.3607112874,     -1.768642232,     -0.202649463,     -0.01653546255,
          -0.1315184346,    -0.1010429096,    0.1637281089,     -0.8013080665,    -1.915594805,     -3.652633751e-14,
          -2.106812279,     -0.3174159872,    -2.798953942,     -0.2389235612,    4.203062081,      7.690383574e-08,
          -2.237403378,     -1.095805402e-08, 0.05937990357,    0.2648259866,     -0.0799285254,    0.196742672,
          -0.9570615415,    -1.670885652e-14, 0.8559739634,     5.551115123e-17,  -0.0009312358165, -0.0002137464399,
          -0.00188481782,   -0.00010560363,   -0.1690957318,    4.685551946e-11,  0.2135706414,     4.229011585e-11,
          -0.0003234768609, -0.0008471008737, 0.0004354250735,  -0.000629313951,  0.2632686421,     -3.330669074e-16,
          -0.0009312358165, 0.0002137464399,  0.6904253737,     -5.551115123e-17, -2.29895663e-05,  -6.759071163e-05,
          8.571154293e-05,  -2.029129991e-05, -0.0002755957962, 6.15929991e-05,   0.01997752724,    -0.03959823579,
          -0.03441868427,   -0.03703684606,   -0.000228978625,  5.163016816e-05,  -0.00188481782,   0.00010560363,
          -2.29895663e-05,  6.759071163e-05,  0.6910137193,     5.551115123e-17,  0.0001789239062,  -9.720928471e-06,
          -0.0005431524679, 3.125021981e-05,  -0.0370646689,    -0.03445096274,   0.03964615719,    -0.02000771091,
          -0.0004552815242, 2.596515681e-05,  -0.1690957318,    -4.68551864e-11,  8.571154293e-05,  2.029129991e-05,
          0.0001789239062,  9.720928471e-06,  0.6510614284,     1.665334537e-16,  -0.3264829535,    -1.12007903e-11,
          -1.628936668e-05, -5.104526467e-05, 2.191572411e-05,  -3.7934521e-05,   -0.1113259098,    -1.13588583e-11,
          -0.1398697498,    1.676108652e-08,  0.009839305325,   -0.03110873861,   -0.2743146761,    0.001115813703,
          0.8078232309,     7.208311281e-09,  -0.3939381356,    7.358152143e-10,  0.01322489101,    0.01202117753,
          -0.01780140861,   0.008930678923,   -0.1179882345,    4.657475086e-09,  0.1004197767,     -0.3872746911,
          -0.4631185606,    -0.1267234616,    0.251940383,      0.3432366261,     0.1357909209,     0.1882467853,
          -0.04657813737,   -0.1481318675,    -0.01754424248,   -0.008096072777,  0.04046503369,    0.04321394303,
          0.01594660416,    -0.1297942252,    -0.1351702458,    -0.2877115343,    0.5308175294,     -0.1497099916,
          -0.302802717,     0.3542104112,     -0.1827815947,    0.1398510191,     0.06269657586,    -0.1100491288,
          0.04046503056,    -0.0383308965,    -0.04195025864,   0.008096072041,   -0.02146496291,   -0.09642585851,
          1.583358633,      5.551115123e-15,  0.5400010285,     0.09966977815,    0.8788816574,     0.06123897411,
          -2.364549544,     -2.291840495e-08, 1.339153745,      3.484556643e-09,  -0.02316768222,   -0.0762631529,
          0.03118495745,    -0.05665686084,   0.6842069691,     1.720845688e-15,
      };
    }
    else if constexpr (std::is_floating_point<OneBodyDensityMatrices::Value>::value)
    {
      data = {
          0.4279869817,     -0.0005724911282, -0.0009952107251, -0.08454786586,   0.1067853207,     -0.0005852888673,
          -9.694443874e-05, 0.131634321,      -0.0005724911282, 0.3452101128,     0.07732509099,    5.299452569e-05,
          -0.0001685944122, -0.01345751284,   -0.03081850011,   -0.0001403043966, -0.0009952107251, 0.07732509099,
          0.3455042856,     9.430886761e-05,  -0.0002871986719, -0.03308358791,   0.006219573346,   -0.0002406233405,
          -0.08454786586,   5.299452569e-05,  9.430886761e-05,  0.3255307142,     -0.1632414767,    -3.366691551e-05,
          -8.009937089e-06, -0.0556629549,    0.4974537841,     0.04423363063,    -0.03284416792,   -0.5270959104,
          0.3214980027,     0.001616059574,   -0.005749976824,  0.1941805869,     1.008798472,      0.9822905645,
          1.251738209,      -1.948318603,     1.056353583,      -0.1335433475,    -0.04654405447,   0.4784832548,
          1.017610724,      1.176381561,      1.414874433,      -1.963489414,     1.064676762,      -0.1520816495,
          -0.05908794817,   0.4825180877,     -0.9577974025,    -1.212114133,     -1.518938752,     2.101531079,
          -1.118701694,     0.1621029451,     0.05840707331,    -0.4785307708,    0.4279869817,     -0.0005724911282,
          -0.0009952107251, -0.08454786586,   0.1067853207,     -0.0005852888673, -9.694443874e-05, 0.131634321,
          -0.0005724911282, 0.3452101128,     0.07732509099,    5.299452569e-05,  -0.0001685944122, -0.01345751284,
          -0.03081850011,   -0.0001403043966, -0.0009952107251, 0.07732509099,    0.3455042856,     9.430886761e-05,
          -0.0002871986719, -0.03308358791,   0.006219573346,   -0.0002406233405, -0.08454786586,   5.299452569e-05,
          9.430886761e-05,  0.3255307142,     -0.1632414767,    -3.366691551e-05, -8.00993709e-06,  -0.0556629549,
          -0.06993488329,   -0.01063469842,   -0.1365994325,    0.4039116172,     -0.1969690692,    0.01262303497,
          -0.004435365783,  -0.05899411958,   0.2438472339,     -0.1415980475,    -0.09314308211,   -0.0262279326,
          0.05077686603,    0.01074587779,    0.0101183853,     0.07287041467,    0.07627064426,    0.3044594466,
          -0.2645758779,    -0.1613163063,    0.08637285105,    0.01857458999,    -0.04049312511,   0.0374804478,
          0.7916793164,     0.3198354033,     0.4700603157,     -1.182274783,     0.6695768741,     -0.04971541756,
          -0.01273595169,   0.3421034845,
      };
    }
    break;
  }
  case OBDMI::Integrator::UNIFORM: {
    if constexpr (IsComplex_t<OneBodyDensityMatrices::Value>::value)
    {
      data = {
          0.8261941614,   0,
          -0.04207477301, -0.01042641052,
          -0.09193944978, -0.004771495112,
          -0.1577639198,  2.251655112e-09,
          0.2529218804,   -3.528230208e-10,
          0.01395852521,  -0.01381912806,
          -0.01878889069, -0.0102664046,
          0.2669242402,   -1.665334537e-16,
          -0.04207477301, 0.01042641052,
          0.7237235295,   -2.220446049e-16,
          -0.1362791782,  0.02883929086,
          -0.0345711543,  -0.01578687034,
          0.07089423603,  0.01585967312,
          0.06197413691,  -0.0673375704,
          -0.07846982731, -0.0428311471,
          0.03401975498,  0.01636340115,
          -0.09193944978, 0.004771495112,
          -0.1362791782,  -0.02883929086,
          0.4726910113,   3.330669074e-16,
          0.1392076512,   0.003920553607,
          -0.1398496516,  -0.008039785222,
          0.02312496202,  0.01413794966,
          -0.04626371454, 0.02462349008,
          -0.1442914533,  -0.003858018529,
          -0.1577639198,  -2.251654918e-09,
          -0.0345711543,  0.01578687034,
          0.1392076512,   -0.003920553607,
          0.6406077933,   -2.220446049e-16,
          -0.3821572706,  -6.950154341e-11,
          -0.03309714608, -0.006729454624,
          0.04455046481,  -0.004999419621,
          -0.06308963653, -2.521567177e-09,
          0.9566798997,   3.993940545e-08,
          0.1287007623,   -0.03913606336,
          -0.3450997639,  0.01459536453,
          -1.049841949,   -8.448363198e-08,
          0.7807572133,   5.406439579e-08,
          0.0636014061,   -0.0144126824,
          -0.08561076546, -0.01070736661,
          0.3453233488,   7.91976143e-09,
          -0.1652425535,  -1.673204981,
          -0.2864309551,  -1.827550661,
          0.1092013893,   -0.4672139335,
          0.4181616684,   3.856599084,
          -0.2678120582,  -2.450374785,
          -0.1661115247,  -0.5147466032,
          0.03950582141,  0.7348504126,
          -0.01780951734, -0.3040424203,
          0.2224251057,   -1.243046364,
          0.4872309759,   -1.334778155,
          0.2833545206,   -0.2677130338,
          -0.5628674436,  2.865118995,
          0.3604889296,   -1.820416191,
          0.03950573518,  -0.4243869589,
          -0.1899388751,  0.5147465492,
          0.02397254984,  -0.2258771884,
          -1.317725416,   -2.642330799e-14,
          -1.83826268,    -0.06023003292,
          -0.5311043355,  -0.2084686282,
          3.662896282,    2.914568187e-08,
          -2.270465021,   2.617466888e-08,
          -0.5283222051,  0.08433855146,
          0.7111490664,   0.06265614476,
          -0.1773974531,  -2.775557562e-15,
          0.8261941614,   -5.551115123e-17,
          -0.04207477301, -0.01042641052,
          -0.09193944978, -0.004771495112,
          -0.1577639198,  2.251654183e-09,
          0.2529218804,   -3.528226322e-10,
          0.01395852521,  -0.01381912806,
          -0.01878889069, -0.0102664046,
          0.2669242402,   -1.942890293e-16,
          -0.04207477301, 0.01042641052,
          0.7237235295,   -2.220446049e-16,
          -0.1362791782,  0.02883929086,
          -0.0345711543,  -0.01578687034,
          0.07089423603,  0.01585967312,
          0.06197413691,  -0.0673375704,
          -0.07846982731, -0.0428311471,
          0.03401975498,  0.01636340115,
          -0.09193944978, 0.004771495112,
          -0.1362791782,  -0.02883929086,
          0.4726910113,   5.551115123e-17,
          0.1392076512,   0.003920553607,
          -0.1398496516,  -0.008039785222,
          0.02312496202,  0.01413794966,
          -0.04626371454, 0.02462349008,
          -0.1442914533,  -0.003858018529,
          -0.1577639198,  -2.251654516e-09,
          -0.0345711543,  0.01578687034,
          0.1392076512,   -0.003920553607,
          0.6406077933,   -1.665334537e-16,
          -0.3821572706,  -6.950134912e-11,
          -0.03309714608, -0.006729454624,
          0.04455046481,  -0.004999419621,
          -0.06308963653, -2.521567073e-09,
          -0.09246990704, 1.714692167e-08,
          0.01718168811,  -0.002306490387,
          -0.0203382379,  0.001948481938,
          0.7391743085,   3.440059682e-09,
          -0.4038538237,  3.37227038e-09,
          -0.05196998174, -0.01900034084,
          0.06995427856,  -0.01411562006,
          0.001586062406, 3.942820508e-09,
          0.09257244551,  -0.4150339965,
          -0.543569044,   -0.1978266927,
          0.2876563907,   0.3016008129,
          0.2041868448,   0.2689153029,
          -0.1456849905,  -0.2702324047,
          -0.06303854185, 0.03455118382,
          0.02035291547,  0.04212529665,
          -0.04182725494, -0.2125070964,
          -0.1246073282,  -0.308334295,
          0.6467124191,   -0.2090012216,
          -0.3333005141,  0.3431996583,
          -0.2748460697,  0.1997807403,
          0.1960995565,   -0.2007592471,
          0.02035290751,  -0.06296446567,
          -0.07531415043, -0.03455118719,
          0.05630165983,  -0.15787436,
          1.360996296,    4.440892099e-15,
          0.460662518,    -0.009109012253,
          -0.08032289324, 0.05224157011,
          -2.1641054,     -3.065196319e-09,
          1.447670416,    -1.103484504e-08,
          0.2132345387,   -0.01793768528,
          -0.2870247183,  -0.01332608765,
          0.369682822,    9.436895709e-16,
      };
    }
    else if constexpr (std::is_floating_point<OneBodyDensityMatrices::Value>::value)
    {
      data = {
          0.4130970807,   -0.02625059177,  -0.04835547244,  -0.07888195875,  0.12646094,      6.96985738e-05,
          -0.01452764765, 0.1334621201,    -0.02625059177,  0.3466032331,    -0.001161163685, -0.009392146568,
          0.02751728009,  0.0007624905835, -0.06528656996,  0.008828176914,  -0.04835547244,  -0.001161163685,
          0.2210869708,   0.06764354093,   -0.06590493546,  0.02539089034,   -0.01991866217,  -0.07021671737,
          -0.07888195875, -0.009392146568, 0.06764354093,   0.3203038942,    -0.1910786333,   -0.01991330156,
          0.01977552422,  -0.031544817,    0.4783399299,    0.04478231453,   -0.165252207,    -0.5249209266,
          0.3903785766,   0.02459435045,   -0.04815905069,  0.1726616705,    0.7539812136,    0.7595813826,
          0.2502049139,   -1.719218703,    1.091281367,     0.1944111372,    -0.3327445041,   0.1431164515,
          0.7327357349,   0.9257821561,    0.3266874516,    -1.713993226,    1.090452555,     0.2048993189,
          -0.372436296,   0.1249248691,    -0.6588627081,   -0.9492463565,   -0.3697864819,   1.831448156,
          -1.135232497,   -0.2219918268,   0.3869026056,    -0.08869872655,  0.4130970807,    -0.02625059177,
          -0.04835547244, -0.07888195875,  0.12646094,      6.96985738e-05,  -0.01452764765,  0.1334621201,
          -0.02625059177, 0.3466032331,    -0.001161163685, -0.009392146568, 0.02751728009,   0.0007624905835,
          -0.06528656996, 0.008828176914,  -0.04835547244,  -0.001161163685, 0.2210869708,    0.06764354093,
          -0.06590493546, 0.02539089034,   -0.01991866217,  -0.07021671737,  -0.07888195875,  -0.009392146568,
          0.06764354093,  0.3203038942,    -0.1910786333,   -0.01991330156,  0.01977552422,   -0.031544817,
          -0.0462349621,  0.007437619475,  -0.009194881695, 0.3695871528,    -0.2019269101,   -0.03548515978,
          0.02791932722,  0.0007930292316, 0.253803221,     -0.1431755519,   -0.0640034998,   -0.03236422918,
          0.06227370761,  -0.006365722925, 0.02063497821,   0.08533992072,   0.0918634834,    0.3878849668,
          -0.2614830065,  -0.2373134048,   0.1984294012,    -0.01545314238,  -0.06281061419,  0.1070880099,
          0.6804981478,   0.2257767529,    -0.01404066156,  -1.082052701,    0.7238352026,    0.09764842672,
          -0.150175403,   0.184841411,
      };
    }
    break;
  }
  case OBDMI::Integrator::DENSITY: {
    if constexpr (IsComplex_t<OneBodyDensityMatrices::Value>::value)
    {
      data = {
          0.9924480859,   -1.110223025e-16,
          -0.1885291696,  0.0319276602,
          0.2815361269,   -0.02138019278,
          0.1872112957,   -3.684201194e-09,
          -0.5244043425,  1.095947622e-08,
          -0.03918278828, -0.3023406232,
          0.05274204952,  -0.2246129078,
          0.01383924214,  9.714451465e-17,
          -0.1885291696,  -0.0319276602,
          0.9180621234,   2.775557562e-16,
          0.2402287308,   -0.01308718541,
          0.2057882223,   0.004465819328,
          0.1024448325,   0.02374868021,
          -0.05743147279, -0.07637345821,
          0.03130436103,  -0.08626204101,
          0.1630605924,   -0.02521124324,
          0.2815361269,   0.02138019278,
          0.2402287308,   0.01308718541,
          1.031979911,    3.330669074e-16,
          -0.03937932123, -0.02333745214,
          -0.2094143097,  -0.01161775883,
          -0.1348871245,  -0.190165311,
          0.1601241654,   -0.149665044,
          0.2223111845,   -0.01849191565,
          0.1872112957,   3.684200389e-09,
          0.2057882223,   -0.004465819328,
          -0.03937932123, 0.02333745214,
          1.323421938,    4.440892099e-16,
          0.2905844345,   1.315278547e-08,
          -0.1594450613,  -0.2611690287,
          0.2146213636,   -0.1940259416,
          0.2116494581,   6.064937663e-09,
          0.5632626352,   1.703273084e-08,
          -0.3488988503,  0.02573820756,
          0.2269574287,   -0.03956679718,
          -1.637494571,   -1.812725284e-07,
          -0.8314244571,  -8.494759768e-08,
          0.1951723949,   0.1024546704,
          -0.2627122381,  0.0761148891,
          -0.282499705,   -2.532976356e-09,
          -0.09350615209, -0.756293056,
          0.1463790836,   -2.234477269,
          -0.5935525371,  -5.517868252,
          0.829208332,    7.707654298,
          0.3647990603,   3.119019826,
          0.1734726648,   -0.3225198405,
          0.1526169414,   0.4999332675,
          -0.07866503442, -0.05924377889,
          0.1258641643,   -0.5618602341,
          1.097431538,    -1.4619276,
          1.305834496,    -3.987176014,
          -1.116157862,   5.726119334,
          -0.4910385854,  2.317161517,
          0.1526170644,   -0.3054089574,
          0.08142330397,  0.3225199255,
          0.105887215,    -0.04401298315,
          -0.4339186561,  -1.629252289e-14,
          -2.305548843,   -0.6188863722,
          -5.45730064,    -0.2614613878,
          7.525382195,    1.317113143e-07,
          2.880837963,    2.432125967e-08,
          -0.325402542,   -0.228642532,
          0.4380088526,   -0.1698614745,
          -0.1214395281,  -3.462508058e-15,
          0.9924480859,   5.551115123e-17,
          -0.1885291696,  0.0319276602,
          0.2815361269,   -0.02138019278,
          0.1872112957,   -3.684200667e-09,
          -0.5244043425,  1.095947622e-08,
          -0.03918278828, -0.3023406232,
          0.05274204952,  -0.2246129078,
          0.01383924214,  -2.975050761e-16,
          -0.1885291696,  -0.0319276602,
          0.9180621234,   -1.665334537e-16,
          0.2402287308,   -0.01308718541,
          0.2057882223,   0.004465819328,
          0.1024448325,   0.02374868021,
          -0.05743147279, -0.07637345821,
          0.03130436103,  -0.08626204101,
          0.1630605924,   -0.02521124324,
          0.2815361269,   0.02138019278,
          0.2402287308,   0.01308718541,
          1.031979911,    0,
          -0.03937932123, -0.02333745214,
          -0.2094143097,  -0.01161775883,
          -0.1348871245,  -0.190165311,
          0.1601241654,   -0.149665044,
          0.2223111845,   -0.01849191565,
          0.1872112957,   3.68420057e-09,
          0.2057882223,   -0.004465819328,
          -0.03937932123, 0.02333745214,
          1.323421938,    0,
          0.2905844345,   1.315278556e-08,
          -0.1594450613,  -0.2611690287,
          0.2146213636,   -0.1940259416,
          0.2116494581,   6.064937469e-09,
          0.2051059875,   2.820113419e-08,
          0.1614723089,   -0.04906840403,
          -0.4326817327,  0.01831176767,
          1.706267288,    1.643257963e-08,
          0.4080117321,   4.721028735e-09,
          -0.1551120856,  -0.2788399487,
          0.2087889614,   -0.2071539023,
          0.1814356421,   -1.855115922e-09,
          0.4257831978,   -0.1752088462,
          -0.4998292759,  0.07462493825,
          0.2373005349,   0.3540303805,
          0.2248224859,   0.09288554326,
          -0.1494463458,  0.1278007236,
          -0.05568689306, -0.2034175338,
          0.1004071574,   0.007263799659,
          0.03200635492,  0.1228505216,
          -0.5731263129,  -0.130165,
          0.5767337619,   -0.002042114292,
          -0.330578729,   0.3705292719,
          -0.3026226697,  0.06900584231,
          0.2011625579,   0.09494483239,
          0.1004071511,   0.1154251684,
          -0.1162463317,  0.2034175291,
          -0.04308221061, 0.09126728585,
          0.7806523758,   3.275157923e-15,
          0.1160951465,   0.2239745488,
          1.974993197,    0.01316587446,
          -4.058777904,   -4.039943313e-08,
          -1.782306544,   -2.124291576e-08,
          0.2934590981,   0.1825596902,
          -0.3950112733,  0.1356259601,
          -0.2695462139,  -7.21644966e-16,
      };
    }
    else if constexpr (std::is_floating_point<OneBodyDensityMatrices::Value>::value)
    {
      data = {
          0.9408447786,   -0.3935847362, -0.2468586928,  0.1346424361,  -0.2790550847,  -0.2728613819,  0.118510686,
          0.373636656,    -0.3935847362, 1.465454239,    0.3906355741,  0.3826898337,   -0.03726984102, -0.07963536415,
          -0.4819140243,  -0.2130898108, -0.2468586928,  0.3906355741,  0.8424465027,   0.2589273731,   -0.03575380357,
          -0.253719428,   0.1136907332,  0.09985411935,  0.1346424361,  0.3826898337,   0.2589273731,   1.097297774,
          -0.2458047841,  -0.1441440472, -0.1503264601,  0.09909758853, 0.6175832301,   -0.6985723605,  -0.6243271971,
          -1.376534064,   0.09219592002, -0.02237084205, 0.2231280489,  0.1476634744,   -1.208204987,   1.712369561,
          1.715275159,    -4.409355333,  0.8508272302,   -0.494084464,  0.3242451383,   -0.2368202448,  -1.480471134,
          2.537211763,    2.19567287,    -4.176679345,   0.8267380683,  -0.6223330608,  0.1422109575,   -0.3020474178,
          1.853750075,    -2.613447154,  -2.388944783,   4.650543191,   -0.9996739752,  0.5654611738,   -0.2023390067,
          0.4167421543,   0.9408447786,  -0.3935847362,  -0.2468586928, 0.1346424361,   -0.2790550847,  -0.2728613819,
          0.118510686,    0.373636656,   -0.3935847362,  1.465454239,   0.3906355741,   0.3826898337,   -0.03726984102,
          -0.07963536415, -0.4819140243, -0.2130898108,  -0.2468586928, 0.3906355741,   0.8424465027,   0.2589273731,
          -0.03575380357, -0.253719428,  0.1136907332,   0.09985411935, 0.1346424361,   0.3826898337,   0.2589273731,
          1.097297774,    -0.2458047841, -0.1441440472,  -0.1503264601, 0.09909758853,  0.3273111356,   0.3790594202,
          -0.01488345199, 1.314728254,   -0.321881638,   -0.1063667686, -0.2555566944,  0.1032133755,   0.7484376476,
          -0.8135360105,  -0.4255269663, -0.03246449107, -0.1597379358, -0.09287811809, 0.2144758734,   0.2846013811,
          -0.1713163625,  1.015300497,   -0.563966848,   -0.3464659878, 0.087960779,    0.2142716444,   -0.5686552903,
          -0.3532239607,  0.1204139568,  -0.3049230367,  0.1268687562,  -2.916771087,   0.4142601541,   -0.209253697,
          0.4634440052,   0.1000471674,
      };
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
    data = {
        2.929123464,    -6.661338148e-16, -0.213298936,    0.07456746639,    0.6575312763,    -0.02418921452,
        0.2532299248,   -1.090116444e-08, -1.124057974,    2.768283235e-09,  -0.1397854178,   0.1645671318,
        0.1881583956,   0.1222590747,     -0.001959538554, 1.189152943e-15,  -0.213298936,    -0.07456746639,
        3.143157244,    -4.440892099e-16, 0.4556655323,    -0.003211036867,  0.0525287831,    0.001809268674,
        -0.1435885589,  0.03450143497,    0.1530394508,    0.2388837829,     -0.1875535298,   0.1030355338,
        -0.1497527705,  -0.002470663282,  0.6575312763,    0.02418921452,    0.4556655323,    0.003211036867,
        3.171107807,    2.775557562e-16,  -0.0159540742,   -0.00595703733,   -0.3042315457,   0.01628372095,
        -0.2911526549,  0.06150813339,    0.4400568181,    0.07985836232,    0.0217861975,    0.01698275322,
        0.2532299248,   1.09011636e-08,   0.0525287831,    -0.001809268674,  -0.0159540742,   0.00595703733,
        3.407050358,    -1.33226763e-15,  0.8103626817,    5.595117272e-09,  -0.09528656784,  -0.2711773981,
        0.1282606987,   -0.201461268,     0.4847746576,    -1.082338323e-09, -10.96553943,    4.351111784e-07,
        9.556035377,    -2.515670736,     -22.18302448,    1.083704451,      -5.423715564,    6.291568599e-09,
        -2.752349188,   -1.423614395e-07, 0.5589025108,    -5.208477429,     -0.7523121017,   -3.869447511,
        3.14968853,     -4.90656993e-08,  5.589348626,     -4.521120313,     -4.019367357,    5.455671372,
        3.708424575,    0.886782507,      -2.655112814,    6.627356755,      -2.434203011,    4.52788362,
        0.01134390839,  1.290390669,      0.6530190148,    0.7867560782,     -0.8250466856,   0.1485604797,
        -7.523554275,   -3.358800218,     5.09625517,      3.041012493,      -6.169403958,    1.575377891,
        3.573919835,    4.923552154,      3.27656369,      3.36382536,       0.6530191051,    -1.565040693,
        -0.382517449,   -1.290390527,     1.110555857,     0.1103674375,     1.613012693,     -7.105427358e-15,
        1.024312958,    -0.8097052201,    -7.139930916,    0.1161626301,     -10.54695952,    1.3144171e-07,
        -7.349075649,   7.603942909e-08,  -0.2799244189,   -2.082751122,     0.3767925322,    -1.547303919,
        1.138309381,    6.80705492e-15,   2.929123464,     5.551115123e-17,  -0.213298936,    0.07456746639,
        0.6575312763,   -0.02418921452,   0.2532299248,    -1.090116469e-08, -1.124057974,    2.768283663e-09,
        -0.1397854178,  0.1645671318,     0.1881583956,    0.1222590747,     -0.001959538554, 2.445960101e-16,
        -0.213298936,   -0.07456746639,   3.143157244,     -3.330669074e-16, 0.4556655323,    -0.003211036867,
        0.0525287831,   0.001809268674,   -0.1435885589,   0.03450143497,    0.1530394508,    0.2388837829,
        -0.1875535298,  0.1030355338,     -0.1497527705,   -0.002470663282,  0.6575312763,    0.02418921452,
        0.4556655323,   0.003211036867,   3.171107807,     -6.661338148e-16, -0.0159540742,   -0.00595703733,
        -0.3042315457,  0.01628372095,    -0.2911526549,   0.06150813339,    0.4400568181,    0.07985836232,
        0.0217861975,   0.01698275322,    0.2532299248,    1.090116468e-08,  0.0525287831,    -0.001809268674,
        -0.0159540742,  0.00595703733,    3.407050358,     2.220446049e-16,  0.8103626817,    5.595117161e-09,
        -0.09528656784, -0.2711773981,    0.1282606987,    -0.201461268,     0.4847746576,    -1.082338538e-09,
        1.683965409,    2.80087866e-07,   2.745345297,     0.4187971974,     3.692926248,     0.3113366977,
        -0.7125045768,  8.85922713e-08,   -0.4120597303,   -1.382131761e-07, -0.0155196004,   1.081201079,
        0.02089027859,  0.8032387932,     -1.416803366,    6.419624812e-08,  2.991047084,     -4.036091748,
        0.3674297008,   -4.220747391,     3.079622169,     -3.287678963,     1.094933045,     -4.460680443,
        -1.268015018,   1.325003366,      -1.422085342,    0.05859756034,    0.01385636787,   -1.00610731,
        0.5495667106,   -1.17731452,      -4.026104815,    -2.998466529,     0.2841266859,    -3.743993937,
        -3.076449049,   -2.441195309,     -1.473837065,    -3.313898937,     1.706814235,     0.984362627,
        0.01385651404,  0.9707649019,     -1.430442654,    -0.05859742213,   -0.7397454204,   -0.8746426227,
        2.02228861,     -9.492406861e-15, -0.6157568594,   -0.1308165014,    -1.153531081,    -0.06983017621,
        7.227555974,    2.913955294e-08,  -0.09202041397,  1.006759607e-07,  -0.612848398,    -2.308230891,
        0.8249258024,   -1.714815559,     2.814827606,     -4.829470157e-15,
    };
  }
  else if constexpr (std::is_floating_point<OneBodyDensityMatrices::Value>::value)
  {
    data = {
        2.941887032,    -0.07303676966, 0.2825942788,   -0.04821280862, -0.4076353925,  0.05658482662, 1.181558243,
        0.1397691924,   -0.07303676966, 3.4212182,      0.89089857,     0.2575937162,   -0.1332837838, 0.5354984573,
        0.02230897646,  -0.4943417103,  0.2825942788,   0.89089857,     3.284611478,    0.2100203992,  -0.29693352,
        -0.06688847271, 0.7393865326,   0.07446761567,  -0.04821280862, 0.2575937162,   0.2100203992,  2.780464467,
        0.1479959718,   -0.2108628317,  -0.1058624903,  0.1094655214,   -15.19961226,   3.499431662,   -29.9410813,
        2.456596277,    3.770697701,    -6.175645508,   -8.264208029,   5.538268191,    10.94720948,   -14.05727688,
        -1.51312531,    -9.320123293,   -0.8158545739,  -1.677680626,   6.910273078,    3.656886378,   -4.390200646,
        -0.3436920564,  -8.023924667,   0.3691672347,   0.6188840349,   -1.257921128,   -2.506169099,  1.431044236,
        -0.9737456802,  -5.980733703,   -15.03553614,   -5.358452336,   1.176799479,    -4.97487775,   0.3262480281,
        5.975873531,    2.941887032,    -0.07303676966, 0.2825942788,   -0.04821280862, -0.4076353925, 0.05658482662,
        1.181558243,    0.1397691924,   -0.07303676966, 3.4212182,      0.89089857,     0.2575937162,  -0.1332837838,
        0.5354984573,   0.02230897646,  -0.4943417103,  0.2825942788,   0.89089857,     3.284611478,   0.2100203992,
        -0.29693352,    -0.06688847271, 0.7393865326,   0.07446761567,  -0.04821280862, 0.2575937162,  0.2100203992,
        2.780464467,    0.1479959718,   -0.2108628317,  -0.1058624903,  0.1094655214,   0.3066611943,  3.439932581,
        2.515144799,    -0.8458768467,  1.251686438,    0.2474206699,   0.617613976,    -1.038915172,  5.316217602,
        5.51728856,     5.207727306,    5.671785393,    -2.01767805,    -0.8133173328,  1.442175387,   1.246477853,
        -0.6941323987,  4.061136339,    -0.3928140137,  1.599557998,    -0.1536382161,  0.7947891936,  -1.011399761,
        -0.7830692539,  1.533636075,    -0.8928169035,  -0.8105425087,  6.912706899,    -1.612409448,  -1.464236775,
        -0.5807617769,  2.162165474,
    };
  }
  return data;
}

} // namespace testing
} // namespace qmcplusplus
