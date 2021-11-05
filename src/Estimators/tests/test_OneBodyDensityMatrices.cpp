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
constexpr bool generate_test_data = true;

namespace testing
{
template<typename T>
class OneBodyDensityMatricesTests
{
public:
  using Evaluators  = OneBodyDensityMatricesInput::Evaluator;
  using Integrators = OneBodyDensityMatricesInput::Integrator;
  using Sampling    = OneBodyDensityMatrices::Sampling;
  using MCPWalker   = OneBodyDensityMatrices::MCPWalker;

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

  /** no change test for evaluateMatrix.
   *  I think this should be portable.
   */
  void testEvaluateMatrix(OneBodyDensityMatrices& obdm,
                          ParticleSet& pset,
                          TrialWaveFunction& trial_wavefunction,
                          MCPWalker& walker,
                          StdRandom<T>& rng)
  {
    obdm.evaluateMatrix(pset, trial_wavefunction, walker, rng);
    using Data = OneBodyDensityMatrices::Data::element_type;
    Data data;
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
    auto& returned_data = *(obdm.data_);
    for (size_t id = 0; id < data.size(); ++id)
      CHECK(returned_data[id] == Approx(data[id]));
  }

  void dumpData(OneBodyDensityMatrices& obdm)
  {
    std::cout << "Here is what is in your OneBodyDensityMatrices:\n" << NativePrint(*(obdm.data_)) << '\n';
  }
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

  {
    // Good constructor
    OneBodyDensityMatrices obDenMat(std::move(obdmi), lattice, species_set, wf_factory, pset_target);
    // Good copy constructor
    OneBodyDensityMatrices obDenMat2(obDenMat);
  }
  {
    species_set = testing::makeSpeciesSet(SpeciesCases::NO_MEMBERSIZE);
    CHECK_THROWS_AS(OneBodyDensityMatrices(std::move(obdmi), lattice, species_set, wf_factory, pset_target),
                    UniformCommunicateError);
  }

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

TEST_CASE("OneBodyDensityMatrices::clone()", "[estimators]")
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
  auto clone = original.clone();
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
  auto& species_set = pset_target.getSpeciesSet();
  OneBodyDensityMatrices(std::move(obdmi), pset_target.Lattice, species_set, wf_factory, pset_target);

  std::vector<MCPWalker> walkers;
  int nwalkers = 4;
  for (int iw = 0; iw < nwalkers; ++iw)
    walkers.emplace_back(8);

  std::vector<ParticleSet> psets;
  for (int iw = 0; iw < nwalkers; ++iw)
    psets.emplace_back(pset_target);

  auto& trial_wavefunction = *(wavefunction_pool.getPrimary());
  std::vector<UPtr<TrialWaveFunction>> twfcs(nwalkers);
  for (int iw = 0; iw < nwalkers; ++iw)
    twfcs[iw] = trial_wavefunction.makeClone(psets[iw]);

  // now the framework for testing accumulation is done

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
  if constexpr(generate_test_data)
  {
    std::cout << "Initialize pset_target.R with the following:\n{";
    for(auto r : pset_target.R)
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
  pset_target.R = ParticleSet::ParticlePos_t{{1.751870349, 4.381521229, 2.865202269}, {3.244515371, 4.382273176, 4.21105285},
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
  if constexpr(generate_test_data)
    obdmt.dumpData(obdm);
  outputManager.resume();
}

} // namespace qmcplusplus
