//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2024 QMCPACK developers.
//
// File developed by: Peter W. Doak, doakpw@ornl.gov, Oak Ridge National Lab
//
// File created by: Peter W. Doak, doakpw@ornl.gov, Oak Ridge National Lab
//////////////////////////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_TEST_ENERGY_DENSITY_ESTIMATOR_INTEGRATION
#define QMCPLUSPLUS_TEST_ENERGY_DENSITY_ESTIMATOR_INTEGRATION

#include "Listener.hpp"

namespace qmcplusplus
{
namespace testing
{

double cannedSum()
{
  CrowdEnergyValues<double> kinetic_values = {
      {{"Kinetic"},
       {
           {
               0.4529014498,
               0.3156825975,
               1.456884142,
               0.2485722715,
               1.061016946,
               1.175176701,
               0.3810680072,
               0.5948852591,
           },
           {
               1.668881229,
               0.6048614649,
               4.647461,
               1.358203882,
               -0.1954819386,
               0.5372290975,
               0.4891795562,
               1.169825556,
           },
           {
               1.231372469,
               0.1803754377,
               -1.205883308,
               1.410959272,
               0.6216167405,
               -0.4923786961,
               0.240801627,
               1.489066498,
           },
           {
               1.372523683,
               0.7555570136,
               3.226946012,
               -0.7894047702,
               1.609594113,
               0.2679526229,
               0.4790043943,
               1.245813374,
           },
       }},
  };
  CrowdEnergyValues<double> local_potential_values = {{{"ElecIon"},
                                                       {
                                                           {
                                                               0.1902921988,
                                                               0.305111589,
                                                               -0.4795567012,
                                                               0.4963486847,
                                                               -0.3115299195,
                                                               -0.3619893748,
                                                               0.09889833001,
                                                               0.2794993745,
                                                           },
                                                           {
                                                               -0.5591061161,
                                                               -0.1018774551,
                                                               -1.758163531,
                                                               -0.6303051354,
                                                               0.4878929497,
                                                               0.3358202321,
                                                               0.3278906128,
                                                               -0.3592569675,
                                                           },
                                                           {
                                                               0.2308181408,
                                                               0.3380463882,
                                                               0.4502142537,
                                                               0.09129681658,
                                                               0.5120332566,
                                                               0.4968745527,
                                                               0.5207083896,
                                                               -0.3453223962,
                                                           },
                                                           {
                                                               0.1190383565,
                                                               -0.1566598476,
                                                               -1.061808708,
                                                               0.4489104633,
                                                               -0.4774404553,
                                                               0.3925164497,
                                                               0.1283951707,
                                                               -0.436834872,
                                                           },
                                                       }},
                                                      {{"ElecElec"},
                                                       {{
                                                            -0.6210508854,
                                                            -0.5492735065,
                                                            -0.5151950441,
                                                            -0.5785561754,
                                                            -0.5071021975,
                                                            -0.5744799702,
                                                            -0.5424577691,
                                                            -0.7204598054,
                                                        },
                                                        {
                                                            -0.3894105334,
                                                            -0.5444925249,
                                                            -0.4896973204,
                                                            -0.5645833732,
                                                            -0.4432189127,
                                                            -0.5226105398,
                                                            -0.4715625335,
                                                            -0.528633033,
                                                        },
                                                        {
                                                            -0.04130562213,
                                                            -0.04840737416,
                                                            -0.5634037657,
                                                            -0.6271532836,
                                                            -0.1013881328,
                                                            -0.2272211231,
                                                            -0.2683512191,
                                                            -0.2033476916,
                                                        },
                                                        {
                                                            -0.3962585141,
                                                            -0.5585863635,
                                                            -0.5626605561,
                                                            -0.3180610119,
                                                            -0.3441004857,
                                                            -0.2138720089,
                                                            -0.3116368159,
                                                            -0.5650381055,
                                                        }}}};
  CrowdEnergyValues<double> local_ion_pot_values   = {
      {{"ElecIon"},
         {
           {
               0.04332090412,
               0.1737532775,
           },
           {
               -0.9332901929,
               -1.323815217,
           },
           {
               1.695662968,
               0.5990064344,
           },
           {
               0.1634203786,
               -1.207303821,
           },
       }},
  };


  auto sumOver = [](auto& crowd_energy) -> double {
    double sum{0};
    for (auto iter = crowd_energy.begin(); iter != crowd_energy.end(); ++iter)
    {
      auto& v_walkers = iter->second;
      for (auto& v_particles : v_walkers)
      {
        for (const double& t : v_particles)
        {
          sum += t;
        }
      }
    }
    return sum;
  };
  auto sum_pot = sumOver(local_potential_values);
  auto kinetic = sumOver(kinetic_values);
  auto ion_pot = sumOver(local_ion_pot_values);

  return sum_pot + kinetic + ion_pot;
}
} // namespace testing
} // namespace qmcplusplus

#endif
