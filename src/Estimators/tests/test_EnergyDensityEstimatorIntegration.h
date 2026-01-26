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
  CrowdEnergyValues<double> local_potential_values = {
      {{"ElecIon"},
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
       {
           {
               -0.6207395076,
               -0.5429517933,
               -0.5065277177,
               -0.5783610014,
               -0.4984350535,
               -0.5916988187,
               -0.5487789099,
               -0.7210825667,
           },
           {
               -0.3816392959,
               -0.5430679423,
               -0.482442029,
               -0.5642089765,
               -0.4469150835,
               -0.5218930872,
               -0.4652968702,
               -0.5215640811,
           },
           {
               0.0712246546,
               0.03917229912,
               -0.5630155803,
               -0.6444294957,
               -0.0487069078,
               -0.178371135,
               -0.2205232716,
               -0.1101885998,
           },
           {
               -0.3381089106,
               -0.5573446148,
               -0.5798054791,
               -0.2771065229,
               -0.2847240288,
               -0.1318201812,
               -0.276998751,
               -0.5649121642,
           },
       }},
  };
  CrowdEnergyValues<double> local_ion_pot_values = {
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
