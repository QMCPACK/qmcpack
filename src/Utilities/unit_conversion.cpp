//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//////////////////////////////////////////////////////////////////////////////////////


#include "unit_conversion.h"
#include "Message/Communicate.h"


namespace qmcplusplus
{
namespace Units
{


units energy_unit(const std::string& su)
{
  if (su == "J")
    return J;
  else if (su == "eV")
    return eV;
  else if (su == "Ry")
    return Ry;
  else if (su == "Ha")
    return Ha;
  else if (su == "kJ/mol")
    return kJ_mol;
  else if (su == "K")
    return K;
  else if (su == "joule")
    return J;
  else if (su == "electron_volt")
    return eV;
  else if (su == "rydberg")
    return Ry;
  else if (su == "hartree")
    return Ha;
  else if (su == "kilojoule_per_mole")
    return kJ_mol;
  else if (su == "kelvin")
    return K;
  else
  {
    APP_ABORT("units::energy_unit\n  invalid energy unit: " + su +
              "\n  valid options are: J/joule, eV/electron_volt, Ry/rydberg, Ha/hartree, kJ/mol/kilo_joule_per_mole, "
              "K/kelvin");
    return J;
  }
}


} // namespace Units
} // namespace qmcplusplus
