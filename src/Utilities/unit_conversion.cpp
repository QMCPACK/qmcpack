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
units count_unit(const std::string& su)
{
  if (su == "mol")
    return mol;
  else if (su == "mole")
    return mol;
  else
  {
    APP_ABORT("units::count_unit\n  invalid count unit: " + su + "\n  valid options are: mol");
    return mol;
  }
}


units distance_unit(const std::string& su)
{
  if (su == "m")
    return m;
  else if (su == "A")
    return A;
  else if (su == "B")
    return B;
  else if (su == "nm")
    return nm;
  else if (su == "pm")
    return pm;
  else if (su == "fm")
    return fm;
  else if (su == "meter")
    return m;
  else if (su == "angstrom")
    return A;
  else if (su == "bohr")
    return B;
  else if (su == "nanometer")
    return nm;
  else if (su == "picometer")
    return pm;
  else if (su == "femtometer")
    return fm;
  else
  {
    APP_ABORT("units::distance_unit\n  invalid distance unit: " + su +
              "\n  valid options are: m/meter, A/angstrom, B/bohr, nm/nanometer, pm/picometer, fm/femtometer");
    return m;
  }
}


units time_unit(const std::string& su)
{
  if (su == "s")
    return s;
  else if (su == "ms")
    return ms;
  else if (su == "ns")
    return ns;
  else if (su == "ps")
    return ps;
  else if (su == "fs")
    return fs;
  else if (su == "second")
    return s;
  else if (su == "millisecond")
    return ms;
  else if (su == "nanosecond")
    return ns;
  else if (su == "picosecond")
    return ps;
  else if (su == "femtosecond")
    return fs;
  else
  {
    APP_ABORT("units::time_unit\n  invalid time unit: " + su +
              "\n  valid options are: s/second, ms/millisecond, ns/nanosecond, ps/picosecond, fs/femtosecond");
    return s;
  }
}


units mass_unit(const std::string& su)
{
  if (su == "kg")
    return kg;
  else if (su == "me")
    return me;
  else if (su == "mp")
    return mp;
  else if (su == "amu")
    return amu;
  else if (su == "Da")
    return Da;
  else if (su == "kilogram")
    return kg;
  else if (su == "electron_mass")
    return me;
  else if (su == "proton_mass")
    return mp;
  else if (su == "atomic_mass_unit")
    return amu;
  else if (su == "dalton")
    return Da;
  else
  {
    APP_ABORT("units::mass_unit\n  invalid mass unit: " + su +
              "\n  valid options are: kg/kilogram, me/electron_mass, mp/proton_mass, amu/atomic_mass_unit, Da/dalton");
    return kg;
  }
}


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


units charge_unit(const std::string& su)
{
  if (su == "C")
    return C;
  else if (su == "e")
    return e;
  else if (su == "coulomb")
    return C;
  else if (su == "proton_charge")
    return e;
  else
  {
    APP_ABORT("units::charge_unit\n  invalid charge unit: " + su + "\n  valid options are: C/coulomb, e/proton_charge");
    return C;
  }
}


units pressure_unit(const std::string& su)
{
  if (su == "Pa")
    return Pa;
  else if (su == "bar")
    return bar;
  else if (su == "Mbar")
    return Mbar;
  else if (su == "GPa")
    return GPa;
  else if (su == "atm")
    return atm;
  else if (su == "pascal")
    return Pa;
  else if (su == "megabar")
    return Mbar;
  else if (su == "gigapascal")
    return GPa;
  else if (su == "atmosphere")
    return atm;
  else
  {
    APP_ABORT("units::pressure_unit\n  invalid pressure unit: " + su +
              "\n  valid options are: Pa/pascal, bar/bar, Mbar/megabar, GPa/gigapascal, atm/atmosphere");
    return Pa;
  }
}


units force_unit(const std::string& su)
{
  if (su == "N")
    return N;
  else if (su == "pN")
    return pN;
  else if (su == "newton")
    return N;
  else if (su == "piconewton")
    return pN;
  else
  {
    APP_ABORT("units::force_unit\n  invalid force unit: " + su + "\n  valid options are: N/newton, pN/piconewton");
    return N;
  }
}

} // namespace Units
} // namespace qmcplusplus
