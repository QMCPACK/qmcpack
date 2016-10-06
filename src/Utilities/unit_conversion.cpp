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
    
    

#include <Utilities/unit_conversion.h>
#include <Message/Communicate.h>


namespace qmcplusplus
{
  namespace Units
  {

    units count_unit(const std::string& su)
    {
      units u;
      if(su=="mol")
        u = mol;
      else if(su=="mole")
        u = mol;
      else
        APP_ABORT("units::count_unit\n  invalid count unit: "+su+"\n  valid options are: mol");
      return u;
    }

    
    units distance_unit(const std::string& su)
    {
      units u;
      if(su=="m")
        u = m; 
      else if(su=="A")
        u = A;
      else if(su=="B")
        u = B;
      else if(su=="nm")
        u = nm;
      else if(su=="pm")
        u = pm;
      else if(su=="fm")
        u = fm;
      else if(su=="meter")
        u = m;
      else if(su=="angstrom")
        u = A;
      else if(su=="bohr")
        u = B;
      else if(su=="nanometer")
        u = nm;
      else if(su=="picometer")
        u = pm;
      else if(su=="femtometer")
        u = fm;
      else
        APP_ABORT("units::distance_unit\n  invalid distance unit: "+su+"\n  valid options are: m/meter, A/angstrom, B/bohr, nm/nanometer, pm/picometer, fm/femtometer");
      return u;
    }

    
    units time_unit(const std::string& su)
    {
      units u;
      if(su=="s")
        u = s; 
      else if(su=="ms")
        u = ms;
      else if(su=="ns")
        u = ns;
      else if(su=="ps")
        u = ps;
      else if(su=="fs")
        u = fs;
      else if(su=="second")
        u = s; 
      else if(su=="millisecond")
        u = ms;
      else if(su=="nanosecond")
        u = ns;
      else if(su=="picosecond")
        u = ps;
      else if(su=="femtosecond")
        u = fs;
      else
        APP_ABORT("units::time_unit\n  invalid time unit: "+su+"\n  valid options are: s/second, ms/millisecond, ns/nanosecond, ps/picosecond, fs/femtosecond");
      return u;
    }

    
    units mass_unit(const std::string& su)
    {
      units u;
      if(su=="kg")
        u = kg; 
      else if(su=="me")
        u = me;
      else if(su=="mp")
        u = mp;
      else if(su=="amu")
        u = amu;
      else if(su=="Da")
        u = Da;
      else if(su=="kilogram")
        u = kg; 
      else if(su=="electron_mass")
        u = me;
      else if(su=="proton_mass")
        u = mp;
      else if(su=="atomic_mass_unit")
        u = amu;
      else if(su=="dalton")
        u = Da;
      else
        APP_ABORT("units::mass_unit\n  invalid mass unit: "+su+"\n  valid options are: kg/kilogram, me/electron_mass, mp/proton_mass, amu/atomic_mass_unit, Da/dalton");
      return u;
    }

    
    units energy_unit(const std::string& su)
    {
      units u;
      if(su=="J")
        u = J; 
      else if(su=="eV")
        u = eV;
      else if(su=="Ry")
        u = Ry;
      else if(su=="Ha")
        u = Ha;
      else if(su=="kJ/mol")
        u = kJ_mol;
      else if(su=="K")
        u = K;
      else if(su=="joule")
        u = J; 
      else if(su=="electron_volt")
        u = eV;
      else if(su=="rydberg")
        u = Ry;
      else if(su=="hartree")
        u = Ha;
      else if(su=="kilojoule_per_mole")
        u = kJ_mol;
      else if(su=="kelvin")
        u = K;
      else
        APP_ABORT("units::energy_unit\n  invalid energy unit: "+su+"\n  valid options are: J/joule, eV/electron_volt, Ry/rydberg, Ha/hartree, kJ/mol/kilo_joule_per_mole, K/kelvin");
      return u;
    }

    
    units charge_unit(const std::string& su)
    {
      units u;
      if(su=="C")
        u = C; 
      else if(su=="e")
        u = e;
      else if(su=="coulomb")
        u = C; 
      else if(su=="proton_charge")
        u = e;
      else
        APP_ABORT("units::charge_unit\n  invalid charge unit: "+su+"\n  valid options are: C/coulomb, e/proton_charge");
      return u;
    }

    
    units pressure_unit(const std::string& su)
    {
      units u;
      if(su=="Pa")
        u = Pa; 
      else if(su=="bar")
        u = bar;
      else if(su=="Mbar")
        u = Mbar;
      else if(su=="GPa")
        u = GPa;
      else if(su=="atm")
        u = atm;
      else if(su=="pascal")
        u = Pa; 
      else if(su=="megabar")
        u = Mbar;
      else if(su=="gigapascal")
        u = GPa;
      else if(su=="atmosphere")
        u = atm;
      else
        APP_ABORT("units::pressure_unit\n  invalid pressure unit: "+su+"\n  valid options are: Pa/pascal, bar/bar, Mbar/megabar, GPa/gigapascal, atm/atmosphere");
      return u;
    }

    
    units force_unit(const std::string& su)
    {
      units u;
      if(su=="N")
        u = N; 
      else if(su=="pN")
        u = pN;
      else if(su=="newton")
        u = N; 
      else if(su=="piconewton")
        u = pN;
      else
        APP_ABORT("units::force_unit\n  invalid force unit: "+su+"\n  valid options are: N/newton, pN/piconewton");
      return u;
    }

  }
}
