//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    


#ifndef QMCPLUSPLUS_UNIT_CONVERSION_H
#define QMCPLUSPLUS_UNIT_CONVERSION_H

#include <Configuration.h>

namespace qmcplusplus
{

  namespace Units
  {
    typedef QMCTraits::RealType real;

    namespace constants
    {
      const real kb = 1.3806503e-23; // J/K
    }

    namespace count
    {
      const real mol = 6.0221415e23;
    }

    namespace distance
    {
      const real m  = 1.e0;
      const real A  = 1.e-10*m;
      const real B  = .52917720859e-10*m;
      const real nm = 1.e-9*m;
      const real pm = 1.e-12*m;
      const real fm = 1.e-15*m;
    }

    namespace time
    {
      const real s  = 1.e0;
      const real ms = 1.e-3*s;
      const real ns = 1.e-9*s;
      const real ps = 1.e-12*s;
      const real fs = 1.e-15*s;
    }

    namespace mass
    {    
      const real kg  = 1.e0;
      const real me  = 9.10938291e-31*kg;
      const real mp  = 1.672621777e-27*kg;
      const real amu = 1.660538921e-27*kg;
      const real Da  = amu;
    }

    namespace energy
    {
      using count::mol;
      using constants::kb;
      const real J      = 1.e0;
      const real eV     = 1.60217646e-19*J;
      const real Ry     = 13.6056923*eV;
      const real Ha     = 2*Ry;
      const real kJ_mol = 1000.*J/mol;
      const real K      = J/kb;
    }

    namespace charge
    {
      const real C = 1.e0;
      const real e = 1.60217646e-19*C;
    }

    namespace pressure
    {
      const real Pa   = 1.e0;
      const real bar  = 1.e5*Pa;
      const real Mbar = 1.e6*bar;
      const real GPa  = 1.e9*Pa;
      const real atm  = 1.01325e5*Pa;
    }

    namespace force
    {
      const real N  = 1.e0;
      const real pN = 1e-12*N;
    }










    enum units
    {
      mol=0,
      A,B,m,nm,pm,fm,
      s,ms,ns,ps,fs,
      kg,me,mp,amu,Da,
      J,eV,Ry,Ha,kJ_mol,K,
      C,e,
      Pa,bar,Mbar,GPa,atm,
      N,pN,
      mole,
      angstrom,bohr,meter,nanometer,picometer,femtometer,
      second,millisecond,nanosecond,picosecond,femtosecond,
      kilogram,electron_mass,proton_mass,atomic_mass_unit,dalton,
      joule,electron_volt,rydberg,hartree,kilojoule_per_mole,kelvin,
      coulomb,proton_charge,
      pascal,megabar,gigapascal,atmosphere,
      newton,piconewton,
      nunits
    };


    const real unit_values[nunits] = {
      count::mol,
      distance::A,distance::B,distance::m,distance::nm,distance::pm,distance::fm,
      time::s,time::ms,time::ns,time::ps,time::fs,
      mass::kg,mass::me,mass::mp,mass::amu,mass::Da,
      energy::J,energy::eV,energy::Ry,energy::Ha,energy::kJ_mol,energy::K,
      charge::C,charge::e,
      pressure::Pa,pressure::bar,pressure::Mbar,pressure::GPa,pressure::atm,
      force::N,force::pN,
      count::mol,
      distance::A,distance::B,distance::m,distance::nm,distance::pm,distance::fm,
      time::s,time::ms,time::ns,time::ps,time::fs,
      mass::kg,mass::me,mass::mp,mass::amu,mass::Da,
      energy::J,energy::eV,energy::Ry,energy::Ha,energy::kJ_mol,energy::K,
      charge::C,charge::e,
      pressure::Pa,pressure::Mbar,pressure::GPa,pressure::atm,
      force::N,force::pN
    };



    inline real convert(real value,units units_in,units units_out)
    {
      return value*unit_values[units_in]/unit_values[units_out];
    }


    template<typename array>
    inline void convert_array(array& values,units units_in,units units_out)
    {
      real conv = unit_values[units_in]/unit_values[units_out];
      typename array::iterator v;
      for(v=values.begin();v!=values.end();++v)
        (*v) *= conv;
    }

    
    /// convert from std::string to count unit
    units count_unit(const std::string& su);

    /// convert from std::string to distance unit
    units distance_unit(const std::string& su);

    /// convert from std::string to time unit
    units time_unit(const std::string& su);

    /// convert from std::string to mass unit
    units mass_unit(const std::string& su);

    /// convert from std::string to energy unit
    units energy_unit(const std::string& su);

    /// convert from std::string to charge unit
    units charge_unit(const std::string& su);

    /// convert from std::string to pressure unit
    units pressure_unit(const std::string& su);

    /// convert from std::string to force unit
    units force_unit(const std::string& su);

  }

}


#endif
