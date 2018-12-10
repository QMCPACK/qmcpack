//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Miguel Morales, moralessilva2@llnl.gov, Lawrence Livermore National Laboratory
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Miguel Morales, moralessilva2@llnl.gov, Lawrence Livermore National Laboratory
//////////////////////////////////////////////////////////////////////////////////////
    
    
#ifndef QMCPLUSPLUS_CI_CONFIGURATION_H
#define QMCPLUSPLUS_CI_CONFIGURATION_H
#include <vector>
#include <iostream>

namespace qmcplusplus
{

// Defines a single CI ci_configuration, with respect to the hartree fock ci_configuration.
struct ci_configuration
{
  // vector of bits, each bit determines whether the corresponding state is occupied or not
  std::vector<bool> occup;
  bool taken;
  int nExct; // with respect to base ci_configuration, which we assume is hf

  ci_configuration(): taken(false),nExct(0) {}

  ci_configuration(std::vector<bool> &v, int n): occup(v),taken(false),nExct(n) {}
  ci_configuration(const ci_configuration& c):occup(c.occup),taken(c.taken),nExct(c.nExct) {}

  ~ci_configuration() {}

  bool operator==(const ci_configuration& c) const
  {
    if(nExct!=c.nExct)
    {
      return false;
    }
    if(occup.size() != c.occup.size())
    {
      APP_ABORT("ci_configuration::operator==() - ci_configurations are not compatible.");
    }
    if(count() != c.count())
    {
      app_log() <<"c0: ";
      for(int i=0; i<occup.size(); i++)
        app_log() <<occup[i];
      app_log() << std::endl <<"c1: ";
      for(int i=0; i<c.occup.size(); i++)
        app_log() <<c.occup[i];
      app_log() << std::endl;
      APP_ABORT("ci_configuration::operator==() - ci_configurations are not compatible. Unequal number of occupied states. ");
    }
    for(int i=0; i<occup.size(); i++)
    {
      if(occup[i]^c.occup[i])
      {
        return false;
      }
    }
    return true;
  }

  // this has a very specific use below
  bool isSingle(const ci_configuration &c, int &rem, int &add) const
  {
    if(c.nExct-nExct != 1)
      return false;
    if(occup.size() != c.occup.size())
    {
      APP_ABORT("ci_configuration::isSingle() - ci_configurations are not compatible.");
    }
    if(count() != c.count())
    {
      APP_ABORT("ci_configuration::isSingle() - ci_configurations are not compatible. Unequal number of occupied states. ");
    }
    int nr=0,na=0,r=-1,a=-1;
    for(int i=0; i<occup.size(); i++)
    {
      if(occup[i]^c.occup[i])
      {
        if(occup[i])
        {
          nr++;
          r=i;
        }
        else
        {
          na++;
          a=i;
        }
      }
    }
    if(na == 1 && nr == 1)
    {
      rem=r;
      add=a;
      return true;
    }
    else
      return false;
  }

  int count() const
  {
    int res=0;
    for(int i=0; i<occup.size(); i++)
      if(occup[i])
        res++;
    return res;
  }

  void add_occupation(std::string &str) {
    int str_size=str.size();
    occup.resize(str_size);
    for (int i=0; i<str_size;i++)
      occup[i]=str[i]-'0';
  }


};

inline std::ostream& operator<<(std::ostream& out, const ci_configuration& c)
{
  out<<"ci ci_configuration: ";
  for(int i=0; i<c.occup.size(); i++)
    out <<c.occup[i];
  out<< std::endl;
  return out;
};
}
#endif
