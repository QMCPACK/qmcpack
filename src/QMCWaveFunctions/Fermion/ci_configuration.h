//////////////////////////////////////////////////////////////////
// (c) Copyright 2008-  by Jeongnim Kim
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//   National Center for Supercomputing Applications &
//   Materials Computation Center
//   University of Illinois, Urbana-Champaign
//   Urbana, IL 61801
//   e-mail: jnkim@ncsa.uiuc.edu
//
// Supported by 
//   National Center for Supercomputing Applications, UIUC
//   Materials Computation Center, UIUC
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
#ifndef QMCPLUSPLUS_CI_CONFIGURATION_H
#define QMCPLUSPLUS_CI_CONFIGURATION_H
#include <vector>
#include <iostream>

namespace qmcplusplus
{

  // Defines a single CI configuration, with respect to the hartree fock configuration.
  struct configuration 
  {
    // vector of bits, each bit determines whether the corresponding state is occupied or not
    vector<bool> occup;
    bool taken;
    int nExct; // with respect to base configuration, which we assume is hf

    configuration(): taken(false),nExct(0) {} 

    configuration(vector<bool> &v, int n): occup(v),taken(false),nExct(n) {} 
    configuration(const configuration& c):occup(c.occup),taken(c.taken),nExct(c.nExct) {}  

    ~configuration() {} 

    bool operator==(const configuration& c) const {
      if(nExct!=c.nExct) { return false; }
      if(occup.size() != c.occup.size()) {
       APP_ABORT("configuration::operator==() - configurations are not compatible.");
      } 
      if(count() != c.count()) {
        app_log() <<"c0: ";
        for(int i=0; i<occup.size(); i++)
          app_log() <<occup[i];
        app_log() <<endl <<"c1: ";
        for(int i=0; i<c.occup.size(); i++)
          app_log() <<c.occup[i];
        app_log() <<endl;
        APP_ABORT("configuration::operator==() - configurations are not compatible. Unequal number of occupied states. ");
      }
      for(int i=0; i<occup.size(); i++) {
        if(occup[i]^c.occup[i]) { return false; }
      }
      return true;
    }

    // this has a very specific use below
    bool isSingle(const configuration &c, int &rem, int &add) const
    {
      if(c.nExct-nExct != 1) return false;
      if(occup.size() != c.occup.size()) { 
        APP_ABORT("configuration::isSingle() - configurations are not compatible.");
      }
      if(count() != c.count()) { 
        APP_ABORT("configuration::isSingle() - configurations are not compatible. Unequal number of occupied states. ");
      }
      int nr=0,na=0,r=-1,a=-1;    
      for(int i=0; i<occup.size(); i++)
      {
        if(occup[i]^c.occup[i]) {
          if(occup[i]) {
            nr++;
            r=i;
          } else { 
            na++;
            a=i;
          }
        }
      }
      if(na == 1 && nr == 1) {
        rem=r;
        add=a;
        return true;
      } else
        return false;
    }

    int count() const {
      int res=0;
      for(int i=0; i<occup.size(); i++) if(occup[i]) res++;
      return res;
    }

  }; 

  std::ostream&
  operator<<(std::ostream& out, const configuration& c)
  {
    out<<"ci configuration: ";
    for(int i=0; i<c.occup.size(); i++)
      out <<c.occup[i];
    out<<endl;
    return out;
  }

}
#endif
