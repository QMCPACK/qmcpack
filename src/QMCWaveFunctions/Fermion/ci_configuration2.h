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
    
    
#ifndef QMCPLUSPLUS_CI_CONFIGURATION2_H
#define QMCPLUSPLUS_CI_CONFIGURATION2_H
//#include <vector>
#include <simd/simd.hpp>
#include <algorithm>
#include <iostream>

namespace qmcplusplus
{

// Defines a ci configuration based of the occupied MOs
struct ci_configuration2
{
  // occupied orbitals
  std::vector<size_t> occup;

  ci_configuration2() {}

  ci_configuration2(std::vector<size_t> &v): occup(v.begin(),v.end()) {}
  ci_configuration2(const ci_configuration2& c):occup(c.occup) {}

  ~ci_configuration2() {}

  inline bool operator==(const ci_configuration2& c) const
  {
    if(occup.size() != c.occup.size())
    {
      APP_ABORT("ci_configuration2::operator==() - ci_configuration2s are not compatible.");
    }
    for(int i=0; i<occup.size(); i++)
    {
      if(occup[i] != c.occup[i])
      {
        return false;
      }
    }
    return true;
  }

  /*
  * The routine determines the number excitations that
  * transform the given configuration (argument) into (*this).
  */
  size_t calculateNumOfExcitations(const ci_configuration2& c) const
  {
    if(occup.size() != c.occup.size())
    {
      APP_ABORT("ci_configuration2::operator==() - ci_configuration2s are not compatible.");
    }
    size_t n=0;
    for(size_t i=0; i<occup.size(); i++)
    {
      if(std::find(c.occup.begin(),c.occup.end(),occup[i]) == c.occup.end())
        n++;
    }
    return n;
  }

  /*
   * The routine determines the excitations that
   * transform the given configuration (argument) into (*this).
   * The sign of the permutation is returned.
   *    - n: number of excitations
   *    - pos: positions of MOs in c that need to be changed (not the label of the MO)
   *    - ocp: label of the MO being removed (needed to calculate matrix elements)
   *    - uno: label of the MO that replaces ocp[i] (should be the same as the position in array, as opposed to ocp)
   *
   */
  double calculateExcitations(const ci_configuration2& c, size_t &n, 
      std::vector<size_t>& pos, std::vector<size_t>& ocp, std::vector<size_t>& uno) const
  {
    if(occup.size() != c.occup.size())
    {
      APP_ABORT("ci_configuration2::operator==() - ci_configuration2s are not compatible.");
    }
    n=0;
    for(size_t i=0; i<occup.size(); i++)
    {
      if(std::find(c.occup.begin(),c.occup.end(),occup[i]) == c.occup.end())
      {
        pos[n] = i;
        ocp[n++] = occup[i];
      }
    }
    if(n==0)
      return 1.0;
    size_t cnt=0;
    for(size_t i=0; i<c.occup.size(); i++)
    {
      if(std::find(occup.begin(),occup.end(),c.occup[i]) == occup.end())
        uno[cnt++] = c.occup[i];
    }
    if(cnt != n)
    {
      APP_ABORT(" Error #1 in ci_configuration2::calculateExcitations() \n");
    }
    double res=1.0;
    // this is needed because ci coefficients are given wrt standard ordering,
    // but by defining the determinant through excitations from a reference might change
    // the parity
    // inefficient but easy, is there a sort routine in STL that gives me the parity too???
    auto  ref0(occup);
    for(size_t i=0; i<n; i++)
      ref0[pos[i]] = uno[i];
    for(size_t i=0; i<ref0.size(); i++)
      for(size_t k=i+1; k<ref0.size(); k++)
      {
        if(ref0[i] > ref0[k])
        {
          size_t q=ref0[i];
          ref0[i]=ref0[k];
          ref0[k]=q;
          res*=-1.0;
        }
        else
          if(ref0[i] == ref0[k])
          {
            APP_ABORT(" Error #2 in ci_configuration2::calculateExcitations() \n");
          }
      }
    return res;
  }


};

inline std::ostream& operator<<(std::ostream& out, const ci_configuration2& c)
{
  out<<"ci ci_configuration2: ";
  for(size_t i=0; i<c.occup.size(); i++)
    out <<c.occup[i] <<" ";
  out<< std::endl;
  return out;
};

}
#endif
