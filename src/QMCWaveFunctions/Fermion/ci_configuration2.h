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
#ifndef QMCPLUSPLUS_CI_CONFIGURATION2_H
#define QMCPLUSPLUS_CI_CONFIGURATION2_H
#include <vector>
#include <algorithm>
#include <iostream>

namespace qmcplusplus
{

// Defines a ci configuration based of the occupied MOs
struct ci_configuration2
{
  // occupied orbitals
  vector<int> occup;

  ci_configuration2() {}

  ci_configuration2(vector<int> &v): occup(v) {}
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
  int calculateNumOfExcitations(const ci_configuration2& c) const
  {
    if(occup.size() != c.occup.size())
    {
      APP_ABORT("ci_configuration2::operator==() - ci_configuration2s are not compatible.");
    }
    int n=0;
    for(int i=0; i<occup.size(); i++)
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
  double calculateExcitations(const ci_configuration2& c, int &n, vector<int>& pos, vector<int>& ocp, vector<int>& uno) const
  {
    if(occup.size() != c.occup.size())
    {
      APP_ABORT("ci_configuration2::operator==() - ci_configuration2s are not compatible.");
    }
    n=0;
    for(int i=0; i<occup.size(); i++)
    {
      if(std::find(c.occup.begin(),c.occup.end(),occup[i]) == c.occup.end())
      {
        pos[n] = i;
        ocp[n++] = occup[i];
      }
    }
    if(n==0)
      return 1.0;
    int cnt=0;
    for(int i=0; i<c.occup.size(); i++)
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
    vector<int> ref0(occup);
    for(int i=0; i<n; i++)
      ref0[pos[i]] = uno[i];
    for(int i=0; i<ref0.size(); i++)
      for(int k=i+1; k<ref0.size(); k++)
      {
        if(ref0[i] > ref0[k])
        {
          int q=ref0[i];
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
  for(int i=0; i<c.occup.size(); i++)
    out <<c.occup[i] <<" ";
  out<<endl;
  return out;
};

}
#endif
