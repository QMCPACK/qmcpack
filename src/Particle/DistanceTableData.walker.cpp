//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    



#include "Particle/ParticleSet.h"
#include "Particle/DistanceTableData.h"
using namespace qmcplusplus;

void SymmetricDTD::reset(int m, int nactive)
{
  if(m != N[SourceIndex] || nactive != N[WalkerIndex])
  {
    N[SourceIndex]=m;
    N[VisitorIndex]=m;
    resize(m*(m-1)/2,nactive);
    M.resize(m+1);
    J.resize(m*(m-1));
    M[0] = 0;
    int ij = 0;
    for(int i=0; i<m; i++)
    {
      for(int j=i+1; j<m; j++, ij++)
      {
        J[ij] = j;
      }
      M[i+1] = ij;
    }
  }
}

///evaluate the Distance Table using a set of Particle Positions
void SymmetricDTD::evaluate(const PosVector_t& a, int visitors, int copies)
{
  ///number of columns
  reset(visitors,copies);
  int ij=0;
  for(int i=0; i<visitors; i++)
  {
    for(int j=i+1; j<visitors; j++)
    {
      int i0 = i*copies;
      int j0 = j*copies;
      for(int iw=0; iw<copies; iw++, ij++, i0++, j0++)
      {
        SPPosition_t drij = a(j0)-a(i0);
        value_type sep = sqrt(dot(drij,drij));
        r(ij) = sep;
        rinv(ij) = 1.0/sep;
        dr(ij) = drij;
      }
    }
  }
}

void AsymmetricDTD::reset(int n1, int n2, int nactive)
{
  if( n1!=N[SourceIndex] || n2 != N[VisitorIndex] || nactive != N[WalkerIndex])
  {
    N[SourceIndex] = n1;
    N[VisitorIndex] = n2;
    int m = n1*n2;
    if(m)
    {
      resize(m,nactive);
      M.resize(n1+1);
      J.resize(m);
      M[0] = 0;
      int ij = 0;
      for(int i=0; i<n1; i++)
      {
        for(int j=0; j<n2; j++, ij++)
        {
          J[ij] = j;
        }
        M[i+1] = M[i]+n2;
      }
    }
  }
}

void AsymmetricDTD::evaluate(const PosVector_t& a, int visitors, int copies)
{
  ///number of columns
  int ns = Origin.getTotalNum();
  reset(ns,visitors,copies);
  int ij = 0;
  for(int i=0; i<ns; i++)
  {
    SPPosition_t r0 = Origin.R(i);
    int j0=0;
    for(int j=0; j<visitors; j++)
    {
      for(int ia=0; ia<copies; ia++, ij++, j0++)
      {
        SPPosition_t drij = a(j0)-r0;
        value_type sep = sqrt(dot(drij,drij));
        r(ij) = sep;
        rinv(ij) = 1.0/sep;
        dr(ij) = drij;
      }
    }
  }
}

