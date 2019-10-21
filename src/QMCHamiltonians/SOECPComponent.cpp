//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//                    Cody A. Melton, cmelton@sandia.gov, Sandia National Laboratories
//                    Raymond C. Clay, rclay@sandia.gov, Sandia National Laboratories
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


#include "Particle/DistanceTableData.h"
#include "QMCHamiltonians/SOECPComponent.h"
#include "Numerics/Ylm.h"

namespace qmcplusplus
{

void SOECPComponent::add(int l, RadialPotentialType* pp)
{
  angpp_m.push_back(l);
  sopp_m.push_back(pp);
}

int SOECPComponent::kroneckerDelta(int x, int y)
{
  if (x == y)
    return 1;
  else
    return 0;
}


SOECPComponent::ComplexType SOECPComponent::lmMatrixElements(int l, int m1, int m2, int dim)
{
  ComplexType val;
  switch(dim)
  {
  case 0: //x
    val = 0.5*ComplexType(std::sqrt(l*(l+1)-m2*(m2+1))*kroneckerDelta(m1,m2+1) + std::sqrt(l*(l+1)-m2*(m2-1))*kroneckerDelta(m1,m2-1),0.0);
    return val;
    break;
  case 1:
    val = 0.5*ComplexType(0.0,std::sqrt(l*(l+1)-m2*(m2-1))*kroneckerDelta(m1,m2-1)-std::sqrt(l*(l+1)-m2*(m2+1))*kroneckerDelta(m1,m2+1));
    return val;
    break;
  case 2:
    val = ComplexType(m2*kroneckerDelta(m1,m2),0.0);
    return val;
    break;
  default:
      APP_ABORT("Invalid operator dimension for lmMatrixElements\n");
      return 0;
      break;
  }
}

SOECPComponent::RealType SOECPComponent::evaluateOne(ParticleSet& W,
                       int iat,
                       TrialWaveFunction & Psi,
                       int iel,
                       RealType r,
                       const PosType& dr)
{
  int sknot = 2; //number of spin integral points, needs to be set in input or by default
  RealType dS = 2*M_PI/sknot; //step size for spin

  RealType sold = W.spins[iel];

  ComplexType sint(0.0);
  //spin integral
  for (RealType snew=0.0; snew <= TWOPI; snew += dS)
  {
    ComplexType sx(std::cos(sold+snew),0.0);
    ComplexType sy(std::sin(sold+snew),0.0);
    ComplexType sz(0.0,std::sin(sold-snew));

    //quadrature sum for angular integral
    for (int j = 0; j < nknot; j++)
    {
      deltaV[j] = r*rrotsgrid_m[j] - dr;
      W.makeMove(iel,deltaV[j]);
      //W.makeSpinMove(iel,sp);
      psiratio[j] = Psi.calcRatio(W,iel)*sgridweight_m[j];
      W.rejectMove(iel);
      //W.rejectSpinMove(iep);
      Psi.resetPhaseDiff();
    }

    // Compute radial potential, do not multiply by (2*l+1)
    for (int ip = 0; ip<nchannel; ip++)
      vrad[ip] = sopp_m[ip]->splint(r); //*wgt_angpp_m[ip]

    ComplexType angint(0.0);
    for (int j = 0; j < nknot; j++)
    {
      ComplexType lsum(0.0);
      for (int l=0; l < nchannel; l++)
      {
        ComplexType msums(0.0);
        for (int m1 = -l; m1 <= l; m1++)
        {
          for (int m2 = -l; m2 <= l; m2++)
          {
            ComplexType ldots(0.0);
            ldots += lmMatrixElements(l,m1,m2,0)*sx;
            ldots += lmMatrixElements(l,m1,m2,1)*sy;
            ldots += lmMatrixElements(l,m1,m2,2)*sz;
            msums += std::conj(Ylm(l,m1,dr))*Ylm(l,m2,r*rrotsgrid_m[j])*ldots;
          }
        }
        lsum += vrad[l]*msums;
      }
      angint += psiratio[j]*lsum;
    }
    sint += dS*angint;
  }
  return std::real(sint);
}
}
