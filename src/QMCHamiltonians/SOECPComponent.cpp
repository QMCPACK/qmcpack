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

SOECPComponent::SOECPComponent() 
    : lmax(0), nchannel(0), nknot(0), sknot(0), Rmax(-1)
{}

SOECPComponent::~SOECPComponent()
{
  for (int i=0; i<sopp_m.size(); i++)
    delete sopp_m[i];
}

void SOECPComponent::print(std::ostream & os)
{

}

void SOECPComponent::add(int l, RadialPotentialType* pp)
{
  angpp_m.push_back(l);
  sopp_m.push_back(pp);
}

SOECPComponent* SOECPComponent::makeClone(const ParticleSet & qp)
{
  SOECPComponent *myclone = new SOECPComponent(*this);
  for (int i=0; i<sopp_m.size(); i++)
    myclone->sopp_m[i] = sopp_m[i]->makeClone();
  return myclone;
}

void SOECPComponent::resize_warrays(int n, int m, int s)
{
  psiratio.resize(n);
  deltaV.resize(n);
  vrad.resize(m);
  rrotsgrid_m.resize(n);
  nchannel = sopp_m.size();
  nknot = sgridxyz_m.size();
  sknot = s;
  if (m != nchannel) 
  {
    APP_ABORT("SOECPComponent::resize_warrays has incorrect number of radial channels\n");
  }
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
  if (sknot < 1)
    APP_ABORT("Spin knots must be greater than 1\n");
  RealType dS = TWOPI/sknot; //step size for spin

  RealType sold = W.spins[iel];

  ComplexType sint(0.0);
  //spin integral, with trapezoid rule
  RealType snew(0.0);
  for (int is = 0; is<= sknot; is++)
  {
    ComplexType sx(std::cos(sold+snew),0.0);
    ComplexType sy(std::sin(sold+snew),0.0);
    ComplexType sz(0.0,std::sin(sold-snew));

    //quadrature sum for angular integral
    for (int j = 0; j < nknot; j++)
    {
      deltaV[j] = r*rrotsgrid_m[j] - dr;
      W.makeMoveWithSpin(iel,deltaV[j],snew-sold);
      psiratio[j] = Psi.calcRatio(W,iel)*sgridweight_m[j];
      W.rejectMove(iel);
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
    if (is==0 || is==sknot) //trapezoidal rule for endpoints
      sint += dS*0.5*angint;
    else
      sint += dS*angint;

    snew += dS;
  }

  RealType pairpot = std::real(sint);

  return pairpot;
}

void SOECPComponent::randomize_grid(RandomGenerator_t& myRNG)
{
  RealType phi(TWOPI * myRNG()), psi(TWOPI * myRNG()), cth(myRNG() - 0.5);
  RealType sph(std::sin(phi)), cph(std::cos(phi)), sth(std::sqrt(1.0 - cth * cth)), sps(std::sin(psi)),
      cps(std::cos(psi));
  TensorType rmat(cph * cth * cps - sph * sps, sph * cth * cps + cph * sps, -sth * cps, -cph * cth * sps - sph * cps,
                  -sph * cth * sps + cph * cps, sth * sps, cph * sth, sph * sth, cth);
  for (int i = 0; i < sgridxyz_m.size(); i++)
    rrotsgrid_m[i] = dot(rmat, sgridxyz_m[i]);
}

template<typename T>
void SOECPComponent::randomize_grid(std::vector<T>& sphere, RandomGenerator_t& myRNG)
{
  RealType phi(TWOPI * myRNG()), psi(TWOPI * myRNG()), cth(myRNG() - 0.5);
  RealType sph(std::sin(phi)), cph(std::cos(phi)), sth(std::sqrt(1.0 - cth * cth)), sps(std::sin(psi)),
      cps(std::cos(psi));
  TensorType rmat(cph * cth * cps - sph * sps, sph * cth * cps + cph * sps, -sth * cps, -cph * cth * sps - sph * cps,
                  -sph * cth * sps + cph * cps, sth * sps, cph * sth, sph * sth, cth);
  SpherGridType::iterator it(sgridxyz_m.begin());
  SpherGridType::iterator it_end(sgridxyz_m.end());
  SpherGridType::iterator jt(rrotsgrid_m.begin());
  while (it != it_end)
  {
    *jt = dot(rmat, *it);
    ++it;
    ++jt;
  }
  //copy the randomized grid to sphere
  for (int i = 0; i < rrotsgrid_m.size(); i++)
    for (int j = 0; j < OHMMS_DIM; j++)
      sphere[OHMMS_DIM * i + j] = rrotsgrid_m[i][j];
}

template void SOECPComponent::randomize_grid(std::vector<float>& sphere, RandomGenerator_t& myRNG);
template void SOECPComponent::randomize_grid(std::vector<double>& sphere, RandomGenerator_t& myRNG);
}
