//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2020 QMCPACK developers.
//
// File developed by: Raymond C. Clay, rclay@sandia.gov, Sandia National Laboratories
//                    Cody A. Melton, cmelton@sandia.gov, Sandia National Laboratories
//
// File created by: Raymond C. Clay, rclay@sandia.gov, Sandia National Laboratories
//////////////////////////////////////////////////////////////////////////////////////


#include "Particle/DistanceTable.h"
#include "SOECPComponent.h"
#include "Numerics/Ylm.h"

namespace qmcplusplus
{
SOECPComponent::SOECPComponent() : lmax(0), nchannel(0), nknot(0), sknot(0), Rmax(-1) {}

SOECPComponent::~SOECPComponent()
{
  for (int i = 0; i < sopp_m.size(); i++)
    delete sopp_m[i];
}

void SOECPComponent::print(std::ostream& os) {}

void SOECPComponent::add(int l, RadialPotentialType* pp)
{
  angpp_m.push_back(l);
  sopp_m.push_back(pp);
}

SOECPComponent* SOECPComponent::makeClone(const ParticleSet& qp)
{
  SOECPComponent* myclone = new SOECPComponent(*this);
  for (int i = 0; i < sopp_m.size(); i++)
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
  nknot    = sgridxyz_m.size();
  sknot    = s;
  if (m != nchannel)
  {
    APP_ABORT("SOECPComponent::resize_warrays has incorrect number of radial channels\n");
  }
}

int SOECPComponent::kroneckerDelta(int x, int y) { return (x == y) ? 1 : 0; }

SOECPComponent::ComplexType SOECPComponent::sMatrixElements(RealType s1, RealType s2, int dim)
{
  switch (dim)
  {
  case 0:
    return ComplexType(std::cos(s1 + s2), 0.0);
    break;
  case 1:
    return ComplexType(std::sin(s1 + s2), 0.0);
    break;
  case 2:
    return ComplexType(0.0, std::sin(s1 - s2));
    break;
  default:
    APP_ABORT("SOECPComponent::sMatrixElements invalid operator dimension\n");
    return 0;
    break;
  }
}

SOECPComponent::ComplexType SOECPComponent::lmMatrixElements(int l, int m1, int m2, int dim)
{
  ComplexType val;
  RealType onehalf = 0.5;
  RealType zero    = 0.0;
  switch (dim)
  {
  case 0: //x
    val = onehalf *
        ComplexType(std::sqrt(l * (l + 1) - m2 * (m2 + 1)) * kroneckerDelta(m1, m2 + 1) +
                        std::sqrt(l * (l + 1) - m2 * (m2 - 1)) * kroneckerDelta(m1, m2 - 1),
                    zero);
    return val;
    break;
  case 1:
    val = onehalf *
        ComplexType(0.0,
                    std::sqrt(l * (l + 1) - m2 * (m2 - 1)) * kroneckerDelta(m1, m2 - 1) -
                        std::sqrt(l * (l + 1) - m2 * (m2 + 1)) * kroneckerDelta(m1, m2 + 1));
    return val;
    break;
  case 2:
    val = ComplexType(m2 * kroneckerDelta(m1, m2), zero);
    return val;
    break;
  default:
    APP_ABORT("SOECPComponent::lMatrixElements Invalid operator dimension\n");
    return 0;
    break;
  }
}

SOECPComponent::ComplexType SOECPComponent::getAngularIntegral(RealType sold,
                                                               RealType snew,
                                                               ParticleSet& W,
                                                               TrialWaveFunction& Psi,
                                                               int iel,
                                                               RealType r,
                                                               const PosType& dr)
{
  //quadrature sum for angular integral
  constexpr RealType fourpi = 2.0 * TWOPI;
  for (int j = 0; j < nknot; j++)
  {
    deltaV[j] = r * rrotsgrid_m[j] - dr;
    W.makeMoveWithSpin(iel, deltaV[j], snew - sold);
    psiratio[j] = Psi.calcRatio(W, iel) * sgridweight_m[j] * fourpi;
    W.rejectMove(iel);
    Psi.resetPhaseDiff();
  }

  ComplexType angint(0.0);
  for (int j = 0; j < nknot; j++)
  {
    ComplexType lsum(0.0);
    for (int il = 0; il < nchannel; il++)
    {
      int l = il + 1; //nchannels starts at l=1, so 0th element is p not s
      ComplexType msums(0.0);
      for (int m1 = -l; m1 <= l; m1++)
      {
        for (int m2 = -l; m2 <= l; m2++)
        {
          ComplexType ldots(0.0);
          for (int d = 0; d < 3; d++)
            ldots += lmMatrixElements(l, m1, m2, d) * sMatrixElements(sold, snew, d);
          //Seemingly Numerics/Ylm takes unit vector with order z,x,y...why
          RealType rmag = std::sqrt(dr[0] * dr[0] + dr[1] * dr[1] + dr[2] * dr[2]);
          PosType rr    = dr / rmag;
          PosType tmp;
          tmp[0]         = rr[2];
          tmp[1]         = rr[0];
          tmp[2]         = rr[1];
          ComplexType Y  = Ylm(l, m1, tmp);
          tmp[0]         = rrotsgrid_m[j][2];
          tmp[1]         = rrotsgrid_m[j][0];
          tmp[2]         = rrotsgrid_m[j][1];
          ComplexType cY = std::conj(Ylm(l, m2, tmp));
          msums += Y * cY * ldots;
        }
      }
      lsum += vrad[il] * msums;
    }
    angint += psiratio[j] * lsum;
  }
  return angint;
}

SOECPComponent::RealType SOECPComponent::evaluateOne(ParticleSet& W,
                                                     int iat,
                                                     TrialWaveFunction& Psi,
                                                     int iel,
                                                     RealType r,
                                                     const PosType& dr)
{
  if (sknot < 2)
    APP_ABORT("Spin knots must be greater than 2\n");

  if (sknot % 2 != 0)
    APP_ABORT("Spin knots uses Simpson's rule. Must have even number of knots");

  for (int ip = 0; ip < nchannel; ip++)
  {
    vrad[ip] = sopp_m[ip]->splint(r);
  }

  RealType smin(0.0);
  RealType smax(TWOPI);
  RealType dS = (smax - smin) / sknot; //step size for spin

  RealType sold = W.spins[iel];
  ComplexType sint(0.0);

  for (int is = 1; is <= sknot - 1; is += 2)
  {
    RealType snew      = smin + is * dS;
    ComplexType angint = getAngularIntegral(sold, snew, W, Psi, iel, r, dr);
    sint += RealType(4. / 3.) * dS * angint;
  }
  for (int is = 2; is <= sknot - 2; is += 2)
  {
    RealType snew      = smin + is * dS;
    ComplexType angint = getAngularIntegral(sold, snew, W, Psi, iel, r, dr);
    sint += RealType(2. / 3.) * dS * angint;
  }
  sint += RealType(1. / 3.) * dS * getAngularIntegral(sold, smin, W, Psi, iel, r, dr);
  sint += RealType(1. / 3.) * dS * getAngularIntegral(sold, smax, W, Psi, iel, r, dr);

  RealType pairpot = std::real(sint) / TWOPI;
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
} // namespace qmcplusplus
