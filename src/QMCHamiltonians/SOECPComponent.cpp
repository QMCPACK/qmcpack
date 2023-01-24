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
SOECPComponent::SOECPComponent() : lmax(0), nchannel(0), nknot(0), sknot(0), total_knots(0), Rmax(-1), VP(nullptr) {}

SOECPComponent::~SOECPComponent()
{
  for (int i = 0; i < sopp_m.size(); i++)
    delete sopp_m[i];
  if (VP)
    delete VP;
}

void SOECPComponent::print(std::ostream& os) {}

void SOECPComponent::initVirtualParticle(const ParticleSet& qp)
{
  assert(VP == nullptr);
  outputManager.pause();
  VP = new VirtualParticleSet(qp, total_knots);
  outputManager.resume();
}

void SOECPComponent::deleteVirtualParticle()
{
  if (VP)
    delete VP;
  VP = nullptr;
}

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
  if (VP)
    myclone->VP = new VirtualParticleSet(qp, nknot);
  return myclone;
}

void SOECPComponent::resize_warrays(int n, int m, int s)
{
  vrad.resize(m);
  rrotsgrid_m.resize(n);
  nchannel = sopp_m.size();
  nknot    = sgridxyz_m.size();
  sknot    = s;
  //Need +1 for Simpsons rule to include both end points.
  //sknot here refers to the number of subintervals for integration
  total_knots = nknot * (sknot + 1);
  psiratio.resize(total_knots);
  deltaV.resize(total_knots);
  deltaS.resize(total_knots);
  spin_quad_weights.resize(total_knots);
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
  case 1:
    return ComplexType(std::sin(s1 + s2), 0.0);
  case 2:
    return ComplexType(0.0, std::sin(s1 - s2));
  default:
    APP_ABORT("SOECPComponent::sMatrixElements invalid operator dimension\n");
    return 0;
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
  case 1:
    val = onehalf *
        ComplexType(0.0,
                    std::sqrt(l * (l + 1) - m2 * (m2 - 1)) * kroneckerDelta(m1, m2 - 1) -
                        std::sqrt(l * (l + 1) - m2 * (m2 + 1)) * kroneckerDelta(m1, m2 + 1));
    return val;
  case 2:
    val = ComplexType(m2 * kroneckerDelta(m1, m2), zero);
    return val;
  default:
    APP_ABORT("SOECPComponent::lMatrixElements Invalid operator dimension\n");
    return 0;
  }
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

  RealType sold = W.spins[iel];
  buildTotalQuadrature(r, dr, sold);

  if (VP)
  {
    VP->makeMovesWithSpin(W, iel, deltaV, deltaS, true, iat);
    Psi.evaluateRatios(*VP, psiratio);
  }
  else
  {
    for (int iq = 0; iq < total_knots; iq++)
    {
      W.makeMoveWithSpin(iel, deltaV[iq], deltaS[iq]);
      psiratio[iq] = Psi.calcRatio(W, iel);
      W.rejectMove(iel);
      Psi.resetPhaseDiff();
    }
  }

  ComplexType pairpot;
  for (int iq = 0; iq < total_knots; iq++)
  {
    RealType snew = sold + deltaS[iq];
    ComplexType lsum;
    for (int il = 0; il < nchannel; il++)
    {
      int l = il + 1; //nchannels starts at l=1, so 0th element is p not s
      ComplexType msums;
      for (int m1 = -l; m1 <= l; m1++)
      {
        ComplexType Y = sphericalHarmonic(l, m1, dr);
        for (int m2 = -l; m2 <= l; m2++)
        {
          ComplexType ldots;
          for (int id = 0; id < 3; id++)
            ldots += lmMatrixElements(l, m1, m2, id) * sMatrixElements(sold, snew, id);
          ComplexType cY = std::conj(sphericalHarmonic(l, m2, rrotsgrid_m[iq % nknot]));
          msums += Y * cY * ldots;
        }
      }
      lsum += vrad[il] * msums;
    }
    pairpot += psiratio[iq] * lsum * spin_quad_weights[iq];
  }
  return std::real(pairpot);
}

void SOECPComponent::buildTotalQuadrature(const RealType r, const PosType& dr, const RealType sold)
{
  int count = 0;
  RealType smin(0.0);
  RealType smax(TWOPI);
  RealType dS = (smax - smin) / sknot; //step size for spin

  auto addSpatialQuadrature = [&](const int is, const RealType r, const PosType& dr, const RealType ds,
                                  const RealType spin_weight) {
    for (int iq = 0; iq < nknot; iq++)
    {
      int offset     = is * nknot + iq;
      deltaV[offset] = r * rrotsgrid_m[iq] - dr;
      deltaS[offset] = ds;
      //spin integral norm is 1/(2.0 * pi), spatial integral has 4.0 * pi, so overall weight is 2
      spin_quad_weights[offset] = 2.0 * spin_weight * sgridweight_m[iq];
      count++;
    }
  };

  //simpsons 1/3 rule for spin integral

  //odd points
  for (int is = 1; is <= sknot - 1; is += 2)
  {
    RealType snew = smin + is * dS;
    addSpatialQuadrature(is, r, dr, snew - sold, RealType(4. / 3.) * dS);
  }

  //even points
  for (int is = 2; is <= sknot - 2; is += 2)
  {
    RealType snew = smin + is * dS;
    addSpatialQuadrature(is, r, dr, snew - sold, RealType(2. / 3.) * dS);
  }

  //end points
  addSpatialQuadrature(0, r, dr, smin - sold, RealType(1. / 3.) * dS);
  addSpatialQuadrature(sknot, r, dr, smax - sold, RealType(1. / 3.) * dS);

  assert(count == total_knots);
}

void SOECPComponent::rotateQuadratureGrid(const TensorType& rmat)
{
  for (int i = 0; i < sgridxyz_m.size(); i++)
    rrotsgrid_m[i] = dot(rmat, sgridxyz_m[i]);
}

} // namespace qmcplusplus
