//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


#include "QMCWaveFunctions/LatticeGaussianProduct.h"

namespace qmcplusplus
{
typedef LatticeGaussianProduct::ValueType ValueType;
typedef LatticeGaussianProduct::GradType GradType;

LatticeGaussianProduct::LatticeGaussianProduct(ParticleSet& centers, ParticleSet& ptcls) : CenterRef(centers)
{
  Optimizable    = false;
  ClassName      = "LatticeGaussianProduct";
  NumTargetPtcls = ptcls.getTotalNum();
  NumCenters     = centers.getTotalNum();
#ifdef ENABLE_SOA
  myTableID      = ptcls.addTable(CenterRef, DT_SOA);
#else
  myTableID      = ptcls.addTable(CenterRef, DT_AOS);
#endif
  U.resize(NumTargetPtcls);
  dU.resize(NumTargetPtcls);
  d2U.resize(NumTargetPtcls);
  FirstAddressOfdU = &(dU[0][0]);
  LastAddressOfdU  = FirstAddressOfdU + dU.size() * OHMMS_DIM;
}

LatticeGaussianProduct::~LatticeGaussianProduct() {}

//evaluate the distance table with P
void LatticeGaussianProduct::resetTargetParticleSet(ParticleSet& P) {}
void LatticeGaussianProduct::checkInVariables(opt_variables_type& active) {}
void LatticeGaussianProduct::checkOutVariables(const opt_variables_type& active) {}
void LatticeGaussianProduct::resetParameters(const opt_variables_type& active) {}
void LatticeGaussianProduct::reportStatus(std::ostream& os) {}

/**
     *@param P input configuration containing N particles
     *@param G a vector containing N gradients
     *@param L a vector containing N laplacians
     *@return The wavefunction value  \f$exp(-J({\bf R}))\f$
     *
     *Upon exit, the gradient \f$G[i]={\bf \nabla}_i J({\bf R})\f$
     *and the laplacian \f$L[i]=\nabla^2_i J({\bf R})\f$ are accumulated.
     *While evaluating the value of the Jastrow for a set of
     *particles add the gradient and laplacian contribution of the
     *Jastrow to G(radient) and L(aplacian) for local energy calculations
     *such that \f[ G[i]+={\bf \nabla}_i J({\bf R}) \f]
     *and \f[ L[i]+=\nabla^2_i J({\bf R}). \f]
     */
LatticeGaussianProduct::RealType LatticeGaussianProduct::evaluateLog(ParticleSet& P,
                                             ParticleSet::ParticleGradient_t& G,
                                             ParticleSet::ParticleLaplacian_t& L)
{
  const auto& d_table = P.getDistTable(myTableID);
  int icent           = 0;
  LogValue            = 0.0;
  RealType dist       = 0.0;
  PosType  disp       = 0.0;
  for (int iat = 0; iat < NumTargetPtcls; iat++)
  {
    U[iat]     = 0.0;
    dU[iat]    = 0.0;
    d2U[iat]   = 0.0;
    RealType a = ParticleAlpha[iat];
    if (a > 0.0)
    {
#ifdef ENABLE_SOA
      dist = d_table.Distances[iat][icent];
      disp = -1.0*d_table.Displacements[iat][icent];
#else
      int index     = d_table.M[icent] + iat;
      dist = d_table.r(index);
      disp = d_table.dr(index);
#endif
      LogValue -= a * dist * dist;
      U[iat] += a * dist * dist;
      G[iat] -= 2.0 * a * disp;
      L[iat] -= 6.0 * a;
      icent++;
    }
  }
  return LogValue;
}

/** evaluate the ratio \f$exp(U(iat)-U_0(iat))\f$
 * @param P active particle set
 * @param iat particle that has been moved.
 */
ValueType LatticeGaussianProduct::ratio(ParticleSet& P, int iat)
{
  const auto& d_table = P.getDistTable(myTableID);
  int icent           = ParticleCenter[iat];
  if (icent == -1)
    return 1.0;
#ifdef ENABLE_SOA
  RealType newdist = d_table.Temp_r[icent];
#else
  RealType newdist = d_table.Temp[icent].r1;
#endif
  curVal           = ParticleAlpha[iat] * (newdist * newdist);
  return std::exp(U[iat] - curVal);
}


GradType LatticeGaussianProduct::evalGrad(ParticleSet& P, int iat)
{
  const auto& d_table = P.getDistTable(myTableID);
  int icent           = ParticleCenter[iat];
  if (icent == -1)
    return GradType();
  RealType a      = ParticleAlpha[iat];
#ifdef ENABLE_SOA
  PosType newdisp = -1.0*d_table.Temp_dr[icent];
#else
  PosType newdisp = d_table.Temp[icent].dr1;
#endif
  curGrad         = -2.0 * a * newdisp;
  return curGrad;
}


//   GradType
//   LatticeGaussianProduct::evalGradSource
//   (ParticleSet& P, ParticleSet& source, int isrc,
//    TinyVector<ParticleSet::ParticleGradient_t, OHMMS_DIM> &grad_grad,
//    TinyVector<ParticleSet::ParticleLaplacian_t,OHMMS_DIM> &lapl_grad)
//   {
//   }

ValueType LatticeGaussianProduct::ratioGrad(ParticleSet& P, int iat, GradType& grad_iat)
{
  const auto& d_table = P.getDistTable(myTableID);
  int icent           = ParticleCenter[iat];
  if (icent == -1)
    return 1.0;
  RealType a       = ParticleAlpha[iat];
#ifdef ENABLE_SOA
  RealType newdist = d_table.Temp_r[icent];
  PosType newdisp = -1.0*d_table.Temp_dr[icent];
#else
  RealType newdist = d_table.Temp[icent].r1;
  PosType newdisp  = d_table.Temp[icent].dr1;
#endif
  curVal           = a * newdist * newdist;
  curGrad          = -2.0 * a * newdisp;
  grad_iat += curGrad;
  return std::exp(U[iat] - curVal);
}

void LatticeGaussianProduct::restore(int iat) {}

void LatticeGaussianProduct::acceptMove(ParticleSet& P, int iat)
{
  U[iat]   = curVal;
  dU[iat]  = curGrad;
  d2U[iat] = curLap;
}

void LatticeGaussianProduct::evaluateLogAndStore(ParticleSet& P,
                                     ParticleSet::ParticleGradient_t& dG,
                                     ParticleSet::ParticleLaplacian_t& dL)
{
  const auto& d_table = P.getDistTable(myTableID);
  RealType dist       = 0.0;
  PosType  disp       = 0.0;
  int icent           = 0;
  LogValue            = 0.0;
  U                   = 0.0;
  dU                  = 0.0;
  d2U                 = 0.0;
  for (int iat = 0; iat < NumTargetPtcls; iat++)
  {
    RealType a = ParticleAlpha[iat];
    if (a > 0.0)
    {
#ifdef ENABLE_SOA
      dist = d_table.Distances[iat][icent];
      disp = -1.0*d_table.Displacements[iat][icent];
#else
      int index     = d_table.M[icent] + iat;
      dist = d_table.r(index);
      disp = d_table.dr(index);
#endif
      LogValue -= a * dist * dist;
      U[iat] += a * dist * dist;
      dU[iat] -= 2.0 * a * disp;
      d2U[iat] -= 6.0 * a;
      dG[iat] -= 2.0 * a * disp;
      dL[iat] -= 6.0 * a;
      icent++;
    }
  }
}

/** equivalent to evalaute with additional data management */
void LatticeGaussianProduct::registerData(ParticleSet& P, WFBufferType& buf)
{
  evaluateLogAndStore(P, P.G, P.L);
  // Add U, d2U and dU. Keep the order!!!
  buf.add(U.begin(), U.end());
  buf.add(d2U.begin(), d2U.end());
  buf.add(FirstAddressOfdU, LastAddressOfdU);
}

LatticeGaussianProduct::RealType LatticeGaussianProduct::updateBuffer(ParticleSet& P, WFBufferType& buf, bool fromscratch = false)
{
  evaluateLogAndStore(P, P.G, P.L);
  buf.put(U.first_address(), U.last_address());
  buf.put(d2U.first_address(), d2U.last_address());
  buf.put(FirstAddressOfdU, LastAddressOfdU);
  return LogValue;
}

/** copy the current data from a buffer
 *@param P the ParticleSet to operate on
 *@param buf PooledData which stores the data for each walker
 *
 *copyFromBuffer uses the data stored by registerData
 */
void LatticeGaussianProduct::copyFromBuffer(ParticleSet& P, WFBufferType& buf)
{
  buf.get(U.first_address(), U.last_address());
  buf.get(d2U.first_address(), d2U.last_address());
  buf.get(FirstAddressOfdU, LastAddressOfdU);
}

WaveFunctionComponentPtr LatticeGaussianProduct::makeClone(ParticleSet& tqp) const
{
  LatticeGaussianProduct* j1copy     = new LatticeGaussianProduct(CenterRef, tqp);
  j1copy->ParticleAlpha  = ParticleAlpha;
  j1copy->ParticleCenter = ParticleCenter;
  return j1copy;
}

}; // namespace qmcplusplus
