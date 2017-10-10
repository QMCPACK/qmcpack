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
    
    
#include "QMCWaveFunctions/IonOrbital.h"

namespace qmcplusplus
{

typedef IonOrbital::ValueType ValueType;
typedef IonOrbital::GradType  GradType;

IonOrbital::IonOrbital (ParticleSet &centers, ParticleSet &ptcls) :
  CenterRef(centers), PtclRef(ptcls)
{
  Optimizable=false;
  OrbitalName="IonOrbital";
  NumTargetPtcls = ptcls.getTotalNum();
  NumCenters     = centers.getTotalNum();
  d_table = DistanceTable::add(CenterRef, PtclRef, DT_AOS);
  U.resize(NumTargetPtcls);
  dU.resize(NumTargetPtcls);
  d2U.resize(NumTargetPtcls);
  FirstAddressOfdU = &(dU[0][0]);
  LastAddressOfdU = FirstAddressOfdU + dU.size()*OHMMS_DIM;
}

IonOrbital::~IonOrbital() {  }

//evaluate the distance table with P
void
IonOrbital::resetTargetParticleSet(ParticleSet& P)
{
  d_table = DistanceTable::add(CenterRef,P, DT_AOS);
  //if (dPsi) dPsi->resetTargetParticleSet(P);
}

void IonOrbital::checkInVariables(opt_variables_type& active)        { }
void IonOrbital::checkOutVariables(const opt_variables_type& active) { }
void IonOrbital::resetParameters(const opt_variables_type& active)   { }
void IonOrbital::reportStatus(std::ostream& os)                           { }

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
IonOrbital::RealType
IonOrbital::evaluateLog(ParticleSet& P,
                        ParticleSet::ParticleGradient_t& G,
                        ParticleSet::ParticleLaplacian_t& L)
{
  int icent = 0;
  LogValue = 0.0;
  //P.update();
  //CenterRef.update();
  //d_table->evaluate(PtclRef);
  for (int iat=0; iat<NumTargetPtcls; iat++)
  {
    RealType a = ParticleAlpha[iat];
    if (a > 0.0)
    {
      int index = d_table->M[icent] + iat;
      RealType dist = d_table->r(index);
      PosType  disp  = d_table->dr(index);
      LogValue -= a*dist*dist;
      U[iat]   += a*dist*dist;
      G[iat]   -= 2.0*a*disp;
      L[iat]   -= 6.0*a;
      icent++;
    }
  }
  return LogValue;
}

ValueType
IonOrbital::evaluate(ParticleSet& P,
                     ParticleSet::ParticleGradient_t& G,
                     ParticleSet::ParticleLaplacian_t& L)
{
  return std::exp(evaluateLog(P,G,L));
}

/** evaluate the ratio \f$exp(U(iat)-U_0(iat))\f$
 * @param P active particle set
 * @param iat particle that has been moved.
 */
ValueType
IonOrbital::ratio(ParticleSet& P, int iat)
{
  int icent = ParticleCenter[iat];
  if (icent == -1)
    return 1.0;
  int index = d_table->M[icent] + iat;
  RealType newdist = d_table->Temp[icent].r1;
  curVal = ParticleAlpha[iat]*(newdist*newdist);
  return std::exp(U[iat] - curVal);
}



/** evaluate the ratio \f$exp(U(iat)-U_0(iat))\f$
    and fill-in the differential gradients/laplacians
  * @param P active particle set
  * @param iat particle that has been moved.
  * @param dG partial gradients
  * @param dL partial laplacians
  */
ValueType
IonOrbital::ratio(ParticleSet& P, int iat,
                  ParticleSet::ParticleGradient_t& dG,
                  ParticleSet::ParticleLaplacian_t& dL)
{
  return std::exp(logRatio (P, iat, dG, dL));
}



GradType
IonOrbital::evalGrad(ParticleSet& P, int iat)
{
  int icent = ParticleCenter[iat];
  if (icent == -1)
    return GradType();
  RealType a = ParticleAlpha[iat];
  PosType  newdisp = d_table->Temp[icent].dr1;
  curGrad = -2.0*a*newdisp;
  return curGrad;
}


//   GradType
//   IonOrbital::evalGradSource
//   (ParticleSet& P, ParticleSet& source, int isrc,
//    TinyVector<ParticleSet::ParticleGradient_t, OHMMS_DIM> &grad_grad,
//    TinyVector<ParticleSet::ParticleLaplacian_t,OHMMS_DIM> &lapl_grad)
//   {
//   }

ValueType
IonOrbital::ratioGrad(ParticleSet& P, int iat, GradType& grad_iat)
{
  int icent = ParticleCenter[iat];
  if (icent == -1)
    return 1.0;
  RealType a = ParticleAlpha[iat];
  RealType newdist = d_table->Temp[icent].r1;
  PosType  newdisp = d_table->Temp[icent].dr1;
  curVal = a*newdist*newdist;
  curGrad = -2.0*a*newdisp;
  grad_iat += curGrad;
  return std::exp(U[iat]-curVal);
}

ValueType
IonOrbital::logRatio(ParticleSet& P, int iat,
                     ParticleSet::ParticleGradient_t& dG,
                     ParticleSet::ParticleLaplacian_t& dL)
{
  int icent = ParticleCenter[iat];
  if (icent == -1)
    return 0.0;
  int index = d_table->M[icent] + iat;
  RealType a = ParticleAlpha[iat];
  RealType newdist = d_table->Temp[icent].r1;
  PosType  newdisp = d_table->Temp[icent].dr1;
  curVal = a*newdist*newdist;
  curGrad = -2.0*a*newdisp;
  curLap  = -6.0*a;
  dG[iat] += curGrad - dU[iat];
  dL[iat] += curLap  - d2U[iat];
  return U[iat] - curVal;
  // fprintf (stderr, "new dist = %8.5f  old dist=%8.5f\n",
  // 	     newdist, olddist);
  // fprintf (stderr, "iat = %d   icent = %d\n", iat, icent);
  // dG[iat] -= 2.0*a*(newdisp - olddisp);
  // dL[iat] -= 0.0*a;
  // RealType lograt = -a*(newdist*newdist - olddist*olddist);
  // std::cerr << "lograt = " << lograt << std::endl;
  // return lograt;
}

void
IonOrbital::restore(int iat) {}

void
IonOrbital::acceptMove(ParticleSet& P, int iat)
{
  U[iat] = curVal;
  dU[iat]=curGrad;
  d2U[iat]=curLap;
}


void
IonOrbital::update(ParticleSet& P,
                   ParticleSet::ParticleGradient_t& dG,
                   ParticleSet::ParticleLaplacian_t& dL,
                   int iat)
{
  dG[iat] += curGrad-dU[iat];
  dU[iat]=curGrad;
  dL[iat] += curLap-d2U[iat];
  d2U[iat]=curLap;
  U[iat] = curVal;
}

void
IonOrbital::evaluateLogAndStore(ParticleSet& P,
                                ParticleSet::ParticleGradient_t& dG,
                                ParticleSet::ParticleLaplacian_t& dL)
{
  int icent = 0;
  LogValue = 0.0;
  U=0.0;
  dU=0.0;
  d2U=0.0;
  for (int iat=0; iat<NumTargetPtcls; iat++)
  {
    RealType a = ParticleAlpha[iat];
    if (a > 0.0)
    {
      int index = d_table->M[icent] + iat;
      RealType dist = d_table->r(index);
      PosType  disp  = d_table->dr(index);
      LogValue -= a*dist*dist;
      U[iat]   += a*dist*dist;
      dU[iat]  -= 2.0*a*disp;
      d2U[iat] -= 6.0*a;
      dG[iat]  -= 2.0*a*disp;
      dL[iat]  -= 6.0*a;
      icent++;
    }
  }
}

/** equivalent to evalaute with additional data management */
IonOrbital::RealType
IonOrbital::registerData(ParticleSet& P, PooledData<RealType>& buf)
{
  evaluateLogAndStore(P, P.G, P.L);
  // Add U, d2U and dU. Keep the order!!!
  buf.add(U.begin(), U.end());
  buf.add(d2U.begin(), d2U.end());
  buf.add(FirstAddressOfdU,LastAddressOfdU);
  return LogValue;
}

IonOrbital::RealType
IonOrbital::updateBuffer(ParticleSet& P, BufferType& buf,
                         bool fromscratch=false)
{
  evaluateLogAndStore(P,P.G,P.L);
  buf.put(U.first_address(), U.last_address());
  buf.put(d2U.first_address(), d2U.last_address());
  buf.put(FirstAddressOfdU,LastAddressOfdU);
  return LogValue;
}

/** copy the current data from a buffer
 *@param P the ParticleSet to operate on
 *@param buf PooledData which stores the data for each walker
 *
 *copyFromBuffer uses the data stored by registerData or evaluate(P,buf)
 */
void
IonOrbital::copyFromBuffer(ParticleSet& P, PooledData<RealType>& buf)
{
  buf.get(U.first_address(), U.last_address());
  buf.get(d2U.first_address(), d2U.last_address());
  buf.get(FirstAddressOfdU,LastAddressOfdU);
}

/** return the current value and copy the current data to a buffer
 *@param P the ParticleSet to operate on
 *@param buf PooledData which stores the data for each walker
 */
IonOrbital::RealType
IonOrbital::evaluateLog(ParticleSet& P, PooledData<RealType>& buf)
{
  RealType sumu = 0.0;
  for (int i=0; i<U.size(); i++)
    sumu+=U[i];
  buf.put(U.first_address(), U.last_address());
  buf.put(d2U.first_address(), d2U.last_address());
  buf.put(FirstAddressOfdU,LastAddressOfdU);
  return -sumu;
  // int icent = 0;
  // LogValue = 0.0;
  // U = 0.0;
  // for (int iat=0; iat<NumTargetPtcls; iat++) {
  //   RealType a = ParticleAlpha[iat];
  //   if (a > 0.0) {
  // 	int index = d_table->M[icent] + iat;
  // 	RealType dist = d_table->r(index);
  // 	PosType  disp  = d_table->dr(index);
  // 	LogValue -= a*dist*dist;
  // 	icent++;
  //   }
  // }
  // return LogValue;
  // return evaluateLog(P, P.G, P.L);
  // std::cerr << "IonOrbital::evaluateLog(ParticleSet& P, PooledData<RealType>& buf) called.\n";
  // buf.put(U.first_address(), U.last_address());
  // buf.put(d2U.first_address(), d2U.last_address());
  // buf.put(FirstAddressOfdU,LastAddressOfdU);
  // return -sumu;
}

OrbitalBasePtr
IonOrbital::makeClone(ParticleSet& tqp) const
{
  IonOrbital* j1copy=new IonOrbital(CenterRef,tqp);
  j1copy->ParticleAlpha = ParticleAlpha;
  j1copy->ParticleCenter = ParticleCenter;
  return j1copy;
}

};
