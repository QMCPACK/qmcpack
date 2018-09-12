// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Bryan Clark, bclark@Princeton.edu, Princeton University
//                    Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Miguel Morales, moralessilva2@llnl.gov, Lawrence Livermore National Laboratory
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Raymond Clay III, j.k.rofling@gmail.com, Lawrence Livermore National Laboratory
//                    Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//                    Jeongnim Kim, jeongnim.kim@intel.com, Intel Corp.
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    

#include "QMCWaveFunctions/Fermion/DiracDeterminantBase.h"
#include "Numerics/DeterminantOperators.h"
#include "Numerics/OhmmsBlas.h"
#include "Numerics/MatrixOperators.h"
#include "simd/simd.hpp"

namespace qmcplusplus
{

/** constructor
 *@param spos the single-particle orbital set
 *@param first index of the first particle
 */
DiracDeterminantBase::DiracDeterminantBase(int first):
  NP(0), FirstIndex(first)
  ,UpdateTimer("DiracDeterminantBase::update",timer_level_fine)
  ,RatioTimer("DiracDeterminantBase::ratio",timer_level_fine)
  ,InverseTimer("DiracDeterminantBase::inverse",timer_level_fine)
  ,BufferTimer("DiracDeterminantBase::buffer",timer_level_fine)
  ,SPOVTimer("DiracDeterminantBase::spoval",timer_level_fine)
  ,SPOVGLTimer("DiracDeterminantBase::spovgl",timer_level_fine)
{
  OrbitalName="DiracDeterminantBase";
  registerTimers();
}

///default destructor
DiracDeterminantBase::~DiracDeterminantBase() {}

#if 0
DiracDeterminantBase& DiracDeterminantBase::operator=(const DiracDeterminantBase& s)
{
  Bytes_in_WFBuffer=s.Bytes_in_WFBuffer;
  NP=0;
  resize(s.NumPtcls, s.NumOrbitals);
  return *this;
}
#endif

/** set the index of the first particle in the determinant and reset the size of the determinant
 *@param first index of first particle
 *@param nel number of particles in the determinant
 */
void DiracDeterminantBase::set(int first, int nel)
{
  FirstIndex = first;
  resize(nel,nel);
}

void DiracDeterminantBase::invertPsiM(const ValueMatrix_t& logdetT, ValueMatrix_t& invMat)
{
  InverseTimer.start();
#ifdef MIXED_PRECISION
  simd::transpose(logdetT.data(), NumOrbitals, logdetT.cols(), 
      psiM_hp.data(), NumOrbitals, psiM_hp.cols());
  ParticleSet::Scalar_t PhaseValue_hp;
  detEng_hp.invert(psiM_hp,true);
  LogValue = static_cast<RealType>(detEng_hp.LogDet);
  PhaseValue = static_cast<RealType>(detEng_hp.Phase);
  invMat = psiM_hp;
#else
  simd::transpose(logdetT.data(), NumOrbitals, logdetT.cols(), 
      invMat.data(), NumOrbitals, invMat.cols());
  detEng.invert(invMat,true);
  LogValue = detEng.LogDet;
  PhaseValue = detEng.Phase;
#endif
  InverseTimer.stop();
}



///reset the size: with the number of particles and number of orbtials
void DiracDeterminantBase::resize(int nel, int morb)
{
  int norb=morb;
  if(norb <= 0)
    norb = nel; // for morb == -1 (default)
  psiM.resize(nel,norb);
  dpsiM.resize(nel,norb);
  d2psiM.resize(nel,norb);
  psiV.resize(norb);
  memoryPool.resize(nel*norb);
  psiM_temp.attachReference(memoryPool.data(),nel,norb);
#ifdef MIXED_PRECISION
  psiM_hp.resize(nel,norb);
#endif
  LastIndex = FirstIndex + nel;
  NumPtcls=nel;
  NumOrbitals=norb;

  dpsiV.resize(NumOrbitals);
  d2psiV.resize(NumOrbitals);
  FirstAddressOfdV = &(dpsiM(0,0)[0]); //(*dpsiM.begin())[0]);
  LastAddressOfdV = FirstAddressOfdV + NumPtcls*NumOrbitals*DIM;
  // For forces
  /* Ye Luo, Apr 18th 2015
   * To save the memory used by every walker, the resizing the following giant matrices are commented.
   * When ZVZB forces and stresses are ready for deployment, R. Clay will take care of those matrices.
   */
  /*
  grad_source_psiM.resize(nel,norb);
  grad_lapl_source_psiM.resize(nel,norb);
  grad_grad_source_psiM.resize(nel,norb);
  phi_alpha_Minv.resize(nel,norb);
  grad_phi_Minv.resize(nel,norb);
  lapl_phi_Minv.resize(nel,norb);
  grad_phi_alpha_Minv.resize(nel,norb);
  */
}

DiracDeterminantBase::GradType
DiracDeterminantBase::evalGrad(ParticleSet& P, int iat)
{
  WorkingIndex = iat-FirstIndex;
  RatioTimer.start();
  DiracDeterminantBase::GradType g = simd::dot(psiM[WorkingIndex],dpsiM[WorkingIndex],NumOrbitals);
  RatioTimer.stop();
  return g;
}


/** move was accepted, update the real container
*/
void DiracDeterminantBase::acceptMove(ParticleSet& P, int iat)
{
  PhaseValue += evaluatePhase(curRatio);
  LogValue +=std::log(std::abs(curRatio));
  UpdateTimer.start();
  detEng.updateRow(psiM,psiV.data(),WorkingIndex,curRatio);
  if(UpdateMode == ORB_PBYP_PARTIAL)
  {
    simd::copy(dpsiM[WorkingIndex],  dpsiV.data(),  NumOrbitals);
    simd::copy(d2psiM[WorkingIndex], d2psiV.data(), NumOrbitals);
  }
  UpdateTimer.stop();
  curRatio=1.0;
}

/** move was rejected. copy the real container to the temporary to move on
*/
void DiracDeterminantBase::restore(int iat)
{
  curRatio=1.0;
}

void
DiracDeterminantBase::registerData(ParticleSet& P, WFBufferType& buf)
{
  // Ye: no idea about NP.
  if(NP == 0) //first time, allocate once
  {
    NP=P.getTotalNum();
  }

  if ( Bytes_in_WFBuffer == 0 )
  {
    //add the data: inverse, gradient and laplacian
    Bytes_in_WFBuffer = buf.current();
    buf.add(psiM.first_address(),psiM.last_address());
    buf.add(FirstAddressOfdV,LastAddressOfdV);
    buf.add(d2psiM.first_address(),d2psiM.last_address());
    Bytes_in_WFBuffer = buf.current()-Bytes_in_WFBuffer;
    // free local space
    psiM.free();
    dpsiM.free();
    d2psiM.free();
  }
  else
  {
    buf.forward(Bytes_in_WFBuffer);
  }
  buf.add(LogValue);
  buf.add(PhaseValue);
}

DiracDeterminantBase::RealType DiracDeterminantBase::updateBuffer(ParticleSet& P,
    WFBufferType& buf, bool fromscratch)
{
  if(fromscratch)
  {
    LogValue=evaluateLog(P,P.G,P.L);
  }
  else
  {
    updateAfterSweep(P,P.G,P.L);
  }
  BufferTimer.start();
  buf.forward(Bytes_in_WFBuffer);
  buf.put(LogValue);
  buf.put(PhaseValue);
  BufferTimer.stop();
  return LogValue;
}

void DiracDeterminantBase::copyFromBuffer(ParticleSet& P, WFBufferType& buf)
{
  BufferTimer.start();
  psiM.attachReference(buf.lendReference<ValueType>(psiM.size()));
  dpsiM.attachReference(buf.lendReference<GradType>(dpsiM.size()));
  d2psiM.attachReference(buf.lendReference<ValueType>(d2psiM.size()));
  buf.get(LogValue);
  buf.get(PhaseValue);
  BufferTimer.stop();
}



/** Calculate the log value of the Dirac determinant for particles
 *@param P input configuration containing N particles
 *@param G a vector containing N gradients
 *@param L a vector containing N laplacians
 *@return the value of the determinant
 *
 *\f$ (first,first+nel). \f$  Add the gradient and laplacian
 *contribution of the determinant to G(radient) and L(aplacian)
 *for local energy calculations.
 */
DiracDeterminantBase::RealType
DiracDeterminantBase::evaluateLog(ParticleSet& P,
                                  ParticleSet::ParticleGradient_t& G,
                                  ParticleSet::ParticleLaplacian_t& L)
{
  recompute(P);

  if(NumPtcls==1)
  {
    ValueType y=psiM(0,0);
    GradType rv = y*dpsiM(0,0);
    G[FirstIndex] += rv;
    L[FirstIndex] += y*d2psiM(0,0) - dot(rv,rv);
  }
  else
  {
    for(int i=0, iat=FirstIndex; i<NumPtcls; i++, iat++)
    {
      mGradType rv=simd::dot(psiM[i],dpsiM[i],NumOrbitals);
      mValueType lap=simd::dot(psiM[i],d2psiM[i],NumOrbitals);
      G[iat] += rv;
      L[iat] += lap - dot(rv,rv);
    }
  }
  return LogValue;
}

void
DiracDeterminantBase::evaluateDerivatives(ParticleSet& P,
    const opt_variables_type& active,
    std::vector<RealType>& dlogpsi,
    std::vector<RealType>& dhpsioverpsi)
{
}

WaveFunctionComponentPtr DiracDeterminantBase::makeClone(ParticleSet& tqp) const
{
  APP_ABORT(" Illegal action. Cannot use DiracDeterminantBase::makeClone");
  return 0;
}

DiracDeterminantBase::DiracDeterminantBase(const DiracDeterminantBase& s)
  : WaveFunctionComponent(s), NP(0), FirstIndex(s.FirstIndex)
  ,UpdateTimer(s.UpdateTimer)
  ,RatioTimer(s.RatioTimer)
  ,InverseTimer(s.InverseTimer)
  ,BufferTimer(s.BufferTimer)
  ,SPOVTimer(s.SPOVTimer)
  ,SPOVGLTimer(s.SPOVGLTimer)
{
  registerTimers();
  this->resize(s.NumPtcls,s.NumOrbitals);
}

//SPOSetPtr  DiracDeterminantBase::clonePhi() const
//{
//  return Phi->makeClone();
//}

void DiracDeterminantBase::registerTimers()
{
  UpdateTimer.reset();
  RatioTimer.reset();
  TimerManager.addTimer (&UpdateTimer);
  TimerManager.addTimer (&RatioTimer);
  TimerManager.addTimer (&InverseTimer);
  TimerManager.addTimer (&BufferTimer);
  TimerManager.addTimer (&SPOVTimer);
  TimerManager.addTimer (&SPOVGLTimer);
}

}
