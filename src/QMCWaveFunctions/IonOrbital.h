//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    
/** @file IonOrbital.h
 * @brief Simple gaussian functions used for orbitals for ions
 */
#ifndef QMCPLUSPLUS_ION_ORBITAL
#define QMCPLUSPLUS_ION_ORBITAL
#include "QMCWaveFunctions/OrbitalBase.h"
#include "Particle/DistanceTableData.h"
#include "Particle/DistanceTable.h"

namespace qmcplusplus
{

/** A composite Orbital
 */
struct IonOrbital : public OrbitalBase
{
private:
  ParticleAttrib<RealType> U,d2U;
  ParticleAttrib<PosType> dU;
  RealType *FirstAddressOfdU, *LastAddressOfdU;
  DistanceTableData* d_table;
  ParticleSet &CenterRef, &PtclRef;
  int NumTargetPtcls, NumCenters;
  RealType curVal, curLap;
  PosType curGrad;
public:
  std::vector<RealType> ParticleAlpha;
  std::vector<int> ParticleCenter;

  IonOrbital(ParticleSet &centers, ParticleSet &ptcls);

  ~IonOrbital();

  /** check out optimizable variables
   */
  void checkOutVariables(const opt_variables_type& o);

  /** check in an optimizable parameter
   * @param o a super set of optimizable variables
   */
  void checkInVariables(opt_variables_type& o);

  /** print the state, e.g., optimizables */
  void reportStatus(std::ostream& os);

  /** reset the parameters during optimizations
   */
  void resetParameters(const opt_variables_type& active);

  void resetTargetParticleSet(ParticleSet& P);

  ValueType
  evaluate(ParticleSet& P,
           ParticleSet::ParticleGradient_t& G,
           ParticleSet::ParticleLaplacian_t& L);

  RealType evaluateLog(ParticleSet& P,
                       ParticleSet::ParticleGradient_t& G, ParticleSet::ParticleLaplacian_t& L);

  ValueType ratio(ParticleSet& P, int iat,
                  ParticleSet::ParticleGradient_t& dG,
                  ParticleSet::ParticleLaplacian_t& dL);

  ValueType ratio(ParticleSet& P, int iat);

  void acceptMove(ParticleSet& P, int iat);

  void restore(int iat);

  void update(ParticleSet& P,
              ParticleSet::ParticleGradient_t& dG,
              ParticleSet::ParticleLaplacian_t& dL,
              int iat);

  RealType
  registerData(ParticleSet& P, BufferType& buf);

  RealType
  updateBuffer(ParticleSet& P, BufferType& buf, bool fromscratch);

  void
  copyFromBuffer(ParticleSet& P, BufferType& buf);

  RealType
  evaluateLog(ParticleSet& P,BufferType& buf);

  GradType evalGrad(ParticleSet& P, int iat);

  ValueType ratioGrad(ParticleSet& P, int iat, GradType& grad_iat);


  OrbitalBase* makeClone(ParticleSet& tqp) const;

  ValueType
  logRatio(ParticleSet& P, int iat,
           ParticleSet::ParticleGradient_t& dG,
           ParticleSet::ParticleLaplacian_t& dL);

  void evaluateLogAndStore(ParticleSet& P,
                           ParticleSet::ParticleGradient_t& dG,
                           ParticleSet::ParticleLaplacian_t& dL);


};
}
#endif
