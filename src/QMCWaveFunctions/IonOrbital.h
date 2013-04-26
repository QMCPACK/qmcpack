//////////////////////////////////////////////////////////////////
// (c) Copyright 2009-  by Ken Esler
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//   National Center for Supercomputing Applications &
//   Materials Computation Center
//   University of Illinois, Urbana-Champaign
//   Urbana, IL 61801
//   e-mail: esler@uiuc.edu
//   Tel:    217-244-6319 (NCSA) 217-333-3324 (MCC)
//
// Supported by
//   National Center for Supercomputing Applications, UIUC
//   Materials Computation Center, UIUC
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
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
  vector<RealType> ParticleAlpha;
  vector<int> ParticleCenter;

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
  void reportStatus(ostream& os);

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
/***************************************************************************
 * $RCSfile$   $Author: esler $
 * $Revision: 3848 $   $Date: 2009-05-20 12:38:03 -0500 (Wed, 20 May 2009) $
 * $Id: IonOrbital.h 3848 2009-05-20 17:38:03Z jnkim $
 ***************************************************************************/
