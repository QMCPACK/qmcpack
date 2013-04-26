//////////////////////////////////////////////////////////////////
// (c) Copyright 2003-  by Jeongnim Kim
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//   National Center for Supercomputing Applications &
//   Materials Computation Center
//   University of Illinois, Urbana-Champaign
//   Urbana, IL 61801
//   e-mail: jnkim@ncsa.uiuc.edu
//
// Supported by
//   National Center for Supercomputing Applications, UIUC
//   Materials Computation Center, UIUC
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
#ifndef QMCPLUSPLUS_MULTISLATERDETERMINANTWITHBACKFLOW_ORBITAL_H
#define QMCPLUSPLUS_MULTISLATERDETERMINANTWITHBACKFLOW_ORBITAL_H
#include <Configuration.h>
#include <QMCWaveFunctions/FermionBase.h>
#include <QMCWaveFunctions/Fermion/DiracDeterminantBase.h>
#include "QMCWaveFunctions/Fermion/DiracDeterminantWithBackflow.h"
#include "QMCWaveFunctions/Fermion/BackflowTransformation.h"
#include <QMCWaveFunctions/Fermion/SPOSetProxyForMSD.h>
#include "QMCWaveFunctions/Fermion/MultiSlaterDeterminant.h"
#include "Utilities/NewTimer.h"

namespace qmcplusplus
{

/** @ingroup OrbitalComponent
 *  @brief MultiSlaterDeterminantWithBackflow
 */
class MultiSlaterDeterminantWithBackflow: public MultiSlaterDeterminant
{

public:

  ///constructor
  MultiSlaterDeterminantWithBackflow(ParticleSet& targetPtcl, SPOSetProxyPtr upspo, SPOSetProxyPtr dnspo, BackflowTransformation *tr);

  ///destructor
  ~MultiSlaterDeterminantWithBackflow();

  void checkInVariables(opt_variables_type& active);
  void checkOutVariables(const opt_variables_type& active);
  void resetParameters(const opt_variables_type& active);
  void reportStatus(ostream& os);

  ///set BF pointers
  void setBF(BackflowTransformation* bf)
  {
    BFTrans=bf;
    for(int i=0; i<dets_up.size(); i++)
      dets_up[i]->setBF(bf);
    for(int i=0; i<dets_dn.size(); i++)
      dets_dn[i]->setBF(bf);
  }

  ValueType
  evaluate(ParticleSet& P
           ,ParticleSet::ParticleGradient_t& G
           ,ParticleSet::ParticleLaplacian_t& L);

  RealType
  evaluateLog(ParticleSet& P //const DistanceTableData* dtable,
              , ParticleSet::ParticleGradient_t& G
              , ParticleSet::ParticleLaplacian_t& L);

  GradType evalGrad(ParticleSet& P, int iat);
  ValueType ratioGrad(ParticleSet& P, int iat, GradType& grad_iat);
  ValueType ratio(ParticleSet& P, int iat
                  , ParticleSet::ParticleGradient_t& dG,ParticleSet::ParticleLaplacian_t& dL);

  ValueType ratio(ParticleSet& P, int iat);
  void acceptMove(ParticleSet& P, int iat);
  void restore(int iat);

  void update(ParticleSet& P
              , ParticleSet::ParticleGradient_t& dG, ParticleSet::ParticleLaplacian_t& dL
              , int iat);
  RealType evaluateLog(ParticleSet& P,BufferType& buf);
  RealType registerData(ParticleSet& P, BufferType& buf);
  RealType updateBuffer(ParticleSet& P, BufferType& buf, bool fromscratch=false);
  void copyFromBuffer(ParticleSet& P, BufferType& buf);

  OrbitalBasePtr makeClone(ParticleSet& tqp) const;
  void evaluateDerivatives(ParticleSet& P,
                           const opt_variables_type& optvars,
                           vector<RealType>& dlogpsi,
                           vector<RealType>& dhpsioverpsi);

  void resize(int,int);

  // transformation
  BackflowTransformation *BFTrans;

  // temporary storage for evaluateDerivatives
  Matrix<RealType> dpsia_up, dLa_up;
  Matrix<RealType> dpsia_dn, dLa_dn;
  Array<GradType,3> dGa_up, dGa_dn;
};

}
#endif
/***************************************************************************
 * $RCSfile$   $Author: miguel.mmorales $
 * $Revision: 4791 $   $Date: 2010-05-12 12:08:35 -0500 (Wed, 12 May 2010) $
 * $Id: MultiSlaterDeterminantWithBackflow.h 4791 2010-05-12 17:08:35Z miguel.mmorales $
 ay**************************************************************************/
