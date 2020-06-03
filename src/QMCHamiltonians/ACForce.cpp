//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2019 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Raymond Clay, Sandia National Laboratories
//
// File created by: Raymond Clay, rclay@sandia.gov, Sandia National Laboratories
//////////////////////////////////////////////////////////////////////////////////////


/**@file ACForce.cpp
 *@brief Implementation of ACForce, Assaraf-Caffarel ZVZB style force estimation.
 */
#include "QMCHamiltonians/ACForce.h"
#include <sstream>

namespace qmcplusplus
{
ACForce::ACForce(ParticleSet& source, ParticleSet& target, TrialWaveFunction& psi_in, QMCHamiltonian& H)
    : ions(source), elns(target), psi(psi_in), ham(H), FirstForceIndex(-1), Nions(0)
{
  prefix = "ACForce";
  myName = prefix;
  Nions  = ions.getTotalNum();
  hf_force.resize(Nions);
  pulay_force.resize(Nions);
  wf_grad.resize(Nions);
};

OperatorBase* ACForce::makeClone(ParticleSet& qp, TrialWaveFunction& psi)
{
  APP_ABORT("ACForce::makeClone(ParticleSet&,TrialWaveFunction&) shouldn't be called");
  return nullptr;
}

OperatorBase* ACForce::makeClone(ParticleSet& qp, TrialWaveFunction& psi_in, QMCHamiltonian& ham_in)
{
  OperatorBase* myclone = new ACForce(qp, ions, psi_in, ham_in);
  return myclone;
}
void ACForce::add2Hamiltonian(ParticleSet& qp, TrialWaveFunction& psi, QMCHamiltonian& ham_in)
{
  //The following line is modified
  OperatorBase* myclone = makeClone(qp, psi, ham_in);
  if (myclone)
    ham_in.addOperator(myclone, myName, UpdateMode[PHYSICAL]);
}
ACForce::Return_t ACForce::evaluate(ParticleSet& P)
{
  hf_force    = 0;
  pulay_force = 0;
  wf_grad     = 0;

  //This function returns d/dR of the sum of all observables in the physical hamiltonian.
  //Note that the sign will be flipped based on definition of force = -d/dR.
  Value = ham.evaluateIonDerivs(P, ions, psi, hf_force, pulay_force, wf_grad);
  return 0.0;
};

void ACForce::addObservables(PropertySetType& plist, BufferType& collectables)
{
  if (FirstForceIndex < 0)
    FirstForceIndex = plist.size();
  for (int iat = 0; iat < Nions; iat++)
  {
    for (int x = 0; x < OHMMS_DIM; x++)
    {
      std::ostringstream hfname;
      std::ostringstream pulayname;
      std::ostringstream wfgradname1;
      std::ostringstream wfgradname2;
      hfname << prefix << "_hf_" << iat << "_" << x;
      pulayname << prefix << "_pulay_" << iat << "_" << x;
      wfgradname1 << prefix << "_Ewfgrad_" << iat << "_" << x;
      wfgradname2 << prefix << "_wfgrad_" << iat << "_" << x;
      plist.add(hfname.str());
      plist.add(pulayname.str());
      plist.add(wfgradname1.str());
      plist.add(wfgradname2.str());
    }
  }
};
void ACForce::setObservables(PropertySetType& plist)
{
  int myindex = FirstForceIndex;
  for (int iat = 0; iat < Nions; iat++)
  {
    for (int iondim = 0; iondim < OHMMS_DIM; iondim++)
    {
      //Flipping the sign, since these terms currently store d/dR values.
      // add the minus one to be a force.
      plist[myindex++] = -hf_force[iat][iondim];
      plist[myindex++] = -pulay_force[iat][iondim];
      plist[myindex++] = -Value * wf_grad[iat][iondim];
      plist[myindex++] = -wf_grad[iat][iondim];
    }
  }
};
void ACForce::setParticlePropertyList(PropertySetType& plist, int offset)
{
  int myindex = FirstForceIndex + offset;
  for (int iat = 0; iat < Nions; iat++)
  {
    for (int iondim = 0; iondim < OHMMS_DIM; iondim++)
    {
      plist[myindex++] = -hf_force[iat][iondim];
      plist[myindex++] = -pulay_force[iat][iondim];
      plist[myindex++] = -Value * wf_grad[iat][iondim];
      plist[myindex++] = -wf_grad[iat][iondim];
    }
  }
};

} // namespace qmcplusplus
