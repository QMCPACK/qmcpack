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
#include "ACForce.h"
#include <sstream>
#include "OhmmsData/AttributeSet.h"

namespace qmcplusplus
{
ACForce::ACForce(ParticleSet& source, ParticleSet& target, TrialWaveFunction& psi_in, QMCHamiltonian& H)
    : ions(source),
      elns(target),
      psi(psi_in),
      ham(H),
      FirstForceIndex(-1),
      Nions(ions.getTotalNum()),
      useSpaceWarp(false),
      swt(target, source)
{
  prefix = "ACForce";
  myName = prefix;

  hf_force.resize(Nions);
  pulay_force.resize(Nions);
  wf_grad.resize(Nions);
  sw_pulay.resize(Nions);
  sw_grad.resize(Nions);
  delta = 1e-4;
};

OperatorBase* ACForce::makeClone(ParticleSet& qp, TrialWaveFunction& psi)
{
  APP_ABORT("ACForce::makeClone(ParticleSet&,TrialWaveFunction&) shouldn't be called");
  return nullptr;
}

OperatorBase* ACForce::makeClone(ParticleSet& qp, TrialWaveFunction& psi_in, QMCHamiltonian& ham_in)
{
  OperatorBase* myclone = new ACForce(ions, qp, psi_in, ham_in);
  return myclone;
}

bool ACForce::put(xmlNodePtr cur)
{
  std::string useSpaceWarpString("no");
  std::string ionionforce("yes");
  RealType swpow(4);
  OhmmsAttributeSet attr;
  attr.add(useSpaceWarpString, "spacewarp"); //"yes" or "no"
  attr.add(swpow, "swpow");                  //Real number"
  attr.add(delta, "delta");                  //Real number"
  attr.put(cur);

  useSpaceWarp = (useSpaceWarpString == "yes") || (useSpaceWarpString == "true");
  swt.setPow(swpow);

  if (useSpaceWarp)
    app_log() << "ACForce is using space warp with power=" << swpow << std::endl;
  else
    app_log() << "ACForce is not using space warp\n";

  return true;
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
  sw_pulay    = 0;
  sw_grad     = 0;
  //This function returns d/dR of the sum of all observables in the physical hamiltonian.
  //Note that the sign will be flipped based on definition of force = -d/dR.
  Value = ham.evaluateIonDerivs(P, ions, psi, hf_force, pulay_force, wf_grad);

  if (useSpaceWarp)
  {
    Force_t el_grad;
    el_grad.resize(P.getTotalNum());
    el_grad = 0;

    ham.evaluateElecGrad(P, psi, el_grad, delta);
    swt.computeSWT(P, ions, el_grad, P.G, sw_pulay, sw_grad);
  }
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

      //TODO: Remove when ACForce is production ready.
      //      if(useSpaceWarp)
      //      {
      //        std::ostringstream swctname1;
      //        std::ostringstream swctname2;
      //        std::ostringstream swctname3;
      //        swctname1 << prefix << "_swct1_" << iat << "_" << x;
      //        swctname2 << prefix << "_swct2_" << iat << "_" << x;
      //        swctname3 << prefix << "_swct3_" << iat << "_" << x;
      //        plist.add(swctname1.str());
      //        plist.add(swctname2.str());
      //        plist.add(swctname3.str());
      //      }
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
      plist[myindex++] = -(pulay_force[iat][iondim] + sw_pulay[iat][iondim]);
      plist[myindex++] = -Value * (wf_grad[iat][iondim] + sw_grad[iat][iondim]);
      plist[myindex++] = -(wf_grad[iat][iondim] + sw_grad[iat][iondim]);

      //TODO: Remove when ACForce is production ready
      //      if(useSpaceWarp)
      //      {
      //        plist[myindex++] = -sw_pulay[iat][iondim];
      //        plist[myindex++] = -Value*sw_grad[iat][iondim];
      //        plist[myindex++] = -sw_grad[iat][iondim];
      //      }
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
      plist[myindex++] = -(pulay_force[iat][iondim] + sw_pulay[iat][iondim]);
      plist[myindex++] = -Value * (wf_grad[iat][iondim] + sw_grad[iat][iondim]);
      plist[myindex++] = -(wf_grad[iat][iondim] + sw_grad[iat][iondim]);
      //TODO: Remove when ACForce is production ready
      //      if(useSpaceWarp)
      //      {
      //        plist[myindex++] = -sw_pulay[iat][iondim];
      //        plist[myindex++] = -Value*sw_grad[iat][iondim];
      //        plist[myindex++] = -sw_grad[iat][iondim];
      //      }
    }
  }
};

} // namespace qmcplusplus
