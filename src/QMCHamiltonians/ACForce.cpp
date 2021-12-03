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
#include "OhmmsData/AttributeSet.h"

namespace qmcplusplus
{
ACForce::ACForce(ParticleSet& source, ParticleSet& target, TrialWaveFunction& psi_in, QMCHamiltonian& H)
    : delta_(1e-4),
      ions_(source),
      elns_(target),
      psi_(psi_in),
      ham_(H),
      first_force_index_(-1),
      useSpaceWarp_(false),
      swt_(target, source)
{
  setName("ACForce");

  const std::size_t nIons = ions_.getTotalNum();
  hf_force_.resize(nIons);
  pulay_force_.resize(nIons);
  wf_grad_.resize(nIons);
  sw_pulay_.resize(nIons);
  sw_grad_.resize(nIons);
};

std::unique_ptr<OperatorBase> ACForce::makeClone(ParticleSet& qp, TrialWaveFunction& psi)
{
  APP_ABORT("ACForce::makeClone(ParticleSet&,TrialWaveFunction&) shouldn't be called");
  return nullptr;
}

std::unique_ptr<OperatorBase> ACForce::makeClone(ParticleSet& qp, TrialWaveFunction& psi_in, QMCHamiltonian& ham_in)
{
  std::unique_ptr<ACForce> myclone = std::make_unique<ACForce>(ions_, qp, psi_in, ham_in);
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
  attr.add(delta_, "delta");                 //Real number"
  attr.put(cur);

  useSpaceWarp_ = (useSpaceWarpString == "yes") || (useSpaceWarpString == "true");
  swt_.setPow(swpow);

  if (useSpaceWarp_)
    app_log() << "ACForce is using space warp with power=" << swpow << std::endl;
  else
    app_log() << "ACForce is not using space warp\n";

  return true;
}

bool ACForce::get(std::ostream& os) const { return true; }

void ACForce::add2Hamiltonian(ParticleSet& qp, TrialWaveFunction& psi, QMCHamiltonian& ham_in)
{
  //The following line is modified
  std::unique_ptr<OperatorBase> myclone = makeClone(qp, psi, ham_in);
  if (myclone)
  {
    ham_in.addOperator(std::move(myclone), name_, update_mode_[PHYSICAL]);
  }
}
ACForce::Return_t ACForce::evaluate(ParticleSet& P)
{
  hf_force_    = 0;
  pulay_force_ = 0;
  wf_grad_     = 0;
  sw_pulay_    = 0;
  sw_grad_     = 0;
  //This function returns d/dR of the sum of all observables in the physical hamiltonian.
  //Note that the sign will be flipped based on definition of force = -d/dR.
  value_ = ham_.evaluateIonDerivs(P, ions_, psi_, hf_force_, pulay_force_, wf_grad_);

  if (useSpaceWarp_)
  {
    Force_t el_grad;
    el_grad.resize(P.getTotalNum());
    el_grad = 0;

    ham_.evaluateElecGrad(P, psi_, el_grad, delta_);
    swt_.computeSWT(P, ions_, el_grad, P.G, sw_pulay_, sw_grad_);
  }
  return 0.0;
};

void ACForce::resetTargetParticleSet(ParticleSet& P) {}

void ACForce::addObservables(PropertySetType& plist, BufferType& collectables)
{
  if (first_force_index_ < 0)
    first_force_index_ = plist.size();
  for (int iat = 0; iat < ions_.getTotalNum(); iat++)
  {
    const std::string iatStr(std::to_string(iat));

    for (int x = 0; x < OHMMS_DIM; x++)
    {
      const std::string xStr(std::to_string(x));

      const std::string hfname("ACForce_hf_" + iatStr + "_" + xStr);
      const std::string pulayname("ACForce_pulay_" + iatStr + "_" + xStr);
      const std::string wfgradname1("ACForce_Ewfgrad_" + iatStr + "_" + xStr);
      const std::string wfgradname2("ACForce_wfgrad_" + iatStr + "_" + xStr);

      plist.add(hfname);
      plist.add(pulayname);
      plist.add(wfgradname1);
      plist.add(wfgradname2);
    }
  }
};
void ACForce::setObservables(PropertySetType& plist)
{
  // TODO : bounds check for plist

  int myindex = first_force_index_;
  for (int iat = 0; iat < ions_.getTotalNum(); iat++)
  {
    for (int iondim = 0; iondim < OHMMS_DIM; iondim++)
    {
      //Flipping the sign, since these terms currently store d/dR values.
      // add the minus one to be a force.
      plist[myindex++] = -hf_force_[iat][iondim];
      plist[myindex++] = -(pulay_force_[iat][iondim] + sw_pulay_[iat][iondim]);
      plist[myindex++] = -value_ * (wf_grad_[iat][iondim] + sw_grad_[iat][iondim]);
      plist[myindex++] = -(wf_grad_[iat][iondim] + sw_grad_[iat][iondim]);
    }
  }
};
void ACForce::setParticlePropertyList(PropertySetType& plist, int offset)
{
  int myindex = first_force_index_ + offset;
  for (int iat = 0; iat < ions_.getTotalNum(); iat++)
  {
    for (int iondim = 0; iondim < OHMMS_DIM; iondim++)
    {
      plist[myindex++] = -hf_force_[iat][iondim];
      plist[myindex++] = -(pulay_force_[iat][iondim] + sw_pulay_[iat][iondim]);
      plist[myindex++] = -value_ * (wf_grad_[iat][iondim] + sw_grad_[iat][iondim]);
      plist[myindex++] = -(wf_grad_[iat][iondim] + sw_grad_[iat][iondim]);
    }
  }
};

} // namespace qmcplusplus
