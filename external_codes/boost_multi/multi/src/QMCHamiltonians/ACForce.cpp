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
#include "ParticleBase/ParticleAttribOps.h"

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
      fastDerivatives_(false),
      swt_(target, source),
      reg_epsilon_(0.0),
      f_epsilon_(1.0)
{
  setName("ACForce");

  const std::size_t nIons = ions_.getTotalNum();
  hf_force_.resize(nIons);
  pulay_force_.resize(nIons);
  wf_grad_.resize(nIons);
  sw_pulay_.resize(nIons);
  sw_grad_.resize(nIons);
  psi_in.initializeTWFFastDerivWrapper(elns_, psi_wrapper_);
};

std::unique_ptr<OperatorBase> ACForce::makeClone(ParticleSet& qp, TrialWaveFunction& psi)
{
  APP_ABORT("ACForce::makeClone(ParticleSet&,TrialWaveFunction&) shouldn't be called");
  return nullptr;
}

std::unique_ptr<OperatorBase> ACForce::makeClone(ParticleSet& qp, TrialWaveFunction& psi_in, QMCHamiltonian& ham_in)
{
  std::unique_ptr<ACForce> myclone = std::make_unique<ACForce>(ions_, qp, psi_in, ham_in);
  myclone->fastDerivatives_        = fastDerivatives_;
  myclone->useSpaceWarp_           = useSpaceWarp_;
  myclone->first_force_index_      = first_force_index_;
  myclone->reg_epsilon_            = reg_epsilon_;
  myclone->delta_                  = delta_;
  return myclone;
}

bool ACForce::put(xmlNodePtr cur)
{
  std::string ionionforce("yes");
  RealType swpow(4);
  OhmmsAttributeSet attr;
  attr.add(useSpaceWarp_, "spacewarp", {false}); //"yes" or "no"
  attr.add(swpow, "swpow");                      //Real number"
  attr.add(delta_, "delta");                     //Real number"
  attr.add(reg_epsilon_, "epsilon");
  attr.add(fastDerivatives_, "fast_derivatives", {false});
  attr.put(cur);
  if (reg_epsilon_ < 0)
    throw std::runtime_error("ACForce::put(): epsilon<0 not allowed.");
  if (fastDerivatives_)
    app_log() << "ACForce is using the fast force algorithm\n";
  else
    app_log() << "ACForce is using the default algorithm\n";
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
  if (fastDerivatives_)
    value_ = ham_.evaluateIonDerivsDeterministicFast(P, ions_, psi_, psi_wrapper_, hf_force_, wf_grad_);
  else
    value_ = ham_.evaluateIonDerivsDeterministic(P, ions_, psi_, hf_force_, pulay_force_, wf_grad_);

  if (useSpaceWarp_)
  {
    Forces el_grad;
    el_grad.resize(P.getTotalNum());
    el_grad = 0;

    ham_.evaluateElecGrad(P, psi_, el_grad, delta_);
    swt_.computeSWT(P, ions_, el_grad, P.G, sw_pulay_, sw_grad_);
  }

  //Now we compute the regularizer.
  //WE ASSUME THAT psi_.evaluateLog(P) HAS ALREADY BEEN CALLED AND Grad(logPsi)
  //IS ALREADY UP TO DATE FOR THIS CONFIGURATION.

  f_epsilon_ = compute_regularizer_f(psi_.G, reg_epsilon_);


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
      plist[myindex++] = -hf_force_[iat][iondim] * f_epsilon_;
      plist[myindex++] = -(pulay_force_[iat][iondim] + sw_pulay_[iat][iondim]) * f_epsilon_;
      plist[myindex++] = -value_ * (wf_grad_[iat][iondim] + sw_grad_[iat][iondim]) * f_epsilon_;
      plist[myindex++] = -(wf_grad_[iat][iondim] + sw_grad_[iat][iondim]) * f_epsilon_;
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
      plist[myindex++] = -hf_force_[iat][iondim] * f_epsilon_;
      plist[myindex++] = -(pulay_force_[iat][iondim] + sw_pulay_[iat][iondim]) * f_epsilon_;
      plist[myindex++] = -value_ * (wf_grad_[iat][iondim] + sw_grad_[iat][iondim]) * f_epsilon_;
      plist[myindex++] = -(wf_grad_[iat][iondim] + sw_grad_[iat][iondim]) * f_epsilon_;
    }
  }
};

ACForce::RealType ACForce::compute_regularizer_f(const ParticleSet::ParticleGradient& G, const RealType epsilon)
{
  //epsilon=0 corresponds to no regularization.  However, since
  //epsilon ends up in denominators, return 1 here.
  if (std::abs(epsilon) < 1e-6)
    return 1.0;

  RealType gdotg = 0.0;
#if defined(QMC_COMPLEX)
  gdotg = Dot_CC(G, G);
#else
  gdotg = Dot(G, G);
#endif

  RealType gmag = std::sqrt(gdotg);
  RealType x;

  RealType regvalue = 0.0;
  //x = grad(logpsi)/|grad(logpsi)|^2 = 1/|grad(logpsi)|.
  //
  //Argument of polynomial is x/epsilon=1/(epsilon*|grad(logpsi)|)
  double xovereps = 1.0 / (epsilon * gmag);
  if (xovereps >= 1.0)
    regvalue = 1.0;
  else
  {
    //There's a discrepancy between AIP Advances 10, 085213 (2020) and arXiv:2002.01434 for polynomial.
    //We choose the arXiv, because f(x=1)=1, as opposed to f(x=1)=-4.
    regvalue = 7.0 * std::pow(xovereps, 6.0) - 15.0 * std::pow(xovereps, 4.0) + 9.0 * std::pow(xovereps, 2.0);
  }
  return regvalue;
};
} // namespace qmcplusplus
