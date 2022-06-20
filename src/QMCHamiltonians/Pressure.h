//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


#ifndef QMCPLUSPLUS_BAREPRESSURE_H
#define QMCPLUSPLUS_BAREPRESSURE_H
#include "Particle/ParticleSet.h"
#include "QMCDrivers/WalkerProperties.h"
#include "QMCHamiltonians/OperatorBase.h"
#include "ParticleBase/ParticleAttribOps.h"
#include "OhmmsData/ParameterSet.h"


namespace qmcplusplus
{
/** @ingroup hamiltonian
 @brief Evaluate the Bare Pressure.
 P=/frac{2T+V}{d* /Omega}
 where d is the dimension of space and /Omega is the volume.
**/

struct Pressure : public OperatorBase
{
  using WP = WalkerProperties::Indexes;
  double pNorm;
  //     bool ZV;
  //     bool ZB;

  /** constructor
   *
   * Pressure operators need to be re-evaluated during optimization.
   */
  Pressure(ParticleSet& P)
  {
    update_mode_.set(OPTIMIZABLE, 1);
    pNorm = 1.0 / (P.getLattice().DIM * P.getLattice().Volume);
  }
  ///destructor
  ~Pressure() override {}

  void resetTargetParticleSet(ParticleSet& P) override { pNorm = 1.0 / (P.getLattice().DIM * P.getLattice().Volume); }

  inline Return_t evaluate(ParticleSet& P) override
  {
    value_ = 2.0 * P.PropertyList[WP::LOCALENERGY] - P.PropertyList[WP::LOCALPOTENTIAL];
    value_ *= pNorm;
    return 0.0;
  }

  /** implements the virtual function.
   *
   * Nothing is done but should check the mass
   */

  bool put(xmlNodePtr cur) override { return true; }

  //     bool put(xmlNodePtr cur, ParticleSet& P, QMCHamiltonian* H) {
  //       xmlNodePtr tcur = cur->children;
  //
  //       double RPAKCut= -1.0;
  //       std::string RPAPCorr("ZB");
  //       std::string RPAPfunc("RPA_LR");
  //       ParameterSet nattrib;
  //       OhmmsAttributeSet attrib;
  //       attrib.add(RPAPCorr,"etype" );
  //       attrib.add(RPAPfunc,"functor" );
  //       attrib.put(cur);
  //       nattrib.add(RPAKCut,"kc");
  //       nattrib.put(cur);

  //       if (RPAPCorr=="ZB"){
  //         ZB=true;
  //         ZV=false;
  //         bpcorr = new RPAPressureCorrection(P);
  //         bpcorr-> put(cur, P);
  //         H->addOperator(bpcorr,"ZVterm");
  //       }
  //       else if (RPAPCorr=="ZVZB"){
  //         ZB=true;
  //         ZV=true;
  //         bpcorr = new RPAPressureCorrection(P);
  //         wfderivE = new RPADerivEnergy(P);
  //         wfderivE2 = new RPADerivEnergy2(P);
  //         wfEP = new RPAEnergyPressure(P);
  //         bpcorr-> put(cur, P);
  //         wfderiv = new RPADeriv(P);
  //         wfderiv -> put(cur, bpcorr);
  //         wfderivE -> put(cur, bpcorr, H);
  //         wfderivE2 -> put(cur, bpcorr, H);
  //         wfEP -> put(cur, bpcorr, H);
  //         potkin = new RPAPotKin(P);
  //         H->addOperator(potkin,"PotKin");
  //         H->addOperator(wfEP,"EPterm");
  //         H->addOperator(bpcorr,"ZVterm");
  //         H->addOperator(wfderiv,"dpsi");
  //         H->addOperator(wfderivE,"Edpsi");
  //         H->addOperator(wfderivE2,"Tdpsi");
  //       }
  //       else if (RPAPCorr=="ZV"){
  //         ZV=true;
  //         ZB=false;
  //         bpcorr = new RPAPressureCorrection(P);
  //         bpcorr-> put(cur, P);
  //         H->addOperator(bpcorr,"ZVterm");
  //       }
  //       else if (RPAPCorr=="none"){
  //         ZV=false;
  //         ZB=false;
  //         app_log() <<" using bare estimator "<< std::endl;;
  //       }

  //       return true;
  //     }

  bool get(std::ostream& os) const override
  {
    os << "Pressure";
    return true;
  }

  std::unique_ptr<OperatorBase> makeClone(ParticleSet& qp, TrialWaveFunction& psi) final
  {
    return std::make_unique<Pressure>(qp);
  }
};
} // namespace qmcplusplus
#endif
