//////////////////////////////////////////////////////////////////
// (c) Copyright 2003  by Jeongnim Kim
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//   National Center for Supercomputing Applications &
//   Materials Computation Center
//   University of Illinois, Urbana-Champaign
//   Urbana, IL 61801
//   e-mail: jnkim@ncsa.uiuc.edu
//   Tel:    217-244-6319 (NCSA) 217-333-3324 (MCC)
//
// Supported by
//   National Center for Supercomputing Applications, UIUC
//   Materials Computation Center, UIUC
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
#ifndef QMCPLUSPLUS_BAREPRESSURE_H
#define QMCPLUSPLUS_BAREPRESSURE_H
#include "Particle/ParticleSet.h"
#include "Particle/WalkerSetRef.h"
#include "QMCHamiltonians/QMCHamiltonianBase.h"
#include "ParticleBase/ParticleAttribOps.h"
#include "OhmmsData/ParameterSet.h"



namespace qmcplusplus
{

/** @ingroup hamiltonian
 @brief Evaluate the Bare Pressure.
 P=/frac{2T+V}{d* /Omega}
 where d is the dimension of space and /Omega is the volume.
**/

struct Pressure: public QMCHamiltonianBase
{
  double pNorm;
//     bool ZV;
//     bool ZB;

  /** constructor
   *
   * Pressure operators need to be re-evaluated during optimization.
   */
  Pressure(ParticleSet& P)
  {
    UpdateMode.set(OPTIMIZABLE,1);
    pNorm = 1.0/(P.Lattice.DIM*P.Lattice.Volume);
  }
  ///destructor
  ~Pressure() { }

  void resetTargetParticleSet(ParticleSet& P)
  {
    pNorm = 1.0/(P.Lattice.DIM*P.Lattice.Volume);
  }

  inline Return_t
  evaluate(ParticleSet& P)
  {
    Value=2.0*P.PropertyList[LOCALENERGY]-P.PropertyList[LOCALPOTENTIAL];
    Value*=pNorm;
    return 0.0;
  }

  inline Return_t
  evaluate(ParticleSet& P, vector<NonLocalData>& Txy)
  {
    return evaluate(P);
  }

  /** implements the virtual function.
   *
   * Nothing is done but should check the mass
   */

  bool put(xmlNodePtr cur )
  {
    return true;
  }

//     bool put(xmlNodePtr cur, ParticleSet& P, QMCHamiltonian* H) {
//       xmlNodePtr tcur = cur->children;
//
//       double RPAKCut= -1.0;
//       string RPAPCorr("ZB");
//       string RPAPfunc("RPA_LR");
//       ParameterSet nattrib;
//       OhmmsAttributeSet attrib;
//       attrib.add(RPAPCorr,"etype" );
//       attrib.add(RPAPfunc,"functor" );
//       attrib.put(cur);
//       nattrib.add(RPAKCut,"kc","double");
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
//         app_log() <<" using bare estimator "<<endl;;
//       }

//       return true;
//     }

  bool get(std::ostream& os) const
  {
    os << "Pressure";
    return true;
  }

  QMCHamiltonianBase* makeClone(ParticleSet& qp, TrialWaveFunction& psi)
  {
    return new Pressure(qp);
  }

};
}
#endif

/***************************************************************************
 * $RCSfile$   $Author: jnkim $
 * $Revision: 1581 $   $Date: 2007-01-04 10:02:14 -0600 (Thu, 04 Jan 2007) $
 * $Id: Pressure.h 1581 2007-01-04 16:02:14Z jnkim $
 ***************************************************************************/

