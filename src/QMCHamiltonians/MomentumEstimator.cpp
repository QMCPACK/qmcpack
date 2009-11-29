//////////////////////////////////////////////////////////////////
// (c) Copyright 2008-  by Jeongnim Kim
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
#include <QMCHamiltonians/MomentumEstimator.h>
#include <OhmmsData/AttributeSet.h>
#include <Utilities/SimpleParser.h>

namespace qmcplusplus 
{

  MomentumEstimator::MomentumEstimator(ParticleSet& elns, TrialWaveFunction& psi)
    :refPsi(psi)
  {
    UpdateMode.set(COLLECTABLE,1);
    psi_ratios.resize(elns.getTotalNum());
  }

  void MomentumEstimator::resetTargetParticleSet(ParticleSet& P)
  {
  }

  MomentumEstimator::Return_t MomentumEstimator::evaluate(ParticleSet& P)
  {
    PosType ru;
    for(int i=0; i<OHMMS_DIM;++i) ru[i]=myRNG();

    P.makeVirtualMoves(ru);
    refPsi.get_ratios(P,psi_ratios);
    P.rejectMove(0); //restore P.R[0] to the orginal position

    const int nk=kPoints.size();
    const int np=P.getTotalNum();
    for(int ik=0; ik<nk; ++ik)
    {
      RealType kdotp_primed=dot(newpos,kPoints[ik]);
      for(int i=0; i<np; ++i)
        kdotp[i]=kdotp_primed-dot(kPoints[ik],P.R[i]);
      eval_e2iphi(kdotp.data(),phases[ik]);
    }

    //multiple phases by ratios
    for(int ik=0; ik<nk; ++ik) convert(dot(phases[ik],psi_ratios.data()),nofK[ik]);

    //need normalization factor
    for(int ik=0,j=myIndex; ik<nk; ++ik,++j) P.Collectables[j]+=nofK[ik];

    return 0.0;
  }

  void MomentumEstimator::registerCollectables(vector<observable_helper*>& h5list
      , hid_t gid) const
  {

  }


  void MomentumEstimator::addObservables(PropertySetType& plist, BufferType& collectables)
  {
    myIndex=collectables.size();
  }


  void MomentumEstimator::setObservables(PropertySetType& plist)
  {
    //std::copy(gof_r.first_address(),gof_r.last_address(),plist.begin()+myIndex);
  }

  void MomentumEstimator::setParticlePropertyList(PropertySetType& plist
      , int offset)
  {
    //std::copy(gof_r.first_address(),gof_r.last_address(),plist.begin()+myDebugIndex+offset);
  }

  bool MomentumEstimator::put(xmlNodePtr cur)
  {
    return true;
  }

  bool MomentumEstimator::get(std::ostream& os) const
  {
    return true;
  }

  QMCHamiltonianBase* MomentumEstimator::makeClone(ParticleSet& qp
      , TrialWaveFunction& psi)
  {
    MomentumEstimator* myclone=new MomentumEstimator(qp,psi);
    myclone->resize(kPoints);
    return myclone;
  }

  void MomentumEstimator::resize(const vector<PosType>& kin)
  {
    //copy kpoints
    kPoints=kin;
    nofK.resize(kin.size());
    phases.resize(kin.size(),psi_ratios.size());
  }

  void MomentumEstimator::setRandomGenerator(RandomGenerator_t* rng)
  {
    //simply copy it
    myRNG=*rng;
  }
}

/***************************************************************************
 * $RCSfile$   $Author: jnkim $
 * $Revision: 2945 $   $Date: 2008-08-05 10:21:33 -0500 (Tue, 05 Aug 2008) $
 * $Id: ForceBase.h 2945 2008-08-05 15:21:33Z jnkim $ 
 ***************************************************************************/
