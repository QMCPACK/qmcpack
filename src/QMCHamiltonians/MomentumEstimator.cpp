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
#include <QMCWaveFunctions/TrialWaveFunction.h>
#include <Numerics/e2iphi.h>
#include <Numerics/OhmmsBlas.h>
#include <OhmmsData/AttributeSet.h>
#include <Utilities/SimpleParser.h>
#include <Particle/DistanceTableData.h>

namespace qmcplusplus 
{

  MomentumEstimator::MomentumEstimator(ParticleSet& elns, TrialWaveFunction& psi)
    :M(4), refPsi(psi)
  {
    UpdateMode.set(COLLECTABLE,1);
    psi_ratios.resize(elns.getTotalNum());
    kdotp.resize(elns.getTotalNum());
    phases.resize(elns.getTotalNum());
  }

  void MomentumEstimator::resetTargetParticleSet(ParticleSet& P)
  {
  }

  MomentumEstimator::Return_t MomentumEstimator::evaluate(ParticleSet& P)
  {
    const int np=P.getTotalNum();
    nofK=0.0;
    compQ=0.0;

    //will use temp[i].r1 for the Compton profile
    const vector<DistanceTableData::TempDistType>& temp(P.DistTables[0]->Temp);
    for(int s=0; s<M; ++s)
    {
      PosType newpos;
      for(int i=0; i<OHMMS_DIM;++i) newpos[i]=myRNG();

      P.makeVirtualMoves(newpos);
      refPsi.get_ratios(P,psi_ratios);
      P.rejectMove(0); //restore P.R[0] to the orginal position

      ////debug get_ratios with ratio, use it whenever an OrbitalBase implements get_ratios
      //vector<RealType> r_org(np);
      //for(int i=0; i<np; ++i)
      //{
      //  PosType delta=newpos-P.R[i];
      //  P.makeMove(i,delta);
      //  r_org[i]=refPsi.ratio(P,i);
      //  P.rejectMove(i);
      //  cout << "ratio("<<i<<")=" << r_org[i] << " diff=" << r_org[i]-psi_ratios[i] << endl;
      //}
      //APP_ABORT("Done with test");

      for(int ik=0; ik < kPoints.size(); ++ik)
      {
        RealType kdotp_primed=dot(kPoints[ik],newpos);
        for(int i=0; i<np; ++i) kdotp[i]=kdotp_primed-dot(kPoints[ik],P.R[i]);
        eval_e2iphi(np,kdotp.data(),phases.data());
        nofK[ik]+=real(BLAS::dot(np,phases.data(),psi_ratios.data()));
      }

      for(int iq=0; iq < Q.size(); ++iq)
      {
        for(int i=0; i<np; ++i) kdotp[i]=Q[iq]*temp[i].r1;
        eval_e2iphi(np,kdotp.data(),phases.data());
        compQ[iq]+=real(BLAS::dot(np,phases.data(),psi_ratios.data()));
      }
    }

    //need normalization factor
    int j=myIndex;
    for(int ik=0; ik<nofK.size(); ++ik,++j) P.Collectables[j]+=norm_nofK*nofK[ik];
    for(int iq=0; iq<compQ.size(); ++iq,++j) P.Collectables[j]+=norm_compQ*compQ[iq];

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
    //need to build kPoints list and normalization
    //NormFactor=1.0/static_cast<RealType>(elns.getTotalNum());
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
    myclone->resize(kPoints,Q);
    return myclone;
  }

  void MomentumEstimator::resize(const vector<PosType>& kin, const vector<RealType>& qin)
  {
    //copy kpoints
    kPoints=kin;
    nofK.resize(kin.size());
  
    //copy q
    Q=qin;
    compQ.resize(qin.size());
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
