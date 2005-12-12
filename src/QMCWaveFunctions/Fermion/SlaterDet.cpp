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
#include "QMCWaveFunctions/Fermion/SlaterDet.h"
#include "Message/Communicate.h"
#include "Utilities/OhmmsInfo.h"

namespace qmcplusplus {

    ///destructor
  SlaterDet::~SlaterDet() { }

    ///add a new DiracDeterminant to the list of determinants
    void SlaterDet::add(Determinant_t* det) { 
      int last=Dets.size();
      Dets.push_back(det);
      M[last+1]=M[last]+Dets[last]->rows();
      DetID.insert(DetID.end(),det->rows(),last);
    }

    ///reset all the Dirac determinants, Optimizable is true
    void SlaterDet::reset() {  
      if(Optimizable) for(int i=0; i<Dets.size(); i++) Dets[i]->reset();
    }

    void SlaterDet::resetTargetParticleSet(ParticleSet& P) {
      //BasisSet->resetTargetParticleSet(P);
      LOGMSG("\nSlaterDet::resetTargetParticleSet")
      for(int i=0; i<Dets.size(); i++) Dets[i]->resetTargetParticleSet(P);
    }

    SlaterDet::ValueType 
    SlaterDet::evaluate(ParticleSet& P, 
	     ParticleSet::ParticleGradient_t& G, 
	     ParticleSet::ParticleLaplacian_t& L) {
      ValueType psi = 1.0;
      for(int i=0; i<Dets.size(); i++) psi *= Dets[i]->evaluate(P,G,L);
      return psi;
    }

    SlaterDet::ValueType 
    SlaterDet::evaluateLog(ParticleSet& P, 
	        ParticleSet::ParticleGradient_t& G, 
	        ParticleSet::ParticleLaplacian_t& L) {
      //@attention BasisSet::evaluate is to be called but the due to the bugs, it is commented out.
      //if(BasisSet == 0) 
      //{
      //  ERRORMSG("SlaterDet::BasisSet is not assigned")
      //  OHMMS::Controller->abort();
      //}
      //BasisSet->evaluate(P);

      ValueType psi = 1.0;
      for(int i=0; i<Dets.size(); i++) psi *= Dets[i]->evaluate(P,G,L);
      SignValue = (psi<0.0)?-1.0:1.0;
      LogValue = log(abs(psi));
      return LogValue;
    }

    SlaterDet::ValueType SlaterDet::registerData(ParticleSet& P, PooledData<RealType>& buf){

      //BasisSet->evaluate(P);

      ValueType psi = 1.0;
      for(int i=0; i<Dets.size(); i++) 
        psi *= Dets[i]->registerData(P,buf);
      SignValue = (psi<0.0)?-1.0:1.0;
      LogValue = log(abs(psi));
      return LogValue;
    }
    
    SlaterDet::ValueType SlaterDet::updateBuffer(ParticleSet& P, PooledData<RealType>& buf){
      ValueType psi = 1.0;
      for(int i=0; i<Dets.size(); i++) psi *= Dets[i]->updateBuffer(P,buf);
      SignValue = (psi<0.0)?-1.0:1.0;
      LogValue = log(abs(psi));
      return LogValue;
    }

    void SlaterDet::copyFromBuffer(ParticleSet& P, PooledData<RealType>& buf) {
      for(int i=0; i<Dets.size(); i++) 	Dets[i]->copyFromBuffer(P,buf);
    }

    /** reimplements the virtual function
     *
     * The DiractDeterminants of SlaterDet need to save the inverse
     * of the determinant matrix to evaluate ratio
     */
    void SlaterDet::dumpToBuffer(ParticleSet& P, PooledData<RealType>& buf) {
      for(int i=0; i<Dets.size(); i++) 	Dets[i]->dumpToBuffer(P,buf);
    }

    /** reimplements the virtual function
     *
     * Matching function to dumpToBuffer.
     */
    void SlaterDet::dumpFromBuffer(ParticleSet& P, PooledData<RealType>& buf) {
      for(int i=0; i<Dets.size(); i++) 	Dets[i]->dumpFromBuffer(P,buf);
    }

    SlaterDet::ValueType 
      SlaterDet::evaluate(ParticleSet& P, PooledData<RealType>& buf) {

      //BasisSet->evaluate(P);

      ValueType r=1.0;
      for(int i=0; i<Dets.size(); i++) 	r *= Dets[i]->evaluate(P,buf);
      return r;
    }

}
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
