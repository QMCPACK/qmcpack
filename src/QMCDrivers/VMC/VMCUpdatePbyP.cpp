//////////////////////////////////////////////////////////////////
// (c) Copyright 2003- by Jeongnim Kim
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//   Jeongnim Kim
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
//   Department of Physics, Ohio State University
//   Ohio Supercomputer Center
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
#include "QMCDrivers/VMC/VMCUpdatePbyP.h"
#include "QMCDrivers/DriftOperators.h"
#include "Message/OpenMP.h"

namespace qmcplusplus { 

  VMCUpdatePbyP::VMCUpdatePbyP(ParticleSet& w, TrialWaveFunction& psi, 
      QMCHamiltonian& h, RandomGenerator_t& rg):
    QMCUpdateBase(w,psi,h,rg), nSubSteps(1)
    { 
    }

  VMCUpdatePbyP::~VMCUpdatePbyP()
  {
  }

  bool VMCUpdatePbyP::put(xmlNodePtr cur)
  {
    ParameterSet params;
    params.add(nSubSteps,"subSteps","int"); params.add(nSubSteps,"substeps","int");
    return params.put(cur);
  }

  void VMCUpdatePbyP::advanceWalkers(WalkerIter_t it, WalkerIter_t it_end) 
  {

    while(it != it_end) 
    {
      Walker_t& thisWalker(**it);
      Walker_t::Buffer_t& w_buffer(thisWalker.DataSet);

      W.R = thisWalker.R;
      w_buffer.rewind();
      W.copyFromBuffer(w_buffer);
      Psi.copyFromBuffer(W,w_buffer);

      RealType psi_old = thisWalker.Properties(SIGN);
      RealType psi = psi_old;

      for(int iter=0; iter<nSubSteps; iter++) {

        makeGaussRandomWithEngine(deltaR,RandomGen);
        bool stucked=true;
        for(int iat=0; iat<W.getTotalNum(); iat++) {

          PosType dr = m_sqrttau*deltaR[iat];
          PosType newpos = W.makeMove(iat,dr);

          RealType ratio = Psi.ratio(W,iat);
          RealType prob = std::min(1.0e0,ratio*ratio);
          if(RandomGen() < prob) { 
            stucked=false;
            ++nAccept;
            W.acceptMove(iat);
            Psi.acceptMove(W,iat);
          } else {
            ++nReject; 
            W.rejectMove(iat); 
            Psi.rejectMove(iat);
          }
        }
        if(stucked) {
          ++nAllRejected;
        }
      }

      thisWalker.R = W.R;
      w_buffer.rewind();
      W.updateBuffer(w_buffer);
      RealType logpsi = Psi.updateBuffer(W,w_buffer);
      //W.copyToBuffer(w_buffer);
      //RealType logpsi = Psi.evaluate(W,w_buffer);

      RealType eloc=H.evaluate(W);

      thisWalker.resetProperty(logpsi,Psi.getPhase(),eloc);
      H.saveProperty(thisWalker.getPropertyBase());
      ++it;
    }
  }

  
  /// Constructor.
  VMCUpdatePbyPWithDrift::VMCUpdatePbyPWithDrift(ParticleSet& w, TrialWaveFunction& psi, 
      QMCHamiltonian& h, RandomGenerator_t& rg):
    QMCUpdateBase(w,psi,h,rg) 
    { 
    }

  VMCUpdatePbyPWithDrift::~VMCUpdatePbyPWithDrift()
  {
  }

  void VMCUpdatePbyPWithDrift::advanceWalkers(WalkerIter_t it, WalkerIter_t it_end) 
  {

    while(it != it_end) 
    {
      Walker_t& thisWalker(**it);
      Walker_t::Buffer_t& w_buffer(thisWalker.DataSet);

      W.R = thisWalker.R;
      w_buffer.rewind();
      W.copyFromBuffer(w_buffer);
      Psi.copyFromBuffer(W,w_buffer);

      RealType psi_old = thisWalker.Properties(SIGN);
      RealType psi = psi_old;
      //create a 3N-Dimensional Gaussian with variance=1
      makeGaussRandomWithEngine(deltaR,RandomGen);

      bool moved = false;

      for(int iat=0; iat<W.getTotalNum(); iat++) {

        PosType dr = m_sqrttau*deltaR[iat]+thisWalker.Drift[iat];
        PosType newpos = W.makeMove(iat,dr);

        //RealType ratio = Psi.ratio(W,iat);
        RealType ratio = Psi.ratio(W,iat,dG,dL);

        G = W.G+dG;

        //RealType forwardGF = exp(-0.5*dot(deltaR[iat],deltaR[iat]));
        //dr = (*it)->R[iat]-newpos-Tau*G[iat]; 
        //RealType backwardGF = exp(-oneover2tau*dot(dr,dr));
        RealType logGf = -0.5e0*dot(deltaR[iat],deltaR[iat]);

        RealType scale=getDriftScale(Tau,G);
        //COMPLEX WARNING
        //dr = thisWalker.R[iat]-newpos-scale*G[iat];
        dr = thisWalker.R[iat]-newpos-scale*real(G[iat]);

        RealType logGb = -m_oneover2tau*dot(dr,dr);

        RealType prob = std::min(1.0e0,ratio*ratio*exp(logGb-logGf));

        //alternatively
        if(RandomGen() < prob) { 
          moved = true;
          ++nAccept;
          W.acceptMove(iat);
          Psi.acceptMove(W,iat);
          W.G = G;
          W.L += dL;

          //thisWalker.Drift = scale*G;
          assignDrift(scale,G,thisWalker.Drift);

        } else {
          ++nReject; 
          W.rejectMove(iat); Psi.rejectMove(iat);
        }
      }

      if(moved) {
        w_buffer.rewind();
        W.copyToBuffer(w_buffer);
        psi = Psi.evaluate(W,w_buffer);

        thisWalker.R = W.R;
        RealType eloc=H.evaluate(W);
        thisWalker.resetProperty(log(abs(psi)), psi,eloc);
        H.saveProperty(thisWalker.getPropertyBase());
      }
      else {
        ++nAllRejected;
      }
      ++it;
    }
  }
}

/***************************************************************************
 * $RCSfile: VMCUpdatePbyP.cpp,v $   $Author: jnkim $
 * $Revision: 1.25 $   $Date: 2006/10/18 17:03:05 $
 * $Id: VMCUpdatePbyP.cpp,v 1.25 2006/10/18 17:03:05 jnkim Exp $ 
 ***************************************************************************/
