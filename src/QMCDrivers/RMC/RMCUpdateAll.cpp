//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Raymond Clay III, j.k.rofling@gmail.com, Lawrence Livermore National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    

#include "QMCDrivers/RMC/RMCUpdateAll.h"
#include "QMCDrivers/DriftOperators.h"
#include "ParticleBase/ParticleAttribOps.h"
#include "Message/OpenMP.h"
#include "Configuration.h"
#include "Particle/Reptile.h"
#include <cmath>
#include "OhmmsData/ParameterSet.h"

//////////////////////////////////////////////////////////////////////////
//
//  This driver proposes all-electron moves like in DMCUpdateAll.cpp.
//
//  For H4 and H2 systems, it appears as though the Umrigar pbyp scaled drift has similar instabilities
//  as in the RQMCMultiple.cpp driver.  Hence, the bare propagator with filtered local energies is used.
//
//  The symmetric link action (Ceperley and Pierleoni) is hardcoded into this driver.
//
//////////////////////////////////////////////////////////////////////////////

namespace qmcplusplus
{

/// Constructor.
  RMCUpdateAllWithDrift::RMCUpdateAllWithDrift (MCWalkerConfiguration & w,
						TrialWaveFunction & psi,
						QMCHamiltonian & h,
						RandomGenerator_t & rg,
						std::vector < int >act,
						std::vector <
						int >tp):QMCUpdateBase (w,
									psi,
									h,
									rg),
    Action (act), TransProb (tp)
  {
    scaleDrift = false;
    actionType = SYM_ACTION;
    vmcToDoSteps = 0;
    equilToDoSteps = 0;
    vmcSteps = 0;
    equilSteps = 0;
  }

  RMCUpdateAllWithDrift::~RMCUpdateAllWithDrift ()
  {
  }

  bool RMCUpdateAllWithDrift::put (xmlNodePtr cur)
  {

    QMCUpdateBase::put (cur);
    ParameterSet m_param;
    // bool usedrift=false;
    std::string action = "SLA";
    std::string usedrift = "no";
    m_param.add (usedrift, "useScaledDrift", "string");
    m_param.add (action, "Action", "string");
    m_param.add (equilSteps, "equilsteps", "int");
    m_param.add (equilSteps, "equilSteps", "int");

    // m_param.add(scaleDrift,"scaleDrift");
    m_param.put (cur);

    bool driftoption = (usedrift == "yes" || usedrift == "Yes"
			|| usedrift == "True" || usedrift == "true");

    if (driftoption)
      {
	scaleDrift = true;
        if (omp_get_thread_num()==0)
	  app_log () << "  Using Umrigar scaled drift\n";
	// H.rejectedMove(W,thisWalker);
      }
    else
      {
        if (omp_get_thread_num()==0)
	  app_log () << "  Using non-scaled drift\n";
      }

    if (action == "DMC")
      {
	actionType = DMC_ACTION;
	if (omp_get_thread_num()==0)
          app_log () << "  Using DMC link-action\n";
      }
    else
      {
        if (omp_get_thread_num()==0)
	  app_log () << "  Using Symmetrized Link-Action\n";
      }

    return true;
  }
//This performs an all electron VMC step in the current direction of the reptile.  
//Performs all evaluations to update the action
  void RMCUpdateAllWithDrift::advanceWalkersVMC ()
  {

    IndexType direction = W.reptile->direction;
    IndexType forward = (1 - direction) / 2;
    IndexType backward = (1 + direction) / 2;
    Walker_t & curhead = W.reptile->getHead ();
    W.loadWalker (curhead, false);
    RealType nodecorr = 1;
    if (scaleDrift == true)
      RealType nodecorr =
	setScaledDriftPbyPandNodeCorr (m_tauovermass, curhead.G, drift);
    else
    assignDrift (m_tauovermass, curhead.G, drift);
    //app_log()<<"Sign head = "<<curhead.Properties(SIGN)<< std::endl;
    //app_log()<<"Old phase = "<<Psi.getPhase()<< std::endl;
    makeGaussRandomWithEngine (deltaR, RandomGen);
    RealType r2proposed = Dot (deltaR, deltaR);
    RealType r2accept = 0.0;
//      W.reptile->r2prop += r2proposed;
//      W.reptile->r2samp++;
    if (!W.makeMoveWithDrift (curhead, drift, deltaR, m_sqrttau))
      {
	++nReject;
	H.rejectedMove (W, curhead);
	// curhead.Age+=1;
	//W.reptile->flip();
	return;
      }

    RealType logpsi (Psi.evaluateLog (W));
    RealType logGf = -0.5 * Dot (deltaR, deltaR);
    RealType Action_forward = -0.5 * logGf;
    curhead.Properties (W.reptile->TransProb[forward]) =
      -0.5 * Dot (deltaR, deltaR);
    curhead.Properties (W.reptile->Action[forward]) =
      0.5 * 0.5 * Dot (deltaR, deltaR);

    W.reptile->saveTransProb (curhead, +1, logGf);
    W.reptile->saveAction (curhead, +1, Action_forward);

    Walker_t::ParticlePos_t fromdeltaR (deltaR);


    if (scaleDrift == true)
      setScaledDrift (m_tauovermass, W.G, drift);
    else
      assignDrift (m_tauovermass, W.G, drift);
    fromdeltaR = curhead.R - W.R - drift;
    EstimatorRealType *restrict new_headProp (W.getPropertyBase ());

    RealType logGb = -m_oneover2tau * Dot (fromdeltaR, fromdeltaR);
    RealType Action_backward = -0.5 * logGb;

    W.reptile->saveTransProb (W, -1, logGb);
    //W.reptile->saveAction(W,-1,Action_backward);

    W.Properties (W.reptile->TransProb[backward]) =
      -m_oneover2tau * Dot (fromdeltaR, fromdeltaR);
    W.Properties (W.reptile->Action[backward]) =
      0.5 * m_oneover2tau * Dot (fromdeltaR, fromdeltaR);

    Walker_t & lastbead (W.reptile->getTail ()),
      nextlastbead (W.reptile->getNext ());
    //Implementing the fixed-node approximation.  If phase difference is not a multiple of 2pi, bounce away from node.
    RealType newphase = Psi.getPhase ();
    RealType phasediff = newphase - curhead.Properties (SIGN);
    //Reject & bounce if node crossed.
    if (branchEngine->phaseChanged (Psi.getPhaseDiff ()))
      {
	++nReject;
	H.rejectedMove (W, curhead);
	// curhead.Age+=1;
	// W.reptile->flip();
	//app_log()<<"hit a node.  Bouncing...\n";
	return;
      }
    RealType eloc = H.evaluate (W);
    new_headProp[Action[2]] = 0.5 * Tau * eloc;

    ////////////////////////////////////////////////////////////////////////
    ///  Like DMC, this filters the local energy to ignore divergences near pathological points in phase space.
    ////////////////////////////////////////////////////////////////////////
//      RealType eest = W.reptile->eest;

//      RealType fbet = std::max(eest - curhead.Properties(LOCALENERGY), eest - eloc);
    //   app_log()<<"eval = "<<eest<<" estdev="<<stddev<< std::endl;
//      RealType rawcutoff=100*std::sqrt(W.reptile->evar);
//      RealType cutoffmax = 1.5*rawcutoff;
//      RealType cutoff=1;
//      if (fbet > rawcutoff)
//        cutoff = 1-(fbet - rawcutoff)/(rawcutoff*0.5);
//      if( fbet > cutoffmax )
//        cutoff=0;
    //////////////////////////////////////////////////////////////////////////
//      RealType tauscale = W.reptile->tauscale;
//      W.Properties(W.reptile->Action[2])= 0.5*Tau*eloc*cutoff*tauscale;
    RealType dS = 0;
    RealType acceptProb = 1;
    if (actionType == SYM_ACTION)
      {
	RealType oldhead_logpsi = curhead.Properties (LOGPSI);
	RealType oldtail_logpsi = lastbead.Properties (LOGPSI);
	RealType newhead_logpsi = logpsi;
	RealType newtail_logpsi = nextlastbead.Properties (LOGPSI);

	RealType oldhead_e = curhead.Properties (LOCALENERGY);
	RealType oldtail_e = lastbead.Properties (LOCALENERGY);
	RealType newhead_e = W.Properties (LOCALENERGY);
	RealType newtail_e = nextlastbead.Properties (LOCALENERGY);

	RealType head_forward = W.reptile->getTransProb (curhead, +1);
	RealType head_backward = W.reptile->getTransProb (W, -1);
	RealType tail_forward = W.reptile->getTransProb (lastbead, +1);
	RealType tail_backward = W.reptile->getTransProb (nextlastbead, -1);

	RealType dS_head =
	  branchEngine->symLinkAction (head_forward, head_backward, newhead_e,
				       oldhead_e);
	RealType dS_tail =
	  branchEngine->symLinkAction (tail_forward, tail_backward, newtail_e,
				       oldtail_e);

	//   dS=dS_head-dS_tail;

	//   RealType dS_old = +(curhead.Properties(LOGPSI) + lastbead.Properties(LOGPSI) - logpsi - nextlastbead.Properties(LOGPSI))
	//         + curhead.Properties(W.reptile->Action[2]) + W.Properties(W.reptile->Action[2])
	//         + curhead.Properties(W.reptile->Action[forward]) + W.Properties(W.reptile->Action[backward])
	//         - (lastbead.Properties(W.reptile->Action[2]) + nextlastbead.Properties(W.reptile->Action[2]))
	//         - (lastbead.Properties(W.reptile->Action[forward]) + nextlastbead.Properties(W.reptile->Action[backward]));
	//   acceptProb=std::exp(-dS_old + (nextlastbead.Properties(W.reptile->TransProb[backward]) - curhead.Properties(W.reptile->TransProb[forward])));
	// acceptProb=std::min(1.0,std::exp(-dS + -(curhead.Properties(LOGPSI) + lastbead.Properties(LOGPSI) - logpsi - nextlastbead.Properties(LOGPSI)) + tail_backward - head_forward));

      }
    else
      {
	// dS = curhead.Properties(W.reptile->Action[2]) + W.Properties(W.reptile->Action[2])
	//- (lastbead.Properties(W.reptile->Action[2]) + nextlastbead.Properties(W.reptile->Action[2]));
	RealType dS_head =
	  branchEngine->DMCLinkAction (eloc,
				       curhead.Properties (LOCALENERGY));
	RealType dS_tail =
	  branchEngine->DMCLinkAction (lastbead.Properties (LOCALENERGY),
				       nextlastbead.Properties (LOCALENERGY));
	//dS=branchEngine->DMCLinkAction(eloc,curhead.Properties(LOCALENERGY)) - branchEngine->DMCLinkAction(lastbead.Properties(LOCALENERGY),nextlastbead.Properties(LOCALENERGY));
	//  dS=dS_head - dS_tail;
	//  acceptProb=std::min(1.0,std::exp(-dS ));
      }
    acceptProb =
      std::exp (logGb - logGf + 2.0 * (logpsi - curhead.Properties (LOGPSI)));


    /*      RealType dS = curhead.Properties(W.reptile->Action[2]) + W.Properties(W.reptile->Action[2])
       - (lastbead.Properties(W.reptile->Action[2]) + nextlastbead.Properties(W.reptile->Action[2]));
       RealType acceptProb=std::min(1.0,std::exp(-dS ));            */
    if (RandomGen () < acceptProb)
      {

	//Assuming the VMC step is fine, we are forcing the move.
	r2accept = r2proposed;
//        W.reptile->r2accept+=r2accept;
	MCWalkerConfiguration::Walker_t & overwriteWalker (W.reptile->
							   getNewHead ());

	W.saveWalker (overwriteWalker);
	overwriteWalker.Properties (LOCALENERGY) = eloc;
	overwriteWalker.Properties (W.reptile->Action[forward]) = 0;
	overwriteWalker.Properties (W.reptile->Action[backward]) =
	  W.Properties (W.reptile->Action[backward]);
	overwriteWalker.Properties (W.reptile->Action[2]) =
	  W.Properties (W.reptile->Action[2]);
	overwriteWalker.Properties (W.reptile->TransProb[forward]) =
	  W.Properties (W.reptile->TransProb[forward]);
	overwriteWalker.Properties (W.reptile->TransProb[backward]) =
	  W.Properties (W.reptile->TransProb[backward]);
	overwriteWalker.resetProperty (logpsi, Psi.getPhase (), eloc);
	H.auxHevaluate (W, overwriteWalker,true,false); //properties but not collectables.
	H.saveProperty (overwriteWalker.getPropertyBase ());
	overwriteWalker.Age = 0;
	++nAccept;
      }
    else
      {
	++nReject;
	H.rejectedMove (W, curhead);
	//    curhead.Age+=1;
	//    W.reptile->flip();
	return;
      }
  }


  void RMCUpdateAllWithDrift::initWalkers (WalkerIter_t it,
					   WalkerIter_t it_end)
  {
    IndexType initsteps = W.reptile->nbeads * 2;
    vmcSteps = W.reptile->nbeads + 1;

    for (; it != it_end; ++it)
      {
	W.R = (*it)->R;
	W.update ();
	RealType logpsi (Psi.evaluateLog (W));
	(*it)->G = W.G;
	(*it)->L = W.L;
	RealType nodecorr =
	  setScaledDriftPbyPandNodeCorr (Tau, MassInvP, W.G, drift);
	RealType ene = H.evaluate (W);
	H.auxHevaluate (W);
	(*it)->resetProperty (logpsi, Psi.getPhase (), ene, 0.0, 0.0,
			      nodecorr);
	(*it)->Weight = 1;
	H.saveProperty ((*it)->getPropertyBase ());
      }


    for (int n = 0; n < initsteps; n++)
      advanceWalkersVMC ();

  }

  void RMCUpdateAllWithDrift::advanceWalker (Walker_t& thisWalker, bool recompute)
  {
  }

  void RMCUpdateAllWithDrift::advanceWalkers (WalkerIter_t it,
					      WalkerIter_t it_end,
					      bool measure)
  {
    if (vmcToDoSteps > 0)
      {
	advanceWalkersVMC ();
	vmcToDoSteps--;
      }
    else if (vmcToDoSteps == 0 && equilToDoSteps > 0)
      {
	advanceWalkersRMC ();
	equilToDoSteps--;
      }
    else
      {
	advanceWalkersRMC ();
      }

  }

  void RMCUpdateAllWithDrift::advanceWalkersRMC ()
  {
    IndexType direction = W.reptile->direction;
    IndexType forward = (1 - direction) / 2;
    IndexType backward = (1 + direction) / 2;
    Walker_t & curhead = W.reptile->getHead ();
    Walker_t & centerbead = W.reptile->getCenter ();

//  if(centerbead.Age>=MaxAge)
//  {
//     vmcToDoSteps=vmcSteps;
//     equilToDoSteps=equilSteps;
//     app_log()<<"MaxAge for center bead exceeded.  Reequilibrating. "<<vmcSteps<<" "<<equilSteps<< std::endl;
//  }

    //We are going to monitor the center bead's age to determine whether we force
    //moves.  This is because the ends are less likely to get pinned.

//  centerbead.Age+=1;

    W.loadWalker (curhead, false);
    //RealType nodecorr=1;
    if (scaleDrift == true)
      setScaledDrift (m_tauovermass, curhead.G, drift);
    else
      assignDrift (m_tauovermass, curhead.G, drift);
    //app_log()<<"Sign head = "<<curhead.Properties(SIGN)<< std::endl;
    //app_log()<<"Old phase = "<<Psi.getPhase()<< std::endl;
    makeGaussRandomWithEngine (deltaR, RandomGen);
    RealType r2proposed = Dot (deltaR, deltaR);
    RealType r2accept = 0.0;
//  W.reptile->r2prop += r2proposed;
//  W.reptile->r2samp++;
    if (!W.makeMoveWithDrift (curhead, drift, deltaR, m_sqrttau))
      {
	++nReject;
	H.rejectedMove (W, curhead);
	curhead.Age += 1;
	W.reptile->flip ();
	return;
      }
    RealType logpsi (Psi.evaluateLog (W));
    // app_log()<<"Sign newhead = "<<W.Properties(SIGN)<< std::endl;
    //RealType* restrict old_headProp ((*it)->getPropertyBase());
    //old_headProp[TransProb[forward]]= 0.5*Dot(deltaR,deltaR);

    curhead.Properties (W.reptile->TransProb[forward]) =
      -0.5 * Dot (deltaR, deltaR);
    curhead.Properties (W.reptile->Action[forward]) =
      0.5 * 0.5 * Dot (deltaR, deltaR);

    RealType logGf = -0.5 * Dot (deltaR, deltaR);
    //W.reptile->saveTransProb(curhead,+1,logGf);

    Walker_t::ParticlePos_t fromdeltaR (deltaR);


    if (scaleDrift == true)
      setScaledDrift (m_tauovermass, W.G, drift);
    else
      assignDrift (m_tauovermass, W.G, drift);
    fromdeltaR = curhead.R - W.R - drift;
    EstimatorRealType *restrict new_headProp (W.getPropertyBase ());
    W.Properties (W.reptile->TransProb[backward]) =
      -m_oneover2tau * Dot (fromdeltaR, fromdeltaR);
    W.Properties (W.reptile->Action[backward]) =
      0.5 * m_oneover2tau * Dot (fromdeltaR, fromdeltaR);

    RealType logGb = -m_oneover2tau * Dot (fromdeltaR, fromdeltaR);

    // W.reptile->saveTransProb(W,-1, logGb); 

    Walker_t & lastbead (W.reptile->getTail ()),
      nextlastbead (W.reptile->getNext ());
    //Implementing the fixed-node approximation.  If phase difference is not a multiple of 2pi, bounce away from node.
    RealType newphase = Psi.getPhase ();
    RealType phasediff = newphase - curhead.Properties (SIGN);
    //Reject & bounce if node crossed.
    if (branchEngine->phaseChanged (Psi.getPhaseDiff ()))
      {
	++nReject;
	H.rejectedMove (W, curhead);
	curhead.Age += 1;
	W.reptile->flip ();
	//app_log()<<"hit a node.  Bouncing...\n";
	return;
      }
    RealType eloc = H.evaluate (W);
    W.Properties (LOCALENERGY) = eloc;
    //new_headProp[Action[2]]= 0.5*Tau*eloc;
    ////////////////////////////////////////////////////////////////////////
    ///  Like DMC, this filters the local energy to ignore divergences near pathological points in phase space.
    ////////////////////////////////////////////////////////////////////////
///  RealType eest = W.reptile->eest;
///  RealType fbet = std::max(eest - curhead.Properties(LOCALENERGY), eest - eloc);

///  RealType rawcutoff=100*std::sqrt(W.reptile->evar);
///  RealType cutoffmax = 1.5*rawcutoff;
/// RealType cutoff=1;
///  if (fbet > rawcutoff)
///   cutoff = 1-(fbet - rawcutoff)/(rawcutoff*0.5);
///  if( fbet > cutoffmax )
///    cutoff=0;
    //////////////////////////////////////////////////////////////////////////
///  RealType tauscale = W.reptile->tauscale;
///  W.Properties(W.reptile->Action[2])= 0.5*Tau*eloc*cutoff*tauscale;
    RealType dS = 0;
    RealType acceptProb = 0;

    if (actionType == SYM_ACTION)
      {
	RealType oldhead_logpsi = curhead.Properties (LOGPSI);
	RealType oldtail_logpsi = lastbead.Properties (LOGPSI);
	RealType newhead_logpsi = logpsi;
	RealType newtail_logpsi = nextlastbead.Properties (LOGPSI);

	RealType oldhead_e = curhead.Properties (LOCALENERGY);
	RealType oldtail_e = lastbead.Properties (LOCALENERGY);
	RealType newhead_e = W.Properties (LOCALENERGY);
	RealType newtail_e = nextlastbead.Properties (LOCALENERGY);

	RealType head_forward = W.reptile->getTransProb (curhead, +1);
	RealType head_backward = W.reptile->getTransProb (W, -1);
	RealType tail_forward = W.reptile->getTransProb (lastbead, +1);
	RealType tail_backward = W.reptile->getTransProb (nextlastbead, -1);

	//   RealType head_forward=curhead.Properties(W.reptile->TransProb[forward]);
	//  RealType head_backward=W.Properties(W.reptile->TransProb[backward]);
	//  RealType tail_forward=lastbead.Properties(W.reptile->TransProb[forward]);
	//  RealType tail_backward=nextlastbead.Properties(W.reptile->TransProb[backward]);


	// RealType dS_head=branchEngine->symLinkActionBare(head_forward, head_backward, newhead_e, oldhead_e);
	// RealType dS_tail=branchEngine->symLinkActionBare(tail_forward, tail_backward, newtail_e, oldtail_e);
	RealType dS_head =
	  branchEngine->symLinkAction (head_forward, head_backward, newhead_e,
				       oldhead_e);
	RealType dS_tail =
	  branchEngine->symLinkAction (tail_forward, tail_backward, newtail_e,
				       oldtail_e);



	RealType dS_0 =
	  dS_head - dS_tail + (curhead.Properties (LOGPSI) +
			       lastbead.Properties (LOGPSI) - logpsi -
			       nextlastbead.Properties (LOGPSI));

	/// RealType dS_old = +(curhead.Properties(LOGPSI) + lastbead.Properties(LOGPSI) - logpsi - nextlastbead.Properties(LOGPSI))
	///      + curhead.Properties(W.reptile->Action[2]) + W.Properties(W.reptile->Action[2])
	///      + curhead.Properties(W.reptile->Action[forward]) + W.Properties(W.reptile->Action[backward])
	///      - (lastbead.Properties(W.reptile->Action[2]) + nextlastbead.Properties(W.reptile->Action[2]))
	///      - (lastbead.Properties(W.reptile->Action[forward]) + nextlastbead.Properties(W.reptile->Action[backward]));
	///acceptProb=std::exp(-dS_0 + (nextlastbead.Properties(W.reptile->TransProb[backward]) - curhead.Properties(W.reptile->TransProb[forward])));

	//  acceptProb=std::exp(-dS_0 + tail_backward - head_forward);
	// app_log()<<"logGf (calced) = "
	//  app_log()<<"dS_old="<<dS_old<< std::endl;
	//  app_log()<<"dS_head="<<dS_head<< std::endl;
	//  app_log()<<"dS_tail="<<dS_tail<< std::endl;
	//  app_log()<<"dS' = "<<dS_0<< std::endl;
	//        app_log()<<"W.Properties(LOCALENERGY)="<<W.Properties(LOCALENERGY)<< std::endl;

	//  app_log()<<"---------------\n";
	acceptProb = std::exp (-dS_0 + (nextlastbead.Properties (W.reptile->TransProb[backward]) - curhead.Properties (W.reptile->TransProb[forward])));	//tail_backward - head_forward);
	//     acceptProb=std::min(1.0,std::exp(-dS + -(curhead.Properties(LOGPSI) + lastbead.Properties(LOGPSI) - logpsi - nextlastbead.Properties(LOGPSI)) + tail_backward - head_forward));

      }
    else
      {
	// dS = curhead.Properties(W.reptile->Action[2]) + W.Properties(W.reptile->Action[2])
	//- (lastbead.Properties(W.reptile->Action[2]) + nextlastbead.Properties(W.reptile->Action[2]));
	RealType dS_head =
	  branchEngine->DMCLinkAction (eloc,
				       curhead.Properties (LOCALENERGY));
	RealType dS_tail =
	  branchEngine->DMCLinkAction (lastbead.Properties (LOCALENERGY),
				       nextlastbead.Properties (LOCALENERGY));
	//dS=branchEngine->DMCLinkAction(eloc,curhead.Properties(LOCALENERGY)) - branchEngine->DMCLinkAction(lastbead.Properties(LOCALENERGY),nextlastbead.Properties(LOCALENERGY));
	dS = dS_head - dS_tail;
	acceptProb = std::min ((RealType)1.0, std::exp (-dS));

      }

    //app_log()<<acceptProb<< std::endl;
//       app_log()<<"r2proposed.... = "<<r2proposed<< std::endl;
    if ((RandomGen () < acceptProb) || curhead.Age >= MaxAge)
      {

	r2accept = r2proposed;
	//    W.reptile->r2accept+=r2accept;
	MCWalkerConfiguration::Walker_t & overwriteWalker (W.reptile->
							   getNewHead ());
	if (curhead.Age >= MaxAge)
	  {
	    app_log () << "\tForce Acceptance...\n";
	    equilToDoSteps = equilSteps;
	  }
	W.saveWalker (overwriteWalker);
	overwriteWalker.Properties (LOCALENERGY) = eloc;
	overwriteWalker.Properties (W.reptile->Action[forward]) = 0;
	overwriteWalker.Properties (W.reptile->Action[backward]) =
	  W.Properties (W.reptile->Action[backward]);
	overwriteWalker.Properties (W.reptile->Action[2]) =
	  W.Properties (W.reptile->Action[2]);
	//overwriteWalker.Properties(W.reptile->TransProb[forward])=W.Properties(W.reptile->TransProb[forward]);
	W.reptile->saveTransProb (overwriteWalker, +1, 0);
	W.reptile->saveTransProb (overwriteWalker, -1, logGb);
	// overwriteWalker.Properties(W.reptile->TransProb[backward])=W.Properties(W.reptile->TransProb[backward]);
	overwriteWalker.resetProperty (logpsi, Psi.getPhase (), eloc);
	overwriteWalker.Properties (R2ACCEPTED) = r2accept;
	overwriteWalker.Properties (R2PROPOSED) = r2proposed;

	// lastbead.Properties(R2PROPOSED)=lastbead.Properties(R2ACCEPTED)=nextlastbead.Properties(R2PROPOSED);
	H.auxHevaluate (W, overwriteWalker,true,false); //evaluate properties but not collectables.
	H.saveProperty (overwriteWalker.getPropertyBase ());
	overwriteWalker.Age = 0;

	++nAccept;
      }
    else
      {
	//  app_log()<<"Reject\n";
	curhead.Properties (R2ACCEPTED) = 0;
	curhead.Properties (R2PROPOSED) = r2proposed;
	lastbead.Properties (R2ACCEPTED) = 0;
	//lastbead.Properties(R2PROPOSED)=nextlastbead.Properties(R2PROPOSED);
	//curhead.Properties(R2
	++nReject;
	H.rejectedMove (W, curhead);
	curhead.Age += 1;
	W.reptile->flip ();
	// app_log()<<"Reject\n";
	return;
      }
      W.loadWalker(centerbead,true);
      W.update(false);  //skip S(k) evaluation?  False
      H.auxHevaluate(W,centerbead,false, true); //collectables, but not properties
    
  }

  void RMCUpdateAllWithDrift::accumulate (WalkerIter_t it,
					  WalkerIter_t it_end)
  {
    if (vmcToDoSteps == 0 && equilToDoSteps == 0)
      Estimators->accumulate (W, it, it_end);
    else;

  }


}
