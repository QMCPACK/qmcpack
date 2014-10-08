#include "QMCDrivers/RMC/RMCUpdatePbyP.h"
#include "QMCDrivers/DriftOperators.h"
#include "Message/OpenMP.h"
#include "Configuration.h"
#include "Particle/Reptile.h"
#include <cmath>
//////////////////////////////////////////////////////////////////////////
//  This driver is strongly based off the method used by Lucas Wagner in QWalk.
//
//  This driver proposes an "all-electron" configuration by performing N single particle moves,
//  accepting or rejecting at each step.  Umrigar's scaled drift is used to generate the VMC configurations.
//
//  The configuration is then accepted/rejected based only on the symmetrized DMC action, which
//  amounts to Maroni/Baroni's method.  Note that energy filtering (like in DMC) is used here as well.
//
//  Until serious discrepencies are detected, it is strongly advised that this driver be used over the
//  RMCUpdateAll, as the time step errors seem to be really reduced using this method.
//////////////////////////////////////////////////////////////////////////////


namespace qmcplusplus
{


/*  void add_rmc_timers(vector<NewTimer*>& timers)
  {
    timers.push_back(new NewTimer("RMCUpdatePbyP::advance")); //timer for the walker loop
    timers.push_back(new NewTimer("RMCUpdatePbyP::movePbyP")); //timer for MC, ratio etc
   // timers.push_back(new NewTimer("RMCUpdatePbyP::updateMBO")); //timer for measurements
   // timers.push_back(new NewTimer("RMCUpdatePbyP::energy")); //timer for measurements
    for (int i=0; i<timers.size(); ++i) TimerManager.addTimer(timers[i]);
  }*/
/// Constructor.
  RMCUpdatePbyPWithDrift::RMCUpdatePbyPWithDrift (MCWalkerConfiguration & w,
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
    //add_rmc_timers(myTimers);
    scaleDrift = false;
    actionType = SYM_ACTION;
    myTimers.push_back (new NewTimer ("RMCUpdatePbyP::advance"));	//timer for the walker loop
    myTimers.push_back (new NewTimer ("RMCUpdatePbyP::movePbyP"));	//timer for MC, ratio etc
    myTimers.push_back (new NewTimer ("RMCUpdatePbyP::updateMBO"));	//timer for measurements
    myTimers.push_back (new NewTimer ("RMCUpdatePbyP::energy"));	//timer for measurements
    TimerManager.addTimer (myTimers[0]);
    TimerManager.addTimer (myTimers[1]);
    TimerManager.addTimer (myTimers[2]);
    TimerManager.addTimer (myTimers[3]);
  }

  RMCUpdatePbyPWithDrift::~RMCUpdatePbyPWithDrift ()
  {
  }


  void RMCUpdatePbyPWithDrift::initWalkersForPbyP (WalkerIter_t it,
						   WalkerIter_t it_end)
  {
    UpdatePbyP = true;

    for (; it != it_end; ++it)
      {
	Walker_t & awalker = **it;	//W.reptile->getHead();
	W.R = awalker.R;
	W.update ();
	//W.loadWalker(awalker,UpdatePbyP);
	if (awalker.DataSet.size ())
	  awalker.DataSet.clear ();
	awalker.DataSet.rewind ();
	RealType logpsi = Psi.registerData (W, awalker.DataSet);
	RealType logpsi2 = Psi.updateBuffer (W, awalker.DataSet, false);
	awalker.G = W.G;
	awalker.L = W.L;
	//  268   W.saveWalker(awalker);
	RealType eloc = H.evaluate (W);
	//   BadState |= isnan(eloc);
	//thisWalker.resetProperty(std::log(abs(psi)), psi,eloc);
	awalker.resetProperty (logpsi, Psi.getPhase (), eloc);
//      randomize(awalker);
      }


    //  IndexType initsteps = W.reptile->nbeads + 10;

//  for (int n=0; n < initsteps; n++) advanceWalkersVMC();
  }

  bool RMCUpdatePbyPWithDrift::put (xmlNodePtr cur)
  {

    QMCUpdateBase::put (cur);

    ParameterSet m_param;
    bool usedrift = true;
    string action = "SLA";
    m_param.add (usedrift, "useDrift", "bool");
    m_param.add (action, "Action", "string");
    m_param.add (equilSteps, "equilsteps", "int");
    m_param.add (equilSteps, "equilSteps", "int");
    m_param.put (cur);

    if (usedrift == true)
      {
	app_log () << "  Using Umrigar scaled drift\n";
      }
    else
      {
	app_log () << "  Using non-scaled drift\n";
      }

    if (action == "DMC")
      {
	actionType = DMC_ACTION;
	app_log () << "  Using DMC link-action\n";
      }
    else
      {
	app_log () << "  Using Symmetrized Link-Action\n";
      }

    return true;
  }
  void RMCUpdatePbyPWithDrift::advanceWalkersVMC ()
  {
    //  double starttime = cpu_clock ();
    myTimers[0]->start ();
    // app_log () << "advanceWalkersVMC()::CALLED.. " << cpu_clock () -
    //   starttime << endl;
    IndexType direction = W.reptile->direction;
    IndexType forward = (1 - direction) / 2;
    IndexType backward = (1 + direction) / 2;
    Walker_t & curhead = W.reptile->getHead ();
    Walker_t prophead (curhead);
    Walker_t::Buffer_t & w_buffer (prophead.DataSet);
    W.loadWalker (prophead, true);
    W.R = prophead.R;
    //app_log () << "advanceWalkersVMC()::initialized variables... " <<
    //   cpu_clock () - starttime << endl;
    //  starttime = cpu_clock ();
//      W.update();
    //W.loadWalker(awalker,UpdatePbyP);
//      if (prophead.DataSet.size())
//        prophead.DataSet.clear();
//      prophead.DataSet.rewind();
//      RealType logpsi=Psi.registerData(W,prophead.DataSet);
//      RealType logpsi2=Psi.updateBuffer(W,prophead.DataSet,false);
    // curhead.G=W.G;
    //  curhead.L=W.L;
    // Walker_t& thisWalker(**it);
    //Psi.copyFromBuffer (W, w_buffer);
    Psi.copyFromBuffer (W, w_buffer);
    //create a 3N-Dimensional Gaussian with variance=1
    makeGaussRandomWithEngine (deltaR, RandomGen);
    int nAcceptTemp (0);
    int nRejectTemp (0);
    //copy the old energy and scale factor of drift
    RealType eold (prophead.Properties (LOCALENERGY));
    RealType vqold (prophead.Properties (DRIFTSCALE));
    RealType enew (eold);
    RealType rr_proposed = 0.0;
    RealType rr_accepted = 0.0;
    RealType gf_acc = 1.0;
    myTimers[1]->start ();
    for (int ig = 0; ig < W.groups (); ++ig)	//loop over species
      {
	RealType tauovermass = Tau * MassInvS[ig];
	RealType oneover2tau = 0.5 / (tauovermass);
	RealType sqrttau = std::sqrt (tauovermass);
	for (int iat = W.first (ig); iat < W.last (ig); ++iat)
	  {
	    //get the displacement
	    GradType grad_iat = Psi.evalGrad (W, iat);
	    PosType dr;
	    getScaledDrift (tauovermass, grad_iat, dr);
	    dr += sqrttau * deltaR[iat];
	    //RealType rr=dot(dr,dr);
	    RealType rr = tauovermass * dot (deltaR[iat], deltaR[iat]);
	    rr_proposed += rr;
	    if (rr > m_r2max)
	      {
		++nRejectTemp;
		continue;
	      }
	    //PosType newpos(W.makeMove(iat,dr));
	    if (!W.makeMoveAndCheck (iat, dr))
	      continue;
	    PosType newpos (W.R[iat]);
	    RealType ratio = Psi.ratioGrad (W, iat, grad_iat);
	    bool valid_move = false;
	    //node is crossed reject the move
	    //if(Psi.getPhase() > numeric_limits<RealType>::epsilon())
	    //if(branchEngine->phaseChanged(Psi.getPhase(),thisWalker.Properties(SIGN)))
	    if (branchEngine->phaseChanged (Psi.getPhaseDiff ()))
	      {
		++nRejectTemp;
		++nNodeCrossing;
		W.rejectMove (iat);
		Psi.rejectMove (iat);
	      }
	    else
	      {
		RealType logGf = -0.5 * dot (deltaR[iat], deltaR[iat]);
		//Use the force of the particle iat
		//RealType scale=getDriftScale(m_tauovermass,grad_iat);
		//dr = thisWalker.R[iat]-newpos-scale*real(grad_iat);
		getScaledDrift (tauovermass, grad_iat, dr);
		dr = prophead.R[iat] - newpos - dr;
		RealType logGb = -oneover2tau * dot (dr, dr);
		RealType prob = ratio * ratio * std::exp (logGb - logGf);
		if (RandomGen () < prob)
		  {
		    valid_move = true;
		    ++nAcceptTemp;
		    W.acceptMove (iat);
		    Psi.acceptMove (W, iat);
		    rr_accepted += rr;
		    gf_acc *= prob;	//accumulate the ratio
		  }
		else
		  {
		    ++nRejectTemp;
		    W.rejectMove (iat);
		    Psi.rejectMove (iat);
		  }
	      }
	  }
      }
    myTimers[1]->stop ();
    //  if(UseTMove)
    //    nonLocalOps.reset();
    bool advanced = true;

    if (nAcceptTemp > 0)
      {
	//need to overwrite the walker properties
	MCWalkerConfiguration::Walker_t & newhead (W.reptile->getNewHead ());
	myTimers[2]->start ();
	prophead.Age = 0;
	prophead.R = W.R;
	//w_buffer.rewind();
	//W.updateBuffer(w_buffer);
	RealType logpsi = Psi.updateBuffer (W, w_buffer, false);
	W.saveWalker (prophead);
	myTimers[2]->stop ();
	myTimers[3]->start ();
//      if(UseTMove)
//        enew= H.evaluate(W,nonLocalOps.Txy);
//      else
	enew = H.evaluate (W);
	myTimers[3]->stop ();
	//nodecorr=getNodeCorrection(W.G,thisWalker.Drift);
	//thisWalker.resetProperty(logpsi,Psi.getPhase(),enew,rr_accepted,rr_proposed,nodecorr);
	prophead.resetProperty (logpsi, Psi.getPhase (), enew, rr_accepted,
				rr_proposed, 0.0);
	prophead.Weight = 1.0;
	H.auxHevaluate (W, prophead);
	H.saveProperty (prophead.getPropertyBase ());
	newhead = prophead;
	nAccept++;
      }
    else
      {
	//all moves are rejected: does not happen normally with reasonable wavefunctions
	advanced = false;
	curhead.Age++;
	curhead.Properties (R2ACCEPTED) = 0.0;
	//weight is set to 0 for traces
	// consistent w/ no evaluate/auxHevaluate
	RealType wtmp = prophead.Weight;
	curhead.Weight = 0.0;
	H.rejectedMove (W, curhead);
	curhead.Weight = wtmp;
	++nAllRejected;
	gf_acc = 1.0;
	nReject++;
      }
    // Traces->buffer_sample();
  }


  void RMCUpdatePbyPWithDrift::initWalkers (WalkerIter_t it,
					    WalkerIter_t it_end)
  {
    IndexType initsteps = W.reptile->nbeads * 2;

    vmcSteps = W.reptile->nbeads + 1;

    for (int n = 0; n < initsteps; n++)
      advanceWalkersVMC ();
  }

  void RMCUpdatePbyPWithDrift::advanceWalkersRMC ()
  {
    IndexType direction = W.reptile->direction;
    IndexType forward = (1 - direction) / 2;
    IndexType backward = (1 + direction) / 2;
    Walker_t & curhead = W.reptile->getHead ();
    Walker_t prophead (curhead);
    Walker_t::Buffer_t & w_buffer (prophead.DataSet);
    W.loadWalker (prophead, true);


    Psi.copyFromBuffer (W, w_buffer);

    makeGaussRandomWithEngine (deltaR, RandomGen);
    int nAcceptTemp (0);
    int nRejectTemp (0);
    //copy the old energy and scale factor of drift
    RealType eold (prophead.Properties (LOCALENERGY));
    RealType vqold (prophead.Properties (DRIFTSCALE));
    RealType enew (eold);
    RealType rr_proposed = 0.0;
    RealType rr_accepted = 0.0;
    RealType gf_acc = 1.0;
    myTimers[1]->start ();
    for (int ig = 0; ig < W.groups (); ++ig)	//loop over species
      {
	RealType tauovermass = Tau * MassInvS[ig];
	RealType oneover2tau = 0.5 / (tauovermass);
	RealType sqrttau = std::sqrt (tauovermass);
	for (int iat = W.first (ig); iat < W.last (ig); ++iat)
	  {
	    //get the displacement
	    GradType grad_iat = Psi.evalGrad (W, iat);
	    PosType dr;
	    getScaledDrift (tauovermass, grad_iat, dr);
	    dr += sqrttau * deltaR[iat];
	    //RealType rr=dot(dr,dr);
	    RealType rr = tauovermass * dot (deltaR[iat], deltaR[iat]);
	    rr_proposed += rr;
	    if (rr > m_r2max)
	      {
		++nRejectTemp;
		continue;
	      }
	    //PosType newpos(W.makeMove(iat,dr));
	    if (!W.makeMoveAndCheck (iat, dr))
	      continue;
	    PosType newpos (W.R[iat]);
	    RealType ratio = Psi.ratioGrad (W, iat, grad_iat);
	    bool valid_move = false;
	    //node is crossed reject the move
	    //if(Psi.getPhase() > numeric_limits<RealType>::epsilon())
	    //if(branchEngine->phaseChanged(Psi.getPhase(),thisWalker.Properties(SIGN)))
	    if (branchEngine->phaseChanged (Psi.getPhaseDiff ()))
	      {
		++nRejectTemp;
		++nNodeCrossing;
		W.rejectMove (iat);
		Psi.rejectMove (iat);
	      }
	    else
	      {
		RealType logGf = -0.5 * dot (deltaR[iat], deltaR[iat]);
		//Use the force of the particle iat
		//RealType scale=getDriftScale(m_tauovermass,grad_iat);
		//dr = thisWalker.R[iat]-newpos-scale*real(grad_iat);
		getScaledDrift (tauovermass, grad_iat, dr);
		dr = prophead.R[iat] - newpos - dr;
		RealType logGb = -oneover2tau * dot (dr, dr);
		RealType prob = ratio * ratio * std::exp (logGb - logGf);
		if (RandomGen () < prob)
		  {
		    valid_move = true;
		    ++nAcceptTemp;
		    W.acceptMove (iat);
		    Psi.acceptMove (W, iat);
		    rr_accepted += rr;
		    gf_acc *= prob;	//accumulate the ratio
		  }
		else
		  {
		    ++nRejectTemp;
		    W.rejectMove (iat);
		    Psi.rejectMove (iat);
		  }
	      }
	  }
      }
    myTimers[1]->stop ();
    //  if(UseTMove)
/*
  RealType logpsiold = prophead.Properties(LOGPSI);
  makeGaussRandomWithEngine(deltaR,RandomGen);
  int nAcceptTemp(0);
  int nRejectTemp(0);
  //copy the old energy and scale factor of drift
  RealType eold(prophead.Properties(LOCALENERGY));
  RealType vqold(prophead.Properties(DRIFTSCALE));
  RealType enew(eold);
  RealType rr_proposed=0.0;
  RealType rr_accepted=0.0;
  RealType gf_acc=1.0;
  RealType logGf=0;
  RealType logGb=0;
  RealType ratio_acc=1;
  //   myTimers[1]->start();
  for(int iat=0; iat<NumPtcl; ++iat)
    {
      PosType dr;
      //get the displacement
      //RealType sc=getDriftScale(m_tauovermass,W.G[iat]);
      //PosType dr(m_sqrttau*deltaR[iat]+sc*real(W.G[iat]));
      getScaledDrift(m_tauovermass,W.G[iat],dr);
      dr += m_sqrttau*deltaR[iat];
      //RealType rr=dot(dr,dr);
      RealType rr=m_tauovermass*dot(deltaR[iat],deltaR[iat]);
      rr_proposed+=rr;
      // if (W.reptile->r2
      if(rr>m_r2max)
        {
          ++nRejectTemp;
         // W.reptile->flip();
          continue;
        }
      //PosType newpos(W.makeMove(iat,dr));
      if(!W.makeMoveAndCheck(iat,dr))
        {
          ++nRejectTemp;
          continue;
          // W.reptile->flip();
          // return;
        }
      PosType newpos(W.R[iat]);
      RealType ratio=Psi.ratio(W,iat,dG,dL);
      bool valid_move=false;
      //node is crossed reject the move
      //if(Psi.getPhase() > numeric_limits<RealType>::epsilon())
      if(branchEngine->phaseChanged(Psi.getPhaseDiff()))
        {
          ++nRejectTemp;
          ++nNodeCrossing;
          W.rejectMove(iat);
          Psi.rejectMove(iat);
        }
      else
        {
          G = W.G+dG;
          RealType logGf_pbyp = -0.5*dot(deltaR[iat],deltaR[iat]);
          logGf+=logGf_pbyp;
          //Use the force of the particle iat
          //RealType scale=getDriftScale(m_tauovermass,G[iat]);
          //dr = thisWalker.R[iat]-newpos-scale*real(G[iat]);
          getScaledDrift(m_tauovermass, G[iat], dr);
          dr = prophead.R[iat] - newpos - dr;
          RealType logGb_pbyp = -m_oneover2tau*dot(dr,dr);
          logGb+=logGb_pbyp;
          RealType prob = ratio*ratio*std::exp(logGb_pbyp-logGf_pbyp);
          //this is useless
          //RealType prob = std::min(1.0,ratio*ratio*std::exp(logGb-logGf));
          if(RandomGen() < prob)
            {
              valid_move=true;
              ++nAcceptTemp;
              W.acceptMove(iat);
              Psi.acceptMove(W,iat);
              W.G = G;
              W.L += dL;
              rr_accepted+=rr;
              gf_acc *=prob;//accumulate the ratio
              ratio_acc*=ratio;
            }
          else
            {
              ++nRejectTemp;
              W.rejectMove(iat);
              Psi.rejectMove(iat);
              // W.reptile->flip();
              // return;
            }
        }
    }*/
// In the rare case that all proposed moves fail, we bounce.
    if (nAcceptTemp == 0)
      {

	++nReject;
	H.rejectedMove (W, prophead);
	curhead.Age += 1;
	W.reptile->flip ();
      }
    prophead.R = W.R;
    RealType logpsi = Psi.updateBuffer (W, w_buffer, false);
    W.saveWalker (prophead);
    Walker_t & lastbead (W.reptile->getTail ()),
      nextlastbead (W.reptile->getNext ());
    RealType eloc = H.evaluate (W);
    RealType dS =
      branchEngine->DMCLinkAction (eloc,
				   curhead.Properties (LOCALENERGY)) -
      branchEngine->DMCLinkAction (lastbead.Properties (LOCALENERGY),
				   nextlastbead.Properties (LOCALENERGY));
    RealType acceptProb = std::min (1.0, std::exp (-dS));
    if ((RandomGen () <= acceptProb)
	|| (prophead.Age >= MaxAge || lastbead.Age >= MaxAge))
      {

	MCWalkerConfiguration::Walker_t & overwriteWalker (W.reptile->
							   getNewHead ());
	if (curhead.Age >= MaxAge || lastbead.Age >= MaxAge)
	  app_log () << "\tForce Acceptance...\n";
	// RealType logpsi = Psi.updateBuffer(W,w_buffer,false);
	//W.saveWalker(prophead);

	prophead.Properties (LOCALENERGY) = eloc;
	prophead.Properties (R2ACCEPTED) = rr_accepted;
	prophead.Properties (R2PROPOSED) = rr_proposed;
	H.auxHevaluate (W, prophead);
	//H.saveProperty(overwriteWalker.getPropertyBase());
	H.saveProperty (prophead.getPropertyBase ());
	//overwriteWalker.Age=0;
	prophead.Age = 0;

	overwriteWalker = prophead;
	++nAccept;
      }
    else
      {
	++nReject;
	H.rejectedMove (W, prophead);
	curhead.Properties (R2ACCEPTED) = 0;
	curhead.Properties (R2PROPOSED) = rr_proposed;
	curhead.Age += 1;
	W.reptile->flip ();
	// return;
      }
  }

  void RMCUpdatePbyPWithDrift::advanceWalkers (WalkerIter_t it,
					       WalkerIter_t it_end, bool init)
  {
    if (init == true)
      advanceWalkersVMC ();
    else
      advanceWalkersRMC ();
  }

  void RMCUpdatePbyPWithDrift::accumulate (WalkerIter_t it,
					   WalkerIter_t it_end)
  {
    // if (vmcToDoSteps==0 && equilToDoSteps==0) Estimators->accumulate(W,it,it_end);
//  else;       
    Estimators->accumulate (W, it, it_end);
  }

}
