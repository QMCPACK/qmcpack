//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Raymond Clay III, j.k.rofling@gmail.com, Lawrence Livermore National Laboratory
//                    Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    


#include "QMCDrivers/RMC/RMCSingleOMP.h"
#include "QMCDrivers/RMC/RMCUpdatePbyP.h"
#include "QMCDrivers/RMC/RMCUpdateAll.h"
#include "QMCDrivers/DriftOperators.h"
#include "ParticleBase/ParticleAttribOps.h"
#include "OhmmsApp/RandomNumberControl.h"
#include "Message/OpenMP.h"
#include "Message/CommOperators.h"
#include "tau/profiler.h"
#include "Particle/Reptile.h"
#include "Utilities/RunTimeManager.h"
#if !defined(REMOVE_TRACEMANAGER)
#include "Estimators/TraceManager.h"
#else
typedef int TraceManager;
#endif


namespace qmcplusplus
{

/// Constructor.
  RMCSingleOMP::RMCSingleOMP (MCWalkerConfiguration & w,
			      TrialWaveFunction & psi, QMCHamiltonian & h,
			      HamiltonianPool & hpool,
			      WaveFunctionPool & ppool):QMCDriver (w, psi, h,
								   ppool),
    CloneManager (hpool), prestepsVMC (-1), rescaleDrift ("no"), beta (-1),
    beads (-1), fromScratch (true)
  {
    RootName = "rmc";
    QMCType = "RMCSingleOMP";
    QMCDriverMode.set (QMC_UPDATE_MODE, 1);
    QMCDriverMode.set (QMC_WARMUP, 0);
    m_param.add (rescaleDrift, "drift", "string");
    m_param.add (beta, "beta", "double");
    m_param.add (beads, "beads", "int");
    m_param.add (resizeReptile, "resize", "int");
    m_param.add (prestepsVMC, "vmcpresteps", "int");

    Action.resize (3);
    Action[0] = w.addProperty ("ActionBackward");
    Action[1] = w.addProperty ("ActionForward");
    Action[2] = w.addProperty ("ActionLocal");
    TransProb.resize (2);
    TransProb[0] = w.addProperty ("TransProbBackward");
    TransProb[1] = w.addProperty ("TransProbForward");
  }

  bool RMCSingleOMP::run ()
  {
    resetRun ();
    //start the main estimator
    Estimators->start (nBlocks);
    for (int ip = 0; ip < NumThreads; ++ip)
      Movers[ip]->startRun (nBlocks, false);
#if !defined(REMOVE_TRACEMANAGER)
    Traces->startRun (nBlocks, traceClones);
#endif
    const bool has_collectables = W.Collectables.size ();

    LoopTimer rmc_loop;
    RunTimeControl runtimeControl(RunTimeManager, MaxCPUSecs);
    for (int block = 0; block < nBlocks; ++block)
    {
      rmc_loop.start();
#pragma omp parallel
	{
	  int ip = omp_get_thread_num ();
	  IndexType updatePeriod =
	    (QMCDriverMode[QMC_UPDATE_MODE]) ? Period4CheckProperties : 0;
	  //assign the iterators and resuse them
	  MCWalkerConfiguration::iterator wit (W.begin () + wPerNode[ip]),
	    wit_end (W.begin () + wPerNode[ip + 1]);
	  Movers[ip]->startBlock (nSteps);
	  int now_loc = CurrentStep;

	  RealType cnorm = 1.0; //This is because there is only one reptile per walkerset.

	  for (int step = 0; step < nSteps; ++step)
	    {
	      //collectables are reset, it is accumulated while advancing walkers
	      wClones[ip]->resetCollectables ();
	      Movers[ip]->advanceWalkers (wit, wit_end, false);
	      if (has_collectables)
		wClones[ip]->Collectables *= cnorm;
	      Movers[ip]->accumulate (wit, wit_end);

	      ++now_loc;
	      if (Period4WalkerDump && now_loc % myPeriod4WalkerDump == 0)
		wClones[ip]->saveEnsemble (wit, wit_end);

	      branchEngine->collect (CurrentStep, W, branchClones);	//Ray Clay:  For now, collects and syncs based on first reptile.  Need a better way to do this.
	    }
	  Movers[ip]->stopBlock (false);
	}			//end-of-parallel for
	CurrentStep += nSteps;
	Estimators->stopBlock (estimatorClones);
	//why was this commented out? Are checkpoints stored some other way?
	if (storeConfigs)
	  recordBlock (block);

        rmc_loop.stop();
        bool enough_time_for_next_iteration = runtimeControl.enough_time_for_next_iteration(rmc_loop);
        // Rank 0 decides whether the time limit was reached
        myComm->bcast(enough_time_for_next_iteration);

        if (!enough_time_for_next_iteration)
        {
          app_log() << runtimeControl.time_limit_message("RMC", block);
          break;
        }
      }				//block
    Estimators->stop (estimatorClones);
    //copy back the random states
    for (int ip = 0; ip < NumThreads; ++ip)
      *(RandomNumberControl::Children[ip]) = *(Rng[ip]);
    //return nbeads and stuff to its orginal unset state;
    resetVars ();
    return finalize (nBlocks);
  }

  void RMCSingleOMP::resetRun ()
  {
    m_param.put (qmcNode);
    //For now, assume that nReptiles=NumThreads;
    nReptiles = NumThreads;

    if (beads < 1)
      beads = beta / Tau;
    else
      beta = beads * Tau;

    app_log () << "Projection time:  " << beta << " Ha^-1" << std::endl;
    //Calculate the number of VMC presteps if not given:
    if (prestepsVMC == -1 && fromScratch == true)
      prestepsVMC = beads + 2;
    //Check to see if the MCWalkerConfiguration is in a state suitable for reptation
    if (!W.ReptileList.empty ())
      {
	fromScratch = false;

	app_log () << "Previous RMC reptiles detected...\n";
	if (Tau == W.ReptileList[0]->getTau ()
	    && beads == W.ReptileList[0]->size ())
	  app_log () << "  Using current reptiles\n";	//do nothing
	else			//we need to extrapolate off of the current reptile set.  
	  {
	    //pull the reptile configurations out
	    app_log () << "  Previous Tau/Beta:  " << W.ReptileList[0]->
	      getTau () << "/" << W.ReptileList[0]->getTau () *
	      W.ReptileList[0]->size () << std::endl;
	    app_log () << "  New      Tau/Beta: " << Tau << "/" << beta <<
	      std::endl;
	    app_log () << "    Linear interpolation to get new reptile.\n";
	    std::vector< ReptileConfig_t > repSamps (0);
	    for (IndexType sampid = 0;
		 sampid < W.ReptileList.size () && sampid < nReptiles;
		 sampid++)
	      repSamps.push_back (W.ReptileList[sampid]->
				  getReptileSlicePositions (Tau, beta));

	    //In the event of a disparity in the number of requested reptiles and the ones received....  just copy
	    //Copies cyclically.  First iteration copies the first entry, second the second, and so on.  So we don't replicate just one config.
	    for (IndexType copyid = 0; repSamps.size () < nReptiles; copyid++)
	      repSamps.push_back (repSamps[copyid]);


	    resetReptiles (repSamps, Tau);
	  }
      }

    //Previous run was nothing, VMC, or DMC.  No reptiles--so we initialize based on whatever is there.
    else
      {
	//Initialize on whatever walkers are in MCWalkerConfiguration.
	app_log () << "Using walkers from previous non-RMC run.\n";
	std::vector<ParticlePos_t> wSamps (0);
	MCWalkerConfiguration::iterator wit (W.begin ()), wend (W.end ());
	for (IndexType sampid = 0; wit != wend && sampid < nReptiles; wit++)
	  wSamps.push_back ((**wit).R);

	for (IndexType copyid = 0; wSamps.size () < nReptiles; copyid++)
	  wSamps.push_back (wSamps[copyid]);
	resetReptiles (wSamps, beads, Tau);
      }

    //Now that we know if we're starting from scratch... decide whether to force VMC warmup.
    if (prestepsVMC == -1 && fromScratch == true)
      prestepsVMC = beads + 2;
    makeClones (W, Psi, H);
    myPeriod4WalkerDump =
      (Period4WalkerDump > 0) ? Period4WalkerDump : (nBlocks + 1) * nSteps;

    if (Movers.empty ())
      {
	Movers.resize (NumThreads, 0);
	branchClones.resize (NumThreads, 0);
	estimatorClones.resize (NumThreads, 0);
	traceClones.resize (NumThreads, 0);
	Rng.resize (NumThreads, 0);
	branchEngine->initReptile (W);
#pragma omp parallel for
	for (int ip = 0; ip < NumThreads; ++ip)
	  {
	    std::ostringstream os;
	    estimatorClones[ip] = new EstimatorManagerBase (*Estimators);	//,*hClones[ip]);
	    estimatorClones[ip]->resetTargetParticleSet (*wClones[ip]);
	    estimatorClones[ip]->setCollectionMode (false);
	    Rng[ip] =
	      new RandomGenerator_t (*(RandomNumberControl::Children[ip]));
#if !defined(REMOVE_TRACEMANAGER)
	    traceClones[ip] = Traces->makeClone ();
#endif
	    hClones[ip]->setRandomGenerator (Rng[ip]);
	    branchClones[ip] = new BranchEngineType (*branchEngine);
	    if (QMCDriverMode[QMC_UPDATE_MODE])
	      {
		os <<
		  "  PbyP moves with drift, using RMCUpdatePbyPWithDriftFast"
		  << std::endl;
		Movers[ip] =
		  new RMCUpdatePbyPWithDrift (*wClones[ip], *psiClones[ip],
					      *hClones[ip], *Rng[ip], Action,
					      TransProb);

	      }
	    else
	      {
		os <<
		  "  walker moves with drift, using RMCUpdateAllWithDriftFast"
		  << std::endl;
		Movers[ip] =
		  new RMCUpdateAllWithDrift (*wClones[ip], *psiClones[ip],
					     *hClones[ip], *Rng[ip], Action,
					     TransProb);
	      }
	    Movers[ip]->nSubSteps = nSubSteps;
	    if (ip == 0)
	      app_log () << os.str () << std::endl;
	  }
      }
#if !defined(REMOVE_TRACEMANAGER)
    else
      {
#pragma omp parallel for
	for (int ip = 0; ip < NumThreads; ++ip)
	  {
	    traceClones[ip]->transfer_state_from (*Traces);
	  }
      }
#endif
    app_log ().flush ();
#pragma omp parallel for
    for (int ip = 0; ip < NumThreads; ++ip)
      {
	Movers[ip]->put (qmcNode);
	Movers[ip]->resetRun (branchClones[ip], estimatorClones[ip],
			      traceClones[ip]);
	// wClones[ip]->reptile = new Reptile(*wClones[ip], W.begin()+wPerNode[ip],W.begin()+wPerNode[ip+1]);
	wClones[ip]->reptile = W.ReptileList[ip];
	//app_log()<<"Thread # "<<ip<< std::endl;
	// printf(" Thread# %d  WalkerList.size()=%d \n",ip,wClones[ip]->WalkerList.size());

	// wClones[ip]->reptile->printState();
	wClones[ip]->activeBead = 0;
	wClones[ip]->direction = +1;

	if (QMCDriverMode[QMC_UPDATE_MODE])
	  {
	   // app_log () << ip << " initWalkers for pbyp...\n";
	    Movers[ip]->initWalkersForPbyP (W.ReptileList[ip]->repstart,
					    W.ReptileList[ip]->repend);
	  }
	else
	  {
	    Movers[ip]->initWalkers (W.begin () + wPerNode[ip],
				     W.begin () + wPerNode[ip + 1]);
	  }

	//this will "unroll" the reptile according to forced VMC steps (no bounce).  See beginning of function for logic of setting prestepVMC.
	for (IndexType prestep = 0; prestep < prestepsVMC; prestep++)
	  {
	    Movers[ip]->advanceWalkers (W.begin (), W.begin (), true);
	  }

	//set up initial action and transprob.
	MCWalkerConfiguration::iterator wit (W.begin () + wPerNode[ip]),
	  wit_end (W.begin () + wPerNode[ip + 1]);
      }


    app_log()<<"Finished "<<prestepsVMC<<" VMC presteps\n";
    branchEngine->checkParameters (W);

#pragma omp parallel for
    for (int ip = 0; ip < NumThreads; ++ip)
      {
	for (int prestep = 0; prestep < nWarmupSteps; ++prestep)
	  {
	    Movers[ip]->advanceWalkers (W.begin () + wPerNode[ip],
					W.begin () + wPerNode[ip + 1], false);
	    branchEngine->collect (CurrentStep, W, branchClones);
	  }
	Movers[ip]->updateWalkers (W.begin () + wPerNode[ip],
				   W.begin () + wPerNode[ip + 1]);

      }

    fromScratch = false;
  }

  bool RMCSingleOMP::put (xmlNodePtr q)
  {
    m_param.put (q);
    return true;
  }

  //This will resize the MCWalkerConfiguration and initialize the ReptileList.  Does not care for previous runs.  
  void RMCSingleOMP::resetReptiles (int nReptiles_in, int nbeads_in,
				    RealType tau)
  {
    for (MCWalkerConfiguration::ReptileList_t::iterator it =
	 W.ReptileList.begin (); it != W.ReptileList.end (); it++)
      delete *it;
    W.ReptileList.clear ();
    // Maybe we should be more vigorous in cleaning the MCWC WalkerList?
    std::vector<int> repWalkerSlice;
    int nwtot = nbeads_in * nReptiles_in;
    FairDivideLow (nwtot, nReptiles_in, repWalkerSlice);
    if (W.getActiveWalkers () - nwtot != 0)
      addWalkers (nwtot - W.getActiveWalkers ());

    for (int i = 0; i < nReptiles_in; i++)
      {
	W.ReptileList.
	  push_back (new
		     Reptile (W, W.begin () + repWalkerSlice[i],
			      W.begin () + repWalkerSlice[i + 1]));
	W.ReptileList[i]->setTau (tau);
      }


  }
  //This will resize the MCWalkerConfiguration and initialize Reptile list.  It will then reinitialize the MCWC with a list of Reptile coordinates
  void RMCSingleOMP::resetReptiles (std::vector< ReptileConfig_t > &reptile_samps,
				    RealType tau)
  {
    if (reptile_samps.empty ())
      {
	APP_ABORT
	  ("RMCSingleOMP::resetReptiles(std::vector< ReptileConfig_t > reptile_samps):  No samples!\n");
      }
    else
      {
	IndexType nReptiles_in = reptile_samps.size ();
	IndexType nBeads_in = reptile_samps[0].size ();
	resetReptiles (nReptiles_in, nBeads_in, tau);

	for (IndexType i = 0; i < W.ReptileList.size (); i++)
	  {
	    W.ReptileList[i]->setReptileSlicePositions (reptile_samps[i]);
	  }
      }
  }
  //For # of walker samples, create that many reptiles with nbeads each.  Initialize each reptile to have the value of the walker "seed".
  void RMCSingleOMP::resetReptiles (std::vector< ParticlePos_t > &walker_samps,
				    int nBeads_in, RealType tau)
  {
    if (walker_samps.empty ())
      {
	APP_ABORT
	  ("RMCSingleOMP::resetReptiles(std::vector< ParticlePos_t > walker_samps):  No samples!\n");
      }
    else
      {
	IndexType nReptiles_in = walker_samps.size ();
	resetReptiles (nReptiles_in, nBeads_in, tau);

	for (IndexType i = 0; i < W.ReptileList.size (); i++)
	  {
	    W.ReptileList[i]->setReptileSlicePositions (walker_samps[i]);
	  }
      }

  }

};
