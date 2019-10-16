//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Luke Shulenburger, lshulen@sandia.gov, Sandia National Laboratories
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Raymond Clay III, j.k.rofling@gmail.com, Lawrence Livermore National Laboratory
//                    Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


#include "QMCDrivers/SimpleFixedNodeBranch.h"
#include "QMCDrivers/DMC/WalkerControlFactory.h"
#include <numeric>
#include "OhmmsData/FileUtility.h"
#include "Utilities/RandomGenerator.h"
#include "QMCDrivers/WalkerControlBase.h"
#include "Estimators/EstimatorManagerBase.h"
#include "QMCDrivers/BranchIO.h"
#include "Particle/Reptile.h"

namespace qmcplusplus
{
///enum to yes/no options saved in sParam
enum
{
  COMBOPT,
  USETAUOPT,
  MIXDMCOPT,
  DUMMYOPT
};

SimpleFixedNodeBranch::SimpleFixedNodeBranch(RealType tau, int nideal)
    : vParam(1.0), MyEstimator(0) //, PopHist(5), DMCEnergyHist(5)
{
  BranchMode.set(B_DMCSTAGE, 0);     //warmup stage
  BranchMode.set(B_POPCONTROL, 1);   //use standard DMC
  BranchMode.set(B_USETAUEFF, 1);    //use taueff
  BranchMode.set(B_CLEARHISTORY, 0); //clear history and start with the current average
  BranchMode.set(B_KILLNODES, 0);    //when killing walkers at nodes etrial is updated differently
  vParam[B_TAU]         = tau;
  vParam[B_TAUEFF]      = tau;
  vParam[B_FEEDBACK]    = 1.0;
  vParam[B_FILTERSCALE] = 10;
  R2Accepted(1.0e-10);
  R2Proposed(1.0e-10);
  //set the default values for integer parameters
  iParam[B_WARMUPSTEPS]          = 200;
  iParam[B_ENERGYUPDATEINTERVAL] = 1;
  iParam[B_BRANCHINTERVAL]       = 1;
  iParam[B_TARGETWALKERS]        = 0;
  iParam[B_MAXWALKERS]           = nideal;
  iParam[B_MINWALKERS]           = nideal;
  iParam[B_COUNTER]              = -1;
  //default is no
  sParam.resize(DUMMYOPT, "no");
  //default is classic
  branching_cutoff_scheme = "classic";
  registerParameters();
  reset();
}

/** copy constructor
 *
 * Copy only selected data members and WalkerController is never copied.
 */
SimpleFixedNodeBranch::SimpleFixedNodeBranch(const SimpleFixedNodeBranch& abranch)
    : BranchMode(abranch.BranchMode),
      iParam(abranch.iParam),
      vParam(abranch.vParam),
      MyEstimator(0),
      branching_cutoff_scheme(abranch.branching_cutoff_scheme),
      sParam(abranch.sParam)
{
  registerParameters();
  reset();
}

void SimpleFixedNodeBranch::registerParameters()
{
  m_param.add(iParam[B_WARMUPSTEPS], "warmupSteps", "int");
  m_param.add(iParam[B_WARMUPSTEPS], "warmupsteps", "int");
  m_param.add(iParam[B_ENERGYUPDATEINTERVAL], "energyUpdateInterval", "int");
  m_param.add(iParam[B_BRANCHINTERVAL], "branchInterval", "int");
  m_param.add(iParam[B_TARGETWALKERS], "targetWalkers", "int");
  m_param.add(iParam[B_TARGETWALKERS], "targetwalkers", "int");
  m_param.add(iParam[B_TARGETWALKERS], "target_walkers", "int");
  //trial energy
  m_param.add(vParam[B_EREF], "refEnergy", "AU");
  m_param.add(vParam[B_EREF], "ref_energy", "AU");
  m_param.add(vParam[B_EREF], "en_ref", "AU");
  m_param.add(vParam[B_TAU], "tau", "AU");
  m_param.add(vParam[B_TAU], "timestep", "AU");
  m_param.add(vParam[B_TAU], "timeStep", "AU");
  m_param.add(vParam[B_TAU], "TimeStep", "AU");
  //filterscale:  sets the filtercutoff to sigma*filterscale
  m_param.add(vParam[B_FILTERSCALE], "filterscale", "double");
  //feed back parameter for population control
  m_param.add(vParam[B_FEEDBACK], "feedback", "double");
  //turn on/off effective tau onl for time-step error comparisons
  m_param.add(sParam[USETAUOPT], "useBareTau", "option");
  m_param.add(sParam[MIXDMCOPT], "warmupByReconfiguration", "opt");
  m_param.add(branching_cutoff_scheme, "branching_cutoff_scheme", "option");
}

void SimpleFixedNodeBranch::start(const std::string& froot, bool append)
{
  RootName              = froot;
  MyEstimator->RootName = froot;
  MyEstimator->reset();
}

int SimpleFixedNodeBranch::initWalkerController(MCWalkerConfiguration& walkers, bool fixW, bool killwalker)
{
  BranchMode.set(B_DMC, 1);                               //set DMC
  BranchMode.set(B_DMCSTAGE, iParam[B_WARMUPSTEPS] == 0); //use warmup
  //this is not necessary
  //check if tau is different and set the initial values
  //vParam[B_TAU]=tau;
  bool fromscratch     = false;
  FullPrecRealType tau = vParam[B_TAU];

  int nwtot_now = walkers.getGlobalNumWalkers();

  //this is the first time DMC is used
  if (WalkerController == nullptr)
  {
    if (iParam[B_TARGETWALKERS] == 0)
    {
      Communicate* acomm = MyEstimator->getCommunicator();
      int ncontexts      = acomm->size();
      std::vector<int> nw(ncontexts, 0), nwoff(ncontexts + 1, 0);
      nw[acomm->rank()] = walkers.getActiveWalkers();
      acomm->allreduce(nw);
      for (int ip = 0; ip < ncontexts; ++ip)
        nwoff[ip + 1] = nwoff[ip] + nw[ip];
      walkers.setGlobalNumWalkers(nwoff[ncontexts]);
      walkers.setWalkerOffsets(nwoff);
      iParam[B_TARGETWALKERS] = nwoff[ncontexts];
    }
    WalkerController.reset(createWalkerController(iParam[B_TARGETWALKERS], MyEstimator->getCommunicator(), myNode));
    if (!BranchMode[B_RESTART])
    {
      fromscratch = true;
      app_log() << "  START ALL OVER " << std::endl;
      vParam[B_TAUEFF] = tau;
      BranchMode.set(B_POPCONTROL, !fixW); //fixW -> 0
      BranchMode.set(B_KILLNODES, killwalker);
      iParam[B_MAXWALKERS] = WalkerController->get_n_max();
      iParam[B_MINWALKERS] = WalkerController->get_n_min();
      if (!fixW && sParam[MIXDMCOPT] == "yes")
      {
        app_log() << "Warmup DMC is done with a fixed population " << iParam[B_TARGETWALKERS] << std::endl;
        BackupWalkerController = std::move(WalkerController); //save the main controller
        WalkerController.reset(
            createWalkerController(iParam[B_TARGETWALKERS], MyEstimator->getCommunicator(), myNode, true));
        BranchMode.set(B_POPCONTROL, 0);
      }
      //PopHist.clear();
      //PopHist.reserve(std::max(iParam[B_ENERGYUPDATEINTERVAL],5));
    }
    WalkerController->setWalkerID(walkers);
  }
  //else
  //{
  //  BranchMode.set(B_DMCSTAGE,0);//always reset warmup
  //}
  MyEstimator->reset();
  //update the simulation parameters
  WalkerController->put(myNode);
  //assign current Eref and a large number for variance
  WalkerController->setTrialEnergy(vParam[B_ETRIAL]);
  this->reset();
  if (fromscratch)
  {
    //determine the branch cutoff to limit wild weights based on the sigma and sigmaBound
    setBranchCutoff(vParam[B_SIGMA2], WalkerController->get_target_sigma(), 50, walkers.R.size());
    vParam[B_TAUEFF] = tau * R2Accepted.result() / R2Proposed.result();
  }
  //reset controller
  WalkerController->reset();
  if (BackupWalkerController)
    BackupWalkerController->reset();
  app_log() << "  QMC counter      = " << iParam[B_COUNTER] << std::endl;
  app_log() << "  time step        = " << vParam[B_TAU] << std::endl;
  app_log() << "  effective time step = " << vParam[B_TAUEFF] << std::endl;
  app_log() << "  trial energy     = " << vParam[B_ETRIAL] << std::endl;
  app_log() << "  reference energy = " << vParam[B_EREF] << std::endl;
  app_log() << "  Feedback = " << vParam[B_FEEDBACK] << std::endl;
  app_log() << "  reference variance = " << vParam[B_SIGMA2] << std::endl;
  app_log() << "  target walkers = " << iParam[B_TARGETWALKERS] << std::endl;
  app_log() << "  branching cutoff scheme " << branching_cutoff_scheme << std::endl;
  app_log() << "  branch cutoff = " << vParam[B_BRANCHCUTOFF] << " " << vParam[B_BRANCHMAX] << std::endl;
  app_log() << "  Max and minimum walkers per node= " << iParam[B_MAXWALKERS] << " " << iParam[B_MINWALKERS]
            << std::endl;
  app_log() << "  QMC Status (BranchMode) = " << BranchMode << std::endl;

  return int(round(double(iParam[B_TARGETWALKERS]) / double(nwtot_now)));
}

int SimpleFixedNodeBranch::initWalkerController(MCPopulation& population, bool fixW, bool killwalker)
{
  BranchMode.set(B_DMC, 1);                               //set DMC
  BranchMode.set(B_DMCSTAGE, iParam[B_WARMUPSTEPS] == 0); //use warmup
  //this is not necessary
  //check if tau is different and set the initial values
  //vParam[B_TAU]=tau;
  bool fromscratch     = false;
  FullPrecRealType tau = vParam[B_TAU];

  int nwtot_now = population.get_num_global_walkers();

  if (WalkerController != nullptr)
    throw std::runtime_error(
        "Unified Driver initWalkerController called with existing WalkerController,\n this is a violation of the assumed "
        "state module for SimpleFixedNodeBranch in the Unified Driver design");
  if (iParam[B_TARGETWALKERS] == 0)
  {
    // has "important" side effect of updating the walker offsets
    iParam[B_TARGETWALKERS] = population.update_num_global_walkers(MyEstimator->getCommunicator());
  }
  WalkerController.reset(createWalkerController(iParam[B_TARGETWALKERS], MyEstimator->getCommunicator(), myNode));
  if (!BranchMode[B_RESTART])
  {
    fromscratch = true;
    app_log() << "  START ALL OVER " << std::endl;
    vParam[B_TAUEFF] = tau;
    BranchMode.set(B_POPCONTROL, !fixW); //fixW -> 0
    BranchMode.set(B_KILLNODES, killwalker);
    iParam[B_MAXWALKERS] = WalkerController->get_n_max();
    iParam[B_MINWALKERS] = WalkerController->get_n_min();
    if (!fixW && sParam[MIXDMCOPT] == "yes")
    {
      app_log() << "Warmup DMC is done with a fixed population " << iParam[B_TARGETWALKERS] << std::endl;
      BackupWalkerController = std::move(WalkerController); //save the main controller
      WalkerController.reset(
          createWalkerController(iParam[B_TARGETWALKERS], MyEstimator->getCommunicator(), myNode, true));
      BranchMode.set(B_POPCONTROL, 0);
    }
    //PopHist.clear();
    //PopHist.reserve(std::max(iParam[B_ENERGYUPDATEINTERVAL],5));
  }

  // start used to be buried in the WalkerController initializing MCWC's walkers ID's
  WalkerController->start();

  //else
  //{
  //  BranchMode.set(B_DMCSTAGE,0);//always reset warmup
  //}

  // EstimatorManager gets touched from too many different places
  MyEstimator->reset();
  //update the simulation parameters
  WalkerController->put(myNode);
  //assign current Eref and a large number for variance
  WalkerController->setTrialEnergy(vParam[B_ETRIAL]);
  this->reset();
  if (fromscratch)
  {
    //determine the branch cutoff to limit wild weights based on the sigma and sigmaBound
    // Was check MCWC's particle set for number of R which I take to mean number of particles
    // will this assumption change if various spin freedoms also are added to ParticleSet?
    setBranchCutoff(vParam[B_SIGMA2], WalkerController->get_target_sigma(), 50, population.get_num_particles());
    vParam[B_TAUEFF] = tau * R2Accepted.result() / R2Proposed.result();
  }
  //reset controller
  WalkerController->reset();
  if (BackupWalkerController)
    BackupWalkerController->reset();
  app_log() << "  QMC counter      = " << iParam[B_COUNTER] << std::endl;
  app_log() << "  time step        = " << vParam[B_TAU] << std::endl;
  app_log() << "  effective time step = " << vParam[B_TAUEFF] << std::endl;
  app_log() << "  trial energy     = " << vParam[B_ETRIAL] << std::endl;
  app_log() << "  reference energy = " << vParam[B_EREF] << std::endl;
  app_log() << "  Feedback = " << vParam[B_FEEDBACK] << std::endl;
  app_log() << "  reference variance = " << vParam[B_SIGMA2] << std::endl;
  app_log() << "  target walkers = " << iParam[B_TARGETWALKERS] << std::endl;
  app_log() << "  branching cutoff scheme " << branching_cutoff_scheme << std::endl;
  app_log() << "  branch cutoff = " << vParam[B_BRANCHCUTOFF] << " " << vParam[B_BRANCHMAX] << std::endl;
  app_log() << "  Max and minimum walkers per node= " << iParam[B_MAXWALKERS] << " " << iParam[B_MINWALKERS]
            << std::endl;
  app_log() << "  QMC Status (BranchMode) = " << BranchMode << std::endl;

  return int(round(double(iParam[B_TARGETWALKERS]) / double(nwtot_now)));
}

void SimpleFixedNodeBranch::initReptile(MCWalkerConfiguration& W)
{
  RealType allowedFlux = 50.0;
  BranchMode.set(B_RMC, 1);                               //set RMC
  BranchMode.set(B_RMCSTAGE, iParam[B_WARMUPSTEPS] == 0); //use warmup
  //this is not necessary
  //check if tau is different and set the initial values
  //vParam[B_TAU]=tau;
  bool fromscratch     = false;
  FullPrecRealType tau = vParam[B_TAU];
  //this is the first time DMC is used
  if (WalkerController == 0)
  {
    //  if(iParam[B_TARGETWALKERS]==0)
    //  {
    //    Communicate* acomm=MyEstimator->getCommunicator();
    //    int ncontexts=acomm->size();
    //    std::vector<int> nw(ncontexts,0),nwoff(ncontexts+1,0);
    //    nw[acomm->rank()]=W.getActiveWalkers();
    //   acomm->allreduce(nw);
    //    for(int ip=0; ip<ncontexts; ++ip)
    //      nwoff[ip+1]=nwoff[ip]+nw[ip];
    //    W.setGlobalNumWalkers(nwoff[ncontexts]);
    //    W.setWalkerOffsets(nwoff);
    //    iParam[B_TARGETWALKERS]=nwoff[ncontexts];
    //  }
    if (!BranchMode[B_RESTART])
    {
      fromscratch = true;
      app_log() << "  START ALL OVER " << std::endl;
      vParam[B_TAUEFF] = tau;
      //PopHist.clear();
      //PopHist.reserve(std::max(iParam[B_ENERGYUPDATEINTERVAL],5));
    }
  }
  //else
  //{
  //  BranchMode.set(B_DMCSTAGE,0);//always reset warmup
  //}
  MyEstimator->reset();
  this->reset();
  if (fromscratch)
  {
    //determine the branch cutoff to limit wild weights based on the sigma and sigmaBound
    setBranchCutoff(vParam[B_SIGMA2], allowedFlux, 50, W.R.size());
    vParam[B_TAUEFF] = tau * R2Accepted.result() / R2Proposed.result();
  }
  //reset controller
  app_log() << "  QMC counter      = " << iParam[B_COUNTER] << std::endl;
  app_log() << "  time step        = " << vParam[B_TAU] << std::endl;
  app_log() << "  effective time step = " << vParam[B_TAUEFF] << std::endl;
  app_log() << "  reference energy = " << vParam[B_EREF] << std::endl;
  app_log() << "  Feedback = " << vParam[B_FEEDBACK] << std::endl;
  app_log() << "  reference variance = " << vParam[B_SIGMA2] << std::endl;
  app_log() << "  branching cutoff scheme " << branching_cutoff_scheme << std::endl;
  app_log() << "  branch cutoff = " << vParam[B_BRANCHCUTOFF] << " " << vParam[B_BRANCHMAX] << std::endl;
  app_log() << "  QMC Status (BranchMode) = " << BranchMode << std::endl;
}

void SimpleFixedNodeBranch::flush(int counter)
{
  if (counter == 0 && WalkerController)
    WalkerController->reset();
}

void SimpleFixedNodeBranch::branch(int iter, MCWalkerConfiguration& walkers)
{
  //collect the total weights and redistribute the walkers accordingly, using a fixed tolerance
  //RealType pop_now= WalkerController->branch(iter,walkers,0.1);
  FullPrecRealType pop_now;
  if (BranchMode[B_DMCSTAGE] || iter)
    pop_now = WalkerController->branch(iter, walkers, 0.1);
  else
    pop_now = WalkerController->doNotBranch(iter, walkers); //do not branch for the first step of a warmup
  //population for trial energy modification should not include any released node walkers.
  pop_now -= WalkerController->get_ensemble_property().RNSamples;
  //current energy
  vParam[B_ENOW] = WalkerController->get_ensemble_property().Energy;
  VarianceHist(WalkerController->get_ensemble_property().Variance);
  R2Accepted(WalkerController->get_ensemble_property().R2Accepted);
  R2Proposed(WalkerController->get_ensemble_property().R2Proposed);
  //PopHist(pop_now);
  vParam[B_EREF] = EnergyHist.mean(); //current mean
  if (BranchMode[B_USETAUEFF])
    vParam[B_TAUEFF] = vParam[B_TAU] * R2Accepted.result() / R2Proposed.result();
  if (BranchMode[B_KILLNODES])
    EnergyHist(vParam[B_ENOW] - std::log(WalkerController->get_ensemble_property().LivingFraction) / vParam[B_TAUEFF]);
  else
    EnergyHist(vParam[B_ENOW]);
  if (BranchMode[B_DMCSTAGE]) // main stage
  {
    if (BranchMode[B_POPCONTROL])
    {
      if (ToDoSteps > 0)
        --ToDoSteps;
      else
      {
        vParam[B_ETRIAL] = vParam[B_EREF] + vParam[B_FEEDBACK] * (logN - std::log(pop_now));
        ToDoSteps        = iParam[B_ENERGYUPDATEINTERVAL] - 1;
      }
    }
    else
      vParam[B_ETRIAL] = vParam[B_EREF];
  }
  else //warmup
  {
    if (BranchMode[B_USETAUEFF])
      vParam[B_TAUEFF] = vParam[B_TAU] * R2Accepted.result() / R2Proposed.result();
    if (BranchMode[B_POPCONTROL])
    {
      //RealType emix=((iParam[B_WARMUPSTEPS]-ToDoSteps)<100)?(0.25*vParam[B_EREF]+0.75*vParam[B_ENOW]):vParam[B_EREF];
      //vParam[B_ETRIAL]=emix+Feedback*(logN-std::log(pop_now));
      //vParam[B_ETRIAL]=vParam[B_EREF]+Feedback*(logN-std::log(pop_now));
      if (BranchMode[B_KILLNODES])
        vParam[B_ETRIAL] = (0.00 * vParam[B_EREF] + 1.0 * vParam[B_ENOW]) +
            vParam[B_FEEDBACK] * (logN - std::log(pop_now)) -
            std::log(WalkerController->get_ensemble_property().LivingFraction) / vParam[B_TAU];
      else
        vParam[B_ETRIAL] = vParam[B_ENOW] + (logN - std::log(pop_now)) / vParam[B_TAU];
    }
    else
      vParam[B_ETRIAL] = vParam[B_EREF];
    --ToDoSteps;
    if (ToDoSteps == 0) //warmup is done
    {
      vParam[B_SIGMA2] = VarianceHist.mean();
      setBranchCutoff(vParam[B_SIGMA2], WalkerController->get_target_sigma(), 10, walkers.R.size());
      app_log() << "\n Warmup is completed after " << iParam[B_WARMUPSTEPS] << std::endl;
      if (BranchMode[B_USETAUEFF])
        app_log() << "\n  TauEff     = " << vParam[B_TAUEFF] << "\n TauEff/Tau = " << vParam[B_TAUEFF] / vParam[B_TAU];
      else
        app_log() << "\n  TauEff proposed   = " << vParam[B_TAUEFF] * R2Accepted.result() / R2Proposed.result();
      app_log() << "\n  Etrial     = " << vParam[B_ETRIAL] << std::endl;
      app_log() << " Running average of energy = " << EnergyHist.mean() << std::endl;
      app_log() << "                  Variance = " << vParam[B_SIGMA2] << std::endl;
      app_log() << "branch cutoff = " << vParam[B_BRANCHCUTOFF] << " " << vParam[B_BRANCHMAX] << std::endl;
      ToDoSteps             = iParam[B_ENERGYUPDATEINTERVAL] - 1;
      iParam[B_WARMUPSTEPS] = 0;
      BranchMode.set(B_DMCSTAGE, 1); //set BranchModex to main stage
      //reset the histogram
      EnergyHist.clear();
      EnergyHist(vParam[B_ENOW]);
      if (sParam[MIXDMCOPT] == "yes")
      {
        app_log() << "Switching to DMC with fluctuating populations" << std::endl;
        BranchMode.set(B_POPCONTROL, 1); //use standard DMC
        WalkerController       = std::move(BackupWalkerController);
        BackupWalkerController = 0;
        vParam[B_ETRIAL]       = vParam[B_EREF];
        app_log() << "  Etrial     = " << vParam[B_ETRIAL] << std::endl;
        WalkerController->start();
      }
      //This is not necessary
      //EnergyHist(DMCEnergyHist.mean());
    }
  }
  WalkerController->setTrialEnergy(vParam[B_ETRIAL]);
  //accumulate collectables and energies for scalar.dat
  FullPrecRealType wgt_inv = WalkerController->get_num_contexts() / WalkerController->get_ensemble_property().Weight;
  walkers.Collectables *= wgt_inv;
  MyEstimator->accumulate(walkers);
}

void SimpleFixedNodeBranch::branch(int iter, UPtrVector<Crowd>& crowds,  MCPopulation& population)
{
  //collect the total weights and redistribute the walkers accordingly, using a fixed tolerance
  //RealType pop_now= WalkerController->branch(iter,walkers,0.1);
  RefVector<MCPWalker> walkers(convertUPtrToRefVector(population.get_walkers()));

  FullPrecRealType pop_now;
  if (false) // && BranchMode[B_DMCSTAGE] || iter)
    pop_now = WalkerController->branch(iter, population, 0.1);
  else
    pop_now = WalkerController->doNotBranch(iter, population); //do not branch for the first step of a warmup
  //population for trial energy modification should not include any released node walkers.
  MCDataType<FullPrecRealType>& wc_ensemble_prop = WalkerController->get_ensemble_property();
  pop_now -= wc_ensemble_prop.RNSamples;
  //current energy
  vParam[B_ENOW] = wc_ensemble_prop.Energy;
  VarianceHist(wc_ensemble_prop.Variance);
  R2Accepted(wc_ensemble_prop.R2Accepted);
  R2Proposed(wc_ensemble_prop.R2Proposed);
  //PopHist(pop_now);
  vParam[B_EREF] = EnergyHist.mean(); //current mean
  if (BranchMode[B_USETAUEFF])
    vParam[B_TAUEFF] = vParam[B_TAU] * R2Accepted.result() / R2Proposed.result();
  if (BranchMode[B_KILLNODES])
    EnergyHist(vParam[B_ENOW] - std::log(wc_ensemble_prop.LivingFraction) / vParam[B_TAUEFF]);
  else
    EnergyHist(vParam[B_ENOW]);
  if (BranchMode[B_DMCSTAGE]) // main stage
  {
    if (BranchMode[B_POPCONTROL])
    {
      if (ToDoSteps > 0)
        --ToDoSteps;
      else
      {
        vParam[B_ETRIAL] = vParam[B_EREF] + vParam[B_FEEDBACK] * (logN - std::log(pop_now));
        ToDoSteps        = iParam[B_ENERGYUPDATEINTERVAL] - 1;
      }
    }
    else
      vParam[B_ETRIAL] = vParam[B_EREF];
  }
  else //warmup
  {
    if (BranchMode[B_USETAUEFF])
      vParam[B_TAUEFF] = vParam[B_TAU] * R2Accepted.result() / R2Proposed.result();
    if (BranchMode[B_POPCONTROL])
    {
      //RealType emix=((iParam[B_WARMUPSTEPS]-ToDoSteps)<100)?(0.25*vParam[B_EREF]+0.75*vParam[B_ENOW]):vParam[B_EREF];
      //vParam[B_ETRIAL]=emix+Feedback*(logN-std::log(pop_now));
      //vParam[B_ETRIAL]=vParam[B_EREF]+Feedback*(logN-std::log(pop_now));
      if (BranchMode[B_KILLNODES])
        vParam[B_ETRIAL] = (0.00 * vParam[B_EREF] + 1.0 * vParam[B_ENOW]) +
            vParam[B_FEEDBACK] * (logN - std::log(pop_now)) -
            std::log(wc_ensemble_prop.LivingFraction) / vParam[B_TAU];
      else
        vParam[B_ETRIAL] = vParam[B_ENOW] + (logN - std::log(pop_now)) / vParam[B_TAU];
    }
    else
      vParam[B_ETRIAL] = vParam[B_EREF];
    --ToDoSteps;
    if (ToDoSteps == 0) //warmup is done
    {
      vParam[B_SIGMA2] = VarianceHist.mean();
      setBranchCutoff(vParam[B_SIGMA2], WalkerController->get_target_sigma(), 10, population.get_num_particles());
      app_log() << "\n Warmup is completed after " << iParam[B_WARMUPSTEPS] << std::endl;
      if (BranchMode[B_USETAUEFF])
        app_log() << "\n  TauEff     = " << vParam[B_TAUEFF] << "\n TauEff/Tau = " << vParam[B_TAUEFF] / vParam[B_TAU];
      else
        app_log() << "\n  TauEff proposed   = " << vParam[B_TAUEFF] * R2Accepted.result() / R2Proposed.result();
      app_log() << "\n  Etrial     = " << vParam[B_ETRIAL] << std::endl;
      app_log() << " Running average of energy = " << EnergyHist.mean() << std::endl;
      app_log() << "                  Variance = " << vParam[B_SIGMA2] << std::endl;
      app_log() << "branch cutoff = " << vParam[B_BRANCHCUTOFF] << " " << vParam[B_BRANCHMAX] << std::endl;
      ToDoSteps             = iParam[B_ENERGYUPDATEINTERVAL] - 1;
      iParam[B_WARMUPSTEPS] = 0;
      BranchMode.set(B_DMCSTAGE, 1); //set BranchModex to main stage
      //reset the histogram
      EnergyHist.clear();
      EnergyHist(vParam[B_ENOW]);
      if (sParam[MIXDMCOPT] == "yes")
      {
        app_log() << "Switching to DMC with fluctuating populations" << std::endl;
        BranchMode.set(B_POPCONTROL, 1); //use standard DMC
        WalkerController       = std::move(BackupWalkerController);
        BackupWalkerController = 0;
        vParam[B_ETRIAL]       = vParam[B_EREF];
        app_log() << "  Etrial     = " << vParam[B_ETRIAL] << std::endl;
        WalkerController->start();
      }
      //This is not necessary
      //EnergyHist(DMCEnergyHist.mean());
    }
  }
  WalkerController->setTrialEnergy(vParam[B_ETRIAL]);
  //accumulate collectables and energies for scalar.dat
  FullPrecRealType wgt_inv = WalkerController->get_num_contexts() / wc_ensemble_prop.Weight;
  //walkers.Collectables *= wgt_inv;
  for(UPtr<Crowd>& crowd_ptr: crowds)
    crowd_ptr->accumulate(population.get_num_global_walkers());
}

/**
 *
 */
void SimpleFixedNodeBranch::collect(int iter, MCWalkerConfiguration& W)
{
  //Update the current energy and accumulate.
  MCWalkerConfiguration::Walker_t& head = W.reptile->getHead();
  MCWalkerConfiguration::Walker_t& tail = W.reptile->getTail();
  vParam[B_ENOW]                        = 0.5 * (head.Properties(LOCALENERGY) + tail.Properties(LOCALENERGY));
  // app_log()<<"IN SimpleFixedNodeBranch::collect\n";
  // app_log()<<"\tvParam[B_ENOW]="<<vParam[B_ENOW]<< std::endl;
  EnergyHist(vParam[B_ENOW]);
  vParam[B_EREF] = EnergyHist.mean();
  // app_log()<<"\tvParam[B_EREF]="<<vParam[B_EREF]<< std::endl;
  //Update the energy variance and R2 for effective timestep and filtering.
  VarianceHist(std::pow(vParam[B_ENOW] - vParam[B_EREF], 2));
  R2Accepted(head.Properties(R2ACCEPTED));
  R2Proposed(head.Properties(R2PROPOSED));
  // app_log()<<"\thead.Properties(R2ACCEPTED)="<<head.Properties(R2ACCEPTED)<< std::endl;
  // app_log()<<"\thead.Properties(R2PROPOSED)="<<head.Properties(R2PROPOSED)<< std::endl;
  //  app_log()<<"\tR2Accepted="<<R2Accepted.result()<< std::endl;
  // app_log()<<"\tR2Proposed="<<R2Proposed.result()<< std::endl;
  //  app_log()<<"\tR2Accept/R2Prop="<<R2Accepted.result()/R2Proposed.result()<< std::endl;
  // app_log()<<"\t <E^2> = "<<VarianceHist.mean()<< std::endl;
  // app_log()<<"\t <E>   = "<<EnergyHist.mean()<< std::endl;
  //  app_log()<<"\t <E>^2 = "<<std::pow(EnergyHist.mean(),2)<< std::endl;
  // app_log()<<"\t var = "<<VarianceHist.mean()-pow(EnergyHist.mean(),2)<< std::endl;
  // app_log()<<"--------------\n";
  //current mean
  if (BranchMode[B_USETAUEFF])
  {
    //app_log()<<" BRANCHMODE = "<<BranchMode[B_USETAUEFF]<< std::endl;
    vParam[B_TAUEFF] = vParam[B_TAU] * R2Accepted.result() / R2Proposed.result();
    //  app_log()<<"\tvParam[B_TAU]="<<vParam[B_TAU]<<" "<<vParam[B_TAUEFF]<< std::endl;
  }
  /*
  if(BranchMode[B_RMCSTAGE]) // main stage
  {
    if(BranchMode[B_POPCONTROL])
    {
      if(ToDoSteps>0)
        --ToDoSteps;
      else
      {
        vParam[B_ETRIAL]=vParam[B_EREF]+vParam[B_FEEDBACK]*(logN-std::log(pop_now));
        ToDoSteps=iParam[B_ENERGYUPDATEINTERVAL]-1;
      }
    }
    else
      vParam[B_ETRIAL]=vParam[B_EREF];
  }*/
  //app_log()<<"BranchMode[B_RMCSTAGE]="<<BranchMode[B_RMCSTAGE]<< std::endl;
  if (!BranchMode[B_RMCSTAGE]) //warmup
  {
    if (BranchMode[B_USETAUEFF])
      vParam[B_TAUEFF] = vParam[B_TAU] * R2Accepted.result() / R2Proposed.result();
    // app_log()<<"\t <E^2> = "<<VarianceHist.mean()<< std::endl;
    // app_log()<<"\t <E>   = "<<EnergyHist.mean()<< std::endl;
    // app_log()<<"\t <E>^2 = "<<std::pow(EnergyHist.mean(),2)<< std::endl;
    //app_log()<<"\t var = "<<VarianceHist.mean()-std::pow(EnergyHist.mean(),2)<< std::endl;
    // app_log()<<"\t var = "<<VarianceHist.mean()<< std::endl;
    //  app_log()<<"----\n";
    //app_log()<<"ToDoSteps="<<ToDoSteps<< std::endl;
    vParam[B_ETRIAL] = vParam[B_EREF];
    --ToDoSteps;
    if (ToDoSteps == 0) //warmup is done
    {
      vParam[B_TAUEFF] = vParam[B_TAU] * R2Accepted.result() / R2Proposed.result();
      vParam[B_SIGMA2] = VarianceHist.mean();
      setBranchCutoff(vParam[B_SIGMA2], vParam[B_FILTERSCALE], vParam[B_FILTERSCALE], W.R.size());
      app_log() << "\n Warmup is completed after " << iParam[B_WARMUPSTEPS] << " steps." << std::endl;
      if (BranchMode[B_USETAUEFF])
        app_log() << "\n  TauEff     = " << vParam[B_TAUEFF] << "\n TauEff/Tau = " << vParam[B_TAUEFF] / vParam[B_TAU];
      else
        app_log() << "\n  TauEff proposed   = " << vParam[B_TAUEFF] * R2Accepted.result() / R2Proposed.result();
      app_log() << "\n Running average of energy = " << EnergyHist.mean() << std::endl;
      app_log() << "\n                  Variance = " << vParam[B_SIGMA2] << std::endl;
      app_log() << "\nbranch cutoff = " << vParam[B_BRANCHCUTOFF] << " " << vParam[B_BRANCHMAX] << std::endl;
      ToDoSteps             = iParam[B_ENERGYUPDATEINTERVAL] - 1;
      iParam[B_WARMUPSTEPS] = 0;
      BranchMode.set(B_RMCSTAGE, 1); //set BranchModex to main stage
      //reset the histogram
      EnergyHist.clear();
      EnergyHist(vParam[B_ENOW]);
    }
  }
  //accumulate collectables and energies for scalar.dat
  MyEstimator->accumulate(W);
}

/** Calculates and saves various action components, also does necessary updates for running averages.
 *
 */
void SimpleFixedNodeBranch::reset()
{
  //use effective time step of BranchInterval*Tau
  //Feed = 1.0/(static_cast<RealType>(NumGeneration*BranchInterval)*Tau);
  //logN = Feed*std::log(static_cast<RealType>(Nideal));
  //JNKIM passive
  //BranchMode.set(B_DMC,1);//set DMC
  //BranchMode.set(B_DMCSTAGE,0);//set warmup
  if (WalkerController)
  {
    //this is to compare the time step errors
    BranchMode.set(B_USETAUEFF, sParam[USETAUOPT] == "no");
    if (BranchMode[B_DMCSTAGE]) //
      ToDoSteps = iParam[B_ENERGYUPDATEINTERVAL] - 1;
    else
      ToDoSteps = iParam[B_WARMUPSTEPS];
    if (BranchMode[B_POPCONTROL])
    {
      //logN = Feedback*std::log(static_cast<RealType>(iParam[B_TARGETWALKERS]));
      logN = std::log(static_cast<FullPrecRealType>(iParam[B_TARGETWALKERS]));
      if (vParam[B_FEEDBACK] == 0.0)
        vParam[B_FEEDBACK] = 1.0;
    }
    else
    {
      //may set Eref to a safe value
      //if(EnergyHistory.count()<5) Eref -= vParam[EnergyWindowIndex];
      vParam[B_ETRIAL]   = vParam[B_EREF];
      vParam[B_FEEDBACK] = 0.0;
      logN               = 0.0;
    }
    //       vParam(abranch.vParam)
    WalkerController->start();
  }
  if (BranchMode[B_RMC])
  {
    //this is to compare the time step errors
    // BranchMode.set(B_USETAUEFF,sParam[USETAUOPT]=="no");
    if (BranchMode[B_RMCSTAGE]) //
      ToDoSteps = iParam[B_ENERGYUPDATEINTERVAL] - 1;
    else
      ToDoSteps = iParam[B_WARMUPSTEPS];
  }
}

void SimpleFixedNodeBranch::setRN(bool rn)
{
  RN = rn;
  WalkerController->set_write_release_nodes(rn);
  WalkerController->start();
}


int SimpleFixedNodeBranch::resetRun(xmlNodePtr cur)
{
  app_log() << "BRANCH resetRun" << std::endl;
  //estimator is always reset
  MyEstimator->reset();
  MyEstimator->setCollectionMode(true);
  std::bitset<B_MODE_MAX> bmode(BranchMode);
  IParamType iparam_old(iParam);
  VParamType vparam_old(vParam);
  myNode = cur;
  //store old target
  int nw_target = iParam[B_TARGETWALKERS];
  m_param.put(cur);

  int target_min = -1;
  ParameterSet p;
  p.add(target_min, "minimumtargetwalkers", "int"); //p.add(target_min,"minimumTargetWalkers","int");
  p.add(target_min, "minimumsamples", "int");       //p.add(target_min,"minimumSamples","int");
  p.put(cur);

  if (iParam[B_TARGETWALKERS] < target_min)
  {
    iParam[B_TARGETWALKERS] = target_min;
  }

  bool same_wc = true;
  if (BranchMode[B_DMC] && WalkerController)
  {
    std::string reconfig("no");
    // method is actually IndexType so conceivably indicates much more that reconfig="yes" or "no"
    if (WalkerController->get_method())
      reconfig = "yes";
    std::string reconfig_prev(reconfig);
    ParameterSet p;
    p.add(reconfig, "reconfiguration", "string");
    p.put(cur);
    same_wc = (reconfig == reconfig_prev);
  }

  //everything is the same, do nothing
  if (same_wc && bmode == BranchMode && std::equal(iParam.begin(), iParam.end(), iparam_old.begin()) &&
      std::equal(vParam.begin(), vParam.end(), vparam_old.begin()))
  {
    app_log() << "  Continue with the same input as the previous block." << std::endl;
    app_log().flush();
    //return 1;
  }
  app_log() << " SimpleFixedNodeBranch::resetRun detected changes in <parameter>'s " << std::endl;
  app_log() << " BranchMode : " << bmode << " " << BranchMode << std::endl;

  //vmc does not need to do anything with WalkerController
  if (!BranchMode[B_DMC])
  {
    app_log() << " iParam (old): " << iparam_old << std::endl;
    app_log() << " iParam (new): " << iParam << std::endl;
    app_log() << " vParam (old): " << vparam_old << std::endl;
    app_log() << " vParam (new): " << vParam << std::endl;
    app_log().flush();
    return 1;
  }

  if (WalkerController == nullptr)
  {
    APP_ABORT("SimpleFixedNodeBranch::resetRun cannot initialize WalkerController");
  }

  if (!same_wc)
  {
    app_log() << "Destroy WalkerController. Existing method " << WalkerController->get_method() << std::endl;
    ;
    WalkerController.reset(createWalkerController(iParam[B_TARGETWALKERS], MyEstimator->getCommunicator(), myNode));
    app_log().flush();

    BranchMode[B_POPCONTROL] = (WalkerController->get_method() == 0);
    if (BranchMode[B_POPCONTROL])
    {
      vParam[B_ETRIAL] = vParam[B_EREF];
      if (vParam[B_FEEDBACK] == 0.0)
        vParam[B_FEEDBACK] = 1.0;
    }
  }

  //always add a warmup step using default 10 steps
  R2Accepted.clear();
  R2Proposed.clear();
  //R2Accepted(1.0e-12);
  //R2Proposed(1.0e-12);
  BranchMode[B_DMCSTAGE] = 0;
  WalkerController->put(myNode);
  ToDoSteps = iParam[B_WARMUPSTEPS] = (iParam[B_WARMUPSTEPS]) ? iParam[B_WARMUPSTEPS] : 10;
  setBranchCutoff(vParam[B_SIGMA2], WalkerController->get_target_sigma(), 10);
  WalkerController->reset();
#ifdef QMC_CUDA
  reset(); // needed. Ye
#endif
  if (BackupWalkerController)
    BackupWalkerController->reset();

  iParam[B_MAXWALKERS] = WalkerController->get_n_max();
  iParam[B_MINWALKERS] = WalkerController->get_n_min();

  app_log() << " iParam (old): " << iparam_old << std::endl;
  app_log() << " iParam (new): " << iParam << std::endl;
  app_log() << " vParam (old): " << vparam_old << std::endl;
  app_log() << " vParam (new): " << vParam << std::endl;

  app_log() << std::endl << " Using branching cutoff scheme " << branching_cutoff_scheme << std::endl;

  app_log().flush();

  //  return static_cast<int>(iParam[B_TARGETWALKERS]*1.01/static_cast<double>(nw_target));
  return static_cast<int>(round(static_cast<double>(iParam[B_TARGETWALKERS] / static_cast<double>(nw_target))));
}

void SimpleFixedNodeBranch::checkParameters(MCWalkerConfiguration& w)
{
  std::ostringstream o;
  if (!BranchMode[B_DMCSTAGE])
  {
    FullPrecRealType e, sigma2;
    MyEstimator->getCurrentStatistics(w, e, sigma2);
    vParam[B_ETRIAL] = vParam[B_EREF] = e;
    vParam[B_SIGMA2]                  = sigma2;
    EnergyHist.clear();
    VarianceHist.clear();
    //DMCEnergyHist.clear();
    EnergyHist(vParam[B_EREF]);
    VarianceHist(vParam[B_SIGMA2]);
    //DMCEnergyHist(vParam[B_EREF]);
    o << "SimpleFixedNodeBranch::checkParameters " << std::endl;
    o << "  Average Energy of a population  = " << e << std::endl;
    o << "  Energy Variance = " << vParam[B_SIGMA2] << std::endl;
  }
  app_log() << o.str() << std::endl;
  app_log().flush();
}

void SimpleFixedNodeBranch::checkParameters(const int global_walkers, RefVector<MCPWalker>& walkers)
{
  std::ostringstream o;
  if (!BranchMode[B_DMCSTAGE])
  {
    FullPrecRealType e, sigma2;
    MyEstimator->getCurrentStatistics(global_walkers, walkers, e, sigma2);
    vParam[B_ETRIAL] = vParam[B_EREF] = e;
    vParam[B_SIGMA2]                  = sigma2;
    EnergyHist.clear();
    VarianceHist.clear();
    //DMCEnergyHist.clear();
    EnergyHist(vParam[B_EREF]);
    VarianceHist(vParam[B_SIGMA2]);
    //DMCEnergyHist(vParam[B_EREF]);
    o << "SimpleFixedNodeBranch::checkParameters " << std::endl;
    o << "  Average Energy of a population  = " << e << std::endl;
    o << "  Energy Variance = " << vParam[B_SIGMA2] << std::endl;
  }
  app_log() << o.str() << std::endl;
  app_log().flush();
}

void SimpleFixedNodeBranch::finalize(const int global_walkers, RefVector<MCPWalker>& walkers)
{
  std::ostringstream o;
  if (WalkerController)
  {
    o << "====================================================";
    o << "\n  SimpleFixedNodeBranch::finalize after a DMC block";
    o << "\n    QMC counter                   = " << iParam[B_COUNTER];
    o << "\n    time step                     = " << vParam[B_TAU];
    o << "\n    effective time step           = " << vParam[B_TAUEFF];
    o << "\n    trial energy                  = " << vParam[B_ETRIAL];
    o << "\n    reference energy              = " << vParam[B_EREF];
    o << "\n    reference variance            = " << vParam[B_SIGMA2];
    o << "\n    target walkers                = " << iParam[B_TARGETWALKERS];
    o << "\n    branch cutoff                 = " << vParam[B_BRANCHCUTOFF] << " " << vParam[B_BRANCHMAX];
    o << "\n    Max and minimum walkers per node= " << iParam[B_MAXWALKERS] << " " << iParam[B_MINWALKERS];
    o << "\n    Feedback                      = " << vParam[B_FEEDBACK];
    o << "\n    QMC Status (BranchMode)       = " << BranchMode;
    o << "\n====================================================";
  }
  //running RMC
  else if (BranchMode[B_RMC])
  {
    o << "====================================================";
    o << "\n  SimpleFixedNodeBranch::finalize after a RMC block";
    o << "\n    QMC counter                   = " << iParam[B_COUNTER];
    o << "\n    time step                     = " << vParam[B_TAU];
    o << "\n    effective time step           = " << vParam[B_TAUEFF];
    o << "\n    reference energy              = " << vParam[B_EREF];
    o << "\n    reference variance            = " << vParam[B_SIGMA2];
    o << "\n    cutoff energy                 = " << vParam[B_BRANCHCUTOFF] << " " << vParam[B_BRANCHMAX];
    o << "\n    QMC Status (BranchMode)       = " << BranchMode;
    o << "\n====================================================";
  }
  else
  {
    //running VMC
    FullPrecRealType e, sigma2;
    //MyEstimator->getEnergyAndWeight(e,w,sigma2);
    MyEstimator->getCurrentStatistics(global_walkers, walkers , e, sigma2);
    vParam[B_ETRIAL] = vParam[B_EREF] = e;
    vParam[B_SIGMA2]                  = sigma2;
    //this is just to avoid diving by n-1  == 0
    EnergyHist(vParam[B_EREF]);
    //add Eref to the DMCEnergyHistory
    //DMCEnergyHist(vParam[B_EREF]);
    o << "====================================================";
    o << "\n  SimpleFixedNodeBranch::finalize after a VMC block";
    o << "\n    QMC counter        = " << iParam[B_COUNTER];
    o << "\n    time step          = " << vParam[B_TAU];
    o << "\n    reference energy   = " << vParam[B_EREF];
    o << "\n    reference variance = " << vParam[B_SIGMA2];
    o << "\n====================================================";
  }
  app_log() << o.str() << std::endl;
  write(RootName, true);
}

void SimpleFixedNodeBranch::finalize(MCWalkerConfiguration& w)
{
  std::ostringstream o;
  if (WalkerController)
  {
    o << "====================================================";
    o << "\n  SimpleFixedNodeBranch::finalize after a DMC block";
    o << "\n    QMC counter                   = " << iParam[B_COUNTER];
    o << "\n    time step                     = " << vParam[B_TAU];
    o << "\n    effective time step           = " << vParam[B_TAUEFF];
    o << "\n    trial energy                  = " << vParam[B_ETRIAL];
    o << "\n    reference energy              = " << vParam[B_EREF];
    o << "\n    reference variance            = " << vParam[B_SIGMA2];
    o << "\n    target walkers                = " << iParam[B_TARGETWALKERS];
    o << "\n    branch cutoff                 = " << vParam[B_BRANCHCUTOFF] << " " << vParam[B_BRANCHMAX];
    o << "\n    Max and minimum walkers per node= " << iParam[B_MAXWALKERS] << " " << iParam[B_MINWALKERS];
    o << "\n    Feedback                      = " << vParam[B_FEEDBACK];
    o << "\n    QMC Status (BranchMode)       = " << BranchMode;
    o << "\n====================================================";
  }
  //running RMC
  else if (BranchMode[B_RMC])
  {
    o << "====================================================";
    o << "\n  SimpleFixedNodeBranch::finalize after a RMC block";
    o << "\n    QMC counter                   = " << iParam[B_COUNTER];
    o << "\n    time step                     = " << vParam[B_TAU];
    o << "\n    effective time step           = " << vParam[B_TAUEFF];
    o << "\n    reference energy              = " << vParam[B_EREF];
    o << "\n    reference variance            = " << vParam[B_SIGMA2];
    o << "\n    cutoff energy                 = " << vParam[B_BRANCHCUTOFF] << " " << vParam[B_BRANCHMAX];
    o << "\n    QMC Status (BranchMode)       = " << BranchMode;
    o << "\n====================================================";
  }
  else
  {
    //running VMC
    FullPrecRealType e, sigma2;
    //MyEstimator->getEnergyAndWeight(e,w,sigma2);
    MyEstimator->getCurrentStatistics(w, e, sigma2);
    vParam[B_ETRIAL] = vParam[B_EREF] = e;
    vParam[B_SIGMA2]                  = sigma2;
    //this is just to avoid diving by n-1  == 0
    EnergyHist(vParam[B_EREF]);
    //add Eref to the DMCEnergyHistory
    //DMCEnergyHist(vParam[B_EREF]);
    o << "====================================================";
    o << "\n  SimpleFixedNodeBranch::finalize after a VMC block";
    o << "\n    QMC counter        = " << iParam[B_COUNTER];
    o << "\n    time step          = " << vParam[B_TAU];
    o << "\n    reference energy   = " << vParam[B_EREF];
    o << "\n    reference variance = " << vParam[B_SIGMA2];
    o << "\n====================================================";
  }
  app_log() << o.str() << std::endl;
  write(RootName, true);
}

/**  Parse the xml file for parameters
 *@param cur current xmlNode
 *@param LogOut std::ostream to which the run-time report is sent
 *
 * Few important parameters are:
 * <ul>
 * <li> en_ref: a reference energy
 * <li> num_gen: number of generations \f$N_G\f$ to reach  equilibrium, used in the feedback parameter
 * \f$ feed = \frac{1}{N_G \tau} \f$
 * </ul>
 */
bool SimpleFixedNodeBranch::put(xmlNodePtr cur)
{
  //save it
  myNode = cur;
  //check dmc/vmc and decide to create WalkerControllerBase
  m_param.put(cur);
  reset();
  MyEstimator->setCollectionMode(true); //always collect
  return true;
}

void SimpleFixedNodeBranch::write(const std::string& fname, bool overwrite)
{
  RootName = fname;
  if (MyEstimator->is_manager())
  {
    //\since 2008-06-24
    vParam[B_ACC_ENERGY]  = EnergyHist.result();
    vParam[B_ACC_SAMPLES] = EnergyHist.count();
    BranchIO hh(*this, MyEstimator->getCommunicator());
    bool success = hh.write(fname);
  }
}

void SimpleFixedNodeBranch::read(const std::string& fname)
{
  BranchMode.set(B_RESTART, 0);
  if (fname.empty())
    return;
  vParam[B_ACC_ENERGY]  = EnergyHist.result();
  vParam[B_ACC_SAMPLES] = EnergyHist.count();
  BranchIO hh(*this, MyEstimator->getCommunicator());
  BranchModeType bmode(BranchMode);
  bool success = hh.read(fname);
  if (success && R2Proposed.good() && bmode[B_POPCONTROL] == BranchMode[B_POPCONTROL])
  {
    BranchMode.set(B_RESTART, 1);
    app_log() << "    Restarting, cummulative properties:"
              << "\n      energy     = " << EnergyHist.mean() << "\n      variance   = " << VarianceHist.mean()
              << "\n      r2accepted = " << R2Accepted.mean() << "\n      r2proposed = " << R2Proposed.mean()
              << std::endl;
  }
  else
  {
    if (BranchMode[B_POPCONTROL] != bmode[B_POPCONTROL])
    {
      app_log() << "  Population control method has changed from " << BranchMode[B_POPCONTROL] << " to "
                << bmode[B_POPCONTROL] << std::endl;
      BranchMode[B_POPCONTROL] = bmode[B_POPCONTROL];
    }
  }

  app_log().flush();
}

void SimpleFixedNodeBranch::setBranchCutoff(FullPrecRealType variance,
                                            FullPrecRealType targetSigma,
                                            FullPrecRealType maxSigma,
                                            int Nelec)
{
  if (branching_cutoff_scheme == "DRV")
  {
    // eq.(3), J. Chem. Phys. 89, 3629 (1988).
    // eq.(9), J. Chem. Phys. 99, 2865 (1993).
    vParam[B_BRANCHCUTOFF] = 2.0 / std::sqrt(vParam[B_TAU]);
  }
  else if (branching_cutoff_scheme == "ZSGMA")
  {
    // eq.(6), Phys. Rev. B 93, 241118(R) (2016)
    // do nothing if Nelec is not passed in.
    if (Nelec > 0)
      vParam[B_BRANCHCUTOFF] = 0.2 * std::sqrt(Nelec / vParam[B_TAU]);
  }
  else if (branching_cutoff_scheme == "YL")
  {
    // a scheme from Ye Luo.
    vParam[B_BRANCHCUTOFF] = std::sqrt(variance) * std::min(targetSigma, std::sqrt(1.0 / vParam[B_TAU]));
  }
  else if (branching_cutoff_scheme == "classic")
  {
    // default QMCPACK choice which is the same as v3.0.0 and before.
    vParam[B_BRANCHCUTOFF] = std::min(std::max(variance * targetSigma, maxSigma), 2.5 / vParam[B_TAU]);
  }
  else
    APP_ABORT("SimpleFixedNodeBranch::setBranchCutoff unknown branching cutoff scheme " + branching_cutoff_scheme);

  vParam[B_BRANCHMAX]    = vParam[B_BRANCHCUTOFF] * 1.5;
  vParam[B_BRANCHFILTER] = 1.0 / (vParam[B_BRANCHMAX] - vParam[B_BRANCHCUTOFF]);
}

} // namespace qmcplusplus
