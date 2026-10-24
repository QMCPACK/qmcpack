//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2019 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//                    Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Luke Shulenburger, lshulen@sandia.gov, Sandia National Laboratories
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Raymond Clay III, j.k.rofling@gmail.com, Lawrence Livermore National Laboratory
//                    Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


#include "SimpleFixedNodeBranch.h"
#include "QMCDrivers/DMC/WalkerControlFactory.h"
#include <numeric>
#include "OhmmsData/FileUtility.h"
#include "Utilities/RandomGenerator.h"
#include "QMCDrivers/WalkerControlBase.h"
#include "Estimators/EstimatorManagerBase.h"
#include "QMCDrivers/BranchIO.h"
#include "Particle/Reptile.h"
#include "type_traits/template_types.hpp"
#include "Message/UniformCommunicateError.h"

namespace qmcplusplus
{
using WP = WalkerProperties::Indexes;

///enum to yes/no options saved in sParam
enum
{
  COMBOPT,
  USETAUOPT,
  MIXDMCOPT,
  DUMMYOPT
};

SimpleFixedNodeBranch::SimpleFixedNodeBranch(RealType tau, int nideal) //, PopHist(5), DMCEnergyHist(5)
{
  branch_mode_.set(B_DMCSTAGE, 0);     //warmup stage
  branch_mode_.set(B_POPCONTROL, 1);   //use standard DMC
  branch_mode_.set(B_USETAUEFF, 1);    //use taueff
  branch_mode_.set(B_CLEARHISTORY, 0); //clear history and start with the current average
  branch_mode_.set(B_KILLNODES, 0);    //when killing walkers at nodes etrial is updated differently
  v_param_.fill(1.0);
  v_param_[SBVP::TAU]         = tau;
  v_param_[SBVP::TAUEFF]      = tau;
  v_param_[SBVP::FEEDBACK]    = 1.0;
  v_param_[SBVP::FILTERSCALE] = 10;
  r2_accepted_(1.0e-10);
  r2_proposed_(1.0e-10);
  //set the default values for integer parameters
  iParam[B_WARMUPSTEPS]          = 200;
  iParam[B_ENERGYUPDATEINTERVAL] = 1;
  iParam[B_BRANCHINTERVAL]       = 1;
  iParam[B_TARGETWALKERS]        = 0;
  iParam[B_MAXWALKERS]           = nideal;
  iParam[B_MINWALKERS]           = nideal;
  iParam[B_COUNTER]              = -1;
  //default is no
  s_param_.resize(DUMMYOPT, "no");
  //default is classic
  branching_cutoff_scheme_ = "classic";
  registerParameters();
  reset();
}

/** copy constructor
 *
 * Copy only selected data members and WalkerController is never copied.
 */
SimpleFixedNodeBranch::SimpleFixedNodeBranch(const SimpleFixedNodeBranch& abranch)
    : iParam(abranch.iParam),
      branch_mode_(abranch.branch_mode_),
      v_param_(abranch.v_param_),
      branching_cutoff_scheme_(abranch.branching_cutoff_scheme_),
      s_param_(abranch.s_param_)
{
  registerParameters();
  reset();
}

SimpleFixedNodeBranch::~SimpleFixedNodeBranch() = default;

void SimpleFixedNodeBranch::registerParameters()
{
  m_param_.add(iParam[B_WARMUPSTEPS], "warmupSteps");
  m_param_.add(iParam[B_WARMUPSTEPS], "warmupsteps");
  m_param_.add(iParam[B_ENERGYUPDATEINTERVAL], "energyUpdateInterval");
  m_param_.add(iParam[B_BRANCHINTERVAL], "branchInterval");
  m_param_.add(iParam[B_TARGETWALKERS], "targetWalkers");
  m_param_.add(iParam[B_TARGETWALKERS], "targetwalkers");
  m_param_.add(iParam[B_TARGETWALKERS], "target_walkers");
  //trial energy
  m_param_.add(v_param_[SBVP::EREF], "refEnergy");
  m_param_.add(v_param_[SBVP::EREF], "ref_energy");
  m_param_.add(v_param_[SBVP::EREF], "en_ref");
  m_param_.add(v_param_[SBVP::TAU], "tau");
  m_param_.add(v_param_[SBVP::TAU], "timestep");
  m_param_.add(v_param_[SBVP::TAU], "timeStep");
  m_param_.add(v_param_[SBVP::TAU], "TimeStep");
  //filterscale:  sets the filtercutoff to sigma*filterscale
  m_param_.add(v_param_[SBVP::FILTERSCALE], "filterscale");
  //feed back parameter for population control
  m_param_.add(v_param_[SBVP::FEEDBACK], "feedback");
  //turn on/off effective tau onl for time-step error comparisons
  m_param_.add(s_param_[USETAUOPT], "useBareTau");
  m_param_.add(s_param_[MIXDMCOPT], "warmupByReconfiguration");
  m_param_.add(branching_cutoff_scheme_, "branching_cutoff_scheme");
  m_param_.add(debug_disable_branching_, "debug_disable_branching", {"no", "yes"});
}

void SimpleFixedNodeBranch::startWalkerController(const std::string& froot)
{
  root_name_ = froot;
  if (walker_controller_)
    walker_controller_->start(root_name_);
}

int SimpleFixedNodeBranch::initWalkerController(MCWalkerConfiguration& walkers, bool fixW, bool killwalker)
{
  branch_mode_.set(B_DMC, 1);                               //set DMC
  branch_mode_.set(B_DMCSTAGE, iParam[B_WARMUPSTEPS] == 0); //use warmup
  //this is not necessary
  //check if tau is different and set the initial values
  //v_param_[SBVP::TAU]=tau;
  bool fromscratch     = false;
  FullPrecRealType tau = v_param_[SBVP::TAU];

  int nwtot_now = walkers.getGlobalNumWalkers();

  //this is the first time DMC is used
  if (walker_controller_ == nullptr)
  {
    if (iParam[B_TARGETWALKERS] == 0)
    {
      Communicate& acomm = *my_estimator_->getCommunicator();
      int ncontexts      = acomm.size();
      std::vector<int> nw(ncontexts, 0), nwoff(ncontexts + 1, 0);
      nw[acomm.rank()] = walkers.getActiveWalkers();
      acomm.allreduce(nw);
      for (int ip = 0; ip < ncontexts; ++ip)
        nwoff[ip + 1] = nwoff[ip] + nw[ip];
      walkers.setWalkerOffsets(nwoff);
      iParam[B_TARGETWALKERS] = nwoff[ncontexts];
    }
    // \todo reparsing XML nodes is an antipattern, remove.
    walker_controller_.reset(createWalkerController(iParam[B_TARGETWALKERS], *my_estimator_->getCommunicator(), my_node_));
    if (!branch_mode_[B_RESTART])
    {
      fromscratch = true;
      app_log() << "  START ALL OVER " << std::endl;
      v_param_[SBVP::TAUEFF] = tau;
      branch_mode_.set(B_POPCONTROL, !fixW); //fixW -> 0
      branch_mode_.set(B_KILLNODES, killwalker);
      iParam[B_MAXWALKERS] = walker_controller_->get_n_max();
      iParam[B_MINWALKERS] = walker_controller_->get_n_min();
      if (!fixW && s_param_[MIXDMCOPT] == "yes")
      {
        app_log() << "Warmup DMC is done with a fixed population " << iParam[B_TARGETWALKERS] << std::endl;
        backup_walker_controller_ = std::move(walker_controller_); //save the main controller
        // \todo: Not likely to be ok to call reset on a moved from WalkerController!
        walker_controller_.reset(
            createWalkerController(iParam[B_TARGETWALKERS], *my_estimator_->getCommunicator(), my_node_, true));
        branch_mode_.set(B_POPCONTROL, 0);
      }
      //PopHist.clear();
      //PopHist.reserve(std::max(iParam[B_ENERGYUPDATEINTERVAL],5));
    }
    walker_controller_->setWalkerID(walkers);
  }
  //else
  //{
  //  branch_mode_.set(B_DMCSTAGE,0);//always reset warmup
  //}
  //update the simulation parameters
  walker_controller_->put(my_node_);
  //assign current Eref and a large number for variance
  walker_controller_->setTrialEnergy(v_param_[SBVP::ETRIAL]);
  this->reset();
  if (fromscratch)
  {
    //determine the branch cutoff to limit wild weights based on the sigma and sigmaBound
    setBranchCutoff(v_param_[SBVP::SIGMA2], walker_controller_->get_target_sigma(), 50, walkers.R.size());
    v_param_[SBVP::TAUEFF] = tau * r2_accepted_.result() / r2_proposed_.result();
  }
  //reset controller
  walker_controller_->reset();
  if (backup_walker_controller_)
    backup_walker_controller_->reset();
  app_log() << "  QMC counter      = " << iParam[B_COUNTER] << std::endl;
  app_log() << "  time step        = " << v_param_[SBVP::TAU] << std::endl;
  app_log() << "  effective time step = " << v_param_[SBVP::TAUEFF] << std::endl;
  app_log() << "  trial energy     = " << v_param_[SBVP::ETRIAL] << std::endl;
  app_log() << "  reference energy = " << v_param_[SBVP::EREF] << std::endl;
  app_log() << "  Feedback = " << v_param_[SBVP::FEEDBACK] << std::endl;
  app_log() << "  reference variance = " << v_param_[SBVP::SIGMA2] << std::endl;
  app_log() << "  target walkers = " << iParam[B_TARGETWALKERS] << std::endl;
  app_log() << "  branching cutoff scheme " << branching_cutoff_scheme_ << std::endl;
  app_log() << "  branch cutoff = " << v_param_[SBVP::BRANCHCUTOFF] << " " << v_param_[SBVP::BRANCHMAX] << std::endl;
  app_log() << "  Max and minimum walkers per node= " << iParam[B_MAXWALKERS] << " " << iParam[B_MINWALKERS]
            << std::endl;
  app_log() << "  QMC Status (BranchMode) = " << branch_mode_ << std::endl;
  if (debug_disable_branching_ == "yes")
    app_log() << "  Disable branching for debugging as the user input request." << std::endl;

  return int(round(double(iParam[B_TARGETWALKERS]) / double(nwtot_now)));
}

void SimpleFixedNodeBranch::initReptile(MCWalkerConfiguration& W)
{
  RealType allowedFlux = 50.0;
  branch_mode_.set(B_RMC, 1);                               //set RMC
  branch_mode_.set(B_RMCSTAGE, iParam[B_WARMUPSTEPS] == 0); //use warmup
  //this is not necessary
  //check if tau is different and set the initial values
  //v_param_[SBVP::TAU]=tau;
  bool fromscratch     = false;
  FullPrecRealType tau = v_param_[SBVP::TAU];
  //this is the first time DMC is used
  if (walker_controller_ == 0)
  {
    if (!branch_mode_[B_RESTART])
    {
      fromscratch = true;
      app_log() << "  START ALL OVER " << std::endl;
      v_param_[SBVP::TAUEFF] = tau;
      //PopHist.clear();
      //PopHist.reserve(std::max(iParam[B_ENERGYUPDATEINTERVAL],5));
    }
  }
  //else
  //{
  //  branch_mode_.set(B_DMCSTAGE,0);//always reset warmup
  //}
  this->reset();
  if (fromscratch)
  {
    //determine the branch cutoff to limit wild weights based on the sigma and sigmaBound
    setBranchCutoff(v_param_[SBVP::SIGMA2], allowedFlux, 50, W.R.size());
    v_param_[SBVP::TAUEFF] = tau * r2_accepted_.result() / r2_proposed_.result();
  }
  //reset controller
  app_log() << "  QMC counter      = " << iParam[B_COUNTER] << std::endl;
  app_log() << "  time step        = " << v_param_[SBVP::TAU] << std::endl;
  app_log() << "  effective time step = " << v_param_[SBVP::TAUEFF] << std::endl;
  app_log() << "  reference energy = " << v_param_[SBVP::EREF] << std::endl;
  app_log() << "  Feedback = " << v_param_[SBVP::FEEDBACK] << std::endl;
  app_log() << "  reference variance = " << v_param_[SBVP::SIGMA2] << std::endl;
  app_log() << "  branching cutoff scheme " << branching_cutoff_scheme_ << std::endl;
  app_log() << "  branch cutoff = " << v_param_[SBVP::BRANCHCUTOFF] << " " << v_param_[SBVP::BRANCHMAX] << std::endl;
  app_log() << "  QMC Status (BranchMode) = " << branch_mode_ << std::endl;
}

void SimpleFixedNodeBranch::flush(int counter)
{
  if (counter == 0 && walker_controller_)
    walker_controller_->reset();
}

void SimpleFixedNodeBranch::branch(int iter, MCWalkerConfiguration& walkers)
{
  //collect the total weights and redistribute the walkers accordingly, using a fixed tolerance

  FullPrecRealType pop_now;
  if (debug_disable_branching_ == "no" && (branch_mode_[B_DMCSTAGE] || iter))
    pop_now = walker_controller_->branch(iter, walkers, 0.1);
  else
    pop_now = walker_controller_->doNotBranch(iter, walkers); //do not branch for the first step of a warmup

  //population for trial energy modification should not include any released node walkers.
  pop_now -= walker_controller_->get_ensemble_property().RNSamples;
  //current energy
  v_param_[SBVP::ENOW] = walker_controller_->get_ensemble_property().Energy;
  variance_hist_(walker_controller_->get_ensemble_property().Variance);
  r2_accepted_(walker_controller_->get_ensemble_property().R2Accepted);
  r2_proposed_(walker_controller_->get_ensemble_property().R2Proposed);
  //PopHist(pop_now);
  v_param_[SBVP::EREF] = energy_hist_.mean(); //current mean
  if (branch_mode_[B_USETAUEFF])
    v_param_[SBVP::TAUEFF] = v_param_[SBVP::TAU] * r2_accepted_.result() / r2_proposed_.result();

  if (branch_mode_[B_KILLNODES])
    energy_hist_(v_param_[SBVP::ENOW] -
                 std::log(walker_controller_->get_ensemble_property().LivingFraction) / v_param_[SBVP::TAUEFF]);
  else
    energy_hist_(v_param_[SBVP::ENOW]);

  if (branch_mode_[B_DMCSTAGE]) // main stage
  {
    if (branch_mode_[B_POPCONTROL])
    {
      if (todo_steps_ > 0)
        --todo_steps_;
      else
      {
        v_param_[SBVP::ETRIAL] = v_param_[SBVP::EREF] + v_param_[SBVP::FEEDBACK] * (log_n_ - std::log(pop_now));
        todo_steps_            = iParam[B_ENERGYUPDATEINTERVAL] - 1;
      }
    }
    else
      v_param_[SBVP::ETRIAL] = v_param_[SBVP::EREF];
  }
  else //warmup
  {
    if (branch_mode_[B_USETAUEFF])
      v_param_[SBVP::TAUEFF] = v_param_[SBVP::TAU] * r2_accepted_.result() / r2_proposed_.result();
    if (branch_mode_[B_POPCONTROL])
    {
      //RealType emix=((iParam[B_WARMUPSTEPS]-todo_steps_)<100)?(0.25*v_param_[SBVP::EREF]+0.75*v_param_[SBVP::ENOW]):v_param_[SBVP::EREF];
      //v_param_[SBVP::ETRIAL]=emix+Feedback*(log_n_-std::log(pop_now));
      //v_param_[SBVP::ETRIAL]=v_param_[SBVP::EREF]+Feedback*(log_n_-std::log(pop_now));
      if (branch_mode_[B_KILLNODES])
        v_param_[SBVP::ETRIAL] = (0.00 * v_param_[SBVP::EREF] + 1.0 * v_param_[SBVP::ENOW]) +
            v_param_[SBVP::FEEDBACK] * (log_n_ - std::log(pop_now)) -
            std::log(walker_controller_->get_ensemble_property().LivingFraction) / v_param_[SBVP::TAU];
      else
        v_param_[SBVP::ETRIAL] = v_param_[SBVP::ENOW] + (log_n_ - std::log(pop_now)) / v_param_[SBVP::TAU];
    }
    else
      v_param_[SBVP::ETRIAL] = v_param_[SBVP::EREF];
    --todo_steps_;
    if (todo_steps_ == 0) //warmup is done
    {
      v_param_[SBVP::SIGMA2] = variance_hist_.mean();
      setBranchCutoff(v_param_[SBVP::SIGMA2], walker_controller_->get_target_sigma(), 10, walkers.R.size());
      app_log() << "\n Warmup is completed after " << iParam[B_WARMUPSTEPS] << std::endl;
      if (branch_mode_[B_USETAUEFF])
        app_log() << "\n  TauEff     = " << v_param_[SBVP::TAUEFF]
                  << "\n TauEff/Tau = " << v_param_[SBVP::TAUEFF] / v_param_[SBVP::TAU];
      else
        app_log() << "\n  TauEff proposed   = " << v_param_[SBVP::TAUEFF] * r2_accepted_.result() / r2_proposed_.result();
      app_log() << "\n  Etrial     = " << v_param_[SBVP::ETRIAL] << std::endl;
      app_log() << " Running average of energy = " << energy_hist_.mean() << std::endl;
      app_log() << "                  Variance = " << v_param_[SBVP::SIGMA2] << std::endl;
      app_log() << "branch cutoff = " << v_param_[SBVP::BRANCHCUTOFF] << " " << v_param_[SBVP::BRANCHMAX] << std::endl;
      todo_steps_           = iParam[B_ENERGYUPDATEINTERVAL] - 1;
      iParam[B_WARMUPSTEPS] = 0;
      branch_mode_.set(B_DMCSTAGE, 1); //set BranchModex to main stage
      //reset the histogram
      energy_hist_.clear();
      energy_hist_(v_param_[SBVP::ENOW]);
      if (s_param_[MIXDMCOPT] == "yes")
      {
        app_log() << "Switching to DMC with fluctuating populations" << std::endl;
        branch_mode_.set(B_POPCONTROL, 1); //use standard DMC
        walker_controller_        = std::move(backup_walker_controller_);
        backup_walker_controller_ = 0;
        v_param_[SBVP::ETRIAL]    = v_param_[SBVP::EREF];
        app_log() << "  Etrial     = " << v_param_[SBVP::ETRIAL] << std::endl;
        walker_controller_->start(root_name_);
      }
      //This is not necessary
      //energy_hist_(DMCEnergyHist.mean());
    }
  }
  walker_controller_->setTrialEnergy(v_param_[SBVP::ETRIAL]);
  //normalize collectables
  FullPrecRealType wgt_inv = walker_controller_->get_num_contexts() / walker_controller_->get_ensemble_property().Weight;
  walkers.Collectables *= wgt_inv;
}

/**
 *
 */
void SimpleFixedNodeBranch::collect(int iter, MCWalkerConfiguration& W)
{
  //Update the current energy and accumulate.
  MCWalkerConfiguration::Walker_t& head = W.reptile->getHead();
  MCWalkerConfiguration::Walker_t& tail = W.reptile->getTail();
  v_param_[SBVP::ENOW]                  = 0.5 * (head.Properties(WP::LOCALENERGY) + tail.Properties(WP::LOCALENERGY));
  // app_log()<<"IN SimpleFixedNodeBranch::collect\n";
  // app_log()<<"\tv_param_[SBVP::ENOW]="<<v_param_[SBVP::ENOW]<< std::endl;
  energy_hist_(v_param_[SBVP::ENOW]);
  v_param_[SBVP::EREF] = energy_hist_.mean();
  // app_log()<<"\tv_param_[SBVP::EREF]="<<v_param_[SBVP::EREF]<< std::endl;
  //Update the energy variance and R2 for effective timestep and filtering.
  variance_hist_(std::pow(v_param_[SBVP::ENOW] - v_param_[SBVP::EREF], 2));
  r2_accepted_(head.Properties(WP::R2ACCEPTED));
  r2_proposed_(head.Properties(WP::R2PROPOSED));
  // app_log()<<"\thead.Properties(R2ACCEPTED)="<<head.Properties(R2ACCEPTED)<< std::endl;
  // app_log()<<"\thead.Properties(R2PROPOSED)="<<head.Properties(R2PROPOSED)<< std::endl;
  //  app_log()<<"\tr2_accepted_="<<r2_accepted_.result()<< std::endl;
  // app_log()<<"\tr2_proposed_="<<r2_proposed_.result()<< std::endl;
  //  app_log()<<"\tR2Accept/R2Prop="<<r2_accepted_.result()/r2_proposed_.result()<< std::endl;
  // app_log()<<"\t <E^2> = "<<variance_hist_.mean()<< std::endl;
  // app_log()<<"\t <E>   = "<<energy_hist_.mean()<< std::endl;
  //  app_log()<<"\t <E>^2 = "<<std::pow(energy_hist_.mean(),2)<< std::endl;
  // app_log()<<"\t var = "<<variance_hist_.mean()-pow(energy_hist_.mean(),2)<< std::endl;
  // app_log()<<"--------------\n";
  //current mean
  if (branch_mode_[B_USETAUEFF])
  {
    //app_log()<<" BRANCHMODE = "<<branch_mode_[B_USETAUEFF]<< std::endl;
    v_param_[SBVP::TAUEFF] = v_param_[SBVP::TAU] * r2_accepted_.result() / r2_proposed_.result();
    //  app_log()<<"\tv_param_[SBVP::TAU]="<<v_param_[SBVP::TAU]<<" "<<v_param_[SBVP::TAUEFF]<< std::endl;
  }
  /*
  if(branch_mode_[B_RMCSTAGE]) // main stage
  {
    if(branch_mode_[B_POPCONTROL])
    {
      if(todo_steps_>0)
        --todo_steps_;
      else
      {
        v_param_[SBVP::ETRIAL]=v_param_[SBVP::EREF]+v_param_[SBVP::FEEDBACK]*(log_n_-std::log(pop_now));
        todo_steps_=iParam[B_ENERGYUPDATEINTERVAL]-1;
      }
    }
    else
      v_param_[SBVP::ETRIAL]=v_param_[SBVP::EREF];
  }*/
  //app_log()<<"branch_mode_[B_RMCSTAGE]="<<branch_mode_[B_RMCSTAGE]<< std::endl;
  if (!branch_mode_[B_RMCSTAGE]) //warmup
  {
    if (branch_mode_[B_USETAUEFF])
      v_param_[SBVP::TAUEFF] = v_param_[SBVP::TAU] * r2_accepted_.result() / r2_proposed_.result();
    // app_log()<<"\t <E^2> = "<<variance_hist_.mean()<< std::endl;
    // app_log()<<"\t <E>   = "<<energy_hist_.mean()<< std::endl;
    // app_log()<<"\t <E>^2 = "<<std::pow(energy_hist_.mean(),2)<< std::endl;
    //app_log()<<"\t var = "<<variance_hist_.mean()-std::pow(energy_hist_.mean(),2)<< std::endl;
    // app_log()<<"\t var = "<<variance_hist_.mean()<< std::endl;
    //  app_log()<<"----\n";
    //app_log()<<"todo_steps_="<<todo_steps_<< std::endl;
    v_param_[SBVP::ETRIAL] = v_param_[SBVP::EREF];
    --todo_steps_;
    if (todo_steps_ == 0) //warmup is done
    {
      v_param_[SBVP::TAUEFF] = v_param_[SBVP::TAU] * r2_accepted_.result() / r2_proposed_.result();
      v_param_[SBVP::SIGMA2] = variance_hist_.mean();
      setBranchCutoff(v_param_[SBVP::SIGMA2], v_param_[SBVP::FILTERSCALE], v_param_[SBVP::FILTERSCALE], W.R.size());
      app_log() << "\n Warmup is completed after " << iParam[B_WARMUPSTEPS] << " steps." << std::endl;
      if (branch_mode_[B_USETAUEFF])
        app_log() << "\n  TauEff     = " << v_param_[SBVP::TAUEFF]
                  << "\n TauEff/Tau = " << v_param_[SBVP::TAUEFF] / v_param_[SBVP::TAU];
      else
        app_log() << "\n  TauEff proposed   = " << v_param_[SBVP::TAUEFF] * r2_accepted_.result() / r2_proposed_.result();
      app_log() << "\n Running average of energy = " << energy_hist_.mean() << std::endl;
      app_log() << "\n                  Variance = " << v_param_[SBVP::SIGMA2] << std::endl;
      app_log() << "\nbranch cutoff = " << v_param_[SBVP::BRANCHCUTOFF] << " " << v_param_[SBVP::BRANCHMAX] << std::endl;
      todo_steps_           = iParam[B_ENERGYUPDATEINTERVAL] - 1;
      iParam[B_WARMUPSTEPS] = 0;
      branch_mode_.set(B_RMCSTAGE, 1); //set BranchModex to main stage
      //reset the histogram
      energy_hist_.clear();
      energy_hist_(v_param_[SBVP::ENOW]);
    }
  }
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
  //branch_mode_.set(B_DMC,1);//set DMC
  //branch_mode_.set(B_DMCSTAGE,0);//set warmup
  if (walker_controller_)
  {
    //this is to compare the time step errors
    branch_mode_.set(B_USETAUEFF, s_param_[USETAUOPT] == "no");
    if (branch_mode_[B_DMCSTAGE]) //
      todo_steps_ = iParam[B_ENERGYUPDATEINTERVAL] - 1;
    else
      todo_steps_ = iParam[B_WARMUPSTEPS];
    if (branch_mode_[B_POPCONTROL])
    {
      //log_n_ = Feedback*std::log(static_cast<RealType>(iParam[B_TARGETWALKERS]));
      log_n_ = std::log(static_cast<FullPrecRealType>(iParam[B_TARGETWALKERS]));
      if (v_param_[SBVP::FEEDBACK] == 0.0)
        v_param_[SBVP::FEEDBACK] = 1.0;
    }
    else
    {
      //may set Eref to a safe value
      //if(EnergyHistory.count()<5) Eref -= v_param_[EnergyWindowIndex];
      v_param_[SBVP::ETRIAL]   = v_param_[SBVP::EREF];
      v_param_[SBVP::FEEDBACK] = 0.0;
      log_n_                   = 0.0;
    }
  }
  if (branch_mode_[B_RMC])
  {
    //this is to compare the time step errors
    // branch_mode_.set(B_USETAUEFF,s_param_[USETAUOPT]=="no");
    if (branch_mode_[B_RMCSTAGE]) //
      todo_steps_ = iParam[B_ENERGYUPDATEINTERVAL] - 1;
    else
      todo_steps_ = iParam[B_WARMUPSTEPS];
  }
}

int SimpleFixedNodeBranch::resetRun(xmlNodePtr cur)
{
  app_log() << "BRANCH resetRun" << std::endl;
  std::bitset<B_MODE_MAX> bmode(branch_mode_);
  IParamType iparam_old(iParam);
  VParamType vparam_old(v_param_);
  my_node_ = cur;
  //store old target
  int nw_target = iParam[B_TARGETWALKERS];
  m_param_.put(cur);

  int target_min = -1;
  ParameterSet p;
  p.add(target_min, "minimumtargetwalkers"); //p.add(target_min,"minimumTargetWalkers");
  p.add(target_min, "minimumsamples");       //p.add(target_min,"minimumSamples");
  p.put(cur);

  if (iParam[B_TARGETWALKERS] < target_min)
  {
    iParam[B_TARGETWALKERS] = target_min;
  }

  bool same_wc = true;
  if (branch_mode_[B_DMC] && walker_controller_)
  {
    std::string reconfig("no");
    // method is actually IndexType so conceivably indicates much more that reconfig="yes" or "no"
    if (walker_controller_->get_method())
      reconfig = "runwhileincorrect"; // forces SR during warmup?
    std::string reconfig_prev(reconfig);
    ParameterSet p;
    p.add(reconfig, "reconfiguration");
    p.put(cur);
    if (reconfig != "no" && reconfig != "runwhileincorrect")
    {
      // remove this once bug is fixed
      throw std::runtime_error("Reconfiguration is currently broken and gives incorrect results. Use dynamic "
                               "population control by setting reconfiguration=\"no\" or removing the reconfiguration "
                               "option from the DMC input section. If accessing the broken reconfiguration code path "
                               "is still desired, set reconfiguration to \"runwhileincorrect\" instead of \"yes\".");
    }
    same_wc = (reconfig == reconfig_prev);
  }

  //everything is the same, do nothing
  if (same_wc && bmode == branch_mode_ && std::equal(iParam.begin(), iParam.end(), iparam_old.begin()) &&
      std::equal(v_param_.begin(), v_param_.end(), vparam_old.begin()))
  {
    app_log() << "  Continue with the same input as the previous block." << std::endl;
    app_log().flush();
    //return 1;
  }
  app_log() << " SimpleFixedNodeBranch::resetRun detected changes in <parameter>'s " << std::endl;
  app_log() << " BranchMode : " << bmode << " " << branch_mode_ << std::endl;

  //vmc does not need to do anything with WalkerController
  if (!branch_mode_[B_DMC])
  {
    app_log() << " iParam (old): " << iparam_old << std::endl;
    app_log() << " iParam (new): " << iParam << std::endl;
    app_log() << " vParam (old): " << vparam_old << std::endl;
    app_log() << " vParam (new): " << v_param_ << std::endl;
    app_log().flush();
    return 1;
  }

  if (walker_controller_ == nullptr)
  {
    APP_ABORT("SimpleFixedNodeBranch::resetRun cannot initialize WalkerController");
  }

  if (!same_wc)
  {
    app_log() << "Destroy WalkerController. Existing method " << walker_controller_->get_method() << std::endl;
    ;
    walker_controller_.reset(createWalkerController(iParam[B_TARGETWALKERS], *my_estimator_->getCommunicator(), my_node_));
    app_log().flush();

    branch_mode_[B_POPCONTROL] = (walker_controller_->get_method() == 0);
    if (branch_mode_[B_POPCONTROL])
    {
      v_param_[SBVP::ETRIAL] = v_param_[SBVP::EREF];
      if (v_param_[SBVP::FEEDBACK] == 0.0)
        v_param_[SBVP::FEEDBACK] = 1.0;
    }
  }

  //always add a warmup step using default 10 steps
  r2_accepted_.clear();
  r2_proposed_.clear();
  //r2_accepted_(1.0e-12);
  //r2_proposed_(1.0e-12);
  branch_mode_[B_DMCSTAGE] = 0;
  walker_controller_->put(my_node_);
  todo_steps_ = iParam[B_WARMUPSTEPS] = (iParam[B_WARMUPSTEPS]) ? iParam[B_WARMUPSTEPS] : 10;
  setBranchCutoff(v_param_[SBVP::SIGMA2], walker_controller_->get_target_sigma(), 10);
  walker_controller_->reset();
  if (backup_walker_controller_)
    backup_walker_controller_->reset();

  iParam[B_MAXWALKERS] = walker_controller_->get_n_max();
  iParam[B_MINWALKERS] = walker_controller_->get_n_min();

  app_log() << " iParam (old): " << iparam_old << std::endl;
  app_log() << " iParam (new): " << iParam << std::endl;
  app_log() << " vParam (old): " << vparam_old << std::endl;
  app_log() << " vParam (new): " << v_param_ << std::endl;

  app_log() << std::endl << " Using branching cutoff scheme " << branching_cutoff_scheme_ << std::endl;

  app_log().flush();

  //  return static_cast<int>(iParam[B_TARGETWALKERS]*1.01/static_cast<double>(nw_target));
  return static_cast<int>(round(static_cast<double>(iParam[B_TARGETWALKERS] / static_cast<double>(nw_target))));
}

void SimpleFixedNodeBranch::checkParameters(MCWalkerConfiguration& w)
{
  std::ostringstream o;
  if (!branch_mode_[B_DMCSTAGE])
  {
    FullPrecRealType e, sigma2;
    my_estimator_->getApproximateEnergyVariance(e, sigma2);
    v_param_[SBVP::ETRIAL] = v_param_[SBVP::EREF] = e;
    v_param_[SBVP::SIGMA2]                      = sigma2;
    energy_hist_.clear();
    variance_hist_.clear();
    //DMCEnergyHist.clear();
    energy_hist_(v_param_[SBVP::EREF]);
    variance_hist_(v_param_[SBVP::SIGMA2]);
    //DMCEnergyHist(v_param_[SBVP::EREF]);
    o << "SimpleFixedNodeBranch::checkParameters " << std::endl;
    o << "  Average Energy of a population  = " << e << std::endl;
    o << "  Energy Variance = " << v_param_[SBVP::SIGMA2] << std::endl;
  }
  app_log() << o.str() << std::endl;
  app_log().flush();
}

void SimpleFixedNodeBranch::finalize(MCWalkerConfiguration& w)
{
  std::ostringstream o;
  if (walker_controller_)
  {
    o << "====================================================";
    o << "\n  SimpleFixedNodeBranch::finalize after a DMC block";
    o << "\n    QMC counter                   = " << iParam[B_COUNTER];
    o << "\n    time step                     = " << v_param_[SBVP::TAU];
    o << "\n    effective time step           = " << v_param_[SBVP::TAUEFF];
    o << "\n    trial energy                  = " << v_param_[SBVP::ETRIAL];
    o << "\n    reference energy              = " << v_param_[SBVP::EREF];
    o << "\n    reference variance            = " << v_param_[SBVP::SIGMA2];
    o << "\n    target walkers                = " << iParam[B_TARGETWALKERS];
    o << "\n    branch cutoff                 = " << v_param_[SBVP::BRANCHCUTOFF] << " " << v_param_[SBVP::BRANCHMAX];
    o << "\n    Max and minimum walkers per node= " << iParam[B_MAXWALKERS] << " " << iParam[B_MINWALKERS];
    o << "\n    Feedback                      = " << v_param_[SBVP::FEEDBACK];
    o << "\n    QMC Status (BranchMode)       = " << branch_mode_;
    o << "\n====================================================";
  }
  //running RMC
  else if (branch_mode_[B_RMC])
  {
    o << "====================================================";
    o << "\n  SimpleFixedNodeBranch::finalize after a RMC block";
    o << "\n    QMC counter                   = " << iParam[B_COUNTER];
    o << "\n    time step                     = " << v_param_[SBVP::TAU];
    o << "\n    effective time step           = " << v_param_[SBVP::TAUEFF];
    o << "\n    reference energy              = " << v_param_[SBVP::EREF];
    o << "\n    reference variance            = " << v_param_[SBVP::SIGMA2];
    o << "\n    cutoff energy                 = " << v_param_[SBVP::BRANCHCUTOFF] << " " << v_param_[SBVP::BRANCHMAX];
    o << "\n    QMC Status (BranchMode)       = " << branch_mode_;
    o << "\n====================================================";
  }
  else
  {
    //running VMC
    FullPrecRealType e, sigma2;
    my_estimator_->getApproximateEnergyVariance(e, sigma2);
    v_param_[SBVP::ETRIAL] = v_param_[SBVP::EREF] = e;
    v_param_[SBVP::SIGMA2]                      = sigma2;
    //this is just to avoid diving by n-1  == 0
    energy_hist_(v_param_[SBVP::EREF]);
    //add Eref to the DMCEnergyHistory
    //DMCEnergyHist(v_param_[SBVP::EREF]);
    o << "====================================================";
    o << "\n  SimpleFixedNodeBranch::finalize after a VMC block";
    o << "\n    QMC counter        = " << iParam[B_COUNTER];
    o << "\n    time step          = " << v_param_[SBVP::TAU];
    o << "\n    reference energy   = " << v_param_[SBVP::EREF];
    o << "\n    reference variance = " << v_param_[SBVP::SIGMA2];
    o << "\n====================================================";
  }
  app_log() << o.str() << std::endl;
  write(root_name_, true);
}

/**  Parse the xml file for parameters
 *@param cur current xmlNode
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
  my_node_ = cur;
  //check dmc/vmc and decide to create WalkerControllerBase
  m_param_.put(cur);
  reset();
  return true;
}

void SimpleFixedNodeBranch::write(const std::string& fname, bool overwrite)
{
  root_name_ = fname;
  if (my_estimator_->is_manager())
  {
    //\since 2008-06-24
    v_param_[SBVP::ACC_ENERGY]  = energy_hist_.result();
    v_param_[SBVP::ACC_SAMPLES] = energy_hist_.count();
    BranchIO<SimpleFixedNodeBranch> hh(*this, *my_estimator_->getCommunicator());
    bool success = hh.write(fname);
  }
}

void SimpleFixedNodeBranch::read(const std::string& fname)
{
  branch_mode_.set(B_RESTART, 0);
  if (fname.empty())
    return;
  v_param_[SBVP::ACC_ENERGY]  = energy_hist_.result();
  v_param_[SBVP::ACC_SAMPLES] = energy_hist_.count();
  BranchIO<SimpleFixedNodeBranch> hh(*this, *my_estimator_->getCommunicator());
  BranchModeType bmode(branch_mode_);
  bool success = hh.read(fname);
  if (success && r2_proposed_.good() && bmode[B_POPCONTROL] == branch_mode_[B_POPCONTROL])
  {
    branch_mode_.set(B_RESTART, 1);
    app_log() << "    Restarting, cummulative properties:"
              << "\n      energy     = " << energy_hist_.mean() << "\n      variance   = " << variance_hist_.mean()
              << "\n      r2accepted = " << r2_accepted_.mean() << "\n      r2proposed = " << r2_proposed_.mean()
              << std::endl;
  }
  else
  {
    if (branch_mode_[B_POPCONTROL] != bmode[B_POPCONTROL])
    {
      app_log() << "  Population control method has changed from " << branch_mode_[B_POPCONTROL] << " to "
                << bmode[B_POPCONTROL] << std::endl;
      branch_mode_[B_POPCONTROL] = bmode[B_POPCONTROL];
    }
  }

  app_log().flush();
}

void SimpleFixedNodeBranch::setBranchCutoff(FullPrecRealType variance,
                                            FullPrecRealType targetSigma,
                                            FullPrecRealType maxSigma,
                                            int Nelec)
{
  if (branching_cutoff_scheme_ == "DRV")
  {
    // eq.(3), J. Chem. Phys. 89, 3629 (1988).
    // eq.(9), J. Chem. Phys. 99, 2865 (1993).
    v_param_[SBVP::BRANCHCUTOFF] = 2.0 / std::sqrt(v_param_[SBVP::TAU]);
  }
  else if (branching_cutoff_scheme_ == "ZSGMA")
  {
    // eq.(6), Phys. Rev. B 93, 241118(R) (2016)
    // do nothing if Nelec is not passed in.
    if (Nelec > 0)
      v_param_[SBVP::BRANCHCUTOFF] = 0.2 * std::sqrt(Nelec / v_param_[SBVP::TAU]);
  }
  else if (branching_cutoff_scheme_ == "YL")
  {
    // a scheme from Ye Luo.
    v_param_[SBVP::BRANCHCUTOFF] = std::sqrt(variance) * std::min(targetSigma, std::sqrt(1.0 / v_param_[SBVP::TAU]));
  }
  else if (branching_cutoff_scheme_ == "classic")
  {
    // default QMCPACK choice which is the same as v3.0.0 and before.
    v_param_[SBVP::BRANCHCUTOFF] = std::min(std::max(variance * targetSigma, maxSigma), 2.5 / v_param_[SBVP::TAU]);
  }
  else
    APP_ABORT("SimpleFixedNodeBranch::setBranchCutoff unknown branching cutoff scheme " + branching_cutoff_scheme_);

  v_param_[SBVP::BRANCHMAX]    = v_param_[SBVP::BRANCHCUTOFF] * 1.5;
  v_param_[SBVP::BRANCHFILTER] = 1.0 / (v_param_[SBVP::BRANCHMAX] - v_param_[SBVP::BRANCHCUTOFF]);
}

std::ostream& operator<<(std::ostream& os, SimpleFixedNodeBranch::VParamType& rhs)
{
  for (auto value : rhs)
    os << std::setw(18) << std::setprecision(10) << value;
  return os;
}

} // namespace qmcplusplus
