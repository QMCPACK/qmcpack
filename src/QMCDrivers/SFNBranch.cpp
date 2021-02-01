//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2020 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//
// File refactored from: SimpleFixedNodeBranch.cpp
//////////////////////////////////////////////////////////////////////////////////////


#include "SFNBranch.h"
#include <numeric>
#include "OhmmsData/FileUtility.h"
#include "Utilities/RandomGenerator.h"
#include "QMCDrivers/DMC/WalkerControl.h"
#include "Estimators/EstimatorManagerNew.h"
#include "Particle/Reptile.h"
#include "type_traits/template_types.hpp"

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

SFNBranch::SFNBranch(RealType tau, int nideal) : MyEstimator(nullptr), debug_disable_branching_("no")
{
  BranchMode.set(B_DMCSTAGE, 0);     //warmup stage
  BranchMode.set(B_POPCONTROL, 1);   //use standard DMC
  BranchMode.set(B_USETAUEFF, 1);    //use taueff
  BranchMode.set(B_CLEARHISTORY, 0); //clear history and start with the current average
  BranchMode.set(B_KILLNODES, 0);    //when killing walkers at nodes etrial is updated differently
  vParam.fill(1.0);
  vParam[SBVP::TAU]         = tau;
  vParam[SBVP::TAUEFF]      = tau;
  vParam[SBVP::FEEDBACK]    = 1.0;
  vParam[SBVP::FILTERSCALE] = 10;
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
  // No
  //reset();
}

/** copy constructor
 *
 * Copy only selected data members and WalkerController is never copied.
 */
SFNBranch::SFNBranch(const SFNBranch& abranch)
    : BranchMode(abranch.BranchMode),
      iParam(abranch.iParam),
      vParam(abranch.vParam),
      MyEstimator(0),
      branching_cutoff_scheme(abranch.branching_cutoff_scheme),
      sParam(abranch.sParam),
      debug_disable_branching_(abranch.debug_disable_branching_)
{
  registerParameters();
  //would like to remove this
  reset();
}

SFNBranch::~SFNBranch() = default;

void SFNBranch::registerParameters()
{
  m_param.add(iParam[B_WARMUPSTEPS], "warmupSteps", "int");
  m_param.add(iParam[B_WARMUPSTEPS], "warmupsteps", "int");
  m_param.add(iParam[B_ENERGYUPDATEINTERVAL], "energyUpdateInterval", "int");
  m_param.add(iParam[B_BRANCHINTERVAL], "branchInterval", "int");
  m_param.add(iParam[B_TARGETWALKERS], "targetWalkers", "int");
  m_param.add(iParam[B_TARGETWALKERS], "targetwalkers", "int");
  m_param.add(iParam[B_TARGETWALKERS], "target_walkers", "int");
  //trial energy
  m_param.add(vParam[SBVP::EREF], "refEnergy", "AU");
  m_param.add(vParam[SBVP::EREF], "ref_energy", "AU");
  m_param.add(vParam[SBVP::EREF], "en_ref", "AU");
  m_param.add(vParam[SBVP::TAU], "tau", "AU");
  m_param.add(vParam[SBVP::TAU], "timestep", "AU");
  m_param.add(vParam[SBVP::TAU], "timeStep", "AU");
  m_param.add(vParam[SBVP::TAU], "TimeStep", "AU");
  //filterscale:  sets the filtercutoff to sigma*filterscale
  m_param.add(vParam[SBVP::FILTERSCALE], "filterscale", "double");
  //feed back parameter for population control
  m_param.add(vParam[SBVP::FEEDBACK], "feedback", "double");
  //turn on/off effective tau onl for time-step error comparisons
  m_param.add(sParam[USETAUOPT], "useBareTau", "option");
  m_param.add(sParam[MIXDMCOPT], "warmupByReconfiguration", "opt");
  m_param.add(branching_cutoff_scheme, "branching_cutoff_scheme", "option");
  m_param.add(debug_disable_branching_, "debug_disable_branching", "option");
}

int SFNBranch::initWalkerController(MCPopulation& population, bool fixW, bool killwalker)
{
  BranchMode.set(B_DMC, 1);                               //set DMC
  BranchMode.set(B_DMCSTAGE, iParam[B_WARMUPSTEPS] == 0); //use warmup
  bool fromscratch     = false;
  FullPrecRealType tau = vParam[SBVP::TAU];

  int nwtot_now = population.get_num_global_walkers();

  if (WalkerController != nullptr)
    throw std::runtime_error("Unified Driver initWalkerController called with existing WalkerController\n"
                             "this is a violation SFNBranch is created to a valid state only\n"
                             "once in the Unified Driver design.");

  if (iParam[B_TARGETWALKERS] == 0)
  {
    // has "important" side effect of updating the walker offsets
    population.syncWalkersPerNode(MyEstimator->getCommunicator());
    iParam[B_TARGETWALKERS] = population.get_num_global_walkers();
  }
  app_log() << "  Creating WalkerControl" << std::endl;
  WalkerController = std::make_unique<WalkerControl>(MyEstimator->getCommunicator(), Random);
  WalkerController->setMinMax(iParam[B_TARGETWALKERS], 0);
  if (!BranchMode[B_RESTART])
  {
    fromscratch = true;
    app_log() << "  START ALL OVER " << std::endl;
    vParam[SBVP::TAUEFF] = tau;
    BranchMode.set(B_POPCONTROL, !fixW); //fixW -> 0
    BranchMode.set(B_KILLNODES, killwalker);
    iParam[B_MAXWALKERS] = WalkerController->get_n_max();
    iParam[B_MINWALKERS] = WalkerController->get_n_min();
    if (!fixW && sParam[MIXDMCOPT] == "yes")
    {
      app_log() << "Warmup DMC is done with a fixed population " << iParam[B_TARGETWALKERS] << std::endl;
      BackupWalkerController = std::move(WalkerController); //save the main controller
      WalkerController       = std::make_unique<WalkerControl>(MyEstimator->getCommunicator(), Random);
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

  //update the simulation parameters
  WalkerController->put(myNode);
  //assign current Eref and a large number for variance
  WalkerController->setTrialEnergy(vParam[SBVP::ETRIAL]);
  this->reset();
  int allow_flux = 50; // std::min(static_cast<int>(population.get_num_global_walkers() * 0.10), 50);
  if (fromscratch)
  {
    //determine the branch cutoff to limit wild weights based on the sigma and sigmaBound
    // Was check MCWC's particle set for number of R which I take to mean number of particles
    // will this assumption change if various spin freedoms also are added to ParticleSet?
    setBranchCutoff(vParam[SBVP::SIGMA2], WalkerController->get_target_sigma(), allow_flux,
                    population.get_num_particles());
    vParam[SBVP::TAUEFF] = tau * R2Accepted.result() / R2Proposed.result();
  }

  app_log() << "  QMC counter      = " << iParam[B_COUNTER] << std::endl;
  app_log() << "  time step        = " << vParam[SBVP::TAU] << std::endl;
  app_log() << "  effective time step = " << vParam[SBVP::TAUEFF] << std::endl;
  app_log() << "  trial energy     = " << vParam[SBVP::ETRIAL] << std::endl;
  app_log() << "  reference energy = " << vParam[SBVP::EREF] << std::endl;
  app_log() << "  Feedback = " << vParam[SBVP::FEEDBACK] << std::endl;
  app_log() << "  reference variance = " << vParam[SBVP::SIGMA2] << std::endl;
  app_log() << "  target walkers = " << iParam[B_TARGETWALKERS] << std::endl;
  app_log() << "  branching cutoff scheme " << branching_cutoff_scheme << std::endl;
  app_log() << "  branch cutoff = " << vParam[SBVP::BRANCHCUTOFF] << " " << vParam[SBVP::BRANCHMAX] << std::endl;
  app_log() << "  Max and minimum walkers per node= " << iParam[B_MAXWALKERS] << " " << iParam[B_MINWALKERS]
            << std::endl;
  app_log() << "  QMC Status (BranchMode) = " << BranchMode << std::endl;

  return int(round(double(iParam[B_TARGETWALKERS]) / double(nwtot_now)));
}

void SFNBranch::branch(int iter, MCPopulation& population)
{
  //collect the total weights and redistribute the walkers accordingly, using a fixed tolerance
  //RealType pop_now= WalkerController->branch(iter,walkers,0.1);
  RefVector<MCPWalker> walkers(convertUPtrToRefVector(population.get_walkers()));

  FullPrecRealType pop_now;
  pop_now = WalkerController->branch(iter, population, iter == 0);

  //population for trial energy modification should not include any released node walkers.
  MCDataType<FullPrecRealType>& wc_ensemble_prop = WalkerController->get_ensemble_property();
  pop_now -= wc_ensemble_prop.RNSamples;
  //current energy
  vParam[SBVP::ENOW] = wc_ensemble_prop.Energy;
  VarianceHist(wc_ensemble_prop.Variance);
  R2Accepted(wc_ensemble_prop.R2Accepted);
  R2Proposed(wc_ensemble_prop.R2Proposed);
  //PopHist(pop_now);
  vParam[SBVP::EREF] = EnergyHist.mean(); //current mean
  if (BranchMode[B_USETAUEFF])
    vParam[SBVP::TAUEFF] = vParam[SBVP::TAU] * R2Accepted.result() / R2Proposed.result();

  if (BranchMode[B_KILLNODES])
    EnergyHist(vParam[SBVP::ENOW] - std::log(wc_ensemble_prop.LivingFraction) / vParam[SBVP::TAUEFF]);
  else
    EnergyHist(vParam[SBVP::ENOW]);

  if (BranchMode[B_DMCSTAGE]) // main stage
  {
    if (BranchMode[B_POPCONTROL])
    {
      if (ToDoSteps > 0)
        --ToDoSteps;
      else
      {
        vParam[SBVP::ETRIAL] = vParam[SBVP::EREF] + vParam[SBVP::FEEDBACK] * (logN - std::log(pop_now));
        ToDoSteps            = iParam[B_ENERGYUPDATEINTERVAL] - 1;
      }
    }
    else
      vParam[SBVP::ETRIAL] = vParam[SBVP::EREF];
  }
  else //warmup
  {
    if (BranchMode[B_USETAUEFF])
      vParam[SBVP::TAUEFF] = vParam[SBVP::TAU] * R2Accepted.result() / R2Proposed.result();
    if (BranchMode[B_POPCONTROL])
    {
      //RealType emix=((iParam[B_WARMUPSTEPS]-ToDoSteps)<100)?(0.25*vParam[SBVP::EREF]+0.75*vParam[SBVP::ENOW]):vParam[SBVP::EREF];
      //vParam[SBVP::ETRIAL]=emix+Feedback*(logN-std::log(pop_now));
      //vParam[SBVP::ETRIAL]=vParam[SBVP::EREF]+Feedback*(logN-std::log(pop_now));
      if (BranchMode[B_KILLNODES])
        vParam[SBVP::ETRIAL] = (0.00 * vParam[SBVP::EREF] + 1.0 * vParam[SBVP::ENOW]) +
            vParam[SBVP::FEEDBACK] * (logN - std::log(pop_now)) -
            std::log(wc_ensemble_prop.LivingFraction) / vParam[SBVP::TAU];
      else
        vParam[SBVP::ETRIAL] = vParam[SBVP::ENOW] + (logN - std::log(pop_now)) / vParam[SBVP::TAU];
    }
    else
      vParam[SBVP::ETRIAL] = vParam[SBVP::EREF];
    --ToDoSteps;

    if (ToDoSteps == 0) //warmup is done
    {
      vParam[SBVP::SIGMA2] = VarianceHist.mean();
      setBranchCutoff(vParam[SBVP::SIGMA2], WalkerController->get_target_sigma(), 10, population.get_num_particles());
      app_log() << "\n Warmup is completed after " << iParam[B_WARMUPSTEPS] << std::endl;
      if (BranchMode[B_USETAUEFF])
        app_log() << "\n  TauEff     = " << vParam[SBVP::TAUEFF]
                  << "\n TauEff/Tau = " << vParam[SBVP::TAUEFF] / vParam[SBVP::TAU];
      else
        app_log() << "\n  TauEff proposed   = " << vParam[SBVP::TAUEFF] * R2Accepted.result() / R2Proposed.result();
      app_log() << "\n  Etrial     = " << vParam[SBVP::ETRIAL] << std::endl;
      app_log() << " Running average of energy = " << EnergyHist.mean() << std::endl;
      app_log() << "                  Variance = " << vParam[SBVP::SIGMA2] << std::endl;
      app_log() << "branch cutoff = " << vParam[SBVP::BRANCHCUTOFF] << " " << vParam[SBVP::BRANCHMAX] << std::endl;
      ToDoSteps             = iParam[B_ENERGYUPDATEINTERVAL] - 1;
      iParam[B_WARMUPSTEPS] = 0;
      BranchMode.set(B_DMCSTAGE, 1); //set BranchModex to main stage
      //reset the histogram
      EnergyHist.clear();
      EnergyHist(vParam[SBVP::ENOW]);
      if (sParam[MIXDMCOPT] == "yes")
      {
        app_log() << "Switching to DMC with fluctuating populations" << std::endl;
        BranchMode.set(B_POPCONTROL, 1); //use standard DMC
        WalkerController       = std::move(BackupWalkerController);
        BackupWalkerController = 0;
        vParam[SBVP::ETRIAL]   = vParam[SBVP::EREF];
        app_log() << "  Etrial     = " << vParam[SBVP::ETRIAL] << std::endl;
        WalkerController->start();
      }
      //This is not necessary
      //EnergyHist(DMCEnergyHist.mean());
    }
  }
  WalkerController->setTrialEnergy(vParam[SBVP::ETRIAL]);
  //accumulate collectables and energies for scalar.dat
  FullPrecRealType wgt_inv = WalkerController->get_num_contexts() / wc_ensemble_prop.Weight;
  //walkers.Collectables *= wgt_inv;
}

void SFNBranch::printStatus() const
{
  std::ostringstream o;
  if (WalkerController)
  {
    o << "====================================================";
    o << "\n  End of a DMC block";
    o << "\n    QMC counter                   = " << iParam[B_COUNTER];
    o << "\n    time step                     = " << vParam[SBVP::TAU];
    o << "\n    effective time step           = " << vParam[SBVP::TAUEFF];
    o << "\n    trial energy                  = " << vParam[SBVP::ETRIAL];
    o << "\n    reference energy              = " << vParam[SBVP::EREF];
    o << "\n    reference variance            = " << vParam[SBVP::SIGMA2];
    o << "\n    target walkers                = " << iParam[B_TARGETWALKERS];
    o << "\n    branch cutoff                 = " << vParam[SBVP::BRANCHCUTOFF] << " " << vParam[SBVP::BRANCHMAX];
    o << "\n    Max and min walkers per node  = " << iParam[B_MAXWALKERS] << " " << iParam[B_MINWALKERS];
    o << "\n    Feedback                      = " << vParam[SBVP::FEEDBACK];
    o << "\n    QMC Status (BranchMode)       = " << BranchMode;
    o << "\n====================================================";
  }
  //running RMC
  else if (BranchMode[B_RMC])
  {
    o << "====================================================";
    o << "\n  End of a RMC block";
    o << "\n    QMC counter                   = " << iParam[B_COUNTER];
    o << "\n    time step                     = " << vParam[SBVP::TAU];
    o << "\n    effective time step           = " << vParam[SBVP::TAUEFF];
    o << "\n    reference energy              = " << vParam[SBVP::EREF];
    o << "\n    reference variance            = " << vParam[SBVP::SIGMA2];
    o << "\n    cutoff energy                 = " << vParam[SBVP::BRANCHCUTOFF] << " " << vParam[SBVP::BRANCHMAX];
    o << "\n    QMC Status (BranchMode)       = " << BranchMode;
    o << "\n====================================================";
  }
  app_log() << o.str() << std::endl;
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
bool SFNBranch::put(xmlNodePtr cur)
{
  //save it
  myNode = cur;
  //check dmc/vmc and decide to create WalkerControllerBase
  m_param.put(cur);
  reset();
  return true;
}

/** Calculates and saves various action components, also does necessary updates for running averages.
 *
 */
void SFNBranch::reset()
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
      if (vParam[SBVP::FEEDBACK] == 0.0)
        vParam[SBVP::FEEDBACK] = 1.0;
    }
    else
    {
      //may set Eref to a safe value
      //if(EnergyHistory.count()<5) Eref -= vParam[EnergyWindowIndex];
      vParam[SBVP::ETRIAL]   = vParam[SBVP::EREF];
      vParam[SBVP::FEEDBACK] = 0.0;
      logN                   = 0.0;
    }
    //       vParam(abranch.vParam)
    WalkerController->start();
  }
  else
    std::cerr << "Calling reset with no WalkerController and therefore nothing to do. Why?\n";
}

void SFNBranch::setEnergyVariance(FullPrecRealType energy, FullPrecRealType variance)
{
  vParam[SBVP::ETRIAL] = vParam[SBVP::EREF] = energy;
  vParam[SBVP::SIGMA2]                      = variance;
}

void SFNBranch::setBranchCutoff(FullPrecRealType variance,
                                FullPrecRealType targetSigma,
                                FullPrecRealType maxSigma,
                                int Nelec)
{
  if (branching_cutoff_scheme == "DRV")
  {
    // eq.(3), J. Chem. Phys. 89, 3629 (1988).
    // eq.(9), J. Chem. Phys. 99, 2865 (1993).
    vParam[SBVP::BRANCHCUTOFF] = 2.0 / std::sqrt(vParam[SBVP::TAU]);
  }
  else if (branching_cutoff_scheme == "ZSGMA")
  {
    // eq.(6), Phys. Rev. B 93, 241118(R) (2016)
    // do nothing if Nelec is not passed in.
    if (Nelec > 0)
      vParam[SBVP::BRANCHCUTOFF] = 0.2 * std::sqrt(Nelec / vParam[SBVP::TAU]);
  }
  else if (branching_cutoff_scheme == "YL")
  {
    // a scheme from Ye Luo.
    vParam[SBVP::BRANCHCUTOFF] = std::sqrt(variance) * std::min(targetSigma, std::sqrt(1.0 / vParam[SBVP::TAU]));
  }
  else if (branching_cutoff_scheme == "classic")
  {
    // default QMCPACK choice which is the same as v3.0.0 and before.
    vParam[SBVP::BRANCHCUTOFF] = std::min(std::max(variance * targetSigma, maxSigma), 2.5 / vParam[SBVP::TAU]);
  }
  else
    APP_ABORT("SFNBranch::setBranchCutoff unknown branching cutoff scheme " + branching_cutoff_scheme);

  vParam[SBVP::BRANCHMAX]    = vParam[SBVP::BRANCHCUTOFF] * 1.5;
  vParam[SBVP::BRANCHFILTER] = 1.0 / (vParam[SBVP::BRANCHMAX] - vParam[SBVP::BRANCHCUTOFF]);
}

std::ostream& operator<<(std::ostream& os, SFNBranch::VParamType& rhs)
{
  for (auto value : rhs)
    os << std::setw(18) << std::setprecision(10) << value;
  return os;
}

} // namespace qmcplusplus
