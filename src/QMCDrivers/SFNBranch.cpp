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

SFNBranch::SFNBranch(RealType tau, RealType feedback, DMCRefEnergyScheme refenergy_update_scheme)
    : WarmUpToDoSteps(0),
      EtrialUpdateToDoSteps(0),
      myNode(NULL),
      ref_energy_collector(refenergy_update_scheme, std::max(1, static_cast<int>(1.0 / (feedback * tau))))
{
  BranchMode.set(B_DMCSTAGE, 0);     //warmup stage
  BranchMode.set(B_POPCONTROL, 1);   //use standard DMC
  BranchMode.set(B_USETAUEFF, 1);    //use taueff
  BranchMode.set(B_CLEARHISTORY, 0); //clear history and start with the current average
  BranchMode.set(B_KILLNODES, 0);    //when killing walkers at nodes etrial is updated differently
  vParam.fill(1.0);
  vParam[SBVP::TAU]         = tau;
  vParam[SBVP::TAUEFF]      = tau;
  vParam[SBVP::FEEDBACK]    = feedback;
  vParam[SBVP::FILTERSCALE] = 10;
  vParam[SBVP::SIGMA_BOUND] = 10;
  R2Accepted(1.0e-10);
  R2Proposed(1.0e-10);
  //set the default values for integer parameters
  iParam[B_WARMUPSTEPS]          = 200;
  iParam[B_ENERGYUPDATEINTERVAL] = 1;
  iParam[B_BRANCHINTERVAL]       = 1;
  iParam[B_TARGETWALKERS]        = 0;
  iParam[B_COUNTER]              = -1;
  //default is no
  sParam.resize(DUMMYOPT, "no");
  //default is classic
  branching_cutoff_scheme = "classic";
  registerParameters();
}

SFNBranch::~SFNBranch() = default;

void SFNBranch::registerParameters()
{
  m_param.add(iParam[B_WARMUPSTEPS], "warmupSteps");
  m_param.add(iParam[B_WARMUPSTEPS], "warmupsteps");
  m_param.add(iParam[B_ENERGYUPDATEINTERVAL], "energyUpdateInterval");
  m_param.add(iParam[B_BRANCHINTERVAL], "branchInterval");
  m_param.add(iParam[B_TARGETWALKERS], "targetWalkers");
  m_param.add(iParam[B_TARGETWALKERS], "targetwalkers");
  m_param.add(iParam[B_TARGETWALKERS], "target_walkers");
  //trial energy
  m_param.add(vParam[SBVP::EREF], "refEnergy");
  m_param.add(vParam[SBVP::EREF], "ref_energy");
  m_param.add(vParam[SBVP::EREF], "en_ref");
  m_param.add(vParam[SBVP::TAU], "tau");
  m_param.add(vParam[SBVP::TAU], "timestep");
  m_param.add(vParam[SBVP::TAU], "timeStep");
  m_param.add(vParam[SBVP::TAU], "TimeStep");
  //filterscale:  sets the filtercutoff to sigma*filterscale
  m_param.add(vParam[SBVP::FILTERSCALE], "filterscale");
  m_param.add(vParam[SBVP::SIGMA_BOUND], "sigmaBound");
  //turn on/off effective tau onl for time-step error comparisons
  m_param.add(sParam[USETAUOPT], "useBareTau");
  m_param.add(sParam[MIXDMCOPT], "warmupByReconfiguration");
  m_param.add(branching_cutoff_scheme, "branching_cutoff_scheme");
}

int SFNBranch::initParam(const MCPopulation& population,
                         FullPrecRealType ene,
                         FullPrecRealType var,
                         bool fixW,
                         bool killwalker)
{
  app_log() << "  START ALL OVER " << std::endl;
  BranchMode.set(B_POPCONTROL, !fixW); //fixW -> 0
  BranchMode.set(B_KILLNODES, killwalker);

  BranchMode.set(B_DMC, 1);                               //set DMC
  BranchMode.set(B_DMCSTAGE, iParam[B_WARMUPSTEPS] == 0); //use warmup
  vParam[SBVP::ETRIAL] = ene;
  vParam[SBVP::EREF]   = ene;
  vParam[SBVP::SIGMA2] = var;
  vParam[SBVP::TAUEFF] = vParam[SBVP::TAU] * R2Accepted.result() / R2Proposed.result();
  /// FIXME, magic number 50
  setBranchCutoff(vParam[SBVP::SIGMA2], vParam[SBVP::SIGMA_BOUND], 50, population.get_golden_electrons().getTotalNum());

  int nwtot_now = population.get_num_global_walkers();
  if (iParam[B_TARGETWALKERS] == 0)
    iParam[B_TARGETWALKERS] = nwtot_now;

  app_log() << "  QMC counter             = " << iParam[B_COUNTER] << std::endl;
  app_log() << "  time step               = " << vParam[SBVP::TAU] << std::endl;
  app_log() << "  effective time step     = " << vParam[SBVP::TAUEFF] << std::endl;
  app_log() << "  trial energy            = " << vParam[SBVP::ETRIAL] << std::endl;
  app_log() << "  reference energy        = " << vParam[SBVP::EREF] << std::endl;
  app_log() << "  reference variance      = " << vParam[SBVP::SIGMA2] << std::endl;
  app_log() << "  Feedback                = " << vParam[SBVP::FEEDBACK] << std::endl;
  app_log() << "  target walkers          = " << iParam[B_TARGETWALKERS] << std::endl;
  app_log() << "  branching cutoff scheme = " << branching_cutoff_scheme << std::endl;
  app_log() << "  branch cutoff, max      = " << vParam[SBVP::BRANCHCUTOFF] << " " << vParam[SBVP::BRANCHMAX]
            << std::endl;
  app_log() << "  QMC Status (BranchMode) = " << BranchMode << std::endl;

  return int(round(double(iParam[B_TARGETWALKERS]) / double(nwtot_now)));
}

void SFNBranch::updateParamAfterPopControl(const MCDataType<FullPrecRealType>& wc_ensemble_prop, int Nelec)
{
  //target weight
  const auto logN = std::log(static_cast<FullPrecRealType>(iParam[B_TARGETWALKERS]));
  //population weight before branching
  const FullPrecRealType pop_weight = wc_ensemble_prop.Weight;
  //current energy
  vParam[SBVP::ENOW] = wc_ensemble_prop.Energy;

  R2Accepted(wc_ensemble_prop.R2Accepted);
  R2Proposed(wc_ensemble_prop.R2Proposed);
  if (BranchMode[B_USETAUEFF])
    vParam[SBVP::TAUEFF] = vParam[SBVP::TAU] * R2Accepted.result() / R2Proposed.result();

  if (BranchMode[B_DMCSTAGE]) // main stage after warmup
  {
    if (WarmUpToDoSteps != 0)
      throw UniformCommunicateError("Bug: WarmUpToDoSteps should be 0 after warmup.");

    // assuming ENOW only fluctuates around the mean (EREF) once warmup completes.
    const auto ene = BranchMode[B_KILLNODES]
        ? vParam[SBVP::ENOW] - std::log(wc_ensemble_prop.LivingFraction) / vParam[SBVP::TAUEFF]
        : vParam[SBVP::ENOW];
    ref_energy_collector.pushWeightEnergyVariance(wc_ensemble_prop.Weight, ene, wc_ensemble_prop.Variance);
    // update the reference energy
    auto [ene_avg, var_avg] = ref_energy_collector.getEnergyVariance();
    vParam[SBVP::EREF]      = ene_avg;

    // update Etrial based on EREF
    if (BranchMode[B_POPCONTROL])
    {
      --EtrialUpdateToDoSteps;
      if (EtrialUpdateToDoSteps == 0)
      {
        vParam[SBVP::ETRIAL]  = vParam[SBVP::EREF] + vParam[SBVP::FEEDBACK] * (logN - std::log(pop_weight));
        EtrialUpdateToDoSteps = iParam[B_ENERGYUPDATEINTERVAL];
      }
    }
    else
    {
      throw UniformCommunicateError("Bug: FIXME SBVP::EREF should be calculated based on weights");
      /// FIXME vParam[SBVP::ETRIAL] = vParam[SBVP::EREF];
    }
  }
  else //warmup
  {
    if (WarmUpToDoSteps == 0)
      throw UniformCommunicateError("Bug: WarmUpToDoSteps should be larger than 0 during warmup.");

    // update Etrial based on ENOW as ENOW is not yet converged in warmup stage
    if (BranchMode[B_POPCONTROL])
    {
      if (BranchMode[B_KILLNODES])
        vParam[SBVP::ETRIAL] = (0.00 * vParam[SBVP::EREF] + 1.0 * vParam[SBVP::ENOW]) +
            vParam[SBVP::FEEDBACK] * (logN - std::log(pop_weight)) -
            std::log(wc_ensemble_prop.LivingFraction) / vParam[SBVP::TAU];
      else
        vParam[SBVP::ETRIAL] = vParam[SBVP::ENOW] + (logN - std::log(pop_weight)) / vParam[SBVP::TAU];
    }
    else
    {
      throw UniformCommunicateError("Bug: FIXME SBVP::EREF should be calculated based on weights");
      /// FIXME vParam[SBVP::ETRIAL] = vParam[SBVP::ENOW];
    }

    --WarmUpToDoSteps;
    if (WarmUpToDoSteps == 0) //warmup is done
    {
      if (ref_energy_collector.count())
        throw UniformCommunicateError("Bug: ref_energy_collector should not have been used during warmup.");

      vParam[SBVP::SIGMA2] = wc_ensemble_prop.Variance;
      setBranchCutoff(vParam[SBVP::SIGMA2], vParam[SBVP::SIGMA_BOUND], 10, Nelec);
      app_log() << "\n Warmup is completed after " << iParam[B_WARMUPSTEPS] << std::endl;
      if (BranchMode[B_USETAUEFF])
        app_log() << "\n  TauEff     = " << vParam[SBVP::TAUEFF]
                  << "\n TauEff/Tau = " << vParam[SBVP::TAUEFF] / vParam[SBVP::TAU];
      else
        app_log() << "\n  TauEff proposed   = " << vParam[SBVP::TAUEFF] * R2Accepted.result() / R2Proposed.result();
      app_log() << "\n  Etrial     = " << vParam[SBVP::ETRIAL] << std::endl;
      app_log() << " Population average of energy = " << vParam[SBVP::ENOW] << std::endl;
      app_log() << "                  Variance = " << vParam[SBVP::SIGMA2] << std::endl;
      app_log() << "branch cutoff = " << vParam[SBVP::BRANCHCUTOFF] << " " << vParam[SBVP::BRANCHMAX] << std::endl;

      BranchMode.set(B_DMCSTAGE, 1); //set BranchModex to main stage
      if (sParam[MIXDMCOPT] == "yes")
      {
        app_log() << "Switching to DMC with fluctuating populations" << std::endl;
        BranchMode.set(B_POPCONTROL, 1); //use standard DMC
        vParam[SBVP::ETRIAL] = vParam[SBVP::EREF];
        app_log() << "  Etrial     = " << vParam[SBVP::ETRIAL] << std::endl;
      }
    }
  }
}

void SFNBranch::printStatus() const
{
  std::ostringstream o;
  //running RMC
  if (BranchMode[B_RMC])
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
  else // running DMC
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
    o << "\n    Feedback                      = " << vParam[SBVP::FEEDBACK];
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
  m_param.put(cur);
  WarmUpToDoSteps       = iParam[B_WARMUPSTEPS];
  EtrialUpdateToDoSteps = iParam[B_ENERGYUPDATEINTERVAL];
  return true;
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
