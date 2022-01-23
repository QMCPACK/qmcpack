//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Miguel Morales, moralessilva2@llnl.gov, Lawrence Livermore National Laboratory
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//                    Raymond Clay III, j.k.rofling@gmail.com, Lawrence Livermore National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


#include "QMCCostFunctionBatched.h"
#include "QMCDrivers/WFOpt/CostFunctionCrowdData.h"
#include "Particle/MCWalkerConfiguration.h"
#include "QMCWaveFunctions/TrialWaveFunction.h"
#include "Message/CommOperators.h"
#include "QMCDrivers/Optimizers/DescentEngine.h"
#include "Concurrency/ParallelExecutor.hpp"
//#define QMCCOSTFUNCTION_DEBUG

namespace qmcplusplus
{
QMCCostFunctionBatched::QMCCostFunctionBatched(MCWalkerConfiguration& w,
                                               TrialWaveFunction& psi,
                                               QMCHamiltonian& h,
                                               SampleStack& samples,
                                               int num_opt_crowds,
                                               int crowd_size,
                                               Communicate* comm)
    : QMCCostFunctionBase(w, psi, h, comm),
      samples_(samples),
      opt_batch_size_(crowd_size),
      opt_num_crowds_(num_opt_crowds),
      check_config_timer_(
          *timer_manager.createTimer("QMCCostFunctionBatched::checkConfigurations", timer_level_medium)),
      corr_sampling_timer_(
          *timer_manager.createTimer("QMCCostFunctionBatched::correlatedSampling", timer_level_medium)),
      fill_timer_(
          *timer_manager.createTimer("QMCCostFunctionBatched::fillOverlapHamiltonianMatrices", timer_level_medium))

{
  app_log() << " Using QMCCostFunctionBatched::QMCCostFunctionBatched" << std::endl;
}


/** Clean up the vector */
QMCCostFunctionBatched::~QMCCostFunctionBatched() = default;

void QMCCostFunctionBatched::GradCost(std::vector<Return_rt>& PGradient,
                                      const std::vector<Return_rt>& PM,
                                      Return_rt FiniteDiff)
{
  if (FiniteDiff > 0)
  {
    QMCTraits::RealType dh = 1.0 / (2.0 * FiniteDiff);
    for (int i = 0; i < NumOptimizables; i++)
    {
      for (int j = 0; j < NumOptimizables; j++)
        OptVariables[j] = PM[j];
      OptVariables[i]               = PM[i] + FiniteDiff;
      QMCTraits::RealType CostPlus  = this->Cost();
      OptVariables[i]               = PM[i] - FiniteDiff;
      QMCTraits::RealType CostMinus = this->Cost();
      PGradient[i]                  = (CostPlus - CostMinus) * dh;
    }
  }
  else
  {
    for (int j = 0; j < NumOptimizables; j++)
      OptVariables[j] = PM[j];
    resetPsi();
    //evaluate new local energies and derivatives
    NumWalkersEff = correlatedSampling(true);
    //Estimators::accumulate has been called by correlatedSampling
    curAvg_w = SumValue[SUM_E_WGT] / SumValue[SUM_WGT];
    //    Return_t curAvg2_w = curAvg_w*curAvg_w;
    curVar_w = SumValue[SUM_ESQ_WGT] / SumValue[SUM_WGT] - curAvg_w * curAvg_w;
    std::vector<Return_rt> EDtotals(NumOptimizables, 0.0);
    std::vector<Return_rt> EDtotals_w(NumOptimizables, 0.0);
    std::vector<Return_rt> E2Dtotals_w(NumOptimizables, 0.0);
    std::vector<Return_rt> URV(NumOptimizables, 0.0);
    std::vector<Return_rt> HD_avg(NumOptimizables, 0.0);
    Return_rt wgtinv   = 1.0 / SumValue[SUM_WGT];
    Return_rt delE_bar = 0;
    {
      for (int iw = 0; iw < rank_local_num_samples_; iw++)
      {
        const Return_rt* restrict saved = RecordsOnNode_[iw];
        Return_rt weight                = saved[REWEIGHT] * wgtinv;
        Return_rt eloc_new              = saved[ENERGY_NEW];
        delE_bar += weight * std::pow(std::abs(eloc_new - EtargetEff), PowerE);
        const Return_rt* HDsaved = HDerivRecords_[iw];
        for (int pm = 0; pm < NumOptimizables; pm++)
          HD_avg[pm] += HDsaved[pm];
      }
    }
    myComm->allreduce(HD_avg);
    myComm->allreduce(delE_bar);
    for (int pm = 0; pm < NumOptimizables; pm++)
      HD_avg[pm] *= 1.0 / static_cast<Return_rt>(NumSamples);
    {
      for (int iw = 0; iw < rank_local_num_samples_; iw++)
      {
        const Return_rt* restrict saved = RecordsOnNode_[iw];
        Return_rt weight                = saved[REWEIGHT] * wgtinv;
        Return_rt eloc_new              = saved[ENERGY_NEW];
        Return_rt delta_l               = (eloc_new - curAvg_w);
        bool ltz(true);
        if (eloc_new - EtargetEff < 0)
          ltz = false;
        Return_rt delE           = std::pow(std::abs(eloc_new - EtargetEff), PowerE);
        Return_rt ddelE          = PowerE * std::pow(std::abs(eloc_new - EtargetEff), PowerE - 1);
        const Return_rt* Dsaved  = DerivRecords_[iw];
        const Return_rt* HDsaved = HDerivRecords_[iw];
        for (int pm = 0; pm < NumOptimizables; pm++)
        {
          EDtotals_w[pm] += weight * (HDsaved[pm] + 2.0 * Dsaved[pm] * delta_l);
          URV[pm] += 2.0 * (eloc_new * HDsaved[pm] - curAvg * HD_avg[pm]);
          if (ltz)
            EDtotals[pm] += weight * (2.0 * Dsaved[pm] * (delE - delE_bar) + ddelE * HDsaved[pm]);
          else
            EDtotals[pm] += weight * (2.0 * Dsaved[pm] * (delE - delE_bar) - ddelE * HDsaved[pm]);
        }
      }
    }
    myComm->allreduce(EDtotals);
    myComm->allreduce(EDtotals_w);
    myComm->allreduce(URV);
    Return_rt smpinv = 1.0 / static_cast<Return_rt>(NumSamples);
    {
      for (int iw = 0; iw < rank_local_num_samples_; iw++)
      {
        const Return_rt* restrict saved = RecordsOnNode_[iw];
        Return_rt weight                = saved[REWEIGHT] * wgtinv;
        Return_rt eloc_new              = saved[ENERGY_NEW];
        Return_rt delta_l               = (eloc_new - curAvg_w);
        Return_rt sigma_l               = delta_l * delta_l;
        const Return_rt* Dsaved         = DerivRecords_[iw];
        const Return_rt* HDsaved        = HDerivRecords_[iw];
        for (int pm = 0; pm < NumOptimizables; pm++)
        {
          E2Dtotals_w[pm] +=
              weight * 2.0 * (Dsaved[pm] * (sigma_l - curVar_w) + delta_l * (HDsaved[pm] - EDtotals_w[pm]));
        }
      }
    }
    myComm->allreduce(E2Dtotals_w);
    for (int pm = 0; pm < NumOptimizables; pm++)
      URV[pm] *= smpinv;
    for (int j = 0; j < NumOptimizables; j++)
    {
      PGradient[j] = 0.0;
      if (std::abs(w_var) > 1.0e-10)
        PGradient[j] += w_var * E2Dtotals_w[j];
      if (std::abs(w_en) > 1.0e-10)
        PGradient[j] += w_en * EDtotals_w[j];
      if (std::abs(w_w) > 1.0e-10)
        PGradient[j] += w_w * URV[j];
      if (std::abs(w_abs) > 1.0e-10)
        PGradient[j] += w_abs * EDtotals[j];
    }
    IsValid = true;
    if (NumWalkersEff < MinNumWalkers * NumSamples)
    {
      WARNMSG("CostFunction-> Number of Effective Walkers is too small " << NumWalkersEff << "Minimum required"
                                                                         << MinNumWalkers * NumSamples)
      IsValid = false;
    }
  }
}

void QMCCostFunctionBatched::getConfigurations(const std::string& aroot)
{
  app_log() << "  Using Nonlocal PP in Opt: " << includeNonlocalH << std::endl;
  outputManager.pause();

  if (H_KE_node_names_.size() == 0)
  {
    H_KE_node_names_.reserve(2);
    H_KE_node_names_.emplace_back("Kinetic");
    if (includeNonlocalH != "no")
    {
      OperatorBase* a(H.getHamiltonian(includeNonlocalH));
      outputManager.resume();
      if (a)
      {
        app_log() << " Found non-local Hamiltonian element named " << includeNonlocalH << std::endl;
        H_KE_node_names_.emplace_back(includeNonlocalH);
      }
      else
        app_log() << " Did not find non-local Hamiltonian element named " << includeNonlocalH << std::endl;
    }
  }
  outputManager.resume();

  rank_local_num_samples_ = samples_.getNumSamples();

  if (dLogPsi.size() != rank_local_num_samples_)
  {
    delete_iter(dLogPsi.begin(), dLogPsi.end());
    delete_iter(d2LogPsi.begin(), d2LogPsi.end());
    int nptcl = W.getTotalNum();
    dLogPsi.resize(rank_local_num_samples_);
    d2LogPsi.resize(rank_local_num_samples_);
    for (int i = 0; i < rank_local_num_samples_; ++i)
      dLogPsi[i] = new ParticleGradient(nptcl);
    for (int i = 0; i < rank_local_num_samples_; ++i)
      d2LogPsi[i] = new ParticleLaplacian(nptcl);
  }
}

// Input - sample_size  - number of samples to process
//       - batch_size  -  process samples in batch_size at a time
// Output - num_batches - number of batches to use
//        - final_batch_size - the last batch size.  May be smaller than batch_size
//                             if the number of samples is not a multiple of the batch size
void compute_batch_parameters(int sample_size, int batch_size, int& num_batches, int& final_batch_size)
{
  num_batches      = sample_size / batch_size;
  final_batch_size = batch_size;
  if (sample_size % batch_size != 0)
  {
    num_batches += 1;
    final_batch_size = sample_size % batch_size;
  }
}

/** evaluate everything before optimization */
void QMCCostFunctionBatched::checkConfigurations()
{
  ScopedTimer tmp_timer(check_config_timer_);

  RealType et_tot = 0.0;
  RealType e2_tot = 0.0;

  // Ensure number of samples did not change after getConfigurations
  assert(rank_local_num_samples_ == samples_.getNumSamples());

  if (RecordsOnNode_.size1() == 0)
  {
    RecordsOnNode_.resize(rank_local_num_samples_, SUM_INDEX_SIZE);
    if (needGrads)
    {
      DerivRecords_.resize(rank_local_num_samples_, NumOptimizables);
      HDerivRecords_.resize(rank_local_num_samples_, NumOptimizables);
    }
  }
  else if (RecordsOnNode_.size1() != rank_local_num_samples_)
  {
    RecordsOnNode_.resize(rank_local_num_samples_, SUM_INDEX_SIZE);
    if (needGrads)
    {
      DerivRecords_.resize(rank_local_num_samples_, NumOptimizables);
      HDerivRecords_.resize(rank_local_num_samples_, NumOptimizables);
    }
  }
  OperatorBase* nlpp = (includeNonlocalH == "no") ? nullptr : H.getHamiltonian(includeNonlocalH);
  bool compute_nlpp  = useNLPPDeriv && nlpp;
  //set the optimization mode for the trial wavefunction
  Psi.startOptimization();
  //    synchronize the random number generator with the node
  (*MoverRng[0]) = (*RngSaved[0]);
  H.setRandomGenerator(MoverRng[0]);


  // Create crowd-local storage for evaluation
  outputManager.pause();
  opt_eval_.resize(opt_num_crowds_);
  for (int i = 0; i < opt_num_crowds_; i++)
    opt_eval_[i] = std::make_unique<CostFunctionCrowdData>(opt_batch_size_, W, Psi, H, H_KE_node_names_, *MoverRng[0]);
  outputManager.resume();


  // Divide samples among the crowds
  std::vector<int> samples_per_crowd(opt_num_crowds_ + 1);
  FairDivide(rank_local_num_samples_, opt_num_crowds_, samples_per_crowd);

  // lambda to execute on each crowd
  auto evalOptConfig = [](int crowd_id, UPtrVector<CostFunctionCrowdData>& opt_crowds,
                          std::vector<int>& samples_per_crowd, int crowd_size, std::vector<ParticleGradient*>& gradPsi,
                          std::vector<ParticleLaplacian*>& lapPsi, Matrix<Return_rt>& RecordsOnNode,
                          Matrix<Return_rt>& DerivRecords, Matrix<Return_rt>& HDerivRecords, const SampleStack& samples,
                          opt_variables_type& optVars, bool needGrads, bool compute_nlpp,
                          const std::string& includeNonlocalH) {
    CostFunctionCrowdData& opt_data = *opt_crowds[crowd_id];

    int local_samples = samples_per_crowd[crowd_id + 1] - samples_per_crowd[crowd_id];
    int num_batches;
    int final_batch_size;

    compute_batch_parameters(local_samples, crowd_size, num_batches, final_batch_size);

    for (int inb = 0; inb < num_batches; inb++)
    {
      int curr_crowd_size = crowd_size;
      if (inb == num_batches - 1)
        curr_crowd_size = final_batch_size;

      int base_sample_index = inb * crowd_size + samples_per_crowd[crowd_id];

      auto wf_list_no_leader = opt_data.get_wf_list(curr_crowd_size);
      auto p_list_no_leader  = opt_data.get_p_list(curr_crowd_size);
      auto h_list_no_leader  = opt_data.get_h_list(curr_crowd_size);
      const RefVectorWithLeader<ParticleSet> p_list(p_list_no_leader[0], p_list_no_leader);
      const RefVectorWithLeader<TrialWaveFunction> wf_list(wf_list_no_leader[0], wf_list_no_leader);
      const RefVectorWithLeader<QMCHamiltonian> h_list(h_list_no_leader[0], h_list_no_leader);

      ResourceCollectionTeamLock<ParticleSet> mw_pset_lock(opt_data.getSharedResource().pset_res, p_list);
      ResourceCollectionTeamLock<TrialWaveFunction> twfs_res_lock(opt_data.getSharedResource().twf_res, wf_list);
      ResourceCollectionTeamLock<QMCHamiltonian> hams_res_lock(opt_data.getSharedResource().ham_res, h_list);

      auto ref_dLogPsi  = convertPtrToRefVectorSubset(gradPsi, base_sample_index, curr_crowd_size);
      auto ref_d2LogPsi = convertPtrToRefVectorSubset(lapPsi, base_sample_index, curr_crowd_size);

      // Load samples into the crowd data
      for (int ib = 0; ib < curr_crowd_size; ib++)
      {
        samples.loadSample(p_list[ib], base_sample_index + ib);

        // Set the RNG used in QMCHamiltonian.  This is used to offset the grid
        // during spherical integration in the non-local pseudopotential.
        // The RNG state gets reset to the same starting point in correlatedSampling
        // to use the same grid offsets in the correlated sampling values.
        // Currently this code sets the RNG to the same state for every configuration
        // on this node.  Every configuration of electrons is different, and so in
        // theory using the same spherical integration grid should not be a problem.
        // If this needs to be changed, one possibility is to advance the RNG state
        // differently for each configuration.  Make sure the same initialization is
        // performed in correlatedSampling.
        *opt_data.get_rng_ptr_list()[ib] = opt_data.get_rng_save();
        h_list[ib].setRandomGenerator(opt_data.get_rng_ptr_list()[ib].get());
      }

      // Compute distance tables.
      ParticleSet::mw_update(p_list);

      // Log psi and prepare for difference the log psi
      opt_data.zero_log_psi();

      TrialWaveFunction::mw_evaluateDeltaLogSetup(wf_list, p_list, opt_data.get_log_psi_fixed(),
                                                  opt_data.get_log_psi_opt(), ref_dLogPsi, ref_d2LogPsi);

      if (needGrads)
      {
        // Compute parameter derivatives of the wavefunction
        int nparam = optVars.size();
        RecordArray<Return_t> dlogpsi_array(nparam, curr_crowd_size);
        RecordArray<Return_t> dhpsioverpsi_array(nparam, curr_crowd_size);
        TrialWaveFunction::mw_evaluateParameterDerivatives(wf_list, p_list, optVars, dlogpsi_array, dhpsioverpsi_array);


        auto energy_list =
            QMCHamiltonian::mw_evaluateValueAndDerivatives(h_list, wf_list, p_list, optVars, dlogpsi_array,
                                                           dhpsioverpsi_array, compute_nlpp);

        for (int ib = 0; ib < curr_crowd_size; ib++)
        {
          int is = base_sample_index + ib;
          for (int j = 0; j < nparam; j++)
          {
            DerivRecords[is][j]  = std::real(dlogpsi_array.getValue(j, ib));
            HDerivRecords[is][j] = std::real(dhpsioverpsi_array.getValue(j, ib));
          }
          RecordsOnNode[is][LOGPSI_FIXED] = opt_data.get_log_psi_fixed()[ib];
          RecordsOnNode[is][LOGPSI_FREE]  = opt_data.get_log_psi_opt()[ib];
        }

        for (int ib = 0; ib < curr_crowd_size; ib++)
        {
          int is    = base_sample_index + ib;
          auto etmp = energy_list[ib];
          opt_data.get_e0() += etmp;
          opt_data.get_e2() += etmp * etmp;

          RecordsOnNode[is][ENERGY_NEW]   = etmp;
          RecordsOnNode[is][ENERGY_TOT]   = etmp;
          RecordsOnNode[is][REWEIGHT]     = 1.0;
          RecordsOnNode[is][ENERGY_FIXED] = h_list[ib].getLocalPotential();

          if (includeNonlocalH != "no")
          {
            OperatorBase* nlpp = h_list[ib].getHamiltonian(includeNonlocalH);
            if (nlpp)
              RecordsOnNode[is][ENERGY_FIXED] -= nlpp->getValue();
          }
        }
      }
      else
      {
        // Energy
        auto energy_list = QMCHamiltonian::mw_evaluate(h_list, wf_list, p_list);

        for (int ib = 0; ib < curr_crowd_size; ib++)
        {
          int is    = base_sample_index + ib;
          auto etmp = energy_list[ib];
          opt_data.get_e0() += etmp;
          opt_data.get_e2() += etmp * etmp;

          RecordsOnNode[is][ENERGY_NEW]   = etmp;
          RecordsOnNode[is][ENERGY_TOT]   = etmp;
          RecordsOnNode[is][ENERGY_FIXED] = h_list[ib].getLocalPotential();
          RecordsOnNode[is][REWEIGHT]     = 1.0;
        }
      }
    }
  };

  ParallelExecutor<> crowd_tasks;
  crowd_tasks(opt_num_crowds_, evalOptConfig, opt_eval_, samples_per_crowd, opt_batch_size_, dLogPsi, d2LogPsi,
              RecordsOnNode_, DerivRecords_, HDerivRecords_, samples_, OptVariablesForPsi, needGrads, compute_nlpp,
              includeNonlocalH);
  // Sum energy values over crowds
  for (int i = 0; i < opt_eval_.size(); i++)
  {
    et_tot += opt_eval_[i]->get_e0();
    e2_tot += opt_eval_[i]->get_e2();
  }

  OptVariablesForPsi.setComputed();
  //     app_log() << "  VMC Efavg = " << eft_tot/static_cast<Return_t>(wPerNode[NumThreads]) << std::endl;
  //Need to sum over the processors
  std::vector<Return_rt> etemp(3);
  etemp[0] = et_tot;
  etemp[1] = static_cast<Return_rt>(rank_local_num_samples_);
  etemp[2] = e2_tot;
  // Sum energy values over nodes
  myComm->allreduce(etemp);
  Etarget    = static_cast<Return_rt>(etemp[0] / etemp[1]);
  NumSamples = static_cast<int>(etemp[1]);
  app_log() << "  VMC Eavg = " << Etarget << std::endl;
  app_log() << "  VMC Evar = " << etemp[2] / etemp[1] - Etarget * Etarget << std::endl;
  app_log() << "  Total weights = " << etemp[1] << std::endl;
  app_log().flush();
  setTargetEnergy(Etarget);
  ReportCounter = 0;

  //collect SumValue for computedCost
  NumWalkersEff           = etemp[1];
  SumValue[SUM_WGT]       = etemp[1];
  SumValue[SUM_WGTSQ]     = etemp[1];
  SumValue[SUM_E_WGT]     = etemp[0];
  SumValue[SUM_ESQ_WGT]   = etemp[2];
  SumValue[SUM_E_BARE]    = etemp[0];
  SumValue[SUM_ESQ_BARE]  = etemp[2];
  SumValue[SUM_ABSE_BARE] = 0.0;
}

#ifdef HAVE_LMY_ENGINE
void QMCCostFunctionBatched::engine_checkConfigurations(cqmc::engine::LMYEngine<Return_t>* EngineObj,
                                                        DescentEngine& descentEngineObj,
                                                        const std::string& MinMethod)
{
  APP_ABORT("LMYEngine not implemented with batch optimization");
}
#endif


void QMCCostFunctionBatched::resetPsi(bool final_reset)
{
  if (OptVariables.size() < OptVariablesForPsi.size())
    for (int i = 0; i < equalVarMap.size(); ++i)
      OptVariablesForPsi[equalVarMap[i][0]] = OptVariables[equalVarMap[i][1]];
  else
    for (int i = 0; i < OptVariables.size(); ++i)
      OptVariablesForPsi[i] = OptVariables[i];
  if (final_reset)
  {
    Psi.stopOptimization();
  }
  //cout << "######### QMCCostFunctionBatched::resetPsi " << std::endl;
  //OptVariablesForPsi.print(std::cout);
  //cout << "-------------------------------------- " << std::endl;
  Psi.resetParameters(OptVariablesForPsi);
  for (int i = 0; i < opt_eval_.size(); i++)
  {
    for (int j = 0; j < opt_eval_[i]->get_wf_ptr_list().size(); j++)
    {
      opt_eval_[i]->get_wf_ptr_list()[j]->resetParameters(OptVariablesForPsi);
    }
  }
}

QMCCostFunctionBatched::Return_rt QMCCostFunctionBatched::correlatedSampling(bool needGrad)
{
  ScopedTimer tmp_timer(corr_sampling_timer_);

  {
    //    synchronize the random number generator with the node
    (*MoverRng[0]) = (*RngSaved[0]);
    H.setRandomGenerator(MoverRng[0]);
  }

  //Return_rt wgt_node = 0.0, wgt_node2 = 0.0;
  const bool nlpp    = (includeNonlocalH != "no");
  Return_rt wgt_tot  = 0.0;
  Return_rt wgt_tot2 = 0.0;

  // Ensure number of samples did not change after getConfiguration
  assert(rank_local_num_samples_ == samples_.getNumSamples());

  Return_rt inv_n_samples = 1.0 / samples_.getGlobalNumSamples();

  bool compute_nlpp             = useNLPPDeriv && (includeNonlocalH != "no");
  bool compute_all_from_scratch = (includeNonlocalH != "no"); //true if we have nlpp

  // Divide samples among crowds
  std::vector<int> samples_per_crowd(opt_num_crowds_ + 1);
  FairDivide(rank_local_num_samples_, opt_num_crowds_, samples_per_crowd);

  // lambda to execute on each crowd
  auto evalOptCorrelated = [](int crowd_id, UPtrVector<CostFunctionCrowdData>& opt_crowds,
                              const std::vector<int>& samples_per_crowd, int crowd_size,
                              std::vector<ParticleGradient*>& gradPsi, std::vector<ParticleLaplacian*>& lapPsi,
                              Matrix<Return_rt>& RecordsOnNode, Matrix<Return_rt>& DerivRecords,
                              Matrix<Return_rt>& HDerivRecords, const SampleStack& samples,
                              const opt_variables_type& optVars, bool compute_all_from_scratch, Return_rt vmc_or_dmc,
                              bool needGrad, bool compute_nlpp) {
    CostFunctionCrowdData& opt_data = *opt_crowds[crowd_id];


    int local_samples = samples_per_crowd[crowd_id + 1] - samples_per_crowd[crowd_id];

    int num_batches;
    int final_batch_size;
    compute_batch_parameters(local_samples, crowd_size, num_batches, final_batch_size);

    for (int inb = 0; inb < num_batches; inb++)
    {
      int curr_crowd_size = crowd_size;
      if (inb == num_batches - 1)
      {
        curr_crowd_size = final_batch_size;
      }

      int base_sample_index = inb * crowd_size + samples_per_crowd[crowd_id];

      auto p_list_no_leader  = opt_data.get_p_list(curr_crowd_size);
      auto wf_list_no_leader = opt_data.get_wf_list(curr_crowd_size);
      auto h0_list_no_leader = opt_data.get_h0_list(curr_crowd_size);
      const RefVectorWithLeader<ParticleSet> p_list(p_list_no_leader[0], p_list_no_leader);
      const RefVectorWithLeader<TrialWaveFunction> wf_list(wf_list_no_leader[0], wf_list_no_leader);
      const RefVectorWithLeader<QMCHamiltonian> h0_list(h0_list_no_leader[0], h0_list_no_leader);

      ResourceCollectionTeamLock<ParticleSet> mw_pset_lock(opt_data.getSharedResource().pset_res, p_list);
      ResourceCollectionTeamLock<TrialWaveFunction> twfs_res_lock(opt_data.getSharedResource().twf_res, wf_list);
      ResourceCollectionTeamLock<QMCHamiltonian> hams_res_lock(opt_data.getSharedResource().ham_res, h0_list);

      // Load this batch of samples into the crowd data
      for (int ib = 0; ib < curr_crowd_size; ib++)
      {
        samples.loadSample(p_list[ib], base_sample_index + ib);
        // Copy the saved RNG state
        *opt_data.get_rng_ptr_list()[ib] = opt_data.get_rng_save();
        h0_list[ib].setRandomGenerator(opt_data.get_rng_ptr_list()[ib].get());
      }

      // Update distance tables, etc for the loaded sample positions
      ParticleSet::mw_update(p_list, true);

      // Evaluate difference in log psi

      std::vector<std::unique_ptr<ParticleSet::ParticleGradient>> dummyG_ptr_list;
      std::vector<std::unique_ptr<ParticleSet::ParticleLaplacian>> dummyL_ptr_list;
      RefVector<ParticleSet::ParticleGradient> dummyG_list;
      RefVector<ParticleSet::ParticleLaplacian> dummyL_list;
      if (compute_all_from_scratch)
      {
        int nptcl = gradPsi[0]->size();
        dummyG_ptr_list.reserve(curr_crowd_size);
        dummyL_ptr_list.reserve(curr_crowd_size);
        for (int i = 0; i < curr_crowd_size; i++)
        {
          dummyG_ptr_list.emplace_back(std::make_unique<ParticleGradient>(nptcl));
          dummyL_ptr_list.emplace_back(std::make_unique<ParticleLaplacian>(nptcl));
        }
        dummyG_list = convertUPtrToRefVector(dummyG_ptr_list);
        dummyL_list = convertUPtrToRefVector(dummyL_ptr_list);
      }
      opt_data.zero_log_psi();

      TrialWaveFunction::mw_evaluateDeltaLog(wf_list, p_list, opt_data.get_log_psi_opt(), dummyG_list, dummyL_list,
                                             compute_all_from_scratch);

      Return_rt inv_n_samples = 1.0 / samples.getGlobalNumSamples();

      for (int ib = 0; ib < curr_crowd_size; ib++)
      {
        int is = base_sample_index + ib;
        wf_list[ib].G += *gradPsi[is];
        wf_list[ib].L += *lapPsi[is];
        // This is needed to get the KE correct in QMCHamiltonian::mw_evaluate below
        p_list[ib].G += *gradPsi[is];
        p_list[ib].L += *lapPsi[is];
        Return_rt weight            = vmc_or_dmc * (opt_data.get_log_psi_opt()[ib] - RecordsOnNode[is][LOGPSI_FREE]);
        RecordsOnNode[is][REWEIGHT] = weight;
        // move to opt_data
        opt_data.get_wgt() += inv_n_samples * weight;
        opt_data.get_wgt2() += inv_n_samples * weight * weight;
      }

      if (needGrad)
      {
        // Parameter derivatives
        int nparam = optVars.size();
        RecordArray<Return_t> dlogpsi_array(nparam, curr_crowd_size);
        RecordArray<Return_t> dhpsioverpsi_array(nparam, curr_crowd_size);

        TrialWaveFunction::mw_evaluateParameterDerivatives(wf_list, p_list, optVars, dlogpsi_array, dhpsioverpsi_array);


        // Energy
        auto energy_list =
            QMCHamiltonian::mw_evaluateValueAndDerivatives(h0_list, wf_list, p_list, optVars, dlogpsi_array,
                                                           dhpsioverpsi_array, compute_nlpp);

        for (int ib = 0; ib < curr_crowd_size; ib++)
        {
          int is                        = base_sample_index + ib;
          auto etmp                     = energy_list[ib];
          RecordsOnNode[is][ENERGY_NEW] = etmp + RecordsOnNode[is][ENERGY_FIXED];
          for (int j = 0; j < nparam; j++)
          {
            if (optVars.recompute(j))
            {
              DerivRecords[is][j]  = std::real(dlogpsi_array.getValue(j, ib));
              HDerivRecords[is][j] = std::real(dhpsioverpsi_array.getValue(j, ib));
            }
          }
        }
      }
      else
      {
        // Just energy needed if no gradients
        auto energy_list = QMCHamiltonian::mw_evaluate(h0_list, wf_list, p_list);
        for (int ib = 0; ib < curr_crowd_size; ib++)
        {
          int is                        = base_sample_index + ib;
          auto etmp                     = energy_list[ib];
          RecordsOnNode[is][ENERGY_NEW] = etmp + RecordsOnNode[is][ENERGY_FIXED];
        }
      }
    }
  };

  ParallelExecutor<> crowd_tasks;
  crowd_tasks(opt_num_crowds_, evalOptCorrelated, opt_eval_, samples_per_crowd, opt_batch_size_, dLogPsi, d2LogPsi,
              RecordsOnNode_, DerivRecords_, HDerivRecords_, samples_, OptVariablesForPsi, compute_all_from_scratch,
              vmc_or_dmc, needGrad, compute_nlpp);
  // Sum weights over crowds
  for (int i = 0; i < opt_eval_.size(); i++)
  {
    wgt_tot += opt_eval_[i]->get_wgt();
    wgt_tot2 += opt_eval_[i]->get_wgt2();
  }

  //this is MPI barrier
  OHMMS::Controller->barrier();
  //collect the total weight for normalization and apply maximum weight
  myComm->allreduce(wgt_tot);
  myComm->allreduce(wgt_tot2);
  //    app_log()<<"Before Purge"<<wgt_tot<<" "<<wgt_tot2<< std::endl;
  Return_rt wgtnorm = (wgt_tot == 0) ? 0 : wgt_tot;
  wgt_tot           = 0.0;
  {
    for (int iw = 0; iw < rank_local_num_samples_; iw++)
    {
      Return_rt* restrict saved = RecordsOnNode_[iw];
      saved[REWEIGHT] =
          std::min(std::exp(saved[REWEIGHT] - wgtnorm), std::numeric_limits<Return_rt>::max() * (RealType)0.1);
      wgt_tot += inv_n_samples * saved[REWEIGHT];
    }
  }
  myComm->allreduce(wgt_tot);
  //    app_log()<<"During Purge"<<wgt_tot<<" "<< std::endl;
  wgtnorm = (wgt_tot == 0) ? 1 : 1.0 / wgt_tot;
  wgt_tot = 0.0;
  {
    for (int iw = 0; iw < rank_local_num_samples_; iw++)
    {
      Return_rt* restrict saved = RecordsOnNode_[iw];
      saved[REWEIGHT]           = std::min(saved[REWEIGHT] * wgtnorm, MaxWeight);
      wgt_tot += inv_n_samples * saved[REWEIGHT];
    }
  }
  myComm->allreduce(wgt_tot);
  //    app_log()<<"After Purge"<<wgt_tot<<" "<< std::endl;
  for (int i = 0; i < SumValue.size(); i++)
    SumValue[i] = 0.0;
  {
    for (int iw = 0; iw < rank_local_num_samples_; iw++)
    {
      const Return_rt* restrict saved = RecordsOnNode_[iw];
      //      Return_t weight=saved[REWEIGHT]*wgt_tot;
      Return_rt eloc_new = saved[ENERGY_NEW];
      Return_rt delE     = std::pow(std::abs(eloc_new - EtargetEff), PowerE);
      SumValue[SUM_E_BARE] += eloc_new;
      SumValue[SUM_ESQ_BARE] += eloc_new * eloc_new;
      SumValue[SUM_ABSE_BARE] += delE;
      SumValue[SUM_E_WGT] += eloc_new * saved[REWEIGHT];
      SumValue[SUM_ESQ_WGT] += eloc_new * eloc_new * saved[REWEIGHT];
      SumValue[SUM_ABSE_WGT] += delE * saved[REWEIGHT];
      SumValue[SUM_WGT] += saved[REWEIGHT];
      SumValue[SUM_WGTSQ] += saved[REWEIGHT] * saved[REWEIGHT];
    }
  }
  //collect everything
  myComm->allreduce(SumValue);
  //     for (int i=0; i<SumValue.size(); i++) std::cerr <<SumValue[i]<<"  ";
  //     std::cerr << std::endl;
  //     app_log()<<"After purge Energy Variance Weight "
  //      << SumValue[SUM_E_WGT]/SumValue[SUM_WGT] << " "
  //      << SumValue[SUM_ESQ_WGT]/SumValue[SUM_WGT] -(SumValue[SUM_E_WGT]/SumValue[SUM_WGT])*(SumValue[SUM_E_WGT]/SumValue[SUM_WGT]) << " "
  //      << SumValue[SUM_WGT]*SumValue[SUM_WGT]/SumValue[SUM_WGTSQ] << std::endl;
  return SumValue[SUM_WGT] * SumValue[SUM_WGT] / SumValue[SUM_WGTSQ];
}


// Construct the overlap and Hamiltonian matrices for the linear method
// A sum over samples.  Inputs are
//   DerivRecords - derivative of log psi ( d ln (psi) / dp = 1/psi * d psi / dp )
//   HDerivRecords - derivative of Hamiltonian
//   RecordsOnNode - energies and weights (for reweighting)
//   SumValue - sums of energies and weights
// Outputs
//   Left - Hamiltonian matrix
//   Right - overlap matrix
//

QMCCostFunctionBatched::Return_rt QMCCostFunctionBatched::fillOverlapHamiltonianMatrices(Matrix<Return_rt>& Left,
                                                                                         Matrix<Return_rt>& Right)
{
  ScopedTimer tmp_timer(fill_timer_);

  RealType b1, b2;
  if (GEVType == "H2")
  {
    b1 = w_beta;
    b2 = 0;
  }
  else
  {
    b2 = w_beta;
    b1 = 0;
  }

  Right = 0.0;
  Left  = 0.0;

  //     resetPsi();
  //     Return_t NWE = NumWalkersEff=correlatedSampling(true);
  curAvg_w            = SumValue[SUM_E_WGT] / SumValue[SUM_WGT];
  Return_rt curAvg2_w = SumValue[SUM_ESQ_WGT] / SumValue[SUM_WGT];
  //    RealType H2_avg = 1.0/curAvg2_w;
  RealType H2_avg = 1.0 / (curAvg_w * curAvg_w);
  //    RealType H2_avg = 1.0/std::sqrt(curAvg_w*curAvg_w*curAvg2_w);
  RealType V_avg = curAvg2_w - curAvg_w * curAvg_w;
  std::vector<Return_rt> D_avg(getNumParams(), 0.0);
  Return_rt wgtinv = 1.0 / SumValue[SUM_WGT];

  for (int iw = 0; iw < rank_local_num_samples_; iw++)
  {
    const Return_rt* restrict saved = RecordsOnNode_[iw];
    Return_rt weight                = saved[REWEIGHT] * wgtinv;
    const Return_rt* Dsaved         = DerivRecords_[iw];
    for (int pm = 0; pm < getNumParams(); pm++)
    {
      D_avg[pm] += Dsaved[pm] * weight;
    }
  }

  myComm->allreduce(D_avg);

  for (int iw = 0; iw < rank_local_num_samples_; iw++)
  {
    const Return_rt* restrict saved = RecordsOnNode_[iw];
    Return_rt weight                = saved[REWEIGHT] * wgtinv;
    Return_rt eloc_new              = saved[ENERGY_NEW];
    const Return_rt* Dsaved         = DerivRecords_[iw];
    const Return_rt* HDsaved        = HDerivRecords_[iw];

    std::vector<int> params_per_crowd(opt_num_crowds_ + 1);
    FairDivide(getNumParams(), opt_num_crowds_, params_per_crowd);


    auto constructMatrices = [](int crowd_id, std::vector<int>& crowd_ranges, int numParams, const Return_rt* Dsaved,
                                const Return_rt* HDsaved, Return_rt weight, Return_rt eloc_new, RealType H2_avg,
                                RealType V_avg, std::vector<Return_rt>& D_avg, RealType b1, RealType b2,
                                RealType curAvg_w, Matrix<Return_rt>& Left, Matrix<Return_rt>& Right) {
      int local_pm_start = crowd_ranges[crowd_id];
      int local_pm_end   = crowd_ranges[crowd_id + 1];

      for (int pm = local_pm_start; pm < local_pm_end; pm++)
      {
        Return_rt wfe = (HDsaved[pm] + (Dsaved[pm] - D_avg[pm]) * eloc_new) * weight;
        Return_rt wfd = (Dsaved[pm] - D_avg[pm]) * weight;
        Return_rt vterm =
            HDsaved[pm] * (eloc_new - curAvg_w) + (Dsaved[pm] - D_avg[pm]) * eloc_new * (eloc_new - 2.0 * curAvg_w);
        //                 H2
        Right(0, pm + 1) += b1 * H2_avg * vterm * weight;
        Right(pm + 1, 0) += b1 * H2_avg * vterm * weight;
        //                 Variance
        Left(0, pm + 1) += b2 * vterm * weight;
        Left(pm + 1, 0) += b2 * vterm * weight;
        //                 Hamiltonian
        Left(0, pm + 1) += (1 - b2) * wfe;
        Left(pm + 1, 0) += (1 - b2) * wfd * eloc_new;
        for (int pm2 = 0; pm2 < numParams; pm2++)
        {
          //                Hamiltonian
          Left(pm + 1, pm2 + 1) += (1 - b2) * wfd * (HDsaved[pm2] + (Dsaved[pm2] - D_avg[pm2]) * eloc_new);
          //                Overlap
          RealType ovlij = wfd * (Dsaved[pm2] - D_avg[pm2]);
          Right(pm + 1, pm2 + 1) += ovlij;
          //                Variance
          RealType varij = weight * (HDsaved[pm] - 2.0 * (Dsaved[pm] - D_avg[pm]) * eloc_new) *
              (HDsaved[pm2] - 2.0 * (Dsaved[pm2] - D_avg[pm2]) * eloc_new);
          Left(pm + 1, pm2 + 1) += b2 * (varij + V_avg * ovlij);
          //                H2
          Right(pm + 1, pm2 + 1) += b1 * H2_avg * varij;
        }
      }
    };

    ParallelExecutor<> crowd_tasks;
    crowd_tasks(opt_num_crowds_, constructMatrices, params_per_crowd, getNumParams(), Dsaved, HDsaved, weight, eloc_new,
                H2_avg, V_avg, D_avg, b1, b2, curAvg_w, Left, Right);
  }
  myComm->allreduce(Right);
  myComm->allreduce(Left);
  Left(0, 0)  = (1 - b2) * curAvg_w + b2 * V_avg;
  Right(0, 0) = 1.0 + b1 * H2_avg * V_avg;
  if (GEVType == "H2")
    return H2_avg;

  return 1.0;
}
} // namespace qmcplusplus
