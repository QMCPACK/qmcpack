//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    


/** @file QMCFixedSampleLinearOptimize.h
 * @brief Definition of QMCDriver which performs VMC and optimization.
 */
#ifndef QMCPLUSPLUS_QMCFSLINEAROPTIMIZATION_VMCSINGLE_H
#define QMCPLUSPLUS_QMCFSLINEAROPTIMIZATION_VMCSINGLE_H

#include "QMCDrivers/QMCLinearOptimize.h"
#include "Optimize/NRCOptimization.h"
#ifdef HAVE_LMY_ENGINE
#include "formic/utils/matrix.h"
#include "formic/utils/lmyengine/engine.h"
#endif

namespace qmcplusplus
{
/** @ingroup QMCDrivers
 * @brief Implements wave-function optimization
 *
 * Optimization by correlated sampling method with configurations
 * generated from VMC.
 */

class QMCFixedSampleLinearOptimize: public QMCLinearOptimize, private NRCOptimization<QMCTraits::RealType>
{
public:

  ///Constructor.
  QMCFixedSampleLinearOptimize(MCWalkerConfiguration& w, TrialWaveFunction& psi,
                               QMCHamiltonian& h, HamiltonianPool& hpool, WaveFunctionPool& ppool);

  ///Destructor
  ~QMCFixedSampleLinearOptimize();

  ///Run the Optimization algorithm.
  bool run();
  ///process xml node
  bool put(xmlNodePtr cur);
  RealType Func(Return_t dl);

private:
  inline bool ValidCostFunction(bool valid)
  {
    if (!valid)
      app_log()<<" Cost Function is Invalid. If this frequently, try reducing the step size of the line minimization or reduce the number of cycles. " << std::endl;
    return valid;
  }

  // check if the proposed new cost function value is the best available
  bool is_best_cost(const int ii, const std::vector<RealType> & cv, const RealType ic) const;

  // perform the adaptive three-shift update
  bool adaptive_three_shift_run();

  // perform the single-shift update, no sample regeneration
  bool one_shift_run();

  void solveShiftsWithoutLMYEngine(const std::vector<double> & shifts_i,
                                   const std::vector<double> & shiffts_s,
                                   std::vector<std::vector<RealType> > & parameterDirections);

  #ifdef HAVE_LMY_ENGINE
  formic::VarDeps vdeps;
  cqmc::engine::LMYEngine * EngineObj;
  #endif

  // prepare a vector of shifts to try
  std::vector<double> prepare_shifts(const double central_shift) const;

  // previous update
#ifdef HAVE_LMY_ENGINE
  std::vector<formic::ColVec<double> > previous_update;
#endif

  void print_cost_summary_header();
  void print_cost_summary(const double si, const double ss, const RealType mc, const RealType cv, const int ind, const int bi, const bool gu);

  int NumOfVMCWalkers;
  ///Number of iterations maximum before generating new configurations.
  int Max_iterations;
  int nstabilizers;
  RealType stabilizerScale, bigChange, exp0, exp1, stepsize, savedQuadstep;
  std::string MinMethod, GEVtype, StabilizerMethod, GEVSplit;
  RealType w_beta;
  /// number of previous steps to orthogonalize to.
  int eigCG;
  /// total number of cg steps per iterations
  int  TotalCGSteps;
  /// whether to use the adaptive three-shift scheme
  bool doAdaptiveThreeShift;
  /// whether to use the single-shift scheme
  bool doOneShiftOnly;
  /// the previous best identity shift
  RealType bestShift_i;
  /// the previous best overlap shift
  RealType bestShift_s;
  /// current shift_i, shift_s input values
  RealType shift_i_input,shift_s_input;
  /// number of shifts we will try
  int num_shifts;
  /// the maximum relative change in the cost function for the adaptive three-shift scheme
  RealType max_relative_cost_change;
  ///max amount a parameter may change relative to current wave function weight
  RealType max_param_change;
  ///whether we are targeting an excited state
  std::string targetExcitedStr;
  ///whether we are targeting an excited state
  bool targetExcited;
  ///whether we are doing block algorithm
  std::string block_lmStr;
  ///whether we are doing block algorithm
  bool block_lm;
  ///number of blocks used in block algorithm
  int nblocks;
  ///number of old updates kept
  int nolds;
  ///number of directions kept
  int nkept;
  ///number of samples to do in correlated sampling part
  int nsamp_comp; 
  ///the shift to use when targeting an excited state
  RealType omega_shift;
  ///whether to do the first part of block lm
  bool block_first;
  ///whether to do the second part of block lm
  bool block_second;
  ///whethe to do the third part of blocl lm
  bool block_third;
};
}
#endif
