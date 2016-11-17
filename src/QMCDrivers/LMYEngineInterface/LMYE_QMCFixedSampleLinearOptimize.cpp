#include "QMCDrivers/QMCFixedSampleLinearOptimize.h"
#include "Particle/HDFWalkerIO.h"
#include "Particle/DistanceTable.h"
#include "OhmmsData/AttributeSet.h"
#include "Message/CommOperators.h"
#if defined(ENABLE_OPENMP)
#include "QMCDrivers/VMC/VMCSingleOMP.h"
#include "QMCDrivers/QMCCostFunctionOMP.h"
#endif
//#include "QMCDrivers/VMC/VMCSingle.h"
//#include "QMCDrivers/QMCCostFunctionSingle.h"
#include "QMCApp/HamiltonianPool.h"
#include "Numerics/Blasf.h"
#include "Numerics/MatrixOperators.h"
#include <cassert>
#if defined(QMC_CUDA)
#include "QMCDrivers/VMC/VMC_CUDA.h"
#include "QMCDrivers/QMCCostFunctionCUDA.h"
#endif
#include <iostream>
#include <fstream>
#include <utility>

//#include "Eigen/Dense"
#include "formic/utils/lmyengine/engine.h"
#include "formic/utils/lmyengine/var_dependencies.h"

namespace qmcplusplus {

///////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief  For each set of shifts, solves the linear method eigenproblem by calling the LMYEngine
///
/// \param[in]      eval_target           value to target when choosing which linear method
///                                       eigenvector to use as the solution
/// \param[in]      shfits_i              vector of identity shifts
/// \param[in]      shfits_s              vector of overlap shifts
/// \param[out]     parameterDirections   on exit, the update directions for the different shifts
///
///////////////////////////////////////////////////////////////////////////////////////////////////
void QMCFixedSampleLinearOptimize::solveShiftsWithLMYEngine(const RealType eval_target,
                                                            const std::vector<double> & shifts_i,
                                                            const std::vector<double> & shifts_s,
                                                            std::vector<std::vector<RealType> > & parameterDirections) {

  //// get number of shifts to solve
  //const int nshifts = shifts_i.size();

  //// get number of optimizable parameters
  //numParams = optTarget->NumParams();

  //// get dimension of the linear method matrices
  //N = numParams + 1;

  //// declare the vectors and matrices we need
  //std::vector<RealType> wgt_vec;
  //std::vector<RealType> vgs_vec;
  ////Eigen::MatrixXd der_rat_wf;
  ////Eigen::MatrixXd der_rat_en;
  //formic::Matrix<double> der_rat_wf;
  //formic::Matrix<double> der_rat_en;

  //// prepare the vectors and matrices
  //optTarget->prepareDerRatDataMatrices(wgt_vec, vgs_vec, der_rat_wf, der_rat_en);

  //// prepare a variable dependency object with no dependencies
  //formic::VarDeps vdeps(numParams, std::vector<double>());

  //// prepare vector that will hold the wave function updates for each set of shifts
  //std::vector<double> solved_updates;

  //// prepare shift information in format the engine needs
  //const double base_shift_i = shifts_i.at(1);
  //const double base_shift_s = shifts_s.at(1);
  //std::vector<double> shift_scales(shifts_i.size(), 1.0);
  //for (int i = 0; i < shift_scales.size(); i++)
  //  shift_scales.at(i) = shifts_i.at(i) / base_shift_i;

  //// hack for strange engine requirements on shift_scale
  ////shift_scales.push_back( 0.5 * ( shift_scales.at(1) + shift_scales.at(2) ) );

  //// prepare vectors to tell which shifts solved correctly
  //std::vector<bool> good_solve(shift_scales.size(), false);
  //std::vector<int> shift_solved(shift_scales.size(), 0);

  //// call the LMYEngine to solve the eigenproblem
  //cqmc::engine::call_engine(&vdeps, // dep_ptr
  //                          !targetExcited, // ground_state
  //                          true, // build_lm_matrix
  //                          false, // variance_correct
  //                          false, // print_matrix
  //                          numParams, // lm_krylov_iter,
  //                          1, // lm_spam_inner_iter,
  //                          0, // appro_degree
  //                          1.0e-6, // lm_eigen_thresh
  //                          1.0e-3, // lm_min_S_eval
  //                          false, // spam_use
  //                          false, // var_deps_use 
  //                          true, // chase_lowest
  //                          false, // chase_closest
  //                          eval_target, // init_energy
  //                          0.0, // init_var 
  //                          100.0, // lm_max_e_change 
  //                          base_shift_i, // lm_ham_shift_i
  //                          base_shift_s, // lm_ham_shift_s
  //                          omega_shift, // omega
  //                          0.0, // var_weight
  //                          10.0, // lm_max_update_abs
  //                          der_rat_wf, // der_rat
  //                          der_rat_en, // le_dr
  //                          vgs_vec, // vg
  //                          wgt_vec, // weight
  //                          solved_updates, // vf_var
  //                          good_solve, // good_solve
  //                          shift_solved, //solved_shift
  //                          shift_scales, // shift_scale
  //                          app_log()); // output

  //// extract the update for each shift
  //parameterDirections.resize(nshifts);
  //auto pdit = parameterDirections.begin();
  //for (int i = 0; i < nshifts; i++, pdit++) {
  //  pdit->assign(N, 0.0);
  //  if ( good_solve.at(i) )
  //    std::copy(&solved_updates.at(i*N), &solved_updates.at(i*N) + N, pdit->begin());
  //  else
  //    pdit->at(0) = 1.0;
  //}

  //// scale the update relative to the current wave function (this use of a lamda function is egregeous, yes...)
  //std::for_each(parameterDirections.begin(), parameterDirections.end(), [] (std::vector<double> & v) { for (auto it = v.rbegin(); it != v.rend(); it++) *it /= v.at(0); });

}

}
