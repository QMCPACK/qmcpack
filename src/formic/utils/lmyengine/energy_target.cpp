///////////////////////////////////////////////////////////////////////////////////////////////////
// \file cqmc/engine/energy_target.cpp
//
//
//  \brief implementation file for harmonic davidson energy and target function 
//
///////////////////////////////////////////////////////////////////////////////////////////////////

#include<vector>

#include<boost/shared_ptr.hpp>

#include<formic/utils/lmyengine/energy_target_accu.h>
#include<formic/utils/lmyengine/energy_target.h>



//////////////////////////////////////////////////////////////////////////////////////////////////
// \brief  compute energy and target function using input data
// 
// \param[in]   exact_sampling    whether to use exact sampling or not
// \param[in]   ground_state      whether to do ground state calculation
// \param[in]   variance_correct  whether to correct the "finite variance" problem
// \param[in]   print             whether to print out energy statistics
// \param[in]   hd_lm_shift       harmonic davidson shift
// \param[in]   var_weight        weight of the variance correction term
// \param[in]   le_history        local energy history
// \param[in]   vg_history        |value/guiding|^2 history
// \param[in]   w_history         weight history
// \param[out]  energy            average of local energy
// \param[out]  esdev             energy stanard deviation
// \param[out]  eserr             energy statistical error
// \param[out]  target            target function value
// \param[out]  tserr             target function statistical error
//
//////////////////////////////////////////////////////////////////////////////////////////////////
void cqmc::engine::et(const bool exact_sampling,
                      const bool ground_state,
                      const bool variance_correct,
                      const bool print,
                      const double hd_lm_shift,
                      const double var_weight,
                      const std::vector<double> & le_history,
                      const std::vector<double> & vg_history,
                      const std::vector<double> & w_history,
                      double & energy,
                      double & esdev,
                      double & eserr,
                      double & target,
                      double & tserr,
                      std::ostream & output)
{

  // create a computer for energy and target function 
  cqmc::engine::ETCompute etcal(le_history, vg_history, w_history, exact_sampling, ground_state, variance_correct, hd_lm_shift, var_weight);

  // compute energy, target function value and relevent parameters
  etcal.en_tar();
  
  // correct finite variance issue if request
  if ( variance_correct ) 
    etcal.correct_finite_variance();

  // get the energy and relevent quantities
  energy = etcal.energy();

  esdev = std::sqrt(etcal.variance());

  eserr = etcal.eserr(output);

  // get the target function value and relevent quantities 
  if ( !ground_state ) {
    target = etcal.tar_fn_val();

    tserr = etcal.tserr(output);
  }

  // print statistical data if requested
  if ( print )
    etcal.print_statistics(output);

}
