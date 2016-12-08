/////////////////////////////////////////////////////////////////////////////////////////////////
// \file cqmc/engine/energy_target.h
//
//
// \brief  header file for the harmonic davidson energy function 
//
/////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef EIGINE_ENERGY_TARGET_FUNCTION_HEADER
#define EIGINE_ENERGY_TARGET_FUNCTION_HEADER

#include<vector>

#include<boost/shared_ptr.hpp>

namespace cqmc {

  namespace engine {
    
    void et(const bool exact_sampling,
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
            std::ostream & output);



  }

}

#endif
