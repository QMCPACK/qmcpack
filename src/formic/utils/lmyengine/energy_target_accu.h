///////////////////////////////////////////////////////////////////////////////////////////////
// \file energy_target_accu.h
//
//
// \brief  header file for harmonic davidson eigine energy and target function calculate
//
///////////////////////////////////////////////////////////////////////////////////////////////

#ifndef EIGINE_ENERGY_TARGET_ACCU_HEADER
#define EIGINE_ENERGY_TARGET_ACCU_HEADER

#include<complex>
#include<vector>
#include<string>
#include<numeric>
#include<cassert>
#include<algorithm>

#include<boost/format.hpp>
#include<boost/shared_ptr.hpp>

//#include<mpi.h>


namespace cqmc {
  
  namespace engine {

////////////////////////////////////////////////////////////////////////////////////////////////
// \brief a class to calculate energy, harmonic davidson target function value and relevant 
//        quantities
//
////////////////////////////////////////////////////////////////////////////////////////////////

class ETCompute {

  private:

    /// \brief [in] history of sampled local energies 
    std::vector<double> _le_history;

    /// \brief [in] history of sampled |value/guiding|^2 ratios 
    std::vector<double> _vg_history;

    /// \brief [in] history of sampled configuration weight 
    std::vector<double> _w_history;

    /// \brief flag to tell whether exact sampling is being used
    bool _exact_sampling;

    /// \brief flag to tell whether to do ground state calculation
    bool _ground_state;

    /// \brief flag to tell whwther to correct the "finite variance" problem
    bool _variance_correct;

    /// \brief energy;
    double _energy;

    /// \brief energy^2
    double _energy_s;

    /// \brief target function value
    double _target_fn_val;

    /// \brief target function variance
    double _target_fn_var;

    /// \brief variance of energy 
    double _variance;

    /// \brief variance of square of energy
    double _variance_s;

    /// \brief estimated statistical error of energy 
    double _serr;

    /// \brief estimated target function numerator statistical error 
    double _tnserr;

    /// \brief estimated target function denominator statistical error 
    double _tdserr;

    /// \brief estimated target function statistical error 
    double _tserr; 

    /// \brief harmonic davidson shift 
    double _hd_shift;

    /// \brief weight of the variance correction term
    double _var_weight;

    /// \brief a history of sampled local energies times the |value/guiding|^2 raitos
    std::vector<double> _lev_history;

    /// \brief a history of sampled target function numerator
    std::vector<double> _tn_history;

    /// \brief a history of target function numerator times the |value/guiding|^2 ratios
    std::vector<double> _tnv_history;

    /// \brief a history of target function denominator
    std::vector<double> _td_history;

    /// \brief a history of target function denominator times the |value/guiding|^2 ratios
    std::vector<double> _tdv_history;

    /// \brief a history of sampled local energy square 
    std::vector<double> _les_history;

    /// \brief a history of sampled local energy square times |value/guiding|^2 ratios 
    std::vector<double> _lesv_history;
    
    /// \brief maximum length for which to compute autocorrelation
    //int _max_ac_length;

  public:

    //////////////////////////////////////////////////////////////////////////////////////////////////
    // \brief  constructs the computer
    //
    //
    //////////////////////////////////////////////////////////////////////////////////////////////////
    ETCompute(const std::vector<double> & le_history,
              const std::vector<double> & vg_history,
              const std::vector<double> & w_history,
              const bool exact_sampling,
              const bool ground_state,
              const bool variance_correct,
              const double hd_shift,
              const double var_weight);


    ///////////////////////////////////////////////////////////////////////////////////////////////////
    // \brief function that initializes all history list
    //
    //
    ///////////////////////////////////////////////////////////////////////////////////////////////////
    void history_initialize();
    
    ///////////////////////////////////////////////////////////////////////////////////////////////////
    // \brief  returns the variance of the block average for the specific block length
    //
    //
    //
    ///////////////////////////////////////////////////////////////////////////////////////////////////
    double bvar(const int nblocks) const ;
    

    ///////////////////////////////////////////////////////////////////////////////////////////////////
    // \brief performs a recursive blocking analysis of the statistical energy error
    //
    //
    //
    ///////////////////////////////////////////////////////////////////////////////////////////////////
    double recursive_blocking(std::ostream & fout, const bool print = true);
    
    //////////////////////////////////////////////////////////////////////////////////////////////////
    // \brief  performs a recursive blocking analysis of the target function numerator statistical error
    //
    //
    //
    //////////////////////////////////////////////////////////////////////////////////////////////////
    double target_fn_nuerr(std::ostream & fout, const bool print = false);
        
    //////////////////////////////////////////////////////////////////////////////////////////////////
    // \brief performs a recursive blocking analysis of the target function denominator statistical error
    //
    //
    //
    //////////////////////////////////////////////////////////////////////////////////////////////////
    double target_fn_dnerr(std::ostream & fout, const bool print = true);

    //////////////////////////////////////////////////////////////////////////////////////////////////
    // \brief calculate the statistical error of target function value
    //
    //
    //
    //////////////////////////////////////////////////////////////////////////////////////////////////
    double target_fn_serr(std::ostream & fout, const bool print = false);
    

    ///////////////////////////////////////////////////////////////////////////////////////////////////
    // \brief funciton that calculates average energy and target function value
    //
    //
    //
    ///////////////////////////////////////////////////////////////////////////////////////////////////
    void en_tar();

    ///////////////////////////////////////////////////////////////////////////////////////////////////
    // \brief function that add the "variance correct" term to target function 
    //
    //
    //
    ///////////////////////////////////////////////////////////////////////////////////////////////////
    void correct_finite_variance();
    
    ///////////////////////////////////////////////////////////////////////////////////////////////////
    // \brief  prints some statistics about the most recent sample 
    //
    //
    //
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    void print_statistics(std::ostream & fout);
    
    /// \brief function that returns average of local energy 
    double energy() const;

    /// \brief function that returns target function value 
    double tar_fn_val() const;

    /// \brief function that returns variance 
    double variance() const;

    /// \brief function that returns statistical error of local energy 
    double eserr(std::ostream & fout);

    /// \brief function that returns statistical error of target function value 
    double tserr(std::ostream & fout);   

    };
 }
}

#endif

