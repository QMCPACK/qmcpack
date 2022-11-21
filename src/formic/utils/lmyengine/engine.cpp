//////////////////////////////////////////////////////////////////////////////////////////////
/// \file cqmc/engine/engine.cpp
///
///
/// \brief implementation of harmonic davidson engine function
///
///
/////////////////////////////////////////////////////////////////////////////////////////////

#include <vector>
#include <string>
#include <algorithm>
#include <utility>
#include <complex>

#include <boost/shared_ptr.hpp>

#include "formic/utils/openmp.h"

#include "formic/utils/zero_one.h"
#include "engine.h"
#include "formic/utils/lmyengine/updater.h"
#include "formic/utils/lmyengine/energy_target.h"
#include "formic/utils/lmyengine/matrix_builder.h"
#include "formic/utils/lmyengine/eom.h"
#include "formic/utils/lmyengine/var_dependencies.h"
#include "formic/utils/lmyengine/engine_timing.h"
#include "formic/utils/exception.h"
#include "formic/utils/lapack_interface.h"
#include "formic/utils/mpi_interface.h"

/////////////////////////////////////////////////////////////////////////////////
/// \brief constructor with given parameters
///
///
/////////////////////////////////////////////////////////////////////////////////
template<typename S>
cqmc::engine::LMYEngine<S>::LMYEngine(const formic::VarDeps* dep_ptr,
                                   const bool exact_sampling,
                                   const bool ground,
                                   const bool variance_correction,
                                   const bool energy_print,
                                   const bool matrix_print,
                                   const bool build_lm_matrix,
                                   const bool spam,
                                   const bool var_deps_use,
                                   const bool chase_lowest,
                                   const bool chase_closest,
                                   const bool eom,
                                   const bool ssquare,
                                   const bool pm_ortho_first,
                                   const bool jas_fixed,
                                   const bool block_lm,
                                   const int num_samp,
                                   const int num_params,
                                   const int lm_krylov_iter,
                                   const int lm_spam_inner_iter,
                                   const int appro_degree,
                                   const int n_sites,
                                   const int n_pm,
                                   const int n_jas,
                                   const double hd_lm_shift,
                                   const double var_weight,
                                   const double lm_eigen_thresh,
                                   const double lm_min_S_eval,
                                   const double init_cost,
                                   const double init_var,
                                   const double lm_max_e_change,
                                   const double lm_ham_shift_i,
                                   const double lm_ham_shift_s,
                                   const double lm_max_update_abs,
                                   std::vector<double> shift_scale,
                                   std::ostream& output)
    : _exact_sampling(exact_sampling),
      _ground(ground),
      _variance_correction(variance_correction),
      _1rdm(false),
      _energy_print(energy_print),
      _matrix_print(matrix_print),
      _build_lm_matrix(build_lm_matrix),
      _spam(spam),
      _var_deps_use(var_deps_use),
      _chase_lowest(chase_lowest),
      _chase_closest(chase_closest),
      _eom(eom),
      _ssquare(ssquare),
      _pm_ortho_first(pm_ortho_first),
      _jas_fixed(jas_fixed),
      _block_lm(block_lm),
      _num_samp(num_samp),
      _num_params(num_params),
      on_hybrid_(false),
      _lm_krylov_iter(lm_krylov_iter),
      _lm_spam_inner_iter(lm_spam_inner_iter),
      _appro_degree(appro_degree),
      _n_sites(n_sites),
      _n_pm(n_pm),
      _n_jas(n_jas),
      _hd_lm_shift(hd_lm_shift),
      _var_weight(var_weight),
      _lm_eigen_thresh(lm_eigen_thresh),
      _lm_min_S_eval(lm_min_S_eval),
      _init_cost(init_cost),
      _init_var(init_var),
      _lm_max_e_change(lm_max_e_change),
      _lm_ham_shift_i(lm_ham_shift_i),
      _lm_ham_shift_s(lm_ham_shift_s),
      _lm_max_update_abs(lm_max_update_abs),
      _shift_scale(shift_scale),
      _dep_ptr(dep_ptr),
      _mbuilder(_der_rat,
                _le_der,
                _spin_der,
                _vg[0],
                _weight[0],
                hd_lm_shift,
                num_params,
                appro_degree,
                spam,
                ground,
                variance_correction,
                build_lm_matrix,
                eom,
                matrix_print),
      output(output)
{
  // get the number of threads being used 
  const int NumThreads = omp_get_max_threads();

  // check the size of lists, resize it if not correct
  if ( _le_list.size() != NumThreads ) 
    _le_list.resize(NumThreads);
  if ( _vg.size() != NumThreads ) 
    _vg.resize(NumThreads);
  if ( _weight.size() != NumThreads ) 
    _weight.resize(NumThreads);

  // initialize the output quantities
  _wfn_update = false;

  // put number of blocks to be one to avoid divide by zero
  _nblocks = 1;

  /// energy and target function statistics
  _energy  = 0.0;
  _esdev   = 0.0;
  _eserr   = 0.0;
  _target  = 0.0;
  _tserr   = 0.0;


  // solve reuslts
  _good_solve.resize(_shift_scale.size()); 
  _solved_shifts.resize(_shift_scale.size());
  std::fill(_good_solve.begin(), _good_solve.end(), false);
  std::fill(_solved_shifts.begin(), _solved_shifts.end(), false);

  // since we have not initialized the sample data, set it to be false
  _sample_initialized = false;

  // set sample count to be zero
  _samp_count = 0;

  // set the store derivative vector flag to be false
  _store_der_vec = false;

  _fully_initialized = false;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////
// \brief  Function that initialize all parameters in the engine 
// 
// 
////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename S>
void cqmc::engine::LMYEngine<S>::get_param(const formic::VarDeps * dep_ptr,
                                           const bool exact_sampling,
                                           const bool ground, 
                                           const bool var_deps_use,
                                           const bool eom,
                                           const bool ssquare,
                                           const bool block_lm,
                                           const int num_samp,
                                           const int num_params,
                                           const double hd_lm_shift,
                                           const double lm_max_e_change,
                                           const double lm_ham_shift_i,
                                           const double lm_ham_shift_s,
                                           const double lm_max_update_abs,
                                           std::vector<double> shift_scale,
                                           const bool one_rdm,
                                           const bool chase_lowest,
                                           const bool chase_closest,
                                           const bool variance_correction,
                                           const bool energy_print,
                                           const bool matrix_print,
                                           const bool build_lm_matrix,
                                           const bool spam,
                                           const bool pm_ortho_first,
                                           const bool jas_fixed,
                                           const int lm_krylov_iter,
                                           const int lm_spam_inner_iter,
                                           const int appro_degree,
                                           const int n_sites,
                                           const int n_pm,
                                           const int n_jas,
                                           const double var_weight,
                                           const double lm_eigen_thresh, 
                                           const double lm_min_S_eval,
                                           const double init_cost,
                                           const double init_var) {
  
  // set all parameters
  _exact_sampling=exact_sampling, _ground=ground, _variance_correction=variance_correction, _energy_print=energy_print, _matrix_print=matrix_print, _build_lm_matrix=build_lm_matrix;
  _spam=spam, _var_deps_use=var_deps_use, _chase_lowest=chase_lowest, _chase_closest=chase_closest, _eom=eom, _ssquare=ssquare, _pm_ortho_first=pm_ortho_first, _jas_fixed=jas_fixed;
  _sample_initialized=false, _block_lm=block_lm, _store_der_vec=false, _num_samp=num_samp, _samp_count=0, _lm_krylov_iter=lm_krylov_iter, _lm_spam_inner_iter=lm_spam_inner_iter;
  _appro_degree=appro_degree, _n_sites=n_sites, _n_pm=n_pm, _n_jas=n_jas, _hd_lm_shift=hd_lm_shift, _var_weight=var_weight, _lm_eigen_thresh=lm_eigen_thresh, _lm_min_S_eval=lm_min_S_eval;
  _init_cost=init_cost, _init_var=init_var, _lm_max_e_change=lm_max_e_change, _lm_ham_shift_i=lm_ham_shift_i, _lm_ham_shift_s=lm_ham_shift_s, _lm_max_update_abs=lm_max_update_abs;
  _wfn_update=false, _nblocks=1, _energy=0.0, _esdev=0.0, _eserr=0.0, _target=0.0, _tserr=0.0, _shift_scale=shift_scale, _dep_ptr=dep_ptr, _num_params=num_params;
  _1rdm = one_rdm;

  // solve results
  _good_solve.resize(_shift_scale.size());
  _solved_shifts.resize(_shift_scale.size());
  std::fill(_good_solve.begin(), _good_solve.end(), false);
  std::fill(_solved_shifts.begin(), _solved_shifts.end(), false);

  _mbuilder.get_param(hd_lm_shift, num_params, appro_degree, spam, ground, variance_correction, build_lm_matrix, eom, matrix_print, one_rdm);

  _fully_initialized = true;


}


////////////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief  Initialize the Engine
///
/// \param[in]    old_updates   A matrix storing old update coefficients
/// \param[in]    nblocks       Number of blocks that we will divide our total variables into
/// \param[in]    nkeps         Number of search directions per shift to kee[ in each block
/// \param[in]    nrand         Number of random vectors to add to the vector of old updates
///
////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename S>
void cqmc::engine::LMYEngine<S>::initialize(const int nblock,
                                            const int nrand,
                                            const int nkeps,
                                            const std::vector<formic::ColVec<double> > & old_updates,
                                            const bool iterative) {

  // set the block first sample finish flag to be false
  _block_first_sample_finished = false;

  // size the sample data matrices and vectors correctly if doing non-block stochastic sampling
  if ( (!_exact_sampling) && (!_block_lm) && _store_der_vec ) {
      
    // if we build matrix
    if ( _build_lm_matrix ) {

      // if we don't use variable dependency system
      if ( ! _var_deps_use ) {
        _der_rat.reset(_num_samp, _dep_ptr->n_tot()+1);
        _le_der.reset(_num_samp, _dep_ptr->n_tot()+1);
        if ( _eom ) 
          _spin_der.reset(_num_samp, _dep_ptr->n_tot()+1);
      }

      // if we don't use variable dependency system
      else {
        _der_rat.reset(_num_samp, _dep_ptr->n_ind()+1);
        _le_der.reset(_num_samp, _dep_ptr->n_ind()+1);
        if ( _eom )
          _spin_der.reset(_num_samp, _dep_ptr->n_ind()+1);
      }
    }

    // if we don't build matrix
    else {
        
      // if we don't use variable dependency system
      if ( ! _var_deps_use ) {
        _der_rat.reset(_num_samp, _dep_ptr->n_tot()+1);
        _le_der.reset(_num_samp, _dep_ptr->n_tot()+1);
        if ( _eom )
          throw formic::Exception("Must build matrix when doing EOM calculation");
      }

      // if we use variable dependency system
      else {
        _der_rat.reset(_num_samp, _dep_ptr->n_ind()+1);
        _le_der.reset(_num_samp, _dep_ptr->n_ind()+1);
        if ( _eom )
          throw formic::Exception("Must build matrix when doing EOM calculation");
      }
    }
  }

  // if we are doing block algorithm
  if ( _block_lm ) {

    _nblocks = nblock;
    _nkeps = nkeps;

    // convert old update directions into independent variable form
    _ou_ind.resize(old_updates.size());
    for (int i = 0; i < old_updates.size(); i++) {
      _ou_ind.at(i) = formic::ColVec<double>(_dep_ptr->n_ind());
      if ( old_updates.at(i).size() != _dep_ptr->n_tot() ) 
        throw formic::Exception("old_updates.at(%i).size() = %i, but should have been %i in cqmc::engine::get_brlm_update")
                      % i % old_updates.at(i).size() % _dep_ptr->n_tot();
      _dep_ptr->extract_ind_vals(old_updates.at(i).begin(), _ou_ind.at(i).begin());
      //output << "old update " << i << "   " << _ou_ind.at(i).print("%12.2e", "old update") << std::endl;
    }

    // add random vectors to the old update list
    for (int i = 0; i < nrand; i++) 
      _ou_ind.push_back(formic::randColVec<double>(_dep_ptr->n_ind()));

    // get start/stop positions and lengths for each block
    std::vector<int> block_beg, block_end, block_len;
    cqmc::engine::brlm_get_block_info(_dep_ptr->n_ind(), nblock, block_beg, block_end, block_len);

    // reset block algorithm object
    _lmb.reset(_dep_ptr->n_ind(), nblock, _ou_ind, _ground, iterative);

    // size the block updates vector correctly
    _block_ups.resize(nblock);
  }

}

////////////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief  Function that Take Sample Data from the Host Code
/// 
/// \param[in]  der_rat_samp   <n|Psi_i>/<n|Psi> (i = 0 (|Psi>), 1, ... N_var )
/// \param[in]  le_der_samp    <n|H|Psi_i>/<n|Psi> (i = 0 (|Psi>), 1, ... N_var )
/// \param[in]  ls_der_samp    <|S^2|Psi_i>/<n|Psi> (i = 0 (|Psi>), 1, ... N_var )
/// \param[in]  vgs_samp       |<n|value_fn>/<n|guiding_fn>|^2
/// \param[in]  weight_samp    weight for this sample
///
////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename S>
void cqmc::engine::LMYEngine<S>::take_sample(std::vector<S> & der_rat_samp,
                                             std::vector<S> & le_der_samp,
                                             std::vector<S> & ls_der_samp,
                                             double vgs_samp,
                                             double weight_samp) {

    
  // get the thread number 
  int myThread = omp_get_thread_num();

  // if we are doing non-block stochastic sampling, make sure that the length of input vector is the same of derivative matrix column number
  if ( ! _exact_sampling && !_block_lm && _store_der_vec ) {
    bool good_size = (der_rat_samp.size() == _der_rat.cols() && le_der_samp.size() == _le_der.cols());

    // throw out an error if not
    if ( !good_size ) 
      throw formic::Exception("der_rat_samp size is %i but cols in _der_rat is %i in cqmc::engine::LMYEngine::take_sample") % der_rat_samp.size() % _der_rat.cols();
  }

  // local energy
  _le_list[myThread].push_back(le_der_samp.at(0));
    
  // |value/guiding|^2
  _vg[myThread].push_back(vgs_samp);

  // weight
  _weight[myThread].push_back(weight_samp);

  //output << boost::format("entering take_sample function2") << std::endl;

  // if we want to do wfn update
  if ( _wfn_update ) {
    //output << "wfn_update true" << std::endl;

    // for non-block 
    if ( !_block_lm ) { 
      
      if ( !_store_der_vec ) {

        // accumulate data for the matrix builder 
        _mbuilder.take_sample(der_rat_samp,
                              le_der_samp,
                              ls_der_samp,
                              vgs_samp,
                              weight_samp);

        //output << boost::format("entering take_sample function3") << std::endl;
        return;
      }

      // if we are doing stochastic sampling 
      if ( !_exact_sampling ) {
    
        // loop over variables
        for (int i = 0; i < _der_rat.cols(); i++) {
      
          // derivative ratios
          _der_rat.at(_samp_count, i) = der_rat_samp.at(i);

          // energy derivatives
          _le_der.at(_samp_count, i)  = le_der_samp.at(i);

          // spin derivatives
          if ( _eom && _build_lm_matrix ) 
            _spin_der.at(_samp_count, i) = ls_der_samp.at(i);
        }
      }

      // if we are doing exact sampling
      else {
    
        // number of independent variables
        int n_var = _var_deps_use ? _dep_ptr->n_ind() : _dep_ptr->n_tot(); 

        // loop over variables
        for (int i = 0; i < n_var + 1; i ++) {

          // derivative ratio
          _der_rat_vec.push_back(der_rat_samp.at(i));

          // energy derivatives
          _le_der_vec.push_back(le_der_samp.at(i));

          // spin derivatives
          if ( _eom && _build_lm_matrix ) 
            _spin_der_vec.push_back(ls_der_samp.at(i));
        }
      }
    }

    // for block linear method
    else { 

      std::vector<double> der_rat_samp_real;
      std::vector<double> le_der_samp_real;
      
      for (auto it = der_rat_samp.begin(); it != der_rat_samp.end(); it++) {
        der_rat_samp_real.push_back(formic::real(*it));
      }
      for (auto it = le_der_samp.begin(); it != le_der_samp.end(); it++) {
         le_der_samp_real.push_back(formic::real(*it));
      }

      bool first_samp = !_block_first_sample_finished;

      // if it is the first sampling
      if ( first_samp ) {
        // combine vgs and weight to get effective weight
        const double d = vgs_samp * weight_samp;

        // accumulate data for the block algorithm object
        // for ground state calculation
        if ( _ground ) 
          _lmb.acc(d, der_rat_samp_real, le_der_samp_real, true);
        // for excited state calculation
        else {
          // get <n|(w-H)|Psi_i>/<n|Psi>
          std::vector<double> mle_der_samp_real(le_der_samp_real.size(), 0.0);
          
          for (int i = 0; i < mle_der_samp_real.size(); i++) 
            mle_der_samp_real.at(i) = _hd_lm_shift * der_rat_samp_real.at(i) - le_der_samp_real.at(i);

          // accumulate data
          //_lmb.acc(d, mle_der_samp, der_rat_samp);
          _lmb.acc(d, der_rat_samp_real, mle_der_samp_real, false);
        }
      }

      // if it is the second sampling
      else { 

        //if ( _samp_count % 10 == 0 ) {
        //  for (int i = 0; i < le_der_samp.size(); i++) {
        //    std::cout << le_der_samp.at(i) << "  "; 
        //  }
        //  std::cout << std::endl;
        //}
        //output << boost::format("entering second part of sampling") << std::endl;
        //output << boost::format("size of drat_cmpct is %i") % drat_cmpct.size();
        //output << boost::format("size of deng_cmpct is %i") % deng_cmpct.size();
        // get the number of special vectors(initial wfn + old updates)
        const int nsv = 1 + _lmb.ou_mat().cols();

        // get compact der rat and der eng for the current wavefunction
        if ( _ground ) {
          drat_cmpct[myThread].at(0) = der_rat_samp_real.at(0);
          deng_cmpct[myThread].at(0) = le_der_samp_real.at(0);
        }
        else {
          drat_cmpct[myThread].at(0) = der_rat_samp_real.at(0);
          deng_cmpct[myThread].at(0) = _hd_lm_shift * der_rat_samp_real.at(0) - le_der_samp_real.at(0);
        }

        // get compact der rat and der eng for old updates
        for (int k = 0; k < _lmb.ou_mat().cols(); k++) {
          drat_cmpct[myThread].at(1+k) = 0.0;
          deng_cmpct[myThread].at(1+k) = 0.0;
          for (int i = 0; i < _dep_ptr->n_ind(); i++) {
            if ( _ground ) {
              drat_cmpct[myThread].at(1+k) += _lmb.ou_mat().at(i, k) * der_rat_samp_real.at(1+i);
              deng_cmpct[myThread].at(1+k) += _lmb.ou_mat().at(i, k) * le_der_samp_real.at(1+i);
            }
            else {
              drat_cmpct[myThread].at(1+k) += _lmb.ou_mat().at(i, k) * der_rat_samp_real.at(1+i);
              deng_cmpct[myThread].at(1+k) += _lmb.ou_mat().at(i, k) * (_hd_lm_shift * der_rat_samp_real.at(1+i) - le_der_samp_real.at(1+i));
            }
          }
        }

        // get compact der rat and der eng for block directions
        // loop over blocks
        for (int b = 0, q = 0; b < _nblocks; b++) {
          // loop over shifts
          for (int s = 0; s < _block_ups.at(b).size(); s++) {
            // loop over directions for this block and shift
            for (int n = 0; n < _block_ups.at(b).at(s).cols(); n++, q++) {
              drat_cmpct[myThread].at(nsv+q) = 0.0;
              // loop over the variables within this block
              for (int i = 0; i < _lmb.bl(b); i++) {
                if ( _ground ) 
                  drat_cmpct[myThread].at(nsv+q) += _block_ups.at(b).at(s).at(i, n) * der_rat_samp_real.at(1+(_lmb.bb(b)+i));
                else 
                  drat_cmpct[myThread].at(nsv+q) += _block_ups.at(b).at(s).at(i, n) * der_rat_samp_real.at(1+(_lmb.bb(b)+i));
              }
              deng_cmpct[myThread].at(nsv+q) = 0.0;
              for (int i = 0; i < _lmb.bl(b); i++) {
                if ( _ground ) 
                  deng_cmpct[myThread].at(nsv+q) += _block_ups.at(b).at(s).at(i, n) * le_der_samp_real.at(1+(_lmb.bb(b)+i));
                else 
                  deng_cmpct[myThread].at(nsv+q) += _block_ups.at(b).at(s).at(i, n) * (_hd_lm_shift * der_rat_samp_real.at(1+(_lmb.bb(b)+i)) - le_der_samp_real.at(1+(_lmb.bb(b)+i)));
              }
            }
          }
        }

        // get contributions to hamiltonian and overlap matrices
        const double d = vgs_samp * weight_samp;
        for (int j = 0; j < hh_block[myThread].cols(); j++) {
          for (int i = 0; i < hh_block[myThread].rows(); i++) {
            hh_block[myThread].at(i,j) += d * drat_cmpct[myThread].at(i) * deng_cmpct[myThread].at(j);
            if ( _ground ) 
              ss_block[myThread].at(i,j) += d * drat_cmpct[myThread].at(i) * drat_cmpct[myThread].at(j);
            else 
              ss_block[myThread].at(i,j) += d * deng_cmpct[myThread].at(i) * deng_cmpct[myThread].at(j);
          }
        }
      }
    }
  }
  
  // increase sample count by 1
  # pragma omp critical
  _samp_count++; 

}

////////////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief  Function that Take Sample Data from the Host Code
/// 
/// \param[in]  local_en       local energy
/// \param[in]  vgs_samp       |<n|value_fn>/<n|guiding_fn>|^2
/// \param[in]  weight_samp    weight for this sample
///
////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename S>
void cqmc::engine::LMYEngine<S>::take_sample(S local_en,
                                             double vgs_samp,
                                             double weight_samp) 
{

  // get thread number 
  int myThread = omp_get_thread_num();

  // local energy 
  _le_list[myThread].push_back(local_en);

  // |value/guiding|^2
  _vg[myThread].push_back(vgs_samp);

  // weight 
  _weight[myThread].push_back(weight_samp);
  
  // increase sample count by 1 
  # pragma omp critical
  _samp_count++;

}

////////////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief  Function that Take Sample Data from the Host Code
/// 
/// \param[in]  local_en       local energy
/// \param[in]  vgs_samp       |<n|value_fn>/<n|guiding_fn>|^2
/// \param[in]  weight_samp    weight for this sample
///
////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename S>
void cqmc::engine::LMYEngine<S>::take_sample(int nbasis, std::vector<S> & one_rdm_samp)
{

  // get thread number 
  int myThread = omp_get_thread_num();

  // accumulate data for the matrix builder
  _mbuilder.take_sample(nbasis, one_rdm_samp);

}

////////////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief  Function that reduces all block matrices information from all processors to the root
///         processor
///
////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename S>
void cqmc::engine::LMYEngine<S>::sample_finish() {
 
  int my_rank = formic::mpi::rank();

  // get total number of threads
  int NumThreads = omp_get_max_threads();

  // evaluate total weight 
  double _tw = 0.0;
  for (int ip = 0; ip < NumThreads; ip++) {
    for (int i = 0; i < _weight[ip].size(); i++) {
      _tw += _weight[ip].at(i) * _vg[ip].at(i);
    }
  }

  double total_weight = 0.0;
  formic::mpi::allreduce(&_tw, &total_weight, 1, MPI_SUM);

  // for energy evaluation only calculation, do nothing
  if ( !_wfn_update )
    return;

  // if we are doing non-block
  if ( !_block_lm ) {

    if ( !_store_der_vec ) {

      // finish sample for the matrix builder 
      _mbuilder.finish_sample(total_weight);

    }

    // if we are doing exact sampling, build der_rat and le_der matrices out of vectors
    if ( _exact_sampling ) {
      // number of variables
      const int n_var = _var_deps_use ? _dep_ptr->n_ind() : _dep_ptr->n_tot();

      // size the matrix correctly
      _der_rat.reset(_samp_count, n_var+1);
      _le_der.reset(_samp_count, n_var+1);
      if ( _eom && _build_lm_matrix ) 
        _spin_der.reset(_samp_count, n_var+1);

      // build the matrix
      for (int i = 0; i < _der_rat.rows(); i++) {
        for (int j = 0; j < _der_rat.cols(); j++) {
          _der_rat.at(i,j) = _der_rat_vec.at(i * _der_rat.cols() + j);
          _le_der.at(i,j)  = _le_der_vec.at(i * _der_rat.cols() + j);
          if ( _eom && _build_lm_matrix ) 
            _spin_der.at(i,j) = _spin_der_vec.at(i * _der_rat.cols() + j);
        }
      }

      // clear old data
      _der_rat_vec.clear();
      _le_der_vec.clear();
      _spin_der_vec.clear();
    }
  }

  // if block algorithm
  else {

    const bool first = !_block_first_sample_finished;
    
    if ( first ) {
      
      // get total weight on this sample
      double total_weight = _lmb.total_weight();

      // the total weight through all samples
      double all_samp_weight = 0.0;

      // mpi all reduce
      formic::mpi::allreduce(&total_weight, &all_samp_weight, 1, MPI_SUM);

      // call the finalize function for the block algorithm object
      _lmb.mpi_finalize(all_samp_weight);

      // say that we finished the sample
      _block_first_sample_finished = true;
    }

    else {

      // the total weight on this sample 
      double tot_weight = 0.0;

      for (int ip = 0; ip < NumThreads; ip++) {
        for (int i = 0; i < _weight[ip].size(); i++) {
          tot_weight += _vg[ip].at(i) * _weight[ip].at(i);
        }
      }

      // the total weight through all samples
      double all_samp_weight = 0.0;

      // mpi all reduce
      formic::mpi::allreduce(&tot_weight, &all_samp_weight, 1, MPI_SUM);
      
      // sum over threads for block matrices 
      for (int ip = 1; ip < NumThreads; ip++) {
        hh_block[0] += hh_block[ip];
        ss_block[0] += ss_block[ip];
      }

      // get space for mpi reduce
      formic::Matrix<double> hh_block_tot(hh_block[0].rows(), hh_block[0].cols());
      formic::Matrix<double> ss_block_tot(ss_block[0].rows(), ss_block[0].cols());

      // compute average of matrices
      formic::mpi::reduce(&hh_block[0].at(0,0), &hh_block_tot.at(0,0), hh_block[0].size(), MPI_SUM);
      formic::mpi::reduce(&ss_block[0].at(0,0), &ss_block_tot.at(0,0), ss_block[0].size(), MPI_SUM); 

      // compute average on root process
      if ( my_rank == 0 ) {
        hh_block[0] = hh_block_tot / all_samp_weight;
        ss_block[0] = ss_block_tot / all_samp_weight;
      }
      //output << hh_block[0].print("%12.6f", "hh_block");
      //output << ss_block[0].print("%12.6f", "ss_block");

      // say that we need to redo the sample
      _block_first_sample_finished = false;
    }
  }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief  Function that updates variable dependency pointer 
///
////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename S>
void cqmc::engine::LMYEngine<S>::var_deps_ptr_update(const formic::VarDeps * new_dep_ptr) {
  
  // update pointer 
  _dep_ptr = new_dep_ptr;

  // call initialize function to size matrices correctly
  if ( !_block_lm ) 
    this->initialize();
}

////////////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief  Energy and Target Function Calculation
///
////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename S>
void cqmc::engine::LMYEngine<S>::energy_target_compute() { 

  // get the number of threads 
  int NumThreads = omp_get_max_threads();

  for (int ip = 1; ip < NumThreads; ip++) {
    for (int i = 0; i < _le_list[ip].size(); i++) {
      _le_list[0].push_back(_le_list[ip].at(i));
      _vg[0].push_back(_vg[ip].at(i));
      _weight[0].push_back(_weight[ip].at(i));
    }
  }
  // call the corresponding "call engine" function
  this->call_engine(_exact_sampling,
                    _ground,
                    _variance_correction,
                    _energy_print,
                    _hd_lm_shift,
                    _var_weight,
                    _le_list[0],
                    _vg[0],
                    _weight[0],
                    _energy,
                    _esdev,
                    _eserr,
                    _target,
                    _tserr,
                    output);
  
  // set initial energy and variance
  if ( _ground )
    _init_cost = _energy;
  else
    _init_cost = _target;

  _init_var    = _esdev * _esdev;
  //std::cout << "engine energy target compute finish" << std::endl;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief  Do necessary preparations for the wavefunction update 
///
////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename S>
void cqmc::engine::LMYEngine<S>::wfn_update_prep() {

  // if we are not doing block algorithm, do nothing
  if ( !_block_lm ) 
    return;
  else 
    this->get_brlm_update_alg_part_one(_dep_ptr,
                                       _nkeps,
                                       _lm_ham_shift_i,
                                       _lm_ham_shift_s,
                                       _shift_scale,
                                       _lm_max_update_abs,
                                       _good_solve,
                                       _solved_shifts,
                                       output);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief  Wave Function Update Calculation
///
////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename S>
void cqmc::engine::LMYEngine<S>::wfn_update_compute() {

  // if we are doing non-block, call the corresponding "call_engine" function
  if ( !_block_lm ) {
    this->call_engine(_dep_ptr, _ground, _build_lm_matrix, _variance_correction, _matrix_print, 
                      _lm_krylov_iter, _lm_spam_inner_iter, _appro_degree, _lm_eigen_thresh,
                      _lm_min_S_eval, _spam, _var_deps_use, _chase_lowest, _chase_closest,
                      _init_cost, _init_var, _lm_max_e_change, _lm_ham_shift_i, 
                      _lm_ham_shift_s, _hd_lm_shift, _var_weight, _lm_max_update_abs,
                      _der_rat, _le_der, _vg[0], _weight[0], _vf_var, _good_solve, _solved_shifts, 
                      _shift_scale, output);
  }

  // if we are doing the block version
  else {
    // call the second part of block algorithm, this assumes that the first part of the algorithm has beeb called 
    std::vector<double> _vf_var_real(_vf_var.size(), 0.0);
    this -> get_brlm_update_alg_part_two(_dep_ptr,
                                         _lmb.iterative(),
                                         _nkeps,
                                         _lm_ham_shift_i,
                                         _lm_ham_shift_s,
                                         _shift_scale,
                                         _lm_max_update_abs,
                                         _vf_var_real,
                                         _good_solve,
                                         _solved_shifts,
                                         output);
    for (std::vector<double>::iterator it=_vf_var_real.begin(); it!=_vf_var_real.end(); it++) {
      _vf_var.push_back(formic::unity(S()) * (*it));
    }
  }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief  Equation-of-Motion Calculation
///
////////////////////////////////////////////////////////////////////////////////////////////////////////           
template<typename S>
void cqmc::engine::LMYEngine<S>::eom_compute() {

  // call the corresponding "call engine" function
  this->call_engine(_matrix_print,
                    _ground, 
                    _eom,
                    _ssquare,
                    _pm_ortho_first,
                    _jas_fixed,
                    1.0e-6,
                    _n_sites,
                    _n_pm,
                    _n_jas,
                    _der_rat, 
                    _le_der,
                    _spin_der,
                    _evecs,
                    _energy_index,
                    _vg[0],
                    _weight[0],
                    output);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief  Reset the Object
///
////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename S>
void cqmc::engine::LMYEngine<S>::reset() {

  // clear energy and target statistics
  //_energy = 0.0;
  //_esdev  = 0.0;
  //_eserr  = 0.0;
  //_target = 0.0;
  //_tserr  = 0.0;
  
  // number of threads 
  int NumThreads = omp_get_max_threads();

  // clear wavefunction update vector
  _vf_var.clear();

  // set solve results to false
  std::fill(_good_solve.begin(), _good_solve.end(), false);
  std::fill(_solved_shifts.begin(), _solved_shifts.end(), false);

  // clear eom energy and eigenvectors
  _evecs.clear();
  _energy_index.clear();

  // set the sample count to be zero
  _samp_count = 0;

  // clear local energy, vgs and weight list
  for (int ip = 0; ip < NumThreads; ip++) {
    _le_list[ip].clear();
    _vg[ip].clear();
    _weight[ip].clear();
  }

  // reset the block algorithm object
  if ( _block_lm && _wfn_update ) 
    _lmb.reset(_dep_ptr->n_ind(), _nblocks, _ou_ind, _ground); 

  // reset the matrix builder 
  _mbuilder.reset();

}

template<typename S>
void cqmc::engine::LMYEngine<S>::shift_update(std::vector<double> & new_shift) {

  // update shift
  for (int i = 0; i < new_shift.size(); i++) 
    _shift_scale.at(i) = new_shift.at(i);

}




//Transfers the vectors from descent to the engine's LMBlocker object during the hybrid method.
template<typename S>
void cqmc::engine::LMYEngine<S>::setHybridBLM_Input(std::vector< std::vector<double>>& from_descent) 
{

    //Change the LMBlocker object's hybrid variable to true so input vectors will be used later on
    _lmb.setHybrid(true);

    //Clear to avoid retaining old sets of vectors
    _lmb.getInputVector().clear();
    
    for (std::vector<double> v : from_descent)
    {

	_lmb.getInputVector().push_back(v);
    }
}

/////////////////////////////////////////////////////////////////////////////////////////////
/// \brief harmonic davidson engine call to calculate energy and target function 
///
///
///
///
///
/////////////////////////////////////////////////////////////////////////////////////////////
template<typename S>
void cqmc::engine::LMYEngine<S>::call_engine(const bool exact_sampling,
                                             const bool ground_state,
                                             const bool variance_correct,
                                             const bool print,
                                             const double hd_lm_shift,
                                             const double var_weight,
                                             std::vector<S> & le_list,
                                             std::vector<double> & vg,
                                             std::vector<double> & weight,
                                             double & energy,
                                             double & esdev, 
                                             double & eserr,
                                             double & target, 
                                             double & tserr,
                                             std::ostream & output) {

  // simply, just call energy and target function calculation function
  cqmc::engine::et(exact_sampling,
                   ground_state,
                   variance_correct,
                   print,
                   hd_lm_shift,
                   var_weight,
                   le_list,
                   vg,
                   weight,
                   energy,
                   esdev,
                   eserr,
                   target,
                   tserr,
                   output);

}


////////////////////////////////////////////////////////////////////////////////////////////
// \brief harmonic davidson engine call to update wave function
//
//
//
//
//
////////////////////////////////////////////////////////////////////////////////////////////
template<typename S>
void cqmc::engine::LMYEngine<S>::call_engine(const formic::VarDeps * dep_ptr, 
                                             const bool ground_state,
                                             const bool build_lm_matrix,
                                             const bool variance_correct,
                                             const bool print_matrix,
                                             const int lm_krylov_iter,
                                             const int lm_spam_inner_iter,
                                             const int appro_degree,
                                             const double lm_eigen_thresh,
                                             const double lm_min_S_eval,
                                             const bool spam_use,
                                             const bool var_deps_use,
                                             const bool chase_lowest,
                                             const bool chase_closest,
                                             const double init_cost,
                                             const double init_var,
                                             const double lm_max_e_change,
                                             const double lm_ham_shift_i,
                                             const double lm_ham_shift_s,
                                             const double omega,
                                             const double var_weight,
                                             const double lm_max_update_abs,
                                             formic::Matrix<S> & der_rat,
                                             formic::Matrix<S> & le_der,
                                             std::vector<double> & vg, 
                                             std::vector<double> & weight,
                                             std::vector<S> & vf_var,
                                             std::vector<bool> & good_solve,
                                             std::vector<int> & shift_solved,
                                             const std::vector<double> & shift_scale, 
                                             std::ostream & output)
{
  // start timer 
  cqmc::start_timer("eigen solver");

  // get rank number
  int my_rank = formic::mpi::rank();
  //MPI_Comm_rank(MPI_COMM_WORLD, & my_rank);

  // get the number of input shift
  const int num_shift = shift_scale.size();

  // size vector storing wavefunction update variable correctly
  const int correct_vf_var_size = (dep_ptr->n_tot()+1) * num_shift;
  if (vf_var.size() != correct_vf_var_size) 
    vf_var.resize(correct_vf_var_size);
  //if ( var_deps_use ) {
  //  const int correct_vf_var_size = (dep_ptr->n_tot()+1) *  num_shift;
  //  if (vf_var.size() != correct_vf_var_size) 
  //    vf_var.resize(correct_vf_var_size);
  //}

  //else {
  //  const int correct_vf_var_size = der_rat.cols() * num_shift;
  //  if (vf_var.size() != correct_vf_var_size) 
  //    vf_var.resize(correct_vf_var_size);
  //}

  // size vector storing solve results correctly
  if (good_solve.size() != num_shift)
    good_solve.resize(num_shift);

  // create harmoinc davidson updater
  boost::shared_ptr< cqmc::engine::HDLinMethodUpdater<S> > hdupdater ( new cqmc::engine::HDLinMethodUpdater<S>(der_rat, 
                                                                                                               le_der,
                                                                                                               vg, 
                                                                                                               weight,
                                                                                                               shift_scale,
                                                                                                               omega,
                                                                                                               var_weight,
                                                                                                               ground_state,
                                                                                                               variance_correct,
                                                                                                               build_lm_matrix));


   if ( !_store_der_vec ) {
     
     // if we don't use unbiased matrix builder 
     hdupdater -> engine_update_build_matrix(dep_ptr,
                                             lm_krylov_iter, 
                                             lm_spam_inner_iter,
                                             lm_eigen_thresh,
                                             lm_min_S_eval,
                                             spam_use,
                                             var_deps_use,
                                             chase_lowest,
                                             chase_closest,
                                             print_matrix,
                                             init_cost,
                                             init_var,
                                             lm_max_e_change,
                                             lm_ham_shift_i,
                                             lm_ham_shift_s,
                                             omega,
                                             lm_max_update_abs,
                                             _mbuilder.ham(),
                                             _mbuilder.ovl(),
                                             _mbuilder.ssquare(),
                                             vf_var,
                                             good_solve,
                                             shift_solved,
                                             output); 

    // end timer 
    cqmc::stop_timer("eigen solver");

    // print timer
    if ( my_rank == 0 ) 
      output << cqmc::print_timers();
  
    return; 
  }

  // call update function 
  if ( build_lm_matrix && _store_der_vec ) 
    hdupdater -> engine_update_build_matrix(dep_ptr,
                                            lm_krylov_iter, 
                                            lm_spam_inner_iter,
                                            lm_eigen_thresh,
                                            lm_min_S_eval,
                                            spam_use,
                                            var_deps_use,
                                            chase_lowest,
                                            chase_closest,
                                            print_matrix,
                                            init_cost,
                                            init_var,
                                            lm_max_e_change,
                                            lm_ham_shift_i,
                                            lm_ham_shift_s,
                                            omega,
                                            lm_max_update_abs,
                                            vf_var,
                                            good_solve,
                                            shift_solved,
                                            output);
  
  else if ( !build_lm_matrix && !spam_use )
    hdupdater -> engine_update_no_matrix(dep_ptr, 
                                         lm_krylov_iter,
                                         lm_spam_inner_iter,
                                         lm_eigen_thresh,
                                         lm_min_S_eval,
                                         spam_use,
                                         var_deps_use,
                                         chase_lowest,
                                         chase_closest,
                                         print_matrix,
                                         init_cost,
                                         lm_max_e_change,
                                         lm_ham_shift_i,
                                         lm_ham_shift_s,
                                         omega,
                                         lm_max_update_abs,
                                         vf_var,
                                         good_solve,
                                         shift_solved,
                                         output);

  else if ( !build_lm_matrix && spam_use ) {
    std::vector<double> vf_var_real(vf_var.size(), 0.0);
    hdupdater -> engine_update_spam(dep_ptr, 
                                    lm_krylov_iter,
                                    lm_spam_inner_iter,
                                    appro_degree,
                                    lm_eigen_thresh,
                                    lm_min_S_eval,
                                    spam_use,
                                    var_deps_use,
                                    chase_lowest,
                                    chase_closest,
                                    print_matrix,
                                    init_cost,
                                    lm_max_e_change,
                                    lm_ham_shift_i,
                                    lm_ham_shift_s,
                                    omega,
                                    lm_max_update_abs,
                                    vf_var_real,
                                    good_solve,
                                    shift_solved,
                                    output);
    vf_var.clear();
    for (auto it=vf_var_real.begin(); it != vf_var_real.end(); it++)
      vf_var.push_back(formic::unity(S())*(*it));
  }



}


////////////////////////////////////////////////////////////////////////////////////////////////////////////
// \brief  equation of motion JAGP function 
//
//
//
//
//
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename S>
void cqmc::engine::LMYEngine<S>::call_engine(const bool print_matrix, 
                                             const bool ground_state,
                                             const bool eom,
                                             const bool ssquare,
                                             const bool pm_ortho_first,
                                             const bool jas_fixed,
                                             const double eom_inv_thresh,
                                             int n_sites,
                                             int n_pm,
                                             int n_jas,
                                             formic::Matrix<S> & der_rat,
                                             formic::Matrix<S> & le_der,
                                             formic::Matrix<S> & spin_der,
                                             formic::Matrix<std::complex<double> > & evecs,
                                             std::vector<std::pair<std::complex<double>, int> > & energy_index,
                                             std::vector<double> & vg,
                                             std::vector<double> & weight,
                                             std::ostream & output)
{
  
  // if this calculation is not ground state calculation, EOM will not be performed
  if ( !ground_state )
    return;

  // if eom flag is set to be false, EOM will not be performed 
  if ( !eom ) 
    return;

  // creates matrix builder 
  boost::shared_ptr< cqmc::engine::HamOvlpBuilderHD<S> > mbuilder( new cqmc::engine::HamOvlpBuilderHD<S>(der_rat, 
                                                                                                         le_der, 
                                                                                                         spin_der, 
                                                                                                         vg, 
                                                                                                         weight, 
                                                                                                         0.0, 
                                                                                                         _num_params,
                                                                                                         100, 
                                                                                                         false, 
                                                                                                         true, 
                                                                                                         false,
                                                                                                         true,
                                                                                                         ssquare, 
                                                                                                         print_matrix));

  // build the matrix 
  mbuilder -> MatrixBuild(output);

  // test the diagonalization
  //eric_test_diag(mbuilder -> ham(), mbuilder -> ovl(), output);

  // analyze derivative vectors 
  //mbuilder -> derivative_analyze(); 

  // creates eom calculator 
  boost::shared_ptr< cqmc::engine::EOM<S> > eom_calculator( new cqmc::engine::EOM<S>(mbuilder -> ham(), 
                                                                                     mbuilder -> ovl(), 
                                                                                     mbuilder -> ssquare(), 
                                                                                     ssquare, 
                                                                                     print_matrix, 
                                                                                     pm_ortho_first, 
                                                                                     jas_fixed, 
                                                                                     n_sites, 
                                                                                     n_pm, 
                                                                                     n_jas, 
                                                                                     eom_inv_thresh));

  // perform eom calculation 
  //eom_calculator -> eom_calculation(output);
  eom_calculator -> eom_calculation_simple(output);

  // get output data
  energy_index = eom_calculator -> energy_list();
  evecs = eom_calculator -> evecs();

  // print eom statistics 
  eom_calculator -> eom_print(output);
  
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief  Computes a wave function update using block recursive linear method algorithm, 
///         which solves each block many times, each time including some information from a second
///         block in the form of previous update components. 
///         These variale solves give a set of desirable vectors within the block,
///         An SVD is then used to pick a small number of vectors that are best able to recreate
///         all of these desirable vectors
///         Note: This function assumes that all block matrices data has been accumulated 
///               and correctly reduced to root process!!!
///
/// \param[in]      deps               object keeping track of variable dependencies
/// \param[in]      nblock             number of blocks to divide the variables into
/// \param[in]      nkps               number of search directions per shift to keep in each block
/// \param[in]      nrand              number of random vectors to add to the vector of old updates
/// \param[in]      shift_i            strength of the identity shift
/// \param[in]      shift_s            strength of the overlap shift
/// \param[in]      shift_scale        scale factors telling which shifts to use
/// \param[in]      max_update_abs     reject
///
/////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename S>
void cqmc::engine::LMYEngine<S>::get_brlm_update_alg_part_one(const formic::VarDeps * deps,
                                                              const int nkps,
                                                              const double shift_i,
                                                              const double shift_s,
                                                              const std::vector<double> & shift_scale,
                                                              const double max_update_abs,
                                                              std::vector<bool> & good_solve,
                                                              std::vector<int> & shift_solved, 
                                                              std::ostream & output) {

  // get the number of threads
  int NumThreads = omp_get_max_threads();

  // prepare vectors telling which shifts are solved healthily
  while ( good_solve.size() < shift_scale.size() )
    good_solve.push_back(false);
  while ( shift_solved.size() < shift_scale.size() ) 
    shift_solved.push_back(0);

  // get a convenient variable for the number of special vectors (i.e. the initial wavefunction plus the number of old updates directions)
  const int nsv = 1 + _lmb.ou_mat().cols();

  // choose how many updates to take from each block for each shift
  const int n_up_per_shift = 1;

  // choose whether to make use of previoud updates
  const bool use_prev_ups = true;

  // prepare an object to hold update directions for each block
  // block_ups[i][j](p,q) refers to the pth element of the qth update vector for the jth shift for the ith block
  //std::vector<std::vector<formic::Matrix<double> > > block_ups(_nblock);

 
  // compute the useful update directions on root process
  _lmb.solve_for_block_dirs(_dep_ptr, nkps, shift_i, shift_s, shift_scale, _block_ups, output, _hd_lm_shift);

  // get the total number of directions involved in the final basis
  const int n_dir = nsv + nkps * shift_scale.size() * _nblocks;

  // build Hamiltonian and overlap matrix in the basis of good update directions
  // first size the matrices and vectors correctly
  hh_block.resize(NumThreads);
  ss_block.resize(NumThreads);
  drat_cmpct.resize(NumThreads);
  deng_cmpct.resize(NumThreads);
  for (int ip = 0; ip < NumThreads; ip++) {
    hh_block[ip].reset(n_dir, n_dir, 0.0);
    ss_block[ip].reset(n_dir, n_dir, 0.0);
    drat_cmpct[ip].reset(n_dir, 0.0);
    deng_cmpct[ip].reset(n_dir, 0.0);
  }
 
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief  Computes a wave function update using block recursive linear method algorithm, 
///         this function solves the eigenvalue problem within the final basis 
///
/// \param[in]      deps               object keeping track of variable dependencies
/// \param[in]      iterative          whether to use iterative method to solve 
/// \param[in]      nblock             number of blocks to divide the variables into
/// \param[in]      nkps               number of search directions per shift to keep in each block
/// \param[in]      nrand              number of random vectors to add to the vector of old updates
/// \param[in]      shift_i            strength of the identity shift
/// \param[in]      shift_s            strength of the overlap shift
/// \param[in]      shift_scale        scale factors telling which shifts to use
/// \param[in]      max_update_abs     reject
///
/////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename S>
void cqmc::engine::LMYEngine<S>::get_brlm_update_alg_part_two(const formic::VarDeps * deps,
                                                              const bool iterative,
                                                              const int nkps,
                                                              const double shift_i,
                                                              const double shift_s,
                                                              const std::vector<double> & shift_scale,
                                                              const double max_update_abs,
                                                              std::vector<double> & updates,
                                                              std::vector<bool> & good_solve,
                                                              std::vector<int> & shift_solved, 
                                                              std::ostream & output) {

  int my_rank = formic::mpi::rank();

  // size update vector correctly
  updates.assign(shift_scale.size() * ( 1 + _dep_ptr->n_tot()), 0.0);

  // get the number of special vectors(initial wfn + old updates)
  const int nsv = 1 + _lmb.ou_mat().cols();

  // on root process
  if ( my_rank == 0 ) {
    
    // get the matrix for the identity shift in this basis 
    formic::Matrix<double> dd(hh_block[0].rows(), hh_block[0].cols(), 0.0);
    for (int bl = 0, ql = 0; bl < _nblocks; bl++) { // loop over blocks for left vector
      for (int sl = 0; sl < _block_ups.at(bl).size(); sl++) { // loop over shifts within a block for left vector
        for (int nl = 0; nl < _block_ups.at(bl).at(sl).cols(); nl++, ql++) { // loop over updates within a shift for left vector
          for (int br = 0, qr = 0; br < _nblocks; br++) { // loop over blocks for right vector
            for (int sr = 0; sr < _block_ups.at(br).size(); sr++) { // loop over shifts within a block for right vector
              for (int nr = 0; nr < _block_ups.at(br).at(sr).cols(); nr++, qr++) { // loop over updates within a shift for right vector
                for (int i = 0; i < _lmb.bl(bl) && bl == br; i++) 
                  dd.at(nsv+ql, nsv+qr) += _block_ups.at(bl).at(sl).at(i, nl) * _block_ups.at(br).at(sr).at(i,nr);
              }
            }
          }
        }
      }
    }
    for (int br = 0, qr = 0; br < _nblocks; br++) { // loop over blocks for right vector
      for (int sr = 0; sr < _block_ups.at(br).size(); sr++) { // loop over shifts within a block for right vector
        for (int nr = 0; nr < _block_ups.at(br).at(sr).cols(); nr++, qr++) { // loop over updates within a shift for right vector
          for (int j = 0; j < _lmb.ou_mat().cols(); j++) {
            for (int i = 0; i < _lmb.bl(br); i++) {
              dd.at(1+j,nsv+qr) += _lmb.ou_mat().at((_lmb.bb(br)+i), j) * _block_ups.at(br).at(sr).at(i,nr);
              dd.at(nsv+qr,1+j) += _lmb.ou_mat().at((_lmb.bb(br)+i), j) * _block_ups.at(br).at(sr).at(i,nr);
            }
          }
        }
      }
    }

    for (int k = 0; k < _lmb.ou_mat().cols(); k++)
    for (int l = 0; l < _lmb.ou_mat().cols(); l++) {
      dd.at(1+k,1+l) = 0.0;
      for (int i = 0; i < _dep_ptr->n_ind(); i++) 
        dd.at(1+k,1+l) += _lmb.ou_mat().at(i, k) * _lmb.ou_mat().at(i, l);
    }

    //std::cout << "dd constructed" << std::endl;

    // for each shift, we solve the linear method eigenvalue problem for the update in this basis
    //updates.assign(shift_scale.size() * ( 1 + _dep_ptr->n_tot()), 0.0);
    //std::cout << "before matrix wrapper constructed" << updates.size() << std::endl;
    formic::Matrix<double> updates_matrix( 1 + _dep_ptr->n_tot(), shift_scale.size(), &updates.at(0)); // matrix wrapper
    //std::cout << "matrix wrapper constructed" << std::endl;
    // loop over shifts
    for (int shift_p = 0; shift_p < shift_scale.size(); shift_p++) {
     
      //std::cout << "entering loop over shifts" << std::endl;
      // compute the update in the compact basis 
      formic::Matrix<double> update_dir;
      if ( !_lmb.iterative() ) 
        update_dir = cqmc::engine::get_important_brlm_dirs(1, 
                                                           formic::real(this->target_value()),
                                                           shift_i * shift_scale.at(shift_p),
                                                           shift_s * shift_scale.at(shift_p),
                                                           1.0e-4,
                                                           hh_block[0],
                                                           ss_block[0],
                                                           dd,
                                                           output);
      else 
        update_dir = cqmc::engine::get_important_brlm_dirs_davidson(_dep_ptr,
                                                                    1, 
                                                                    _hd_lm_shift,
                                                                    formic::real(_energy),
                                                                    shift_i * shift_scale.at(shift_p),
                                                                    shift_s * shift_scale.at(shift_p),
                                                                    1.0e-6,
                                                                    hh_block[0],
                                                                    ss_block[0],
                                                                    dd,
                                                                    output);

      //output << update_dir.print("%12.2e", "final compact basis update for this shift") << std::endl;

      // prepare vector to hold the update direction in independent variable form
      formic::ColVec<double> iv_up(_dep_ptr->n_ind(), 0.0);

      // expand the part of the update that comes from the different variable blocks
      for (int b = 0, q = 0; b < _nblocks; b++) {
        for (int s = 0; s < _block_ups.at(b).size(); s++) {
          for (int n = 0; n < _block_ups.at(b).at(s).cols(); n++, q++) {
            for (int i = 0; i < _lmb.bl(b); i++) 
               iv_up.at(_lmb.bb(b)+i) += _block_ups.at(b).at(s).at(i, n) * update_dir.at(nsv+q,0);
          }
        }
      }

      // expand the part of the update that comes from previous update directions
      for (int k = 0; k < _lmb.ou_mat().cols(); k++) 
      for (int i = 0; i < _dep_ptr->n_ind(); i++) 
        iv_up.at(i) += _lmb.ou_mat().at(i,k) * update_dir.at(1+k,0);

      // expand the update into the basis of dependent variables
      updates_matrix.at(0,shift_p) = 1.0;
      _dep_ptr->expand_ind_to_all(iv_up.begin(), updates_matrix.col_begin(shift_p)+1);

      // say that the solve for this shift worked
      good_solve.at(shift_p) = true;
      shift_solved.at(shift_p) = 1;
    }

    output << std::endl;

    //output << updates_matrix.print("%12.2e", "final full basis matrix of updates") << std::endl;
  }

  // broadcast solve results and wavefunction updates to all processes
  formic::mpi::bcast(&updates.at(0), updates.size());
  //for (int i = 0; i < updates.size(); i++) 
  //  output << boost::format("%12.6f ") % updates.at(i);

}

//Function for setting the size of the derivative ratio histories according to the total number of optimizable parameters and the number of samples per process
 template<typename S>
 void cqmc::engine::LMYEngine<S>::setUpStorage(int numParams, int numSamples)
 {

     le_der_rat_history.reset(numSamples, numParams+1);
     der_rat_history.reset(numSamples, numParams+1);
     vgs_history.resize(numSamples);
     weight_history.resize(numSamples);
     
 }

//Function for storing a sample (local energy and the sets of derivative ratios)
 template<typename S>
 void cqmc::engine::LMYEngine<S>::store_sample(std::vector<S> & der_rat_samp,std::vector<S> & le_der_samp,std::vector<S> & ls_der_samp,double vgs_samp,double weight_samp,int sample_index)
 {
     for(int i = 0; i < le_der_samp.size();i++)
     {
         le_der_rat_history.at(sample_index,i) = le_der_samp[i];
         der_rat_history.at(sample_index,i) = der_rat_samp[i];
     }
    
    vgs_history[sample_index] = vgs_samp;
    weight_history[sample_index] = weight_samp;
 }

//Function for constructing matrices form stored samples
 //It uses the engine's existing take_sample function, which is where the matrices are actually built
 //This function only passes along the derivative ratios that are stored or a subset of them if parameters are being filtered
 template<typename S>
 void cqmc::engine::LMYEngine<S>::buildMatricesFromDerivatives()
 {
     //Number of samples on a process
     int num_samples = der_rat_history.rows();
     int der_vec_len = der_rat_history.cols();


     for(int i = 0; i < num_samples; i++)
     {
         std::vector<S> der_rat_samp;
         std::vector<S> le_der_samp;

         for(int j = 0;j < der_vec_len;j++)
         {
             der_rat_samp.push_back(der_rat_history.at(i,j));
             le_der_samp.push_back(le_der_rat_history.at(i,j));
         }

         double vgs = vgs_history.at(i);
         double weight = weight_history.at(i);

         if(filter_param_)
         {
             std::vector<S> reduced_der_rat_samp;
             std::vector<S> reduced_le_der_samp;

             reduced_der_rat_samp.push_back(der_rat_samp[0]);
             reduced_le_der_samp.push_back(le_der_samp[0]);

             for(int i = 1; i< der_rat_samp.size();i++)
             {
                 if(parameterSettings[i-1] == true)
                 {

                     reduced_der_rat_samp.push_back(der_rat_samp[i]);
                     reduced_le_der_samp.push_back(le_der_samp[i]);
                 }
             }

             this->take_sample(reduced_der_rat_samp,reduced_le_der_samp,reduced_le_der_samp,vgs,weight);

         }
         else
          {
             this->take_sample(der_rat_samp,le_der_samp,le_der_samp,vgs,weight);
             }


     }
     this->sample_finish();

 }

//Function for clearing stored samples
 template<typename S>
 void cqmc::engine::LMYEngine<S>::clear_histories()
 {
     int my_rank = formic::mpi::rank();

         der_rat_history.clear();
         le_der_rat_history.clear();
         vgs_history.clear();
         weight_history.clear();
         if(my_rank == 0)
         {
             std::cout << "Matrices have been built, clearing stored samples for this iteration" << std::endl;
         }


 }

//Function for changing the parameter number used by the engine
 template<typename S>
 void cqmc::engine::LMYEngine<S>::resetParamNumber(int new_num)
 {

     int my_rank = formic::mpi::rank();

 _num_params=new_num;

 _mbuilder.resetParamNumber(new_num);


 int num_shift = _shift_scale.size();

   // size vector storing wavefunction update variable correctly
   int correct_vf_var_size = (new_num+1) * num_shift;
   if (_vf_var.size() != correct_vf_var_size)
   {  
      if(my_rank ==0)
      {  
         std::cout << "Should be resizing _vf_var vector for storing updates for different shifts, new size is: " << correct_vf_var_size << std::endl;
      } 
     _vf_var.resize(correct_vf_var_size);

   }
 }

//Function for storing blocked LM parameter inputs
template<typename S>
 void cqmc::engine::LMYEngine<S>::store_blocked_lm_info(int nblock,int nkeps)
 {
 _nblocks = nblock;
 _nkeps = nkeps;

 }


//Function for filtering parameters when samples are stored in the engine
 template<typename S>
 void cqmc::engine::LMYEngine<S>::selectParameters()
 {
 int my_rank = formic::mpi::rank();

 if(my_rank == 0)
 {
     std::cout << "Filtering parameters based on the noise of their gradients." << std::endl;
 }

 int num_samples = der_rat_history.rows();
 int num_all_params = der_rat_history.cols()-1;

 parameterSettings.resize(num_all_params);

 int new_param_num = 0;

 std::vector<double> grad_mean_list;
 std::vector<double> grad_sigma_list;
 std::vector<double> grad_ratio_list;

 std::vector<double> temp_e_list;
 std::vector<double> numerHistory;
 std::vector<double> denomHistory;
 std::vector<double> full_weight_history;
 std::vector<double> full_vgs_history;

 double numerVal=0.0;
 double denomVal=0.0;

 double numerSigma=0.0;
 double denomSigma=0.0;


//loop over all stored local energies and associated weights
 for(int j = 0; j < num_samples; j++)
 {
     double etmp = formic::real(le_der_rat_history.at(j,0));

     double vgs = vgs_history.at(j);
     double weight = weight_history.at(j);

     temp_e_list.push_back(etmp*vgs);

   //Full histories should be 1D with values from each thread's vector in the original histories  
     full_vgs_history.push_back(vgs);
     full_weight_history.push_back(weight);

 //For excited state, compute corresponding contributations to numerator and denominator of target function
 if(!_ground)
 {
     double numer_val = (_hd_lm_shift - etmp)*vgs;
     double denom_val = (_hd_lm_shift * _hd_lm_shift - 2 * _hd_lm_shift * etmp + etmp * etmp)*vgs;

     numerHistory.push_back(numer_val);
     denomHistory.push_back(denom_val);
     
 }

 }

 std::vector<double> energy_results = computeSigma_helper(full_weight_history, temp_e_list, full_vgs_history);

 double mean_energy = energy_results[0];

 if(!_ground)
 {
     std::vector<double> numerResults = computeSigma_helper(full_weight_history,numerHistory,full_vgs_history);
     std::vector<double> denomResults = computeSigma_helper(full_weight_history,denomHistory,full_vgs_history);
     std::vector<double> targetResults = computeSigma_helper(full_weight_history,numerHistory,denomHistory);

     numerVal = numerResults[0];
     denomVal = denomResults[0];

     numerSigma = numerResults[1];
     denomSigma = denomResults[1];

 }

 //Now loop over parameter derivative ratios
 
 for(int i = 1; i< num_all_params+1; i++)
 {

    std::vector<double> param_numer_history;
    std::vector<double> param_denom_history;



     for(int j = 0; j < num_samples; j++)
     {
         double param_le_der = formic::real(le_der_rat_history.at(j,i));
         double param_der_rat = formic::real(der_rat_history.at(j,i));
         double vgs = vgs_history.at(j);
         double etmp = formic::real(le_der_rat_history.at(j,0));

         if(_ground)
         {
             param_numer_history.push_back(param_le_der*vgs);
             param_denom_history.push_back(param_der_rat*vgs);



           }
         else
         {
            double numer_der_term = 2*(_hd_lm_shift*param_der_rat - param_le_der);
            double denom_der_term = 2*((_hd_lm_shift*param_der_rat - param_le_der)* (_hd_lm_shift - etmp));

             param_numer_history.push_back(numer_der_term);
             param_denom_history.push_back(denom_der_term);

         }


     }
     //After loop over samples, compute average mean and sigma for the parameter derivative
            double mean_deriv;
            double deriv_sigma;
            double ratio;
         if(_ground)
         {

             std::vector<double> first_term_results = computeSigma_helper(full_weight_history, param_numer_history, full_vgs_history);
             std::vector<double> second_term_results = computeSigma_helper(full_weight_history, param_denom_history, full_vgs_history);

             mean_deriv = first_term_results[0] - mean_energy*second_term_results[0];


             deriv_sigma = std::sqrt(first_term_results[1]*first_term_results[1] + second_term_results[1]*second_term_results[1]);

             grad_mean_list.push_back(mean_deriv);
             grad_sigma_list.push_back(deriv_sigma);

             ratio = mean_deriv/deriv_sigma;
             grad_ratio_list.push_back(std::abs(ratio));



         }
 //Excited State Case
         else
         {
             std::vector<double> results = computeSigma_helper(full_weight_history, param_numer_history, full_vgs_history);
             std::vector<double> secondResults = computeSigma_helper(full_weight_history, param_denom_history, full_vgs_history);


             double derivNumer = results[0]*denomVal - secondResults[0]*numerVal;
             double derivDenom = denomVal*denomVal;
             mean_deriv = derivNumer/derivDenom;

             double numerVar_term1 = results[0]*denomVal*results[0]*denomVal* ( (denomSigma/denomVal)*(denomSigma/denomVal) + (results[1]/results[0])*(results[1]/results[0]) );
             double numerVar_term2 = secondResults[0]*numerVal*secondResults[0]*numerVal* ( (numerSigma/numerVal)*(numerSigma/numerVal) + (secondResults[1]/secondResults[0])*(secondResults[1]/secondResults[0]) );
             double numerVar_tot = numerVar_term1 + numerVar_term2;


             double total_var = mean_deriv*mean_deriv * ( numerVar_tot/(derivNumer*derivNumer) + (2*denomSigma*denomSigma)/(derivDenom*derivDenom));

              deriv_sigma = std::sqrt(total_var);

             grad_mean_list.push_back(mean_deriv);
             grad_sigma_list.push_back(deriv_sigma);


             ratio = mean_deriv/deriv_sigma;
             grad_ratio_list.push_back(std::abs(ratio));

         }

     if(filter_info_ && my_rank == 0)
     {
         std::cout << "Mean for Parameter #" << i-1 << " : " << mean_deriv;
         std::cout << "   Sigma for Parameter #" << i-1 << " : " << deriv_sigma;
         std::cout << "   Ratio for Parameter #" << i-1 << " : " << std::abs(ratio) << std::endl;
     }

     if(std::abs(ratio) < ratio_threshold_)
     {
         parameterSettings[i-1] = false;

     }
     else
     {

         parameterSettings[i-1] = true;
         new_param_num++;
     }



 }

 if(my_rank == 0)
 {
     std::cout << "Number of Parameters left on: " << new_param_num << std::endl;
 }

 //In the unlikely event all parameters were turned off, turn some back on at random to avoid crashing the engine.
 if(new_param_num == 0)
 {
     if(my_rank==0)
     {
         std::cout << "All parameters were turned off by filtration. Turning some on at random." << std::endl;
     }

     //Try to turn on 20% of the total parameter number
     new_param_num = num_all_params/5;

     int i = 0;
     while(i< new_param_num)
     {
         int idx = std::rand() % num_all_params;
         parameterSettings[idx] = true;
         i++;
     }

 }

 if(_block_lm)
 {
     int tot_block_lm_param = _nblocks*_nkeps;

    //If the blocked LM is used, it is possible to have fewer parameters left than what the input number of 
    //blocks and kept directions require. 
     if(new_param_num < tot_block_lm_param)
     {
         if(my_rank == 0)
         {
             std::cout << "Remaining number of parameters less than that expected by blocked LM. Turning on more parameters." << std::endl;

         }

         std::vector<std::tuple <double,int,bool> > combined_ratio_param_list;
         for(int i = 0;i < num_all_params; i++)
         {
             std::tuple<double,int,bool> param_entry = std::make_tuple(grad_ratio_list[i],i,parameterSettings[i]);
             combined_ratio_param_list.push_back(param_entry);

         }

         //Identify the parameters with the largest ratios and turn on 
         //as many as the blocked LM's minimum requirement.
         std::sort(combined_ratio_param_list.begin(),combined_ratio_param_list.end());
         for(int i = num_all_params-1; i > num_all_params - tot_block_lm_param-1; i--)
         {
             int idx = std::get<1>(combined_ratio_param_list[i]);
             double ratio = std::get<0>(combined_ratio_param_list[i]);
             parameterSettings[idx] = true;
         }
         new_param_num = tot_block_lm_param; 
     }

 }


 resetParamNumber(new_param_num);



 }

//Helper function for computing mean derivatives and standard deviation for parameter filtration
 template<typename S>
 std::vector<double> cqmc::engine::LMYEngine<S>::computeSigma_helper(std::vector<double> weights, std::vector<double> numerSamples, std::vector<double> denomSamples)
 {
 double numSamples = weights.size();

 double y[7];
   y[0] = 0.0;                   // normalization constant
   y[1] = 0.0;                   // mean of numerator
   y[2] = 0.0;                   // mean of denominator
   y[3] = 0.0;                   // mean of the square of the numerator terms
   y[4] = 0.0;                   // mean of the square of the denominator terms
   y[5] = 0.0;                   // mean of the product of numerator times denominator
   y[6] = double(numSamples); // number of samples

   for (int i = 0; i < numSamples; i++)
   { 
     double n      = numerSamples[i];
     double d      = denomSamples[i];
     double weight = weights[i];

     y[0] += weight;
     y[1] += weight * n;
     y[2] += d;
     y[3] += weight * n * n;
     y[4] += weight * d * d;
     y[5] += weight * n * d;
   }


     double z[7];
   formic::mpi::allreduce(&y[0], &z[0], 7, MPI_SUM);

   double mf = z[1] / z[0]; // mean of numerator
   double mg = z[2] / z[0]; // mean of denominator
   double sf = z[3] / z[0]; // mean of the square of the numerator terms
   double sg = z[4] / z[0]; // mean of the square of the denominator terms
   double mp = z[5] / z[0]; // mean of the product of numerator times denominator
   double ns = z[6];        // number of samples


   double vf = (sf - mf * mf) * ns / (ns - static_cast<double>(1.0));
   double vg = (sg - mg * mg) * ns / (ns - static_cast<double>(1.0));
   double cv = (mp - mf * mg) * ns / (ns - static_cast<double>(1.0));

   double w_sum_   = z[0];
   double mean     = (mf / mg) / (static_cast<double>(1.0) + (vg / mg / mg - cv / mf / mg) / ns);
   double variance = (mf * mf / mg / mg) * (vf / mf / mf + vg / mg / mg - static_cast<double>(2.0) * cv / mf / mg);
   double stdErr   = std::sqrt(variance/z[0]); 

   std::vector<double> results;
   results.push_back(mean);
   results.push_back(stdErr);

 return results;

 }



#ifndef QMC_COMPLEX
template class cqmc::engine::LMYEngine<double>;
#else
template class cqmc::engine::LMYEngine<std::complex<double> >;
#endif
