//////////////////////////////////////////////////////////////////////////////////////////////
/// \file cqmc/engine/engine.cpp
///
///
/// \brief implementation of harmonic davidson engine function
///
///
/////////////////////////////////////////////////////////////////////////////////////////////

#include<vector>
#include<string>
#include<algorithm>
#include<utility>

#include<boost/shared_ptr.hpp>

#include<formic/utils/openmp.h>

#include <formic/utils/lmyengine/engine.h>
#include <formic/utils/lmyengine/updater.h>
#include <formic/utils/lmyengine/energy_target.h>
#include <formic/utils/lmyengine/matrix_builder.h>
#include <formic/utils/lmyengine/eom.h>
#include <formic/utils/lmyengine/var_dependencies.h>
#include <formic/utils/lmyengine/engine_timing.h>
#include <formic/utils/exception.h>
#include <formic/utils/lapack_interface.h>
#include <formic/utils/mpi_interface.h>

/////////////////////////////////////////////////////////////////////////////////
/// \brief constructor with given parameters
///
///
/////////////////////////////////////////////////////////////////////////////////
cqmc::engine::LMYEngine::LMYEngine(const formic::VarDeps * dep_ptr,
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
                                   std::ostream & output)
:_mbuilder(_der_rat, _le_der, _spin_der, _vg[0], _weight[0], hd_lm_shift, num_params, appro_degree, spam, ground, variance_correction, build_lm_matrix, eom, matrix_print),
_dep_ptr(dep_ptr),
_exact_sampling(exact_sampling),
_ground(ground),
_variance_correction(variance_correction),
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
output(output)
{
  // get the number of threads being used 
  int NumThreads = omp_get_max_threads();

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
void cqmc::engine::LMYEngine::get_param(const formic::VarDeps * dep_ptr,
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
            const bool variance_correction,
            const bool energy_print,
            const bool matrix_print,
            const bool build_lm_matrix,
            const bool spam,
            const bool chase_lowest,
            const bool chase_closest,
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

  // solve results
  _good_solve.resize(_shift_scale.size());
  _solved_shifts.resize(_shift_scale.size());
  std::fill(_good_solve.begin(), _good_solve.end(), false);
  std::fill(_solved_shifts.begin(), _solved_shifts.end(), false);

  _mbuilder.get_param(hd_lm_shift, num_params, appro_degree, spam, ground, variance_correction, build_lm_matrix, eom, matrix_print);

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
void cqmc::engine::LMYEngine::initialize(const int nblock,
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
void cqmc::engine::LMYEngine::take_sample(std::vector<double> & der_rat_samp,
                                          std::vector<double> & le_der_samp,
                                          std::vector<double> & ls_der_samp,
                                          double vgs_samp,
                                          double weight_samp) {

    
  // get the number of threads being used
  int NumThreads = omp_get_num_threads();

  //output << boost::format("entering take_sample function") << std::endl;

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

      bool first_samp = !_block_first_sample_finished;

      // if it is the first sampling
      if ( first_samp ) {
        // combine vgs and weight to get effective weight
        const double d = vgs_samp * weight_samp;

        // accumulate data for the block algorithm object
        // for ground state calculation
        if ( _ground ) 
          _lmb.acc(d, der_rat_samp, le_der_samp, true);
        // for excited state calculation
        else {
          // get <n|(w-H)|Psi_i>/<n|Psi>
          std::vector<double> mle_der_samp(le_der_samp.size(), 0.0);
          for (int i = 0; i < mle_der_samp.size(); i++) 
            mle_der_samp.at(i) = _hd_lm_shift * der_rat_samp.at(i) - le_der_samp.at(i);

          // accumulate data
          //_lmb.acc(d, mle_der_samp, der_rat_samp);
          _lmb.acc(d, der_rat_samp, mle_der_samp, false);
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
          drat_cmpct[myThread].at(0) = der_rat_samp.at(0);
          deng_cmpct[myThread].at(0) = le_der_samp.at(0);
        }
        else {
          drat_cmpct[myThread].at(0) = der_rat_samp.at(0);
          deng_cmpct[myThread].at(0) = _hd_lm_shift * der_rat_samp.at(0) - le_der_samp.at(0);
        }

        // get compact der rat and der eng for old updates
        for (int k = 0; k < _lmb.ou_mat().cols(); k++) {
          drat_cmpct[myThread].at(1+k) = 0.0;
          deng_cmpct[myThread].at(1+k) = 0.0;
          for (int i = 0; i < _dep_ptr->n_ind(); i++) {
            if ( _ground ) {
              drat_cmpct[myThread].at(1+k) += _lmb.ou_mat().at(i, k) * der_rat_samp.at(1+i);
              deng_cmpct[myThread].at(1+k) += _lmb.ou_mat().at(i, k) * le_der_samp.at(1+i);
            }
            else {
              drat_cmpct[myThread].at(1+k) += _lmb.ou_mat().at(i, k) * der_rat_samp.at(1+i);
              deng_cmpct[myThread].at(1+k) += _lmb.ou_mat().at(i, k) * (_hd_lm_shift * der_rat_samp.at(1+i) - le_der_samp.at(1+i));
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
                  drat_cmpct[myThread].at(nsv+q) += _block_ups.at(b).at(s).at(i, n) * der_rat_samp.at(1+(_lmb.bb(b)+i));
                else 
                  drat_cmpct[myThread].at(nsv+q) += _block_ups.at(b).at(s).at(i, n) * der_rat_samp.at(1+(_lmb.bb(b)+i));
              }
              deng_cmpct[myThread].at(nsv+q) = 0.0;
              for (int i = 0; i < _lmb.bl(b); i++) {
                if ( _ground ) 
                  deng_cmpct[myThread].at(nsv+q) += _block_ups.at(b).at(s).at(i, n) * le_der_samp.at(1+(_lmb.bb(b)+i));
                else 
                  deng_cmpct[myThread].at(nsv+q) += _block_ups.at(b).at(s).at(i, n) * (_hd_lm_shift * der_rat_samp.at(1+(_lmb.bb(b)+i)) - le_der_samp.at(1+(_lmb.bb(b)+i)));
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
void cqmc::engine::LMYEngine::take_sample(double local_en,
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
/// \brief  Function that reduces all block matrices information from all processors to the root
///         processor
///
////////////////////////////////////////////////////////////////////////////////////////////////////////
void cqmc::engine::LMYEngine::sample_finish() {
  
  // get rank number and number of ranks
  int my_rank = formic::mpi::rank();
  int num_rank = formic::mpi::size();

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
      return;

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
void cqmc::engine::LMYEngine::var_deps_ptr_update(const formic::VarDeps * new_dep_ptr) {
  
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
void cqmc::engine::LMYEngine::energy_target_compute() { 

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
void cqmc::engine::LMYEngine::wfn_update_prep() {

   // get rank number and number of ranks
  int my_rank = formic::mpi::rank();
  int num_rank = formic::mpi::size();
  //MPI_Comm_rank(MPI_COMM_WORLD, & my_rank);
  //MPI_Comm_size(MPI_COMM_WORLD, & num_rank);
  //std::cout << "entering wfn_update_prep on rank " << my_rank << std::endl;
 
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
void cqmc::engine::LMYEngine::wfn_update_compute() {

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
  this -> get_brlm_update_alg_part_two(_dep_ptr,
                                       _lmb.iterative(),
                                       _nkeps,
                                       _lm_ham_shift_i,
                                       _lm_ham_shift_s,
                                       _shift_scale,
                                       _lm_max_update_abs,
                                       _vf_var,
                                       _good_solve,
                                       _solved_shifts,
                                       output);
  }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief  Equation-of-Motion Calculation
///
////////////////////////////////////////////////////////////////////////////////////////////////////////           
void cqmc::engine::LMYEngine::eom_compute() {

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
void cqmc::engine::LMYEngine::reset() {

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

void cqmc::engine::LMYEngine::shift_update(std::vector<double> & new_shift) {

  // update shift
  for (int i = 0; i < new_shift.size(); i++) 
    _shift_scale.at(i) = new_shift.at(i);

}

/////////////////////////////////////////////////////////////////////////////////////////////
/// \brief harmonic davidson engine call to calculate energy and target function 
///
///
///
///
///
/////////////////////////////////////////////////////////////////////////////////////////////
void cqmc::engine::LMYEngine::call_engine(const bool exact_sampling,
                                          const bool ground_state,
                                          const bool variance_correct,
                                          const bool print,
                                          const double hd_lm_shift,
                                          const double var_weight,
                                          std::vector<double> & le_list,
                                          std::vector<double> & vg,
                                          std::vector<double> & weight,
                                          double & energy,
                                          double & esdev, 
                                          double & eserr,
                                          double & target, 
                                          double & tserr,
                                          std::ostream & output) {

  // get rank number
  int my_rank = formic::mpi::rank();
  //MPI_Comm_rank(MPI_COMM_WORLD, & my_rank);

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
void cqmc::engine::LMYEngine::call_engine(const formic::VarDeps * dep_ptr, 
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
                                          formic::Matrix<double> & der_rat,
                                          formic::Matrix<double> & le_der,
                                          std::vector<double> & vg, 
                                          std::vector<double> & weight,
                                          std::vector<double> & vf_var,
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
  boost::shared_ptr< cqmc::engine::HDLinMethodUpdater > hdupdater ( new cqmc::engine::HDLinMethodUpdater(der_rat, 
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

  else if ( !build_lm_matrix && spam_use ) 
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
                                    vf_var,
                                    good_solve,
                                    shift_solved,
                                    output);



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
void cqmc::engine::LMYEngine::call_engine(const bool print_matrix, 
                                          const bool ground_state,
                                          const bool eom,
                                          const bool ssquare,
                                          const bool pm_ortho_first,
                                          const bool jas_fixed,
                                          const double eom_inv_thresh,
                                          int n_sites,
                                          int n_pm,
                                          int n_jas,
                                          formic::Matrix<double> & der_rat,
                                          formic::Matrix<double> & le_der,
                                          formic::Matrix<double> & spin_der,
                                          formic::Matrix<std::complex<double> > & evecs,
                                          std::vector<std::pair<std::complex<double>, int> > & energy_index,
                                          std::vector<double> & vg,
                                          std::vector<double> & weight,
                                          std::ostream & output)
{
  
  
  // get rank number and number of ranks 
  int my_rank = formic::mpi::rank();
  int num_rank = formic::mpi::size();
  //MPI_Comm_rank(MPI_COMM_WORLD, & my_rank);
  //MPI_Comm_size(MPI_COMM_WORLD, & num_rank);

  // if this calculation is not ground state calculation, EOM will not be performed
  if ( !ground_state )
    return;

  // if eom flag is set to be false, EOM will not be performed 
  if ( !eom ) 
    return;

  // creates matrix builder 
  boost::shared_ptr< cqmc::engine::HamOvlpBuilderHD > mbuilder( new cqmc::engine::HamOvlpBuilderHD(der_rat, 
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
  boost::shared_ptr< cqmc::engine::EOM > eom_calculator( new cqmc::engine::EOM(mbuilder -> ham(), mbuilder -> ovl(), mbuilder -> ssquare(), ssquare, print_matrix, pm_ortho_first, jas_fixed, n_sites, n_pm, n_jas, eom_inv_thresh));

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
void cqmc::engine::LMYEngine::get_brlm_update_alg_part_one(const formic::VarDeps * deps,
                                                           const int nkps,
                                                           const double shift_i,
                                                           const double shift_s,
                                                           const std::vector<double> & shift_scale,
                                                           const double max_update_abs,
                                                           std::vector<bool> & good_solve,
                                                           std::vector<int> & shift_solved, 
                                                           std::ostream & output) {


  // get rank number of number of ranks
  int my_rank = formic::mpi::rank();
  int num_rank = formic::mpi::size();

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
void cqmc::engine::LMYEngine::get_brlm_update_alg_part_two(const formic::VarDeps * deps,
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

  // get rank number and number of ranks
  int my_rank = formic::mpi::rank();
  int num_rank = formic::mpi::size();
  //MPI_Comm_rank(MPI_COMM_WORLD, & my_rank);
  //MPI_Comm_size(MPI_COMM_WORLD, & num_rank);

  // size update vector correctly
  //std::cout << shift_scale.size() << "  " << 1 + _dep_ptr->n_tot() << std::endl;
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
                                                           this->target_value(),
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
                                                                    _energy,
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
