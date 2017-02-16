//////////////////////////////////////////////////////////////////////////////////////////////////////////////
///  \file cqmc/engine/engine.h
///
/// 
///  \brief header file for harmonic davidson engine 
///
//////////////////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef CQMC_ENGINE_ENGINE_HEADER
#define CQMC_ENGINE_ENGINE_HEADER

#include<vector>
#include<sstream>

#include<formic/utils/matrix.h>
#include<formic/utils/lmyengine/block_mat.h>
#include<formic/utils/lmyengine/block_alg.h>
#include<formic/utils/lmyengine/block_detail.h>
#include<formic/utils/lmyengine/matrix_builder.h>

namespace formic {
  class VarDeps;
}

namespace cqmc {
  
  namespace engine {

  class LMYEngine {

private:
  
  /// \brief [in]flag to tell whether we are doing exact sampling
  bool _exact_sampling;

  /// \brief [in]flag to tell whether we are doing ground state simulation
  bool _ground;

  /// \brief [in]flag to tell whether we are doing variance correction
  bool _variance_correction;

  /// \brief [in]flag to tell whether to print out energy statistics
  bool _energy_print;

  /// \brief [in]flag to tell whether to print out linear method matrix
  bool _matrix_print;

  /// \brief [in]flag to tell whether to build linear method matrices
  bool _build_lm_matrix;

  /// \brief [in]flag to tell whether to use SPAM algorithm
  bool _spam;
  
  /// \brief [in]flag to tell whether to use variable dependency system
  bool _var_deps_use;

  /// \brief [in]flag to tell whether to chase the lowest eigenvalue
  bool _chase_lowest;

  /// \brief [in]flag to tell whether to chase the closest eigenvalue
  bool _chase_closest;

  /// \brief [in]flag to tell whether to do eom calculation
  bool _eom;

  /// \brief [in]flag to tell whether to compute S^2 value along with energy
  bool _ssquare;

  /// \brief [in]flag to tell whether to orthonormalize pairing matrix subspace first
  bool _pm_ortho_first;
  
  /// \brief [in]flag to tell whether jastrow factor is fixed
  bool _jas_fixed;

  /// \brief [in]flag to tell whether the sample data is initialized
  bool _sample_initialized;

  /// \brief [in]flag to tell whether to perform block algorithm or not
  bool _block_lm;

  /// \brief [in]flag to tell whether the block algorithm's first sample has finished
  bool _block_first_sample_finished;

  /// \brief [in]flag to tell whether we also want wfn update besides the energy
  bool _wfn_update;

  /// \brief [out]flag to tell whether to store derivative vectors 
  bool _store_der_vec;

  /// \brief [out]flag to well whethet engine has fully initialized
  bool _fully_initialized;

  /// \brief [in]total number of samples 
  int _num_samp;

  /// \brief sample count 
  int _samp_count;

  /// \brief number of optimizable parameters
  int _num_params;

  /// \brief [in]maximum number of iteration allowed in the davidson solver
  int _lm_krylov_iter;

  /// \brief [in]maximum number of iteration allowed in the spam solver's inner loop
  int _lm_spam_inner_iter;

  /// \brief [in]approximation degree in spam
  int _appro_degree;

  /// \brief [in]number of sites
  int _n_sites;

  /// \brief [in]number of pairing matrix elements
  int _n_pm;

  /// \brief [in]number of jastrow factor elements
  int _n_jas;

  /// \brief [in]number of blocks that we will divide total variables into
  int _nblocks;

  /// \brief [in]number of search directions per shift to keep in each block
  int _nkeps;

  /// \brief [in]harmonic davidson shift
  double _hd_lm_shift;

  /// \brief [in]weight of variance in variance correction
  double _var_weight;

  /// \brief [in]resudual threshold below which we say eigen solver has converged
  double _lm_eigen_thresh;

  /// \brief [in]minimum singular value of overlap matrix allowed
  double _lm_min_S_eval;

  /// \brief [in]initial energy/target function
  double _init_cost;

  /// \brief [in]initial variance
  double _init_var;

  /// \brief [in]maximum energy change allowed in eigensolver
  double _lm_max_e_change;

  /// \brief [in]identity shift
  double _lm_ham_shift_i;

  /// \brief [in]overlap shift
  double _lm_ham_shift_s;

  /// \brief [in]maximum wfn update value allowed
  double _lm_max_update_abs;

  /// \brief [out]energy
  double _energy;

  /// \brief [out]standard deviation of energy
  double _esdev;

  /// \brief [out]statistical error of energy
  double _eserr;

  /// \brief [out]target function value
  double _target;

  /// \brief [out]target function statistical error
  double _tserr;

  /// \brief [in]vector that stores history of local energy
  std::vector<std::vector<double> > _le_list;

  /// \brief [in]vector that stores history of |value/guiding|^2
  std::vector<std::vector<double> > _vg;

  /// \brief [in]vector that stores history of weight
  std::vector<std::vector<double> > _weight;

  /// \brief [in]vector that stores the shift that we will be solving
  std::vector<double>  _shift_scale;

  /// \brief [in]vector that stores energy and its index
  std::vector<std::pair<std::complex<double>, int> >  _energy_index;

  /// \brief [out]vector that stores wavefunction updates
  std::vector<double> _vf_var;

  /// \brief [out]vector that stores the solve results of each shift
  std::vector<bool> _good_solve;

  /// \brief [out]vector that stores the whether each shift has been solved
  std::vector<int> _solved_shifts;

  /// \brief [in]matrix that stores derivative ratio vectors
  formic::Matrix<double> _der_rat;

  /// \brief vector that stores derivative ratio vectors (only used in exact sampling)
  std::vector<double> _der_rat_vec;

  /// \brief [in]matrix that stores energy derivative vectors
  formic::Matrix<double> _le_der;

  /// \brief vector that stores local energy derivatives (only used in exact sampling)
  std::vector<double> _le_der_vec;

  /// \brief [in]matrix that stores local spin derivatives
  formic::Matrix<double> _spin_der;

  /// \brief vector that stores local spin derivatives (only used in exact sampling)
  std::vector<double> _spin_der_vec;

  /// \brief [out]matrix that stores eigenvectors from eom calculation
  formic::Matrix<std::complex<double> > _evecs;

  /// \brief [out]update directions from each block for each shift
  /// block_ups[i][j](p,q) refers to the pth element of the qth update for the jth shift for the ith block
  std::vector<std::vector<formic::Matrix<double> > > _block_ups;

  /// \brief [in]vector of vector storing old updates in independent variable form
  std::vector<formic::ColVec<double> > _ou_ind;

  /// \brief [in]pointer to formic's variable dependency system
  const formic::VarDeps * _dep_ptr;

  /// \brief object of block linear method class
  cqmc::engine::LMBlocker _lmb;

  /// \brief object of matrix builder class 
  cqmc::engine::HamOvlpBuilderHD _mbuilder;

  /// \brief [out]Hamiltonian matrix constructed by block linear method's final basis
  std::vector<formic::Matrix<double> > hh_block;

  /// \brief [out]overlap matrix constructed by block linear method's final basis
  std::vector<formic::Matrix<double> > ss_block;

  /// \brief [out]compact der rat
  std::vector<formic::ColVec<double> > drat_cmpct;

  /// \brief [out]compact eng der
  std::vector<formic::ColVec<double> > deng_cmpct;

  /// \brief [in]output stream
  std::ostream & output;

public:
  
  /////////////////////////////////////////////////////////////////////////////////
  /// \brief constructor with given parameters
  ///
  ///
  /////////////////////////////////////////////////////////////////////////////////
  LMYEngine(const formic::VarDeps * dep_ptr,
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
            std::ostream & output);

  ////////////////////////////////////////////////////////////////////////////////////////////////////////
  // \brief  Function that initialize all parameters in the engine 
  // 
  // 
  ////////////////////////////////////////////////////////////////////////////////////////////////////////
  void get_param(const formic::VarDeps * dep_ptr,
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
            const bool variance_correction=false,
            const bool energy_print=true,
            const bool matrix_print=false,
            const bool build_lm_matrix=true,
            const bool spam=false,
            const bool chase_lowest=true,
            const bool chase_closest=false,
            const bool pm_ortho_first=false,
            const bool jas_fixed=false,
            const int lm_krylov_iter=60,
            const int lm_spam_inner_iter=1,
            const int appro_degree=1,
            const int n_sites=0,
            const int n_pm=0,
            const int n_jas=0,
            const double var_weight=0.0,
            const double lm_eigen_thresh=1.0e-8, 
            const double lm_min_S_eval=0.99,
            const double init_cost=0.0,
            const double init_var=0.0);

  ////////////////////////////////////////////////////////////////////////////////////////////////////////
  /// \brief  Initialize the Engine
  ///
  /// \param[in]    old_updates   A matrix storing old update coefficients
  /// \param[in]    nblocks       Number of blocks that we will divide our total variables into
  /// \param[in]    nkeps         Number of search directions per shift to kee[ in each block
  /// \param[in]    nrand         Number of random vectors to add to the vector of old updates
  /// \param[in]    iterative     Whether to solve eigenvalue problem iteratively
  ///
  ////////////////////////////////////////////////////////////////////////////////////////////////////////
  void initialize(const int nblock = 1,
                  const int nrand = 0,
                  const int nkeps = 0,
                  const std::vector<formic::ColVec<double> > & old_updates = std::vector<formic::ColVec<double> >(),
                  const bool iterative = false);
  
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
  void take_sample(std::vector<double> & der_rat_samp,
                   std::vector<double> & le_der_samp,
                   std::vector<double> & ls_der_samp,
                   double vgs_samp,
                   double weight_samp);

  ////////////////////////////////////////////////////////////////////////////////////////////////////////
  /// \brief  Function that Take Sample Data from the Host Code
  /// 
  /// \param[in]  local_en       local energy
  /// \param[in]  vgs_samp       |<n|value_fn>/<n|guiding_fn>|^2
  /// \param[in]  weight_samp    weight for this sample
  ///
  ////////////////////////////////////////////////////////////////////////////////////////////////////////
  void take_sample(double local_en,
                   double vgs_samp,
                   double weight_samp);

  ////////////////////////////////////////////////////////////////////////////////////////////////////////
  /// \brief  Function that reduces all block matrices information from all processors to the root
  ///         processor
  ///
  ////////////////////////////////////////////////////////////////////////////////////////////////////////
  void sample_finish();

  ////////////////////////////////////////////////////////////////////////////////////////////////////////
  /// \brief  Function that updates variable dependency pointer 
  ///
  ////////////////////////////////////////////////////////////////////////////////////////////////////////
  void var_deps_ptr_update(const formic::VarDeps * new_dep_ptr);

  ////////////////////////////////////////////////////////////////////////////////////////////////////////
  /// \brief  Turn on wfn update
  ///
  ////////////////////////////////////////////////////////////////////////////////////////////////////////
  void turn_on_update() { _wfn_update = true; }

  ////////////////////////////////////////////////////////////////////////////////////////////////////////
  /// \brief  Turn off wfn update
  ///
  ////////////////////////////////////////////////////////////////////////////////////////////////////////
  void turn_off_update() { _wfn_update = false; }

  ////////////////////////////////////////////////////////////////////////////////////////////////////////
  /// \brief  Energy and Target Function Calculation
  ///
  ////////////////////////////////////////////////////////////////////////////////////////////////////////
  void energy_target_compute();

  ////////////////////////////////////////////////////////////////////////////////////////////////////////
  /// \brief  Do necessary preparations for the wavefunction update 
  ///
  ////////////////////////////////////////////////////////////////////////////////////////////////////////
  void wfn_update_prep();

  ////////////////////////////////////////////////////////////////////////////////////////////////////////
  /// \brief  Wave Function Update Calculation
  ///
  ////////////////////////////////////////////////////////////////////////////////////////////////////////
  void wfn_update_compute();

  ////////////////////////////////////////////////////////////////////////////////////////////////////////
  /// \brief  Equation-of-Motion Calculation
  ///
  ////////////////////////////////////////////////////////////////////////////////////////////////////////           
  void eom_compute();

  ////////////////////////////////////////////////////////////////////////////////////////////////////////
  /// \brief  Reset the Object
  ///
  ////////////////////////////////////////////////////////////////////////////////////////////////////////
  void reset();

  void shift_update(std::vector<double> & new_shift);

  /// \brief function that returns whether engine is fully initialized
  bool full_init() { return _fully_initialized; }

  /// \brief function that returns whether engine is going to do first or second block lm sampling
  bool block_first() { return !_block_first_sample_finished; }

  /// \brief function that returns the average of local energy
  double energy_mean() { return _energy; }

  /// \brief function that returns the statistical error of energy
  double energy_statistical_err() { return _eserr; } 

  /// \brief function that returns the standard deviation of energy
  double energy_sdev() { return _esdev; }

  /// \brief function that returns the variance of energy
  double energy_variance() { return _esdev * _esdev; } 

  /// \brief function that returns the target function value
  double target_value() { return (_ground ? _energy : _target); }

  /// \brief function that returns the target function statistical error
  double target_statistical_err() { return _tserr; }

  /// \brief function that returns the wavefunction update vector
  std::vector<double> wfn_update() { return _vf_var; }

  /// \brief function that returns the eigenvectors coming out of eom calculation
  formic::Matrix<std::complex<double> > & eom_evecs() { return _evecs; } 

  /// \brief function that returns the energy vector coming out of eom calculation
  std::vector<std::pair<std::complex<double>, int> > & eom_energy() { return _energy_index; }

  /// \brief function that returns the solve results
  std::vector<bool> & good_solve() { return _good_solve; }

  /// \brief function that returns the derivative ratio matrix
  formic::Matrix<double> & der_rat() { return _der_rat; }

  /// \brief function that returns the energy derivative matrix
  formic::Matrix<double> & le_der() { return _le_der; }

  /// \brief function that returns the weight list 
  const std::vector<double> & weight_list() { return _weight[0]; }

  /// \brief function that returns the |value/guiding|^2 list 
  const std::vector<double> & vgs_list() { return _vg[0]; }

  ////////////////////////////////////////////////////////////////////////////////////////////////////////
  ///
  /// \brief  harmonic davidson energy calculation function
  ///
  /// \param[in]   exacting_sampling     whether the sampling is exact sampling
  /// \param[in]   ground_state          whether to do ground state calculation
  /// \param[in]   variance_correction   whether to correct the finite variance problem
  /// \param[in]   print                 whether to print out energy statistics
  /// \param[in]   hd_lm_shift           harmonic davidson shift used to target excited states(if you are  
  ///                                    calculating ground state, then this is ignored
  /// \param[in]   var_weight            weight of variance correction term
  /// \param[in]   le_list               local energy list
  /// \param[in]   vg                    |value/guiding|^2 list
  /// \param[in]   weight                weight list
  /// \param[out]  energy                average of local energy
  /// \param[out]  esdev                 energy standard deviation
  /// \param[out]  eserr                 energy statistical error
  /// \param[out]  target                target function value(0 if doing ground state)
  /// \param[out]  tserr                 target function value's statistical error(0 if doing ground state)
  /// \param[out]  output                output stream
  ///
  ////////////////////////////////////////////////////////////////////////////////////////////////////////
  void call_engine(const bool exact_sampling,
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
                   std::ostream & output);

  ////////////////////////////////////////////////////////////////////////////////////////////////////////
  ///
  /// \brief  harmonic davidson wavefunction update function
  ///
  /// \param[in]   dep_ptr             pointer to formic's variable dependence object(see README.txt for an example)
  /// \param[in]   ground_state        whether to do ground state calculation
  /// \param[in]   build_lm_matrix     whether to build linear method matrix explicitly or just multiply by it
  /// \param[in]   variance_correct    whether to correct the early shift problem caused by finite variance
  /// \param[in]   print_matrix        whether to print hamiltonian and overlap matrix after built
  /// \param[in]   lm_krylov_iter      maximum number of iterations allowed in eigenvalue solver
  /// \param[in]   lm_spam_inner_iter  maximum number of iterations allowed in spam inner loop
  /// \param[in]   lm_eigen_thresh     resudual threshold below which we say eigen solver has converged
  /// \param[in]   lm_min_S_eval       minimum singular value of overlap matrix allowed
  /// \param[in]   spam_use            whether to use spam algorithm or not
  /// \param[in]   var_dep_use         whether to use variable dependence system or not
  /// \param[in]   chase_lowest        whether to chase the lowest eigenvalue or not
  /// \param[in]   chase_closest       whether to chase the closest eigenvalue or not
  /// \param[in]   init_energy         initial energy
  /// \param[in]   init_var            initial variance
  /// \param[in]   lm_max_e_change     maximum energy change allowed in eigen solver
  /// \param[in]   lm_ham_shift_i      vector of identity shift applied to hamiltonian matrix
  /// \param[in]   lm_ham_shift_s      vector of overlap shift applied to hamiltonian matrix
  /// \param[in]   omega               harmonic davidson shift
  /// \param[in]   var_weight          weight of variance correction term
  /// \param[in]   lm_max_update_abs   the largest allowable absolute value for a derivative function coefficient in the linear method
  /// \param[in]   der_rat             bare derivative ratio vector
  /// \param[in]   le_der              energy derivative ratio vector
  /// \param[in]   vg                  |value/guiding|^2 list
  /// \param[in]   weight              weight list
  /// \param[in]   shift_scale         scaling of penalty shift
  /// \param[out]  vf_var              wave function update vector
  /// \param[out]  good_solve          vector to store the result of each solve of each shift 
  /// \param[out]  solved_shift        vector to tell whether this shift has already been solved
  /// \param[out]  output              output stream
  ///
  ////////////////////////////////////////////////////////////////////////////////////////////////////////
  void call_engine(const formic::VarDeps * dep_ptr,
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
                   const double init_energy,
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
                   std::vector<int> & solved_shift,
                   const std::vector<double> & shift_scale,
                   std::ostream & output);

  /////////////////////////////////////////////////////////////////////////////////////////////////////////
  ///
  /// \brief  EOM-JAGP function 
  ///   
  /// \param[in]  print_matrix        whether to print hamiltonian and overlap matrix after built 
  /// \param[in]  eom_inv_thresh      threshold to use in eom's matrix inversion 
  /// \param[in]  ground_state        whether this calculation is ground state calculation
  /// \param[in]  ssquare             whether to compute S^2 value along with energy
  /// \param[in]  pm_ortho_first      whether to orthonormalize pairing matrix subspace first 
  /// \param[in]  jas_fixed           whether jastrow factor is fixed in wavefunction
  /// \param[in]  n_sites             number of sites(orbitals)
  /// \param[in]  n_pm                number of pairing matrix elements 
  /// \param[in]  n_jas               number of jastrow factor matrix elements(both same and opposite spin)
  /// \param[in]  der_rat             bare derivative ratio vector 
  /// \param[in]  le_der              energy derivative ratio vector 
  /// \param[in]  spin_der            spin derivative ratio vector 
  /// \param[in]  vg                  |value/guiding|^2 list
  /// \param[in]  weight              weight list
  /// \param[in]  output              output stream
  ///
  /////////////////////////////////////////////////////////////////////////////////////////////////////////
  void call_engine(const bool print_matrix, 
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
                   std::ostream & output);


  /////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /// \brief  Computes a wave function update using block recursive linear method algorithm, 
  ///         which solves each block many times, each time including some information from a second
  ///         block in the form of previous update components. 
  ///         These variale solves give a set of desirable vectors within the block,
  ///         An SVD is then used to pick a small number of vectors that are best able to recreate
  ///         all of these desirable vectors
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
  void get_brlm_update_alg_part_one(const formic::VarDeps * deps,
                                    const int nkps,
                                    const double shift_i,
                                    const double shift_s,
                                    const std::vector<double> & shift_scale,
                                    const double max_update_abs,
                                    std::vector<bool> & good_solve,
                                    std::vector<int> & shift_solved, 
                                    std::ostream & output);

  /////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /// \brief  Computes a wave function update using block recursive linear method algorithm, 
  ///         this function solves the eigenvalue problem within the final basis 
  ///
  /// \param[in]      deps               object keeping track of variable dependencies
  /// \param[in]      iterative          whether to solve the eigenvalue problem iteratively
  /// \param[in]      nblock             number of blocks to divide the variables into
  /// \param[in]      nkps               number of search directions per shift to keep in each block
  /// \param[in]      nrand              number of random vectors to add to the vector of old updates
  /// \param[in]      shift_i            strength of the identity shift
  /// \param[in]      shift_s            strength of the overlap shift
  /// \param[in]      shift_scale        scale factors telling which shifts to use
  /// \param[in]      max_update_abs     reject
  ///
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////
  void get_brlm_update_alg_part_two(const formic::VarDeps * deps,
                                    const bool iterative,
                                    const int nkps,
                                    const double shift_i,
                                    const double shift_s,
                                    const std::vector<double> & shift_scale,
                                    const double max_update_abs,
                                    std::vector<double> & updates,
                                    std::vector<bool> & good_solve,
                                    std::vector<int> & shift_solved, 
                                    std::ostream & output);

};

  }

}

#endif
