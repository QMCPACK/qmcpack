////////////////////////////////////////////////////////////////////////////////////////////////////
// \file eigen_solver.h       an abstract class to solver generalized eigenvalue problem iteratively 
//
// 
//
//
//
////////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef ENGINE_EIGEN_SOLVER_HEADER
#define ENGINE_EIGEN_SOLVER_HEADER

#include<vector>
#include<string>
#include<numeric>
#include<cassert>
#include<algorithm>
#include<cmath>
#include<complex>
#include<iostream>

#include<boost/format.hpp>
#include<boost/shared_ptr.hpp>

//#include<mpi.h>

#include<formic/utils/lmyengine/var_dependencies.h>
#include<formic/utils/matrix.h>

namespace cqmc {
  
  namespace engine { 

  class EigenSolver {

protected:

  /// \brief pointer to the variable dependency object 
  const formic::VarDeps * _dep_ptr;

  /// \brief number of wave function derivatives
  int _nder;

  /// \brief dimension of first derivative subspace 
  int _nfds;

  /// \brief current estimate of energy/target function  
  double _cost;

  /// \brief initial wavefunction energy/target function 
  double _init_cost;

  /// \brief initial wavefunction's variance
  double _init_variance;

  /// \brief maximum alowed energy change
  double _max_energy_change;

  /// \brief linear method identity shift applied to the Hamiltonian diagonal 
  double _hshift_i;

  /// \brief linear method overlap shift applied to the Hamiltonian diagonal
  double _hshift_s;

  /// \brief harmonic Davidson shift applied to Hamiltonian and overlap matrix(only apply when doing excited state calculaions)
  double _hd_shift;

  /// \brief weight of variance correction term
  double _var_weight;

  /// \brief residual threshold
  double _residual_threshold;

  /// \brief total weight
  double _total_weight;

  /// \brief average of |value/guiding|^2 function 
  double _vgsa;

  /// \brief value used for variance correct calculation(tau)
  double _tau;

  /// \brief flag to tell whether to use variable dependencies
  bool _var_deps_use;

  /// \brief flag to tell whether the lowest eigenvalue has a complex component 
  bool _eval_was_complex;

  /// \brief flag to tell whether to chase the lowest eigenvalue or not 
  bool _chase_lowest;

  /// \brief flag to tell whether to chase the lowest eigenvalue or not 
  bool _chase_closest;
  
  /// \brief flag to tell whether to do ground state calculation 
  bool _ground;

  /// \brief flag to tell whether to correct the finite variance issue
  bool _variance_correct;

  /// \brief calculated eigenvector on each process 
  formic::ColVec<double> _evecs;

  /// \brief calculated eigenvector for independent variables
  formic::ColVec<double> _ind_evecs;

  /// \brief converted wavefunction coefficients
  formic::ColVec<double> _wf_coefficients;

  /// \brief derivative ratio vectors(merged with |value/guiding|^2 and weights)
  formic::Matrix<double> & _der_rat;

  /// \brief local energy derivative vectors(merged with |value/guiding|^2 and weights) 
  formic::Matrix<double> & _le_der;
 
  
  public:

  /////////////////////////////////////////////////////////////////////////////
  // \brief solve the subspace generalized eigenvalue prblem with
  //        nonsymmetric H and symmetric S
  //
  // \param[in]  outer whether to solve for outer loop or inner loop
  //
  /////////////////////////////////////////////////////////////////////////////
  
  virtual void solve_subspace_nonsymmetric(const bool outer = true) = 0;

  ////////////////////////////////////////////////////////////////////////////
  // \brief solve the subspace eigenvalue problem
  //
  //
  //
  //
  ////////////////////////////////////////////////////////////////////////////

  virtual void solve_subspace(const bool outer = true) = 0;

  ////////////////////////////////////////////////////////////////////////////
  // \brief  constructor 
  //
  //
  //
  //
  ////////////////////////////////////////////////////////////////////////////
  EigenSolver(const formic::VarDeps * dep_ptr,
              const int nfds,
              const double lm_eigen_thresh,
              const bool var_deps_use,
              const bool chase_lowest,
              const bool chase_closest,
              const bool ground,
              const bool variance_correct,
              const std::vector<double> & vf,
              const double init_cost,
              const double init_variance,
              const double hd_shift,
              const double var_weight,
              const double lm_max_e_change,
              const double total_weight,
              const double vgsa,
              formic::Matrix<double> & der_rat,
              formic::Matrix<double> & le_der)
  :_dep_ptr(dep_ptr),
  _nder(nfds-1),
  _nfds(nfds),
  _residual_threshold(lm_eigen_thresh),
  _cost(init_cost),
  _init_variance(init_variance),
  _var_deps_use(var_deps_use),
  _chase_lowest(chase_lowest),
  _chase_closest(chase_closest),
  _ground(ground),
  _variance_correct(variance_correct),
  _init_cost(init_cost),
  _max_energy_change(lm_max_e_change),
  _total_weight(total_weight),
  _vgsa(vgsa),
  _tau(0.0),
  _hshift_i(0.0),
  _hshift_s(0.0),
  _hd_shift(hd_shift),
  _var_weight(var_weight),
  _der_rat(der_rat),
  _le_der(le_der)
  {
    // if we use variable dependence system, change _nder and _nfds to number of 
    // independent variables 
    if ( _var_deps_use ) {
      _nder = dep_ptr -> n_ind();
      _nfds = _nder + 1;
    }
  }

  ////////////////////////////////////////////////////////////////////////////////////
  // \brief evaluate tau and compute the "variance corrected" Hamiltonian
  // 
  // 
  ////////////////////////////////////////////////////////////////////////////////////

  virtual void tau_and_correct_ham() = 0;

  ////////////////////////////////////////////////////////////////////////////////////
  // \brief adds a new krylov vector 
  //
  // \param[in]  v  vector that needs to be added
  //
  ////////////////////////////////////////////////////////////////////////////////////

  virtual void add_krylov_vector(const formic::ColVec<double> & v) = 0;

  ////////////////////////////////////////////////////////////////////////////////
  // \brief adds a new Krylov basis vector for spam inner loop
  //
  //
  //
  ////////////////////////////////////////////////////////////////////////////////

  virtual void add_krylov_vector_inner(const formic::ColVec<double> & v) = 0;

  ////////////////////////////////////////////////////////////////////////////////
  // \brief adds a bunch of new Krylov vectors(which is a matrix) for spam outer 
  //        loop
  //
  //
  ////////////////////////////////////////////////////////////////////////////////

  virtual void add_krylov_vectors_outer(const formic::Matrix<double> & m) = 0;

  ////////////////////////////////////////////////////////////////////////////////////
  // \brief function that perfoms hamiltonian matrix-vector multiplication 
  // 
  // \param[in]   x              input vector 
  // \param[in]   matrix_built   whether we have already built the matrix or not
  // \param[in]   transpose      whether to use transposed hamiltonian or not 
  // \param[in]   approximate    whether to use approximated hamiltonian or not
  // \param[out]  y              result vector
  //
  ////////////////////////////////////////////////////////////////////////////////////

  virtual void HMatVecOp(const formic::ColVec<double> & x, formic::ColVec<double> & y, const bool transpose = false, const bool approximate = false) = 0;  

  ////////////////////////////////////////////////////////////////////////////////////
  // \brief function that performs hamiltonian matrix-matrix multiplication 
  //
  // \param[in]   x              input matrix 
  // \param[in]   matrix_built   whether we have already built the matrix or not
  // \param[in]   transpose      whether to use the transposed hamiltonian or not
  // \param[in]   approximate    whether to use approximated hamiltonian or not
  // \param[out]  y              result matrix
  //
  ////////////////////////////////////////////////////////////////////////////////////

  virtual void HMatMatOp(const formic::Matrix<double> & x, formic::Matrix<double> & y, const bool transpose, const bool approximate = false) = 0;

  ////////////////////////////////////////////////////////////////////////////////////
  // \brief function that performs overlap matrix-vector multiplication 
  //
  // \param[in]   x              input vector 
  // \param[in]   matrix_built   whether we have already built the matrix or not 
  // \param[in]   approximate    whether to use approximated overlap or not 
  // \param[out]  y              result vector 
  // NOTE: Unlike hamiltonian matrix-vector multiplication function, no transpose flag
  //       in this function because overlap matrix is assumed to be symmetric
  //
  ////////////////////////////////////////////////////////////////////////////////////

  virtual void SMatVecOp(const formic::ColVec<double> & x, formic::ColVec<double> & y, const bool approximate = false) = 0;

  ////////////////////////////////////////////////////////////////////////////////////
  // \brief function that performs overlap matrix-matrix multiplication 
  //
  // \param[in]   x              input matrix 
  // \param[in]   matrix_built   whether we have already built the matrix or not 
  // \param[in]   approximate    whether to use approximated overlap or not 
  // \param[out]  y              result matrix 
  // NOTE: Unlike hamiltonian matrix-matrix multiplication function, no transpose flag
  //       in this function because overlap matrix is assumed to be symmetric
  //
  ////////////////////////////////////////////////////////////////////////////////////

  virtual void SMatMatOp(const formic::Matrix<double> & x, formic::Matrix<double> & y, const bool approximate = false) = 0;

  ////////////////////////////////////////////////////////////////////////////////////
  // \brief change linear method identity and overlap shift 
  //
  //
  //
  ////////////////////////////////////////////////////////////////////////////////////
  
  void update_lm_shift(const double new_i_shift, const double new_s_shift)
  {
    // change shift
    _hshift_i = new_i_shift;
    _hshift_s = new_s_shift;

  }

  ////////////////////////////////////////////////////////////////////////////////////
  // \brief updates hamiltonian * krylov vector and hamiltonian projection based on 
  //        new shift 
  //
  //
  ////////////////////////////////////////////////////////////////////////////////////

  virtual void update_hvecs_sub(const double new_i_shift, const double new_s_shift) = 0;

  ////////////////////////////////////////////////////////////////////////////////////
  // \brief reset the eigen solver 
  //
  //
  // \brief clear wfn coefficient vector and call reset function of child class
  //
  ////////////////////////////////////////////////////////////////////////////////////

  void reset()
  {
   
    // clear eigenvectors and wavefunction coefficients
    _evecs.reset(0);
    _ind_evecs.reset(0);
    _wf_coefficients.reset(0);

    // set current energy as initial wave function's energy
    _cost = _init_cost;

    // call child class reset function 
    this -> child_reset();
  }

  ////////////////////////////////////////////////////////////////////////////////////
  // \brief child class's reset function 
  // 
  //
  //
  //
  //
  ////////////////////////////////////////////////////////////////////////////////////

  virtual void child_reset() = 0;

  ////////////////////////////////////////////////////////////////////////////////////
  // \brief iteratively solve the eigenvalue problem 
  //
  //
  //
  //
  ////////////////////////////////////////////////////////////////////////////////////

  virtual bool iterative_solve(double & eval, std::ostream & output) = 0;

  ////////////////////////////////////////////////////////////////////////////////////
  // \brief solve the (generalized) eigenvalue problem, which for now just call 
  //        the iterative_solve function when implemented in child class
  //
  //
  //
  ////////////////////////////////////////////////////////////////////////////////////

  virtual bool solve(double & eval, std::ostream & output) = 0;

  ///////////////////////////////////////////////////////////////////////////////////////////////
  // \brief converts eigenvectors into wave function coefficients
  //        solving this question Hc = Sd for d, S is the ground overlap matrix 
  //
  //
  ///////////////////////////////////////////////////////////////////////////////////////////////

  void convert_to_wf_coeff()
  {
					
    // calculate wave function parameters by scaling the eigenvector and make the first element 1
    const double scaling = _evecs.at(0);
    _wf_coefficients.reset(_nfds);
    _wf_coefficients = _evecs;
    _wf_coefficients /= scaling;

    // children's class's convert
    this -> child_convert_to_wf_coeff();
  }

  //////////////////////////////////////////////////////////////////////////////////////////////
  // \brief child's class's wfn conversion function 
  //
  //
  //
  //
  //////////////////////////////////////////////////////////////////////////////////////////////

  virtual void child_convert_to_wf_coeff() = 0;

  //////////////////////////////////////////////////////////////////////////////////////////////
  // \brief return the wave function coefficients
  //
  //
  //
  //////////////////////////////////////////////////////////////////////////////////////////////

  formic::ColVec<double> wf_coeff()
  {
    return _wf_coefficients;
  }

  //////////////////////////////////////////////////////////////////////////////////////////////
  // \brief build linear matrix(used when doing debug)
  //
  //
  //
  //////////////////////////////////////////////////////////////////////////////////////////////
  
  //virtual const formic::Matrix<double> debug_build() = 0;


};

} // end of namespace engine

} // end of namespace cqmc

#endif
