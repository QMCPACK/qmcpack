//////////////////////////////////////////////////////////////////////////////////////////////////////////
// \file spam_solver.h       a class to solve generalized eigenvalue problem with SPAM method
//
//  SPAM is short for Subspace Projected Approximation Matrix
//  Cite this paper:
//  The Subspace Projected Approximation Matrix(SPAM) Modification of the Davidson Method 
//  R.Shepard, A.Wagner, J.Tilson and Michael Minkoff. Journal of Computational Physics 172, 472-514(2001)
//
//  A short review of SPAM method:
//  H:    full Hamiltonian H_ij = <psi_i|H|psi_j>
//  S:    full overlap     S_ij = <psi_i|psi_j>
//  H(1): approximate hamiltonian
//  S(1): approximate overlap
//  X:    krylov space
//  P:    projector to krylov space, P = X * X^T
//  Q:    projector to the rest of space: Q = 1 - P   
//  H_bar: hybrid hamiltonian H_bar = PHP + PHQ + QHP + QH(1)Q
//  S_bar: hybrid overlap S_bar = PSP + PSQ + QSP + QS(1)Q
//  SPAM Algorithm:
//
//  for iter_outer in 1 : max_iter_outer:
//    project H and S into current krylov space, get sub_H, sub_S ---------------(1)
//    solve for sub_H * x = lambda * sub_S * x for x and lambda
//    generate eigenvector evec by taking linear combination of krylov vectors with coefficient x
//    generate residual vector r = H * evec - lambda * S * evec
//    evaluate residual norm = ||r||
//    if ||r|| < threshold
//      exit
//    else
//      add r to a intermediate space, which contain current krylov vectors
//      orthogonalize r with other vectors in the space
//      for iter_inner in 1 : max_iter_inner:
//        project H_bar and S_bar in to intermediate space, get sub_hy_H, sub_hy_S
//        solve for sub_hy_H * x = lambda * sub_hy_S * x for x and lambda
//        generate eigenvector evec by taking linear combination of krylov vectors with coefficient x
//        generate residual vector r = H_bar * evec - lambda * S_bar * evec
//        evaluate residual norm = ||r||
//        if ||r|| < threshold
//          exit
//        else
//          add r to a intermediate space, which contain current krylov vectors
//      add all vectors in generated in the inner loop into krylov space, go to (1)
//////////////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef ENGINE_SPAM_SOLVER_HEADER
#define ENGINE_SPAM_SOLVER_HEADER

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
#include<formic/utils/lmyengine/eigen_solver.h>
#include<formic/utils/matrix.h>

namespace cqmc {
  
  namespace engine { 
  
  class SpamLMHD : public EigenSolver { 

private: 

  /// \brief total number of krylov vectors
  int _nkry;

  /// \brief number of krylov vectors that have been multiplied by full hamiltonian and overlap matrix 
  int _nkry_full;

  /// \brief maximum number of iterations allowed in SPAM outer loop
  int _n_max_iter;

  /// \brief maximum number of iterations allowed in SPAM inner loop
  int _inner_maxIter;

  /// \brief integer to tell the approximate degree
  int _appro_degree;

  /// \brief flag that tell whether to print SPAM inner iteration results
  bool _inner_print;

  /// \brief the smallest subspace overlap matrix singular value of inner loop
  double _smallest_sin_value_inner;

  /// \brief the smallest subspace overlap matrix singular value of outer loop
  double _smallest_sin_value_outer;

  /// \brief subspace overlap matrix singular value threshold
  double _singular_value_threshold;

  /// \brief the subspace eigenvalue of inner loop 
  double _sub_eval_inner;

  /// \brief the subspace eigenvalue of outer loop
  double _sub_eval_outer;

  /// \brief best residual of outer loop
  double _best_residual_outer;

  /// \brief best residual of inner loop
  double _best_residual_inner;

  /// \brief initial energy used in resetting 
  double _init_energy;

  /// \brief energy of outer loop 
  double _energy_outer;

  /// \brief energy of inner loop
  double _energy_inner;

  /// \brief approximate prefactor
  double _appro_factor;

  /// \brief approximate derivative vectors
  formic::Matrix<double> & _der_rat_appro;

  /// \brief approximate energy derivatives
  formic::Matrix<double> & _le_der_appro;

  /// \brief normalized krylov space basis vectors 
  formic::Matrix<double> _kvecs;

  /// \brief krylov vectors that are about to be added in to outer loop
  formic::Matrix<double> _kvecs_about_to_add;

  /// \brief vectors resulting from Hamiltonian action on the krylov basis 
  formic::Matrix<double> _hvecs; 

  /// \brief vectors resulting from Hamiltonian transpose action on krylov basis 
  formic::Matrix<double> _thvecs;

  /// \brief vectors resulting from approximated Hamiltonian action on krylov basis 
  formic::Matrix<double> _ahvecs;

  /// \brief vectors resulting from approximated Hamiltonian transpose action on krylov basis 
  formic::Matrix<double> _athvecs;

  /// \brief vectors resulting from hybrid Hamiltonian action on krylov basis 
  formic::Matrix<double> _hhvecs;

  /// \brief vectors resulting from overlap matrix action on the krylov basis 
  formic::Matrix<double> _svecs;

  /// \brief vectors resulting from approximated overlap matrix action on the krylov basis 
  formic::Matrix<double> _asvecs;

  /// \brief vectors resulting from hybrid overlap matrix action on the krylov basis 
  formic::Matrix<double> _hsvecs;

  /// \brief subspace Hamiltonian matrix 
  formic::Matrix<double> _subH;

  /// \brief hybrid subspace Hamiltonian matrix 
  formic::Matrix<double> _hy_subH;

  /// \brief subspace overlap matrix 
  formic::Matrix<double> _subS;

  /// \brief hybrid subspace overlap matrix 
  formic::Matrix<double> _hy_subS;

	/// \brief work vector 
  formic::ColVec<double> _wv1;

  /// \brief work vector 
  formic::ColVec<double> _wv2;

  /// \brief work vector 
  formic::RowVec<double> _wv3;

  /// \brief work vector 
  formic::ColVec<std::complex<double> > _wv4;

  /// \brief work vector 
  formic::ColVec<std::complex<double> > _wv5;

  /// \brief work vector
  formic::ColVec<double> _wv6;

  /// \brief work matrix 
  formic::Matrix<double> _wm1;

  /// \brief work matrix 
  formic::Matrix<double> _wm2;

  /// \brief work matrix 
  formic::Matrix<double> _wm3;

  /// \brief work matrix
  formic::Matrix<double> _wm4;
		
  /// \brief subspace eigenvector of outer loop 
  formic::ColVec<double> _sub_evec_outer;

  /// \brief subspace eigenvector of inner loop 
  formic::ColVec<double> _sub_evec_inner;

  /// \brief calculated eigenvector of inner loop
  formic::ColVec<double> _evecs_inner;

  public:

  ///////////////////////////////////////////////////////////////////////////////
  // \brief solve the subspace generalized eigenvalue problem
  //        with nonsymmetric H and S 
  //
  //
  ///////////////////////////////////////////////////////////////////////////////

  void solve_subspace_nonsymmetric(const bool outer);

  ///////////////////////////////////////////////////////////////////////////////
  // \brief solves the subspace eigenvalue problem
  //
  //
  //
  ///////////////////////////////////////////////////////////////////////////////

  void solve_subspace(const bool outer);

public:
  
  ///////////////////////////////////////////////////////////////////////////////
  // \brief constructor with given parameters
  //
  //
  //
  ///////////////////////////////////////////////////////////////////////////////

  SpamLMHD(const formic::VarDeps * dep_ptr, 
           const int nfds,
           const int lm_krylov_iter,
           const int inner_maxIter,
           const int appro_degree,
           const bool inner_print,
           const double lm_eigen_thresh,
           const double lm_min_S_eval,
           const double appro_factor,
           const bool var_deps_use,
           const bool chase_lowest,
           const bool chase_closest,
           const bool ground,
           const std::vector<double> & vf,
           const double init_energy,
           const double hd_shift,
           const double lm_max_e_change,
           const double total_weight,
           const double vgsa,
           formic::Matrix<double> & der_rat,
           formic::Matrix<double> & le_der,
           formic::Matrix<double> & der_rat_appro,
           formic::Matrix<double> & le_der_appro);

  ////////////////////////////////////////////////////////////////////////////////////
  // \brief evaluate tau and compute the "variance corrected" Hamiltonian
  // 
  // 
  ////////////////////////////////////////////////////////////////////////////////////

  void tau_and_correct_ham()
  {}

  ////////////////////////////////////////////////////////////////////////////////////
  // \brief adds a new krylov vector 
  //
  // \param[in]  v  vector that needs to be added
  //
  ////////////////////////////////////////////////////////////////////////////////////

  void add_krylov_vector(const formic::ColVec<double> & v)
  {}

  ////////////////////////////////////////////////////////////////////////////////
  // \brief adds a new Krylov basis vector for spam inner loop
  //
  //
  //
  ////////////////////////////////////////////////////////////////////////////////

  void add_krylov_vector_inner(const formic::ColVec<double> & v);

  ////////////////////////////////////////////////////////////////////////////////
  // \brief adds a bunch of new Krylov vectors(which is a matrix) for spam outer 
  //        loop
  //
  //
  ////////////////////////////////////////////////////////////////////////////////

  void add_krylov_vectors_outer(const formic::Matrix<double> & m);

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

  void HMatVecOp(const formic::ColVec<double> & x, formic::ColVec<double> & y, const bool transpose, const bool approximate);  

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

  void HMatMatOp(const formic::Matrix<double> & x, formic::Matrix<double> & y, const bool transpose, const bool approximate);

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

  void SMatVecOp(const formic::ColVec<double> & x, formic::ColVec<double> & y, const bool approximate);

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

  void SMatMatOp(const formic::Matrix<double> & x, formic::Matrix<double> & y, const bool approximate = false);

  ////////////////////////////////////////////////////////////////////////////////////
  // \brief updates hamiltonian * krylov vector and hamiltonian projection based on 
  //        new shift 
  //
  //
  ////////////////////////////////////////////////////////////////////////////////////

  void update_hvecs_sub(const double new_i_shift, const double new_s_shift);

  ////////////////////////////////////////////////////////////////////////////////////
  // \brief reset the eigen solver 
  // 
  // \brief clear subspace projection of Hamiltonian and overlap matrix, clear Krylov
  //        subspace and action of Hamiltonian and overlap matrix
  ////////////////////////////////////////////////////////////////////////////////////

  void child_reset();
	
  /////////////////////////////////////////////////////////////////////////////////////
  // \brief solves the eigenvalue problem via the normal davidson method
  //
  //
  //
  /////////////////////////////////////////////////////////////////////////////////////
  	
  bool iterative_solve(double & eval, std::ostream & output);
		
  ///////////////////////////////////////////////////////////////////////////////////////////////
  // \brief solves the eigenvalue problem
  //
  //
  //
  ///////////////////////////////////////////////////////////////////////////////////////////////

  bool solve(double & eval, std::ostream & output);
	
  ///////////////////////////////////////////////////////////////////////////////////////////////
  // \brief converts eigenvectors into wave function coefficients
  //        solving this question Hc = Sd for d, S is the ground overlap matrix 
  //
  //
  ///////////////////////////////////////////////////////////////////////////////////////////////

  void child_convert_to_wf_coeff();
};

}

}

#endif
