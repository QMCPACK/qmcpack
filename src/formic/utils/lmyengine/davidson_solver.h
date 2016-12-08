//////////////////////////////////////////////////////////////////////////////////////////////////
// \file davidson_solver.h       a class to perform generalized eigenvalue problem
//
//
//
//  solve with normal davidson method(either build or not build the matrix)
//////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef EIGINE_DAVIDSON_SOLVER_HEADER
#define EIGINE_DAVIDSON_SOLVER_HEADER

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

  class DavidsonLMHD : public EigenSolver{

private:

	/// \brief number of krylov vectors
	int _nkry;
		 
	/// \brief maximum number of iteration allowed in the solver
	int _n_max_iter;

	/// \brief the smallest subspace overlap matrix singular value
	double _smallest_sin_value;

	/// \brief the subspace eigenvalue(in harmonic Davidson this is not the energy)
	double _sub_eval;	
			
	/// \brief subspace overlap matrix singular value threshold
	double _singular_value_threshold;

	// \brief best residual
	double _best_residual;

  /// \brief difference between new and old tau
  double _tau_diff;

  // \brief flag to tell whether we have built the matrix or not 
  bool _build_lm_matrix;

	/// \brief the Hamiltonian matrix
	formic::Matrix<double> & _hmat;

	/// \brief the overlap matrix 
	formic::Matrix<double> & _smat;

  /// \brief the overlap matrix of normal linear method
  formic::Matrix<double> & _lmsmat;

	/// \brief the preconditioning matrix
	//formic::Matrix<double> _mmat;

	/// \brief inverse-square-root of the diagonal of the overlap matrix
  //std::vector<double> _ovl_diag_neg_half;

	/// \brief normalized krylov space basis vectors
	formic::Matrix<double> _kvecs;

	/// \brief vectors resulting from Hamiltonian action on the krylov basis 
	formic::Matrix<double> _hvecs;

	/// \brief vectors resulting from Hamiltonian transpose action on the krylov basis 
	formic::Matrix<double> _htvecs;

	/// \brief vectors resulting from overlap matrix action on the krylov basis 
	formic::Matrix<double> _svecs;

	/// \brief subspace Hamiltonian matrix
	formic::Matrix<double> _subH;

	/// \brief subspace overlap matrix 
	formic::Matrix<double> _subS;

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
		
	/// \brief work vector
  formic::ColVec<double> _wvX;

	/// \brief work vector
  formic::ColVec<double> _wvY;

	/// \brief subspace eigenvector 
  formic::ColVec<double> _sub_evec;


  public:
		
	///////////////////////////////////////////////////////////////////
	// \brief solve the subspace generalized eigenvalue problem 
	//        with nonsymmetric H and S 
	//
	//
	///////////////////////////////////////////////////////////////////

	void solve_subspace_nonsymmetric(const bool outer);
	
	///////////////////////////////////////////////////////////////////////////////
	// \brief solves the subspace eigenvalue problem
	//
	//
	//
	///////////////////////////////////////////////////////////////////////////////

	void solve_subspace(const bool outer); 
	

public:
		
	////////////////////////////////////////////////////////////////////////////////
	// \brief constructor with given parameters
	//
	//
	//
	/////////////////////////////////////////////////////////////////////////////////
	
	DavidsonLMHD(const formic::VarDeps * dep_ptr,
               const int nfds,
               const int lm_krylov_iter,
               const double lm_eigen_thresh,
               const double lm_min_S_eval,
               const bool var_deps_use,
               const bool chase_lowest,
               const bool chase_closest,
               const bool ground,
               const bool variance_correct,
               const bool build_lm_matrix,
               const std::vector<double> & vf,
               const double init_cost,
               const double init_variance,
               const double hd_shift,
               const double var_weight,
               const double lm_max_e_change,
               const double total_weight,
               const double vgsa,
               formic::Matrix<double> & der_rat,
               formic::Matrix<double> & le_der,
               formic::Matrix<double> & hmat,
               formic::Matrix<double> & smat,
               formic::Matrix<double> & lmsmat);

  ////////////////////////////////////////////////////////////////////////////////////
  // \brief evaluate tau and compute the "variance corrected" Hamiltonian
  // 
  // 
  ////////////////////////////////////////////////////////////////////////////////////

  void tau_and_correct_ham();

	////////////////////////////////////////////////////////////////////////////////////
	// \brief adds a new Krylov basis vector for normal davidson method 
	//
	//
	//
	////////////////////////////////////////////////////////////////////////////////////
		
	void add_krylov_vector(const formic::ColVec<double> & v);

  ////////////////////////////////////////////////////////////////////////////////
  // \brief adds a new Krylov basis vector for spam inner loop
  //
  //
  //
  ////////////////////////////////////////////////////////////////////////////////

  void add_krylov_vector_inner(const formic::ColVec<double> & v)
  {}

  ////////////////////////////////////////////////////////////////////////////////
  // \brief adds a bunch of new Krylov vectors(which is a matrix) for spam outer 
  //        loop
  //
  //
  ////////////////////////////////////////////////////////////////////////////////

  void add_krylov_vectors_outer(const formic::Matrix<double> & m)
  {}

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

  void HMatVecOp(const formic::ColVec<double> & x, formic::ColVec<double> & y, const bool transpose = false, const bool approximate = false);  

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

  void HMatMatOp(const formic::Matrix<double> & x, formic::Matrix<double> & y, const bool transpose, const bool approximate = false)
  {}

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

  void SMatVecOp(const formic::ColVec<double> & x, formic::ColVec<double> & y, const bool approximate = false);

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

  void SMatMatOp(const formic::Matrix<double> & x, formic::Matrix<double> & y, const bool approximate = false)
  {}

  //void MMatVecOp(const formic::ColVec<double> & x, formic::ColVec<double> & y, const int mmat_first);

  ////////////////////////////////////////////////////////////////////////////////////
  // \brief updates hamiltonian * krylov vector and hamiltonian projection based on 
  //        new shift 
  //
  //
  ////////////////////////////////////////////////////////////////////////////////////

  void update_hvecs_sub(const double new_i_shift, const double new_s_shift);

	////////////////////////////////////////////////////////////////////////////////////
	// \brief update Hamiltonian and overlap matrix for the solver 
	//
	//
	//
	////////////////////////////////////////////////////////////////////////////////////

	void update_hamovlp(formic::Matrix<double> & hmat, formic::Matrix<double> & smat);
	
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

  //////////////////////////////////////////////////////////////////////////////////////////////
  // \brief build linear matrix(used when doing debug)
  //
  //
  //
  //////////////////////////////////////////////////////////////////////////////////////////////
  
  //const formic::Matrix<double> debug_build()
  //{
  //  formic::Matrix<double> retval(2, 2, 0.0);
  //  retval.at(0,0) = 1.0;
  //  retval.at(1,1) = 1.0;
  //  return retval;
  //}

	//////////////////////////////////////////////////////////////////////////////////////////////
	// \brief return the wave function coefficients
	//
	//
	//
	//////////////////////////////////////////////////////////////////////////////////////////////

	//formic::ColVec<double> wf_coeff();
	

};

}

}

#endif

