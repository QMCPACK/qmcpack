//////////////////////////////////////////////////////////////////////////////////////////////////
// \brief a class to perform generalized eigenvalue problem
//
//
//
//  solve with normal davidson method
//////////////////////////////////////////////////////////////////////////////////////////////////

#include<vector>
#include<string>
#include<numeric>
#include<cassert>
#include<algorithm>
#include<cmath>
#include<complex>
#include<iostream>
//#include<mpi.h>

#include<boost/format.hpp>
#include<boost/shared_ptr.hpp>

#include<formic/utils/exception.h>
#include<formic/utils/matrix.h>
#include<formic/utils/mpi_interface.h>
#include<formic/utils/lapack_interface.h>
#include<formic/utils/lmyengine/eigen_solver.h>
#include<formic/utils/lmyengine/davidson_solver.h>


//////////////////////////////////////////////////////////////////////////////////////////////////
// \brief solve the subspace generalized eigenvalue problem 
//        with nonsymmetric H and symmetric S 
//
//
//////////////////////////////////////////////////////////////////////////////////////////////////

void cqmc::engine::DavidsonLMHD::solve_subspace_nonsymmetric(const bool outer) 
{

  // get rank number and number of ranks 
  int my_rank = formic::mpi::rank(); 
  //int num_rank;
  //MPI_Comm_rank(MPI_COMM_WORLD, & my_rank);
  //MPI_Comm_size(MPI_COMM_WORLD, & num_rank);

  //one and zero in complex form
  const std::complex<double> complex_one(1.0, 0.0);
  const std::complex<double> complex_zero(0.0, 0.0);

  //get subspace dimension
  const int m = _nkry;

  //create vectors and matrices used in svd routine
  formic::Matrix<double> u, v, vt;
  formic::ColVec<double> sin_vals;
  int truncate_index;

  // make sure the subspace matrix is not empty
  if ( _subS.rows() == 0 || _subH.rows() == 0 ) 
    throw formic::Exception("subspace matrix is empty upon solving subspace eigenvalue problem");

  //perform svd to subspace overlap matrix 
  _subS.svd(u, sin_vals, vt);
  v = vt.clone();
  v.tip();

  //record the smallest singular value of subspace overlap matrix
  _smallest_sin_value = std::abs(sin_vals.at(0));
  for(int i = 0; i < m; i++){
    _smallest_sin_value = std::min(_smallest_sin_value, std::abs(sin_vals(i)));
    truncate_index = i;

    //check if singular value is smaller that the singular value threshold
    if( _smallest_sin_value < _singular_value_threshold ){
      break;
    }
  }

  //get the number of colums of new U and V matrix by add 1 to truncate index
  truncate_index ++;

  //throw away those columns in U and V matrix which corresponds to singular values below the threshold
  u.conservativeResize(m, truncate_index);
  v.conservativeResize(m, truncate_index);

  //throw away those small singular values
  sin_vals.conservativeResize(truncate_index);

  //convert the truncated singular value vector to a truncated_index * truncated_index diagonal matrix
  formic::Matrix<double> trun_sin_val_matrix(truncate_index, truncate_index, 0.0);
  for (int i = 0; i < truncate_index; i++)
    trun_sin_val_matrix.at(i,i) = sin_vals.at(i);

  //calculate the inverse of this matrix
  for(int i = 0; i < truncate_index; i++){
    for(int j = 0; j < truncate_index; j++){
      trun_sin_val_matrix.at(i, j) = (trun_sin_val_matrix.at(i, j) == 0.0 ? trun_sin_val_matrix.at(i, j) : 1.0 / trun_sin_val_matrix.at(i, j));
    }
  }

  //calculate matrix S_trun^-1 * U^-1 * H * V
  formic::Matrix<double> new_sub_H(truncate_index, truncate_index, 0.0); 
  new_sub_H = trun_sin_val_matrix * u.t() * _subH * v;

  //solve this standard eigenvalue problem ( new_H * y = lambda * y, y = V^T * x, x is the original eigenvector)
  formic::ColVec<std::complex<double> > e_evals;
  formic::Matrix<std::complex<double> > e_evecs;
  new_sub_H.nonsym_eig(e_evals, e_evecs);

  // set an eigen array to hold energies
  formic::ColVec<std::complex<double> > energy_list = e_evals.clone();

  //if we are doing excited state calculations, convert the resulting eigenvalues to energy( this is important in harmonic davidson)
  //if ( !_ground ) {
  //  energy_list = _hd_shift * complex_one - complex_one / energy_list;
  //}
  
  int selected = 0;
  //if we want to chase the closest, selected the eigenvalue that is most similar to the previous eigenvalue( this is essestially an attempt to stay in the same solution)
  if (_chase_closest) {
    std::complex<double> closest_cost = e_evals.at(0);
    for (int j = 1; j < truncate_index; j++){
      if(std::abs( complex_one * _cost - e_evals.at(j) ) < std::abs(complex_one * _cost - closest_cost )){
        selected = j;
        closest_cost = e_evals.at(j);
      }
    }

    // if the eigenvalue has an imaginary component, we abort 
    _eval_was_complex = false;
    if( std::abs(closest_cost.imag()) > 1.0e-6 ) {
      _eval_was_complex = true;
      return;
    }

    // if the eigenvalue is real, we record it and the corresponding eigenvector 
    // record energy 
    _cost = closest_cost.real();

    // record the eigenvalue 
    _sub_eval = _cost; 

    // record the eigenvector y 
    _wv4 = e_evecs.col_as_vec(selected);

    // convert y to x 
    _wv5.reset(m);
    // convert v to complex form to make sure that all quantities in dgemm call are of the same type
    formic::Matrix<std::complex<double> > v_complex(v.rows(), v.cols(), complex_zero);
    for (int i = 0; i < v_complex.rows(); i++) {
      for (int j = 0; j < v_complex.cols(); j++) {
        v_complex.at(i,j) = std::complex<double>(v.at(i,j), 0.0);
      }
    }
    formic::xgemm('N', 'N', m, 1, truncate_index, complex_one, &v_complex.at(0,0), m, &_wv4.at(0), truncate_index, complex_zero, &_wv5.at(0), m);
			
    // take real part of the vector x, put that into eigenvector 
    _sub_evec.reset(m);
    for (int i = 0; i < m; i++) {
      _sub_evec.at(i) = _wv5.at(i).real();
    }
  }

  // if we want to chase the lowest, selected the lowest eigenvalue
  if (_chase_lowest) {
    std::complex<double> lowest_eval = e_evals.at(0);
    for (int j = 1; j < truncate_index; j++){
      if (e_evals.at(j).real() < lowest_eval.real()) {
        selected = j;
        lowest_eval = e_evals.at(j);
      }
    }

    // if the eigenvalue has an imaginary component, we abort
    _eval_was_complex = false;
    if( std::abs(lowest_eval.imag()) > 1.0e-6 ) {
      _eval_was_complex = true;
      return;
    }

    // if the eigenvalue is ture, we record it and the corresponding eigenvector
    // if we are doing excited state calculation, convert and record energy
    //if ( !_ground )
    //  _energy = _hd_shift - 1.0 / lowest_eval.real();
    //else
    _cost = lowest_eval.real();

    // record the eigenvalue
    _sub_eval = lowest_eval.real();

    // record the eigenvector y
    _wv4 = e_evecs.col_as_vec(selected);

    // convert y to x
    _wv5.reset(m);
    // convert v to complex form to make sure that all quantities in dgemm call are of the same type
    formic::Matrix<std::complex<double> > v_complex(v.rows(), v.cols(), complex_zero);
    for (int i = 0; i < v_complex.rows(); i++) {
      for (int j = 0; j < v_complex.cols(); j++) {
        v_complex.at(i,j) = std::complex<double>(v.at(i,j), 0.0);
      }
    }
    formic::xgemm('N', 'N', m, 1, truncate_index, complex_one, &v_complex.at(0,0), m, &_wv4.at(0), truncate_index, complex_zero, &_wv5.at(0), m);
			
    // take real part of vector x, put that into eigenvector
    _sub_evec.reset(m);
    for (int i = 0; i < m; i++){
      _sub_evec.at(i) = _wv5.at(i).real();
    }
  }
}


///////////////////////////////////////////////////////////////////////////////
// \brief solves the subspace eigenvalue problem
//
//
//
///////////////////////////////////////////////////////////////////////////////

void cqmc::engine::DavidsonLMHD::solve_subspace(const bool outer) 
{
  this -> solve_subspace_nonsymmetric(outer);
  return;

}

////////////////////////////////////////////////////////////////////////////////
// \brief constructor
//
//
//
/////////////////////////////////////////////////////////////////////////////////
	
cqmc::engine::DavidsonLMHD::DavidsonLMHD(const formic::VarDeps * dep_ptr,
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
                                         formic::Matrix<double> & lmsmat)
:EigenSolver(dep_ptr,
             nfds,
             lm_eigen_thresh,
             var_deps_use,
             chase_lowest,
             chase_closest,
             ground,
             variance_correct,
             vf,
             init_cost,
             init_variance,
             hd_shift,
             var_weight,
             lm_max_e_change,
             total_weight,
             vgsa,
             der_rat,
             le_der),
_build_lm_matrix(build_lm_matrix),
_nkry(0),
_n_max_iter(lm_krylov_iter),
_singular_value_threshold(lm_min_S_eval),
_hmat(hmat),
_smat(smat),
_lmsmat(lmsmat),
_smallest_sin_value(0.0)
{

  // get rank number and number of ranks 
  int my_rank = formic::mpi::rank(); 
  //int num_rank;
  //MPI_Comm_rank(MPI_COMM_WORLD, & my_rank);
  //MPI_Comm_size(MPI_COMM_WORLD, & num_rank);

  // compute preconditioning matrix as S^(-1/2)
//  if ( _build_lm_matrix ) {
//
//    // first get the inverse square roots of the overlap's diagonal
//    _ovl_diag_neg_half.resize(_smat.rows());
//    for (int i = 0; i < _smat.rows(); i++)
//      _ovl_diag_neg_half[i] = 1.0 / std::sqrt(std::max(1.0e-10, _smat(i,i)));
//
//    // transform non-psi block of overlap matrix by these inverse square roots
//    formic::Matrix<double> stilde;
//    stilde.reset(_smat.rows()-1, _smat.cols()-1);
//    for (int i = 0; i < stilde.rows(); i++)
//    for (int j = 0; j < stilde.cols(); j++)
//      stilde.at(i,j) = _ovl_diag_neg_half[i+1] * _smat.at(i+1,j+1) * _ovl_diag_neg_half[j+1]; // + ( i == j ? 1.0e-3 : 0.0 );
//
//    // diagonalize this transformed overlap
//    formic::Matrix<std::complex<double> > v;
//    formic::Matrix<double> v_real;
//    formic::ColVec<std::complex<double> > w;
//    formic::ColVec<double> w_real;
//    stilde.nonsym_eig(w, v);
//    w_real.reset(w.size());
//    v_real.reset(v.rows(), v.cols());
//
//    // compute the pseudo-stilde^(-1/2) matrix
//    int nonzero_count = 0;
//    for (int i = 0; i < w.size(); i++) {
//      if ( w.at(i).real() > 1.0e-8 ) {
//        nonzero_count++;
//        w_real.at(i) = 1.0 / std::sqrt(w.at(i).real());
//      } else {
//        w_real.at(i) = 0.0;
//      }
//    }
//
//    for (int i = 0; i < v.rows(); i++)
//    for (int j = 0; j < v.cols(); j++)
//      v_real.at(i,j) = v.at(i,j).real();
//
//    formic::Matrix<double> diagonal_w(v.rows(), v.cols(), 0.0);
//    for (int i = 0; i < diagonal_w.cols(); i++)
//      diagonal_w.at(i,i) = w_real.at(i);
//    formic::Matrix<double> temp = v_real * diagonal_w * v_real.inv();
//
//    // fill in the preconditioning matrix
//    _mmat.reset(_smat.rows(), _smat.cols());
//    for (int i = 0; i < _smat.rows(); i++) {
//      _mmat.at(i,0) = 0.0;
//      _mmat.at(0,i) = 0.0;
//    }
//    _mmat.at(0,0) = 1.0;
//    for (int i = 0; i < stilde.rows(); i++)
//    for (int j = 0; j < stilde.cols(); j++)
//      //_mmat(i+1,j+1) = ( i == j ? 1.0 : 0.0 );
//      _mmat.at(i+1,j+1) = temp.at(i,j);

    //// test that we have preconditioned corrrectly
    //Eigen::MatrixXd test1;
    //test1.resize(_smat.rows(), _smat.cols());
    //for (int i = 0; i < _smat.rows(); i++)
    //for (int j = 0; j < _smat.cols(); j++)
    //  test1(i,j) = _ovl_diag_neg_half[i] * _smat(i,j) * _ovl_diag_neg_half[j];
    //Eigen::MatrixXd test2 = _mmat * test1 * _mmat;
    //formic::of << std::endl;
    //formic::of << "test matrix:" << std::endl;
    //for (int i = 0; i < test2.rows(); i++) {
    //  for (int j = 0; j < test2.cols(); j++) {
    //    formic::of << boost::format(" %12.8f") % test2(i,j);
    //  }
    //  formic::of << std::endl;
    //}
    //formic::of << std::endl;

//    //Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(_smat);
//    //Eigen::VectorXd w = es.eigenvalues();
//    //Eigen::MatrixXd v = es.eigenvectors();
//    Eigen::JacobiSVD<Eigen::MatrixXd> svd(_smat, Eigen::ComputeThinU | Eigen::ComputeThinV);
//    Eigen::VectorXd w = svd.singularValues();
//    formic::of << std::endl;
//    for (int i = 0; i < w.size(); i++)
//      formic::of << boost::format("w( %4i ) = %12.4e") % i % w(i) << std::endl;
//    formic::of << std::endl;
//    for (int i = 0; i < w.size(); i++)
//      w(i) = ( w(i) > 1.0e-8 ? 1.0 / std::max(1.0e-15, w(i)) : 0.0 );
//      //w(i) = 1.0 / std::sqrt(std::max(1.0e-10, w(i)));
//    //_mmat.resize(_smat.rows(), _smat.cols());
//    _mmat = svd.matrixV() * w.asDiagonal() * svd.matrixU().inverse();
//    double max_error = 0.0;
//    double max_elem = 0.0;
//    Eigen::MatrixXd test = _mmat * _smat;
//    formic::of << "test matrix:" << std::endl;
//    for (int i = 0; i < w.size(); i++) {
//      for (int j = 0; j < w.size(); j++) {
//        max_error = std::max(max_error, std::abs(test(i,j) - (i==j?1.0:0.0)));
//        max_elem  = std::max(max_elem , std::abs(_smat(i,j)));
//        formic::of << boost::format(" %12.8f") % test(i,j);
//      }
//      formic::of << std::endl;
//    }
//    formic::of << std::endl;
//    formic::of << boost::format("w.size() = %i") % w.size() << std::endl;
//    formic::of << boost::format("_smat.rows() = %i") % _smat.rows() << std::endl;
//    formic::of << boost::format("inverting _smat with max_elem = %.2e gave an error of %.2e") % max_elem % max_error << std::endl;

//  }

  // initialize the eigenvector
  _evecs.reset( ( _var_deps_use ? 1 + _dep_ptr->n_tot() : _nfds ), 0.0 );
  _evecs.at(0) = 1.0;

  // initialize the independent eigenvector 
  _ind_evecs.reset( ( _var_deps_use ? 1 + _dep_ptr->n_ind() : _nfds), 0.0 );
  _ind_evecs.at(0) = 1.0;
}

////////////////////////////////////////////////////////////////////////////////////
// \brief evaluate tau and compute the "variance corrected" Hamiltonian
// 
// 
////////////////////////////////////////////////////////////////////////////////////

void cqmc::engine::DavidsonLMHD::tau_and_correct_ham()
{
  // get rank number and number of ranks
  int my_rank = formic::mpi::rank();
  //int num_rank;
  //MPI_Comm_rank(MPI_COMM_WORLD, & my_rank);
  //MPI_Comm_rank(MPI_COMM_WORLD, & num_rank);

  if ( my_rank == 0 ) {
    
    // first remember the old tau
    double old_tau = _tau;

    // evaluate new tau
    _tau = dotc(_ind_evecs, _smat * _ind_evecs)  / dotc(_ind_evecs, _lmsmat * _ind_evecs);

    // evaluate the difference between old and new tau
    _tau_diff = std::abs(old_tau - _tau);

    // modify hamiltonian matrix based on the difference between new and old tau
    double square_diff = _tau * _tau - old_tau * old_tau;

    // update variance modified hamiltonian matrix 
    _hmat += -1.0 * _var_weight * 10.0 * square_diff * _lmsmat;

  }
}

////////////////////////////////////////////////////////////////////////////////////
// \brief adds a new Krylov basis vector for normal davidson method 
//
//
//
////////////////////////////////////////////////////////////////////////////////////
		
void cqmc::engine::DavidsonLMHD::add_krylov_vector(const formic::ColVec<double> & v)
{
	
  // get rank number and number of ranks 
  int my_rank = formic::mpi::rank(); 
  //int num_rank;
  //MPI_Comm_rank(MPI_COMM_WORLD, & my_rank);
  //MPI_Comm_size(MPI_COMM_WORLD, & num_rank);

  // check vector length
  if (my_rank == 0 && !_build_lm_matrix && v.size() != _der_rat.cols())  
    throw formic::Exception("bad vector length of %d in DavidsonLMHD::add_krylov_vector: expected length of %d") % v.size() % _der_rat.cols();

  else if (my_rank == 0 && _build_lm_matrix && v.size() != _nfds) 
    throw formic::Exception("bad vector length of %d in DavidsonLMHD::add_krylov_vector: expected length of %d") % v.size() % _nfds;

  // increment krylov subspace size and remember the old size 
  const int nold = _nkry++;

  // copy the new vector to work vector 
  _wv1.reset(v.size());
  _wv1 = v.clone();

  // perform graham schmidt orthogonolization against the existing krylov vectors
  if (my_rank == 0) {
    for (int i = 0; i < nold; i++) {
      _wv1 -= (dotc(_svecs.col_as_vec(i), _wv1) * _kvecs.col_as_vec(i));
    }
  }

  // if we don't build the matrix, we need to broadcast this vector to all processes
  if ( !_build_lm_matrix ) 
    formic::mpi::bcast(&_wv1.at(0), _wv1.size());
    //MPI_Bcast(&_wv1.at(0), _wv1.size(), MPI_DOUBLE, 0, MPI_COMM_WORLD);

  // compute the product of the overlap matrix times the new krylov vector on the root process
  if ( _build_lm_matrix ) {
    if ( my_rank == 0 ) {
      //this -> MMatVecOp(_wv1, _wvX, 1); // preconditioning
      this -> SMatVecOp(_wv1, _wv2);
      //this -> MMatVecOp(_wvY, _wv2, 0); // preconditioning
    }
  }

  else {
    this -> SMatVecOp(_wv1, _wv2);
    // reduce the result vector
    formic::ColVec<double> _wv2_avg(_wv2.size());
    formic::mpi::reduce(&_wv2.at(0), &_wv2_avg.at(0), _wv2.size(), MPI_SUM);
    _wv2 = _wv2_avg.clone();
  }

  // compute the product of the hamiltonian matrix times the new krylov vector
  formic::ColVec<double> hs(_nfds);
  if ( _build_lm_matrix ) {
    if ( my_rank == 0 ) {
      //this -> MMatVecOp(_wv1, _wvX, 1); // preconditioning
      this -> HMatVecOp(_wv1, hs);
      //this -> MMatVecOp(_wvY,   hs, 0); // preconditioning
    }
  }
  else {
    this -> HMatVecOp(_wv1, hs);
    // reduce the result vector 
    formic::ColVec<double> hs_avg(hs.size());
    formic::mpi::reduce(&hs.at(0), &hs_avg.at(0), hs.size(), MPI_SUM);
    hs = hs_avg.clone();
    //if (my_rank == 0) {
    //  for (int i = 0; i < hs_avg.size(); i++) 
    //    std::cout << boost::format("%10.8e ") % h
  }

  // modify hamiltonian product to account for "identity" shift 
  if (my_rank == 0) {
    for (int i = 1; i < hs.size(); i++) {
      hs.at(i) += _hshift_i * _wv1.at(i);
    }
    //this -> MMatVecOp(_wv1, _wvX, 1); // preconditioning
    //_wvY.reset(_wvX.size());
    //_wvY.at(0) = 0.0;
    //for (int i = 1; i < hs.size(); i++)
    //  _wvY.at(i) = _hshift_i * _wvX.at(i);
    //this -> MMatVecOp(_wvY, _wvX, 0); // preconditioning
    //for (int i = 0; i < hs.size(); i++)
    //  hs.at(i) += _wvX.at(i);
  }

  // modify hamiltonian product to account for "overlap" shift
  if (my_rank == 0 && nold > 0) {
    hs += _hshift_s * _wv2;
  }

  // normalize the new krylov vector and save the vector and its operation on overlap matrix 
  if (my_rank == 0) {
    const double norm = std::sqrt(dotc(_wv1,_wv2));
    _wv1 /= norm;
    _wv2 /= norm;
    hs /= norm;
  
    _kvecs.conservativeResize(_nfds, _nkry);
    std::copy(_wv1.begin(), _wv1.end(), _kvecs.col_begin(_nkry-1));
    _svecs.conservativeResize(_nfds, _nkry);
    std::copy(_wv2.begin(), _wv2.end(), _svecs.col_begin(_nkry-1));

    // multiply new vector by the Hamiltonian
    _hvecs.conservativeResize(_nfds, _nkry);
    std::copy(hs.begin(), hs.end(), _hvecs.col_begin(_nkry-1));

    // update subspace projection of the Hamiltonian
    _subH.conservativeResize(_nkry, _nkry);
    _wv3 = _wv1.t() * _hvecs;
    _wv6 = _kvecs.t() * hs;
    std::copy(_wv3.begin(), _wv3.end(), _subH.row_begin(_nkry-1));
    std::copy(_wv6.begin(), _wv6.end(), _subH.col_begin(_nkry-1));

    // update subspace projection of the overlap
    _subS.conservativeResize(_nkry, _nkry);
    _wv3 = _wv1.t() * _svecs;
    _wv6 = _kvecs.t() * _wv2;
    std::copy(_wv3.begin(), _wv3.end(), _subS.row_begin(_nkry-1));
    std::copy(_wv6.begin(), _wv6.end(), _subS.row_begin(_nkry-1));
  }

}

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

void cqmc::engine::DavidsonLMHD::HMatVecOp(const formic::ColVec<double> & x, formic::ColVec<double> & y, const bool transpose, const bool approximate)
{

  int my_rank = formic::mpi::rank(); 
  //int num_rank;
  //MPI_Comm_rank(MPI_COMM_WORLD, & my_rank);
  //MPI_Comm_size(MPI_COMM_WORLD, & num_rank);

  // size the resulting vector correctly
  y.reset(x.size());

  // if we have already built the matrix, we just multiply it 
  if ( _build_lm_matrix ) {
    
      // if the approximate flag is set to be true, then throw out a exception
      if ( approximate ) 
        throw formic::Exception("Davidson solver's matrix-vector multiplication routine doesn't support approximate matrix!"); 

      // call level-2 blas function
      formic::dgemv('N', _hmat.rows(), _hmat.cols(), 1.0, &_hmat.at(0, 0), _hmat.rows(), &x.at(0), 1, 0.0, &y.at(0), 1);
      //for (int i = 0; i < y.size(); i++) 
        //std::cout << boost::format("%10.8e ") % y(i);
      //std::cout << std::endl;

      return;

  }
 
  // if we do not have the matrix, then we need to do der_rat * le_der * x on each process
  // the matrix free implementation is a little bit complicated, I will explain it here 
  // we have two derivative vector matrox, L(i, j) = <i|H|[psi_j>/<i|psi>, D(i, j) = <i|psi_j>/<i|psi> 
  // i denotes configuration(electron position in real space and number vector in Hilbert space) 
  // psi_j denotes wfn derivative w.r.t. jth variable, and j=0 means undifferentiated wfn
  // note that |value/guiding|^2 and weights should be absorbed in L and D matrix
  // in ground state calculation, H = D^T * L, Hx = D^T * Lx
  // in excited state calculation, H = D^T * (omega * D - L), temp1 = omage * D * x
  // temp2 = Lx

  else {

    // the number of samples on each process 
    const int Ns = _le_der.rows();

    // the number of independent variables 
    const int Nind = _le_der.cols();

    // check whether derivative vector matrices have the same size
    if ( _le_der.rows() != _der_rat.rows() && _le_der.cols() != _der_rat.cols())
      throw formic::Exception("the input derivative vector matrices are of different size!");

    // if the approximate flag is set to be ture, then throw out an exception 
    if ( approximate )
      throw formic::Exception("Davidson solver's matrix-vector multiplication routine doesn't support approximate matrix!");

    // if we are doing ground state calculation
    if ( _ground ) {

      // temp vector to store le_der * x
      formic::ColVec<double> temp(_le_der.rows());

      // call blas level-2 function
      formic::dgemv('N', Ns, Nind, 1.0, &_le_der.at(0,0), Ns, &x.at(0), 1, 0.0, &temp.at(0), 1);

      // call blas level-2 function 
      formic::dgemv('T', Ns, Nind, 1.0, &_der_rat.at(0,0), Ns, &temp.at(0), 1, 0.0, &y.at(0), 1);
          
      //if (my_rank == 0) {
      //  for (int i = 0; i < y.size(); i++) 
      //    std::cout << boost::format("%10.8e ") % y(i);
      //  std::cout << std::endl;
      //}

      return;
    }
   
    else if ( !_ground ) {
     
      // temp vectpr to store omega * der_rat * x
      formic::ColVec<double> temp1(Ns);

      // temp vector to store le_der * x
      formic::ColVec<double> temp2(Ns);

      // call blas level-2 function 
      formic::dgemv('N', Ns, Nind, _hd_shift, &_der_rat.at(0,0), Ns, &x.at(0), 1, 0.0, &temp1.at(0), 1);

      // call blas level-2 function 
      formic::dgemv('N', Ns, Nind, 1.0, &_le_der.at(0,0), Ns, &x.at(0), 1, 0.0, &temp2.at(0), 1);

      // combine these two temp vector together 
      temp1 -= temp2;

      // left multiply by _der_rat^T
      formic::dgemv('T', Ns, Nind, 1.0, &_der_rat.at(0,0), Ns, &temp1.at(0), 1, 0.0, &y.at(0), 1);
      //for (int i = 0; i < y.size(); i++) 
        //std::cout << boost::format("%10.8e ") % y(i);
      //std::cout << std::endl;

      return;
    }
  }
}

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

void cqmc::engine::DavidsonLMHD::SMatVecOp(const formic::ColVec<double> & x, formic::ColVec<double> & y, const bool approximate)
{
 
  // size the resulting vector correctly 
  y.reset(x.size());

  // if we have already built the matrix, we just multiply it
  if ( _build_lm_matrix ) {

      // if the approximate flag is set to be true, then throw out a exception
      if ( approximate ) 
        throw formic::Exception("Davidson solver's matrix-vector multiplication routine doesn't support approximate matrix!"); 

      // call level-2 blas function
      formic::dgemv('N', _smat.rows(), _smat.cols(), 1.0, &_smat.at(0,0), _smat.rows(), &x.at(0), 1, 0.0, &y.at(0), 1);

      //for (int i = 0; i < y.size(); i++) 
        //std::cout << boost::format("%10.8e ") % y(i);
      //std::cout << std::endl;

      // return product vector 
      return;
    
  }
 
  // if we do not have the matrix, then we need to do der_rat * le_der * x on each process
  // the matrix free implementation is a little bit complicated, I will explain it here 
  // we have two derivative vector matrox, L(i, j) = <i|H|[psi_j>/<i|psi>, D(i, j) = <i|psi_j>/<i|psi> 
  // i denotes configuration(electron position in real space and number vector in Hilbert space) 
  // psi_j denotes wfn derivative w.r.t. jth variable, and j=0 means undifferentiated wfn
  // note that |value/guiding|^2 and weights should be absorbed in L and D matrix
  // in ground state calculation, H = D^T * D, Hx = D^T * Dx
  // in excited state calculation, H = (omega * D - L)^T * (omega * D - L), temp1 = omage * D * x
  // temp2 = Lx

  else {

    // the number of samples on each process 
    const int Ns = _le_der.rows();

    // the number of independent variables + 1
    const int Nind = _le_der.cols();

    // check to see whether derivative vector matrices have same size 
    if ( _le_der.rows() != _der_rat.rows() && _le_der.cols() != _der_rat.cols() )
      throw formic::Exception("input derivative vectors are of different size");

    // if the approximate flag is set to be ture, then throw out an exception 
    if ( approximate )
      throw formic::Exception("Davidson solver's matrix-vector multiplication routine doesn't support approximate matrix!");

    // if we are doing ground state calculation
    if ( _ground ) {

      // temp vector to store der_rat * x
      formic::ColVec<double> temp(Ns);

      // call blas level-2 function
      formic::dgemv('N', Ns, Nind, 1.0, &_der_rat.at(0,0), Ns, &x.at(0), 1, 0.0, &temp.at(0), 1);

      // call blas level-2 function 
      formic::dgemv('T', Ns, Nind, 1.0, &_der_rat.at(0,0), Ns, &temp.at(0), 1, 0.0, &y.at(0), 1);

      //for (int i = 0; i < y.size(); i++) 
        //std::cout << boost::format("%10.8e ") % y(i);
      //std::cout << std::endl;

      return;
    }
   
    else if ( !_ground ) {
     
      // temp vectpr to store omega * der_rat * x
      formic::ColVec<double> temp1(Ns);

      // temp vector to store le_der * x
      formic::ColVec<double> temp2(Ns);

      // temp vector that store omega * der_rat^T * (omega * der_rat - le_der) * x
      formic::ColVec<double> temp3(x.size());

      // call blas level-2 function 
      formic::dgemv('N', Ns, Nind, _hd_shift, &_der_rat.at(0,0), Ns, &x.at(0), 1, 0.0, &temp1.at(0), 1);

      // call blas level-2 function 
      formic::dgemv('N', Ns, Nind, 1.0, &_le_der.at(0,0), Ns, &x.at(0), 1, 0.0, &temp2.at(0), 1);

      // combine these two temp vector together 
      temp1 -= temp2;

      // omega * D^T * (omega * D - L) * x  
      formic::dgemv('T', Ns, Nind, _hd_shift, &_der_rat.at(0,0), Ns, &temp1.at(0), 1, 0.0, &y.at(0), 1);

      // L^T * (omega * D - L) * x
      formic::dgemv('T', Ns, Nind, 1.0, &_le_der.at(0, 0), Ns, &temp1.at(0), 1, 0.0, &temp3.at(0), 1);

      // (omega * D^T - L^T) * (omega * D - L) * x
      y -= temp3;
      //for (int i = 0; i < y.size(); i++) 
        //std::cout << boost::format("%10.8e ") % y(i);
      //std::cout << std::endl;

      return;
    }
  }
} 

////////////////////////////////////////////////////////////////////////////////////
// \brief function that performs preconditioner matrix-vector multiplication 
//
// \param[in]   x              input vector 
// \param[in]   matrix_built   whether we have already built the matrix or not 
// \param[in]   approximate    whether to use approximated overlap or not 
// \param[out]  y              result vector
// NOTE: Unlike hamiltonian matrix-vector multiplication function, no transpose flag
//       in this function because overlap matrix is assumed to be symmetric
//
////////////////////////////////////////////////////////////////////////////////////

//void cqmc::engine::DavidsonLMHD::MMatVecOp(const formic::ColVec<double> & x, formic::ColVec<double> & y, const int mmat_first)
//{
// 
//  // size the resulting vector correctly 
//  y.reset(x.size());
//
//  // if we have already built the matrix, we just multiply it
//  if ( _build_lm_matrix ) {
//
//      formic::ColVec<double> z;
//      z.reset(x.size());
//
//      if ( mmat_first ) {
//        formic::dgemv('N', _mmat.rows(), _mmat.cols(), 1.0, &_mmat.at(0,0), _mmat.rows(), &x.at(0), 1, 0.0, &z.at(0), 1);
//        for (int i = 0; i < x.size(); i++)
//          y.at(i) = z.at(i) * _ovl_diag_neg_half[i];
//          //y(i) = z(i) / std::sqrt(std::max(1.0e-10, _smat(i,i)));
//      } else {
//        for (int i = 0; i < x.size(); i++)
//          z.at(i) = x.at(i) * _ovl_diag_neg_half[i];
//          //z(i) = x(i) / std::sqrt(std::max(1.0e-10, _smat(i,i)));
//        formic::dgemv('N', _mmat.rows(), _mmat.cols(), 1.0, &_mmat.at(0,0), _mmat.rows(), &z.at(0), 1, 0.0, &y.at(0), 1);
//      }
//
//      //for (int i = 0; i < x.size(); i++)
//      //  y(i) = x(i) / std::sqrt(std::max(1.0e-10, _smat(i,i)));
//
//      // return product vector 
//      return;
//    
//  // currently no preconditioning if we are not building the matrix
//  } else {
//    y = x;
//    return;
//  }
//
//}


//////////////////////////////////////////////////////////////////////////////////////////
// \brief apply linear method identity shift to the diagonal element of the Hamiltonian
//
//
//
//////////////////////////////////////////////////////////////////////////////////////////

//void cqmc::engine::DavidsonLMHD::apply_lm_shift()
//{
//  // get rank number and number of ranks 
//  int my_rank; 
//  int num_rank;
//  MPI_Comm_rank(MPI_COMM_WORLD, & my_rank);
//  MPI_Comm_size(MPI_COMM_WORLD, & num_rank);

//  // apply identity shift in first derivative space
//  if (my_rank == 0) {
//    for (int i = 1; i < _nfds; i++){
//      _hmat(i, i) = _hmat(i, i) + _hshift_i;
//    }
//  }
//}


/////////////////////////////////////////////////////////////////////////////////////
// \brief changes the Linear Method Hamiltonian shift
//
// 
//
//////////////////////////////////////////////////////////////////////////////////////
		
//void cqmc::engine::DavidsonLMHD::update_lm_shift(const double new_i_shift, const double new_s_shift)
//{
//  //remember this new shift
//  _hshift_i = new_i_shift;
//  _hshift_s = new_s_shift;
//			
//}

////////////////////////////////////////////////////////////////////////////////////
// \brief updates hamiltonian * krylov vector and hamiltonian projection based on 
//        new shift 
//
//
////////////////////////////////////////////////////////////////////////////////////

void cqmc::engine::DavidsonLMHD::update_hvecs_sub(const double new_i_shift, const double new_s_shift)
{

  // get rank number and number of ranks 
  int my_rank = formic::mpi::rank();
  //int num_rank;
  //MPI_Comm_rank(MPI_COMM_WORLD, & my_rank);
  //MPI_Comm_size(MPI_COMM_WORLD, & num_rank);

  // get the different between new shift and old shift 
  const double diff_shift_i = new_i_shift - _hshift_i;
  const double diff_shift_s = new_s_shift - _hshift_s;

  if (my_rank == 0) { 
    // update "identity shift" for the hamiltonian product
    for (int j = 0; j < _nkry; j ++) {
      for (int i = 1; i < _nfds; i ++) {
        _hvecs.at(i,j) += diff_shift_i * _kvecs.at(i,j);
      }
    }

    // update "overlap shift" for the hamiltonian product 
    for (int j = 1; j < _nkry; j++) {
      for (int i = 0; i < _nfds; i++) {
        _hvecs.at(i,j) += diff_shift_s * _svecs.at(i,j);
      }
    }

    // update projection of hamiltonian matrix 
    _subH = _kvecs.t() * _hvecs;
  }
}
  

////////////////////////////////////////////////////////////////////////////////////
// \brief update Hamiltonian and overlap matrix for the solver 
//
//
//
////////////////////////////////////////////////////////////////////////////////////

void cqmc::engine::DavidsonLMHD::update_hamovlp(formic::Matrix<double> & hmat, formic::Matrix<double> & smat)
{
	
  // get rank number and number of ranks 
  int my_rank = formic::mpi::rank(); 
  //int num_rank;
  //MPI_Comm_rank(MPI_COMM_WORLD, & my_rank);
  //MPI_Comm_size(MPI_COMM_WORLD, & num_rank);

  // check if the input matrices have the right size 
  bool right_size = (hmat.rows() == _nfds && hmat.cols() == _nfds && smat.rows() == _nfds && smat.cols()); 

  if (my_rank == 0 && !right_size)
    throw formic::Exception("bad matrix size found in update hamiltonian and overlap function");

  else {
    _hmat = hmat;
    _smat = smat;
  }


}

////////////////////////////////////////////////////////////////////////////////////
// \brief reset the eigen solver 
// 
// \brief clear subspace projection of Hamiltonian and overlap matrix, clear Krylov
//        subspace and action of Hamiltonian and overlap matrix
////////////////////////////////////////////////////////////////////////////////////

void cqmc::engine::DavidsonLMHD::child_reset()
{
  // clear subspace projection of Hamiltonian and overlap matrix
  _subH.reset(0, 0);
  _subS.reset(0, 0);

  // clear Krylov subspace 
  _nkry = 0;
  _kvecs.reset(0, 0);

  // clear Hamiltonian and overlap matrix's action on krylov subspace 
  _hvecs.reset(0, 0);
  _svecs.reset(0, 0);
  _htvecs.reset(0, 0);
			
  // clear eigenvectors and wavefunction coefficients
  _sub_evec.reset(0);

  // clear all values calculated from last solve
  _sub_eval = 0.0;
  _smallest_sin_value = 1e10;

  // set the eigenvector to initial guess
  _evecs.reset( ( _var_deps_use ? 1 + _dep_ptr->n_tot() : _nfds ), 0.0 );
  _evecs.at(0) = 1.0;

  // set the independent eigenvector to initial guess
  _ind_evecs.reset( ( _var_deps_use ? 1 + _dep_ptr->n_ind() : _nfds), 0.0 );
  _ind_evecs.at(0) = 1.0;
}


/////////////////////////////////////////////////////////////////////////////////////
// \brief solves the eigenvalue problem via the normal davidson method
//
//
//
/////////////////////////////////////////////////////////////////////////////////////
		
bool cqmc::engine::DavidsonLMHD::iterative_solve(double & eval, std::ostream & output)
{
	
	
  // get rank number and number of ranks 
  int my_rank = formic::mpi::rank(); 
  //int num_rank;
  //MPI_Comm_rank(MPI_COMM_WORLD, & my_rank);
  //MPI_Comm_size(MPI_COMM_WORLD, & num_rank);

  // initialize the solution vector to the unit vector along the first direction
  _evecs.reset( ( _var_deps_use ? 1 + _dep_ptr->n_tot() : _nfds ), 0.0 );
  _evecs.at(0) = 1.0;

  // ensure that at least one vector is in the krylov subspace 
  if( my_rank == 0 && _nkry == 0)
    throw formic::Exception("Empty krylov subspace upon entry to iterative_solve. Did you forget to add the initial vector?");
			
  // return value, whether the solver is successful or not
  bool retval = true;

  // converge flag
  bool converged = false;

  // best residual
  _best_residual = 1.0e100;
			
  // times of iteration 
  int iter = 0;

  // print out that we have started the iteration
  if (my_rank == 0) 
    output << boost::format("iteration solving starts here(engine) \n") << std::endl << std::endl;

  while(true){

    // average of smallest singular value 
    double smallest_sin_val_avg = 0.0;

    // solve subspace eigenvalue problem on root process
    if (my_rank == 0) {
      this -> solve_subspace(true);
    }

    // send resulting eigenvalues to all processes and record it as the new best estimate(if we build the matrix)
    formic::mpi::bcast(&_cost, 1);
    eval = _cost;

    // check if the cost function has imaginary component and stop iterating when it does
    formic::mpi::bcast(&_eval_was_complex, 1);
    if (_eval_was_complex){
      if(my_rank == 0)
        output << boost::format("davidson iteration %4i: stopping due to imaginary component in cost function") % iter << std::endl;
      break;
    }

    // if cost function change is unreasonable, stop iterating and set bad solve flag
    if (std::abs(_cost - _init_cost) > _max_energy_change){
      retval = false; 
      if( my_rank == 0 )
        output << boost::format("davidson iteration %4i stopping due to too large eigenvalue change") % iter << std::endl;
      break;
    }

    // if the overlap matrix becomes singular, stop iterating 
    formic::mpi::bcast(&_smallest_sin_value, 1);
    if (std::abs(_smallest_sin_value) < _singular_value_threshold) {
      if (my_rank == 0)
        output << boost::format("davidson iteration %4i stopping due to small subspace S singular value of %.2e") % iter % _smallest_sin_value << std::endl;
      break;
    }

    // construct new krylov vectors from subspace eigenvector
    if (my_rank == 0){
      _wv6.reset(_nfds);
      _wv6 = _kvecs * _sub_evec;
      // get the normalization factor
      const double temp_norm = 1.0 / _wv6.norm2();

      // normalize this new vector 
      _wv6 /= temp_norm;

      // add up linear combination of Hamiltonian and overlap products to make the new residual vector 
      _wv1.reset(_nfds);
      _wv1 = _hvecs * _sub_evec;
      _wv1 = _wv1 - _sub_eval * _svecs * _sub_evec;
    }

    // if we don't build matrix, then send this new vector to all processes 
    if ( !_build_lm_matrix ) 
      formic::mpi::bcast(&_wv1.at(0), _wv1.size());

    // compute the residual norm and send it to all processes
    double current_residual;
    double current_residual_avg;

    current_residual = _wv1.norm2();
    formic::mpi::bcast(&current_residual, 1);

    // if this is the best residual, save it and save the new eigenvector estimate
    if (my_rank == 0 && current_residual < _best_residual){
      _best_residual = current_residual; 

      // get current eigenvector estimate, which corresponds to the set of independent variables
      //Eigen::VectorXd ind_evecs;
      //this -> MMatVecOp(_wv6, _ind_evecs, 1); // preconditioning
      _ind_evecs = _wv6;

      // if our actual variables are dependent on the set of independent variables worked with here, expand the eigenvector into the full set of variables
      if ( _var_deps_use ) {

        // size the eigenvector correctly
        _evecs.reset(1 + _dep_ptr -> n_tot(), 0.0);
        _evecs.at(0) = _ind_evecs.at(0);

        // get some temporary vectors
        formic::ColVec<double> _evec_temp(_dep_ptr -> n_tot());
        formic::ColVec<double> _ind_temp(_dep_ptr -> n_ind());
        for ( int i = 0; i < _ind_temp.size(); i++ ) {
          _ind_temp.at(i) = _ind_evecs.at(i+1);
        }
        _dep_ptr -> expand_ind_to_all(&_ind_temp.at(0), &_evec_temp.at(0));

        for ( int i = 0; i < _evec_temp.size(); i++ ) {
          _evecs.at(i+1) = _evec_temp.at(i);
        }
      }
      
      // otherwise just copy the eigenvector into output since the independent and total variable sets are the same
      else {
        _evecs = _ind_evecs.clone();
      }
    }

    // print itertion result
    if (my_rank == 0) {
      
      // if we are doing ground state calculation, then print out energy
      if ( _ground ) 
        output << boost::format("davidson iteration %4i:   krylov dim = %3i   energy = %20.12f       residual = %.2e           smallest_sin_value = %.2e")
        % iter
        % _nkry
        % _cost
        % current_residual
        % _smallest_sin_value
        << std::endl;

      // if we are doing excited state calculation, then print out target function value 
      else 
        output << boost::format("davidson iteration %4i:   krylov dim = %3i   tar_fn = %20.12f       residual = %.2e           smallest_sin_value = %.2e")
        % iter
        % _nkry
        % _cost
        % current_residual
        % _smallest_sin_value
        << std::endl;
    }

    // check for convergence
    converged = current_residual < _residual_threshold;

    // if we build matrix and iteration has already converged, we exit iteration 
    if ( converged ) 
      break;

    // if iteration hasn't converged, we increment the iteration count and stop if maximum number of iterations has been reached
    if( iter ++ >= _n_max_iter) 
      break;

    // add the new krylov basis vector 
    this -> add_krylov_vector(_wv1);


  }

  // print iteration information
  if (converged && my_rank == 0)
    output << boost::format("davidson solver converged in %10i iterations") % iter << std::endl << std::endl;

  else if (my_rank == 0)
    output << boost::format("davidson solver did not converge after %10i iterations") % iter << std::endl << std::endl;

  return retval;

}


///////////////////////////////////////////////////////////////////////////////////////////////
// \brief solves the eigenvalue problem
//
//
//
///////////////////////////////////////////////////////////////////////////////////////////////

bool cqmc::engine::DavidsonLMHD::solve(double & eval, std::ostream & output)
{
  //if (_spam_use)
  //	return this -> iterative_solve_outer(eval);

  //if (! _spam_use)

  // get rank number and number of ranks 
  int my_rank = formic::mpi::rank(); 
  //int num_rank;
  //MPI_Comm_rank(MPI_COMM_WORLD, & my_rank);
  //MPI_Comm_size(MPI_COMM_WORLD, & num_rank);

  if ( !_variance_correct ) 
    return this -> iterative_solve(eval, output);

  else {
    
    // first evaluate tau and correct hamiltonian matrix 
    this -> tau_and_correct_ham();

    // bcast tau 
    formic::mpi::bcast(&_tau, 1);
    
    // flag to tell whether tau is converged
    bool tau_converged = false;

    // number of iteration of tau
    int iter_tau = 0;

    // return value
    bool retval = false;

    // reset eigensolver
    this -> child_reset();

    // set the current energy as initial wavefunction's energy
    _cost = _init_cost;
        
    // add initial guess
    // get the dimension of matrix 
    const int N = _hmat.cols();
    formic::ColVec<double> temp(N);
    for (int j = 0; j < temp.size(); j++) 
      temp.at(j) = (j == 0 ? 1.0 : 0.0);
    this -> add_krylov_vector(temp);

    // iteratively solve the generalized eigenvalue problem
    retval = this -> iterative_solve(eval, output);

    // return solve result
    return retval;

    //while ( true ) {

    //  // print iteration information
    //  if ( my_rank == 0 ) 
    //    output << boost::format("tau iteration %4i: tau = %20.12f") % iter_tau % _tau << std::endl;     

    //  // if tau has converged, stop iteration
    //  //if ( tau_converged ) 
    //  //  break;

    //  // if the maximum number of allowed iteration has reached, stop iteration
    //  if ( iter_tau++ > 0 ) 
    //    break;

    //  // if not converged and it is not the first iteration, reset the eigensolver and add initial guess
    //  if ( !tau_converged ) {
    //    
    //    // reset eigensolver
    //    this -> child_reset();

    //    // set the current energy as initial wavefunction's energy
    //    _energy = _init_energy;
    //    
    //    // add initial guess
    //    // get the dimension of matrix 
    //    const int N = _hmat.cols();
    //    Eigen::VectorXd temp(N);
    //    for (int j = 0; j < temp.size(); j++) 
    //      temp(j) = (j == 0 ? 1.0 : 0.0);
    //    this -> add_krylov_vector(temp);
    //  }

    //  // print iteration information
    //  if ( my_rank == 0 ) 
    //    output << boost::format("tau iteration %4i: tau = %20.12f") % iter_tau % _tau << std::endl;

    //  // iteratively solve the generalized eigenvalue problem
    //  retval = this -> iterative_solve(eval, output);

    //  // re-evaluate tau and correct hamiltonian based on new tau
    //  this -> tau_and_correct_ham();
    //  MPI_Bcast(&_tau, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    //  // check if converged
    //  tau_converged = (_tau_diff < 1.0e-3);
    //  MPI_Bcast(&tau_converged, 1, MPI::BOOL, 0, MPI_COMM_WORLD);

    //}

    //// print iteration information
    //if ( tau_converged && my_rank == 0) 
    //  output << boost::format("tau iteration converged in %10i iterations") % iter_tau << std::endl; 

    //else if (my_rank == 0) 
    //  output << boost::format("tau iteration did not converge after %10i iterations") % iter_tau << std::endl;
    //
    //// return iteration results
    //return retval;
  }
}


///////////////////////////////////////////////////////////////////////////////////////////////
// \brief converts eigenvectors into wave function coefficients
//        solving this question Hc = Sd for d, S is the ground overlap matrix 
//
//
///////////////////////////////////////////////////////////////////////////////////////////////

void cqmc::engine::DavidsonLMHD::child_convert_to_wf_coeff()
{
					
  return;

}

//////////////////////////////////////////////////////////////////////////////////////////////
// \brief return the wave function coefficients
//
//
//
//////////////////////////////////////////////////////////////////////////////////////////////

//Eigen::VectorXd cqmc::engine::DavidsonLMHD::wf_coeff()
//{
//  return _wf_coefficients;
//}

