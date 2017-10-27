//////////////////////////////////////////////////////////////////////////////////////////////////
// \brief a class to perform solve generalized eigenvalue problem by SPAM 
//
//
//
//  solve with spam method
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
#include<formic/utils/lapack_interface.h>
#include<formic/utils/mpi_interface.h>
#include<formic/utils/lmyengine/eigen_solver.h>
#include<formic/utils/lmyengine/spam_solver.h>

//////////////////////////////////////////////////////////////////////////////////////////////////
// \brief solve the subspace generalized eigenvalue problem 
//        with nonsymmetric H and symmetric S for outer loop
//
//
//////////////////////////////////////////////////////////////////////////////////////////////////

void cqmc::engine::SpamLMHD::solve_subspace_nonsymmetric(const bool outer) 
{
  
  // get rank number and number of rank
  int my_rank = formic::mpi::rank();
  int num_rank = formic::mpi::rank();
  //MPI_Comm_rank(MPI_COMM_WORLD, & my_rank);
  //MPI_Comm_size(MPI_COMM_WORLD, & num_rank);

  // one in complex form
  const std::complex<double> complex_one(1.0, 0.0);
  const std::complex<double> complex_zero(0.0, 0.0);

  //print subspace matrix to debug 
  //if (my_rank == 0) {
  //  std::cout << boost::format("subspace S matrix is for rank %d \n") % my_rank << std::endl << std::endl; 
  //  for (int i = 0; i < _subS.rows(); i++) {
  //    for (int j = 0; j < _subS.cols(); j++) {
  //      std::cout << boost::format("%10.8e ") % _subS(i, j);
  //    }
  //    std::cout << boost::format("\n");
  //  }
  //}

  ////print subspace matrix to debug 
  //if (my_rank == 0) {
  //  std::cout << boost::format("subspace H matrix is for rank %d \n") % my_rank << std::endl << std::endl; 
  //  for (int i = 0; i < _subH.rows(); i++) {
  //    for (int j = 0; j < _subH.cols(); j++) {
  //      std::cout << boost::format("%10.8e ") % _subH(i, j);
  //    }
  //    std::cout << boost::format("\n");
  //  }
  //}
  // get subspace dimension 
  const int m  = _nkry;

  // create vectors and matrices used in svd routine
  formic::Matrix<double> u;
  formic::Matrix<double> vt;
  formic::ColVec<double> sin_vals;
  int truncate_index;

  // make sure the subspace matrix is not empty
  if ( outer ) {
    if ( _subS.rows() == 0 || _subH.rows() == 0 ) 
      throw formic::Exception("subspace matrix is empty upon solving subspace eigenvalue problem(outer)");
  }

  else {
    if ( _hy_subS.rows() == 0 || _hy_subH.rows() == 0 )
      throw formic::Exception("subspace matrix is empty upon solving subspace eigenvalue problem(inner)");
  }

  // perform svd to subspace overlap matrix 
  if ( outer )
    _subS.svd(u, sin_vals, vt);
  else 
    _hy_subS.svd(u, sin_vals, vt);
  formic::Matrix<double> v = vt.t();

  // record the smallest singular value of subspace overlap matrix 
  if ( outer ) 
    _smallest_sin_value_outer = std::abs(sin_vals.at(0));
  else 
    _smallest_sin_value_inner = std::abs(sin_vals.at(0));
  for (int i = 0; i < m; i++) {
    if ( outer )
      _smallest_sin_value_outer = std::min(_smallest_sin_value_outer, std::abs(sin_vals.at(i)));
    else 
      _smallest_sin_value_inner = std::min(_smallest_sin_value_inner, std::abs(sin_vals.at(i)));
    truncate_index = i;
    
    // check if singular value is smaller than singular value threshold
    if ( outer && _smallest_sin_value_outer < _singular_value_threshold) 
      break;
    else if ( !outer && _smallest_sin_value_inner < _singular_value_threshold)
      break;
  }

  // get the number of colums of new U and V matrix by add 1 to truncate index
  truncate_index ++;

  // throw away those columns in U and V matrix which corresponds to singular values below the threshold
  u.conservativeResize(m, truncate_index);
  v.conservativeResize(m, truncate_index);

  // throw away those small singular values
  sin_vals.conservativeResize(truncate_index);

  // convert the truncated singular value vector to a truncated_index * truncated_index diagonal matrix
  formic::Matrix<double> trun_sin_val_matrix(truncate_index, truncate_index, 0.0);
  for (int i = 0; i < truncate_index; i++)
    trun_sin_val_matrix.at(i,i) = sin_vals.at(i);

  // calculate the inverse of this matrix
  for(int i = 0; i < truncate_index; i++){
    for(int j = 0; j < truncate_index; j++){
      trun_sin_val_matrix.at(i,j) = (trun_sin_val_matrix.at(i,j) == 0 ? trun_sin_val_matrix.at(i,j) : 1 / trun_sin_val_matrix.at(i,j));
    }
  }

  // calculate matrix S_trun^-1 * U^-1 * H * V
  formic::Matrix<double> new_sub_H = trun_sin_val_matrix * u.t() * (outer ? _subH : _hy_subH) * v;

  // set a vector to hold eigenvalues
  formic::ColVec<std::complex<double> > e_evals;

  // set an eigen array to hold energies
  formic::ColVec<std::complex<double> > energy_list;

  // a matrix to hold eigenvectors
  formic::Matrix<std::complex<double> > evecs_list;

  // solve this standard eigenvalue problem ( new_H * y = lambda * y, y = V^T * x, x is the original eigenvector)
  new_sub_H.nonsym_eig(e_evals, evecs_list);

  // if we are doing excited state calculations, convert the resulting eigenvalues to energy( this is important in harmonic davidson)
  if ( !_ground ) {
    energy_list = _hd_shift * complex_one - complex_one / e_evals;
  }

  // if we are doing ground state calculation, then do nothing
  if ( _ground ) {
    energy_list = e_evals.clone();
  }

  int selected = 0;
  // if we want to chase the closest, selected the eigenvalue that is real and most similar to the previous eigenvalue( this is essestially an attempt to stay in the same solution)
  if ( _chase_closest ) { // currently turn on this if statement for spam 
    std::complex<double> closest_energy = energy_list.at(0);
    std::complex<double> corrs_eval = e_evals.at(0);
    for (int j = 1; j < truncate_index; j++){
      // first check if this energy is real 
      if ( std::abs((energy_list.at(j)).imag()) < 1.0e-6 ) {  
        
        // then select the energy that is closest to the previous one
        if(std::abs( complex_one * (outer ? _energy_outer : _energy_inner) - energy_list.at(j) ) < std::abs(complex_one * (outer ? _energy_outer : _energy_inner) - closest_energy )){
          selected = j;
          closest_energy = energy_list.at(j);
          corrs_eval = e_evals.at(j);
        }
      }
    }

    // if the eigenvalue has an imaginary component, we abort 
    _eval_was_complex = false;
    if( std::abs(closest_energy.imag()) > 1.0e-6 ) {
      _eval_was_complex = true;
      return;
    }

    // if the eigenvalue is real, we record it and the corresponding eigenvector 
    if ( outer ) {
      // record energy 
      _energy_outer = closest_energy.real();

      // record the eigenvalue
      _sub_eval_outer = corrs_eval.real();
    }

    else {
      // record energy
      _energy_inner = closest_energy.real();

      // record the eigenvalue
      _sub_eval_inner = corrs_eval.real();
    }

    // record the eigenvector y
    _wv4 = evecs_list.col_as_vec(selected);

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
    //_wv5 = V * _wv4;

    // take real part of the vector x, put that into eigenvector
    (outer ? _sub_evec_outer : _sub_evec_inner).reset(m);
    for (int i = 0; i < m; i++) {
      (outer ? _sub_evec_outer : _sub_evec_inner)(i) = _wv5.at(i).real();
    }
  }

  // if we want to chase the lowest, select the lowest eigenvalue(currently turn off this if statement for spam)
  if ( _chase_lowest ) {
    
    // the vector that stores all lower-than-current energy(target function)
    std::vector<std::complex<double> > lower_than_current_list;
    std::vector<int> lower_than_current_index;
    //Eigen::ArrayXcd eval_list = (es.eigenvalues()).array();
    std::complex<double> lowest_eval = e_evals.at(0);

    double inner_eval = 0.0;
    if ( !_ground ) 
      inner_eval = 1.0 / (_hd_shift - _energy_inner);
    else 
      inner_eval = _energy_inner;

    // if it's outer iteration, we just make sure that we choose the lowest energy(target function)
    //if ( outer ) {
    for (int j = 1; j < truncate_index; j++) {
      if (e_evals.at(j).real() < lowest_eval.real()) {
        selected = j;
        lowest_eval = e_evals.at(j);
      }
    }
    //}

    // if it's inner iteration
    //else {
    //  
    //  // we first need to get all energy(target function) lower than current
    //  for (int j = 0; j < truncate_index; j++) {
    //    if (eval_list(j).real() < inner_eval) {
    //      lower_than_current_list.push_back(eval_list(j));
    //      lower_than_current_index.push_back(j);
    //    }
    //  }

    //  lowest_eval = lower_than_current_list.at(0);
    //  selected = lower_than_current_index.at(0);
    //  // then we select from lower list the closest energy(target function)
    //  for (int i = 1; i < lower_than_current_list.size(); i++) {
    //    
    //    // first check if this energy is real
    //    if ( std::abs((lower_than_current_list.at(i)).imag()) < 1.0e-6 ) {

    //      // then select energy
    //      if ( std::abs( complex_one * inner_eval - lower_than_current_list.at(i) ) < std::abs(complex_one * inner_eval - lowest_eval )) {
    //        selected = lower_than_current_index.at(i);
    //        lowest_eval = lower_than_current_list.at(i);
    //      }
    //    }
    //  }
    //}
          

    // if the eigenvalue has an imaginary component, we abort
    _eval_was_complex = false;
    if ( std::abs(lowest_eval.imag()) > 1.0e-6 ) {
      _eval_was_complex = true;
      return;
    }

    // if the eigenvalue is real, we record it 
    if ( !_ground ) {
      if ( outer ) 
        _energy_outer = _hd_shift - 1.0 / lowest_eval.real();
      else 
        _energy_inner = _hd_shift - 1.0 / lowest_eval.real();
    }

    else {
      if ( outer )
        _energy_outer = lowest_eval.real();
      else 
        _energy_inner = lowest_eval.real();
    }

    // record the eigenvalue
    (outer ? _sub_eval_outer : _sub_eval_inner) = lowest_eval.real();

    // record the eigenvector y
    _wv4 = evecs_list.col_as_vec(selected);

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
    //_wv5 = V * _wv4;

    // take real part of vector x, put that into eigenvector
    (outer ? _sub_evec_outer : _sub_evec_inner).reset(m);
    for (int i = 0; i < m; i++) {
      (outer ? _sub_evec_outer : _sub_evec_inner)(i) = _wv5.at(i).real();
    }
  }
}

///////////////////////////////////////////////////////////////////////////////
// \brief solves the subspace eigenvalue problem
//
//
//
///////////////////////////////////////////////////////////////////////////////

void cqmc::engine::SpamLMHD::solve_subspace(const bool outer) 
{
  this -> solve_subspace_nonsymmetric(outer);
}

/////////////////////////////////////////////////////////////////////////////////
// \brief adds a new krylov vector for inner loop of spam
//
//
//
//
/////////////////////////////////////////////////////////////////////////////////

void cqmc::engine::SpamLMHD::add_krylov_vector_inner(const formic::ColVec<double> & v)
{
  
  // get rank number and number of rank
  int my_rank = formic::mpi::rank();
  int num_rank = formic::mpi::size();
  //MPI_Comm_rank(MPI_COMM_WORLD, & my_rank);
  //MPI_Comm_size(MPI_COMM_WORLD, & num_rank);

  // check vector length
  if (my_rank == 0 && v.size() != _der_rat.cols())
    throw formic::Exception("bad vector length of %d in SpamLMHD::add_krylov_vector_outer: expected length of %d") % v.size() % _der_rat.cols();

  // increment krylov subspace size and remember the old size
  const int nold = _nkry++;

  _wv1.reset(v.size());
  _wv1 = v;

  // perform gram-schmidt orthogonalization against the existing krylov vectors
  if (my_rank == 0) {
    for (int i = 0; i < nold; i++) {
      _wv1 -= dotc(_kvecs.col_as_vec(i), _wv1) * _kvecs.col_as_vec(i);
    }
  }

  // broadcast this vector to all processes
  formic::mpi::bcast(&_wv1.at(0), _wv1.size());

  // compute the product of approximate hamiltonian times the new krylov vector
  formic::ColVec<double> hs(_nfds);
  this -> HMatVecOp(_wv1, hs, false, true);
  formic::ColVec<double> hs_avg(hs.size());
  formic::mpi::reduce(&hs.at(0), &hs_avg.at(0), hs.size(), MPI_SUM);
  hs = hs_avg.clone();

  // compute the product of approximate overlap matrix times this new krylov vector
  this -> SMatVecOp(_wv1, _wv2, true);
  formic::ColVec<double> _wv2_avg(_wv2.size());
  formic::mpi::reduce(&_wv2.at(0), &_wv2_avg.at(0), _wv2.size(), MPI_SUM);
  _wv2 = _wv2_avg.clone();


  // modify the hamiltonian product to account for "identity shift"
  if (my_rank == 0) {
    for (int i = 1; i < hs.size(); i++) {
      hs.at(i) += _hshift_i * _wv1.at(i);
    }
  }

  // modify hamiltonian product to account for "overlap" shift 
  if (my_rank == 0 && nold > 0) {
    hs += _hshift_s * _wv2;
  }


  // normalize the new krylov vector and save the vector and its operation on overlap matrix
  if (my_rank == 0) {
    const double norm = std::sqrt(_wv1.norm2());
    _wv1 /= norm;
    _wv2 /= norm;
    hs /= norm;

    // krylov space
    _kvecs.conservativeResize(_nfds, _nkry);
    std::copy(_wv1.begin(), _wv1.end(), _kvecs.col_begin(_nkry-1));

    // S_hybrid * krylov space
    _hsvecs.conservativeResize(_nfds, _nkry);

    // temp vector that store (S * X)^T * xnew by blas level-2
    formic::ColVec<double> new_col_s_temp1(_nkry_full);
    formic::dgemv('T', _nfds, _nkry_full, 1.0, &_svecs.at(0,0), _nfds, &_wv1.at(0), 1, 0.0, &new_col_s_temp1.at(0), 1); 

    // tenp vector that store X * (S * X)^T * xnew + S(1) * xnew
    formic::ColVec<double> new_col_s_temp2(_nfds);
    formic::dgemv('N', _nfds, _nkry_full, 1.0, &_kvecs.at(0,0), _nfds, &new_col_s_temp1.at(0), 1, 0.0, &new_col_s_temp2.at(0), 1);
    //daxpy(_nfds, 1.0, &_wv2(0), 1, &new_col_s_temp2(0), 1);
    new_col_s_temp2 += _wv2;

    // temp vector that store X * S(1)*xnew
    formic::ColVec<double> new_col_s_temp3(_nkry_full);
    formic::dgemv('T', _nfds, _nkry_full, 1.0, &_kvecs.at(0,0), _nfds, &_wv2.at(0), 1, 0.0, &new_col_s_temp3.at(0), 1);

    formic::ColVec<double> new_col_s(_nfds);
    formic::dgemv('N', _nfds, _nkry_full, 1.0, &_kvecs.at(0,0), _nfds, &new_col_s_temp3.at(0), 1, 0.0, &new_col_s.at(0), 1);
    //daxpy(_nfds, -1.0, &new_col_s(0), 1, &new_col_s_temp2(0), 1);
    new_col_s = -1.0 * new_col_s + new_col_s_temp2;

    std::copy(new_col_s.begin(), new_col_s.end(), _hsvecs.col_begin(_nkry-1));

    // H_hybrid * krylov space
    _hhvecs.conservativeResize(_nfds, _nkry);

    // temp vector that store (H^T * X)^T * xnew
    formic::ColVec<double> new_col_h_temp1(_nkry_full);
    formic::dgemv('T', _nfds, _nkry_full, 1.0, &_thvecs.at(0,0), _nfds, &_wv1.at(0), 1, 0.0, &new_col_h_temp1.at(0), 1);

    // temp vector that store X * (H^T * X)^T * xnew + H(1)*xnew
    formic::ColVec<double> new_col_h_temp2(_nfds);
    formic::dgemv('N', _nfds, _nkry_full, 1.0, &_kvecs.at(0,0), _nfds, &new_col_h_temp1.at(0), 1, 0.0, &new_col_h_temp2.at(0), 1);
    //daxpy(_nfds, 1.0, &hs(0), 1, &new_col_h_temp2(0), 1);
    new_col_h_temp2 += hs;

    // temp vector that store X^T * H(1)*xnew
    formic::ColVec<double> new_col_h_temp3(_nkry_full);
    formic::dgemv('T', _nfds, _nkry_full, 1.0, &_kvecs.at(0,0), _nfds, &hs.at(0), 1, 0.0, &new_col_h_temp3.at(0), 1);

    formic::ColVec<double> new_col_h(_nfds);
    formic::dgemv('N', _nfds, _nkry_full, 1.0, &_kvecs.at(0,0), _nfds, &new_col_h_temp3.at(0), 1, 0.0, &new_col_h.at(0), 1);
    //daxpy(_nfds, -1.0, &new_col_h(0), 1, &new_col_h_temp2(0), 1);
    new_col_h = -1.0 * new_col_h + new_col_h_temp2;

    std::copy(new_col_h.begin(), new_col_h.end(), _hhvecs.col_begin(_nkry-1));

    // update subspace projection of hybrid Hamiltonian 
    _hy_subH.conservativeResize(_nkry, _nkry);
    _wv3 = _wv1.t() * _hhvecs;
    _wv6 = _kvecs.t() * new_col_h;
    std::copy(_wv3.begin(), _wv3.end(), _hy_subH.row_begin(_nkry-1));
    std::copy(_wv6.begin(), _wv6.end(), _hy_subH.col_begin(_nkry-1));
    //_hy_subH.bottomRows(1) = _wv3;
    //_hy_subH.rightCols(1) = _wv6;

    // update subspace projection of hybrid overlap
    _hy_subS.conservativeResize(_nkry, _nkry);
    _wv3 = _wv1.t() * _hsvecs;
    _wv6 = _kvecs.t() * new_col_s;

    std::copy(_wv3.begin(), _wv3.end(), _hy_subS.row_begin(_nkry-1));
    std::copy(_wv6.begin(), _wv6.end(), _hy_subS.col_begin(_nkry-1));
    //_hy_subS.bottomRows(1) = _wv3;
    //_hy_subS.rightCols(1) = _wv6;

    // add this vector to intermediate krylov space
    _kvecs_about_to_add.conservativeResize(_nfds, _nkry - _nkry_full);
     std::copy(_wv1.begin(), _wv1.end(), _kvecs_about_to_add.col_begin(_nkry-_nkry_full-1));   
    //_kvecs_about_to_add.rightCols(1) = _wv1;
  }


} 

////////////////////////////////////////////////////////////////////////////////////
// \brief adds a bunch of new Krylov vectors for spam outer loop
// 
// NOTE: This function assumes that the input vectors are already orthonormal!!
//
//
////////////////////////////////////////////////////////////////////////////////////
void cqmc::engine::SpamLMHD::add_krylov_vectors_outer(const formic::Matrix<double> & m)
{
 
  // get rank number and number of ranks
  int my_rank = formic::mpi::rank();
  int num_rank = formic::mpi::size();

  // check matrix size
  if (my_rank == 0 && m.rows() != _der_rat.cols())
    throw formic::Exception("bad matrix size of %d in SpamLMHD::add_krylov_vector_outer: expected length of %d") % m.rows() % _der_rat.cols();

  // get the number of new krylov vectors
  const int Nnew = m.cols();

  // remember the old size
  const int nold = _nkry;

  // if this is the first krylov vector, we increment the number of krylov vectors by the number of new krylov vectors
  if ( nold == 0 ) 
    _nkry += Nnew;

  // record the number of krylov vectors that have been multiplied by full hamiltonian and overlap, which is different from the total number of krylov vectors
  _nkry_full = _nkry;

  // put input matrix into work matrix 
  _wm1.reset(m.rows(), Nnew);
  _wm1 = m;

  // compute the product of Hamiltonian times these new krylov vectors
  formic::Matrix<double> hs(_nfds, Nnew);
  this -> HMatMatOp(_wm1, hs, false, false);
  formic::Matrix<double> hs_avg(_nfds, Nnew);
  formic::mpi::reduce(&hs.at(0,0), &hs_avg.at(0,0), hs.size(), MPI_SUM);
  hs = hs_avg.clone();

  // compute the product of Hamiltonian transpose times these new krylov vectors
  formic::Matrix<double> ths(_nfds, Nnew);
  this -> HMatMatOp(_wm1, ths, true, false);
  formic::Matrix<double> ths_avg(ths.rows(), ths.cols());
  formic::mpi::reduce(&ths.at(0,0), &ths_avg.at(0,0), ths.size(), MPI_SUM);
  ths = ths_avg.clone();

  // compute the product of the overlap matrix times these new krylov vectors
  this -> SMatMatOp(_wm1, _wm2, false);
  formic::Matrix<double> _wm2_avg(_wm2.rows(), _wm2.cols());
  formic::mpi::reduce(&_wm2.at(0,0), &_wm2_avg.at(0,0), _wm2.size(), MPI_SUM);
  _wm2 = _wm2_avg.clone();

  // modify hamiltonian product to account for "identity" shift
  if (my_rank == 0) {
    for (int i = 1; i < Nnew; i++) {
      for (int j = 0; j < hs.rows(); j++) {
        hs.at(j,i) += _hshift_i * _wm1.at(j,i);
        ths.at(j,i) += _hshift_i * _wm1.at(j,i);
      }
    }
  }

  // modify hamiltonian product to account for "overlap" shift 
  if (my_rank == 0 && nold > 0) {
    hs += _hshift_s * _wm2;
    ths += _hshift_s * _wm2;
  }



  // save these krylov vectors and their operation on matrix
  if ( my_rank == 0 ) {
    
    // krylov space 
    if ( nold == 0 ) {
      _kvecs.conservativeResize(_nfds, _nkry);
      for (int i = 0; i < Nnew; i++) 
        std::copy(_wm1.col_begin(i), _wm1.col_end(i), _kvecs.col_begin(_nkry-Nnew+i));
      //_kvecs.rightCols(Nnew) = _wm1;
    }

    // H * krylov space 
    _hvecs.conservativeResize(_nfds, _nkry);
    for (int i = 0; i < Nnew; i++) 
      std::copy(hs.col_begin(i), hs.col_end(i), _hvecs.col_begin(_nkry-Nnew+i));
    //_hvecs.rightCols(Nnew) = hs;
      
    // H_hybrid * krylov space 
    _hhvecs.reset(_nfds, _nkry);
    _hhvecs = _hvecs.clone();

    // H^T * krylov space 
    _thvecs.conservativeResize(_nfds, _nkry);
    for (int i = 0; i < Nnew; i++) 
      std::copy(ths.col_begin(i), ths.col_end(i), _thvecs.col_begin(_nkry-Nnew+i));
    //_thvecs.rightCols(Nnew) = ths;

    // S * krylov space 
    _svecs.conservativeResize(_nfds, _nkry);
    for (int i = 0; i < Nnew; i++) 
      std::copy(_wm2.col_begin(i), _wm2.col_end(i), _svecs.col_begin(_nkry-Nnew+i));
    //_svecs.rightCols(Nnew) = _wm2;

    // S_hybrid * krylov space 
    _hsvecs.reset(_nfds, _nkry);
    _hsvecs = _svecs.clone();

    // update subspace projection of hamiltonian
    _subH.conservativeResize(_nkry, _nkry);
    _wm3 = _wm1.t() * _hvecs;
    _wm4 = _kvecs.t() * hs;
    for (int i = 0; i < Nnew; i++) 
      std::copy(_wm4.col_begin(i), _wm4.col_end(i), _subH.col_begin(_nkry-Nnew+i));

    for (int i = 0; i < Nnew; i++) 
      std::copy(_wm3.col_begin(i), _wm3.col_end(i), _subH.row_begin(_nkry-Nnew+i));
    //_subH.bottomRows(Nnew) = _wm3;
    //_subH.rightCols(Nnew) = _wm4;

    // set the hybrid Hamiltonian subspace projection matrix same as the full hamiltonian
    _hy_subH.reset(_nkry, _nkry);
    _hy_subH = _subH.clone();

    // update subspace projection of the overlap 
    _subS.conservativeResize(_nkry, _nkry);
    _wm3 = _wm1.t() * _svecs;
    _wm4 = _kvecs.t() * _wm2;
    for (int i = 0; i < Nnew; i++) 
      std::copy(_wm4.col_begin(i), _wm4.col_end(i), _subS.col_begin(_nkry-Nnew+i));

    for (int i = 0; i < Nnew; i++) 
      std::copy(_wm3.col_begin(i), _wm3.col_end(i), _subS.row_begin(_nkry-Nnew+i));
    //_subS.bottomRows(Nnew) = _wm3;
    //_subS.rightCols(Nnew) = _wm4;

    // set the hybrid overlap subspace projection matrix same as the full overlap
    _hy_subS.reset(_nkry, _nkry);
    _hy_subS = _subS.clone();

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

void cqmc::engine::SpamLMHD::HMatVecOp(const formic::ColVec<double> & x, formic::ColVec<double> & y, const bool transpose, const bool approximate)
{


  int my_rank = formic::mpi::rank(); 
  int num_rank = formic::mpi::size();
  //MPI_Comm_rank(MPI_COMM_WORLD, & my_rank);
  //MPI_Comm_size(MPI_COMM_WORLD, & num_rank);

  // size the resulting vector correctly
  y.reset(x.size());

    
  // the number of samples on each process
  int Ns = _le_der.rows();

  // if we multiply by approximated matrix, we need to change Ns by approximate degree
  if ( approximate ) 
    Ns /= _appro_degree;

  // the number of independent variables 
  const int Nind = _le_der.cols();

  // check whether derivative vector matrices have the same size
  if ( _le_der.rows() != _der_rat.rows() && _le_der.cols() != _der_rat.cols())
    throw formic::Exception("the input derivative vector matrices are of different size!");

  // if we are doing ground state calculation
  if ( _ground && !approximate ) {

    // if we do H*x
    if ( !transpose ) {

      // temp vector to store le_der * x
      formic::ColVec<double> temp(Ns);

      // call blas level-2 function
      formic::dgemv('N', Ns, Nind, 1.0, &_le_der.at(0,0), Ns, &x.at(0), 1, 0.0, &temp.at(0), 1);

      // call blas level-2 function 
      formic::dgemv('T', Ns, Nind, 1.0, &_der_rat.at(0,0), Ns, &temp.at(0), 1, 0.0, &y.at(0), 1);
          
      return;
    }

    // if we do H^T * x
    else {
        
      // temp vector that store der_rat * x
      formic::ColVec<double> temp(Ns);

      // call blas level-2 function 
      formic::dgemv('N', Ns, Nind, 1.0, &_der_rat.at(0,0), Ns, &x.at(0), 1, 0.0, &temp.at(0), 1);

      // call blas level-2 function 
      formic::dgemv('T', Ns, Nind, 1.0, &_le_der.at(0,0), Ns, &temp.at(0), 1, 0.0, &y.at(0), 1);

      return;
    }
  }


  // if we are doing ground state calculation and multiply by approximate matrix 
  else if ( _ground && approximate ) {
      
    // if we want H * x
    if ( !transpose ) {
        
      // temp vector that stores le_der_appro * x
      formic::ColVec<double> temp(Ns);
      
      // call blas level-2 function 
      formic::dgemv('N', Ns, Nind, _appro_factor, &_le_der_appro.at(0,0), Ns, &x.at(0), 1, 0.0, &temp.at(0), 1);

      // call blas level-2 function 
      formic::dgemv('T', Ns, Nind, 1.0, &_der_rat_appro.at(0,0), Ns, &temp.at(0), 1, 0.0, &y.at(0), 1);

      return;
    }

    // if we want H^T * x
    else {
        
      // temp vector that stores der_rat * x
      formic::ColVec<double> temp(Ns);

      // call blas level-2 function 
      formic::dgemv('N', Ns, Nind, _appro_factor, &_der_rat_appro.at(0,0), Ns, &x.at(0), 1, 0.0, &temp.at(0), 1);

      // call blas level-2 function 
      formic::dgemv('T', Ns, Nind, 1.0, &_le_der_appro.at(0,0), Ns, &temp.at(0), 1, 0.0, &y.at(0), 1);

      return;
    }
  }

  // if we are doing excited state calculation and full matrix
  else if ( !_ground && !approximate ) {
      
    // if we want H*x
    if ( !transpose ) {
     
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

      return;
    }

    // if we want H^T * x
    else {
        
      // temp vector that store _der_rat * x
      formic::ColVec<double> temp1(Ns);

      // call blas level-2 function 
      formic::dgemv('N', Ns, Nind, 1.0, &_der_rat.at(0,0), Ns, &x.at(0), 1, 0.0, &temp1.at(0), 1);

      // temp vector that store _le_der^T * _der_rat * x
      formic::ColVec<double> temp2(Nind);

      // call blas level-2 function 
      formic::dgemv('T', Ns, Nind, 1.0, &_le_der.at(0,0), Ns, &temp1.at(0), 1, 0.0, &temp2.at(0), 1);

      // call bals level-2 function 
      formic::dgemv('T', Ns, Nind, _hd_shift, &_der_rat.at(0,0), Ns, &temp1.at(0), 1, 0.0, &y.at(0), 1);

      // get the resulting vector 
      y -= temp2;

      return;
    }
  }

  // if we are doing excited state calculation and approximate matrix 
  else if ( !_ground && approximate ) {
    
    // if we want H*x
    if ( !transpose ) {
      
      // temp vector that stores omega * der_rat * x
      formic::ColVec<double> temp1(Ns);

      // temp vector that stores le_der * x
      formic::ColVec<double> temp2(Ns);

      // call blas level-2 function 
      formic::dgemv('N', Ns, Nind, _hd_shift, &_der_rat_appro.at(0, 0), Ns, &x.at(0), 1, 0.0, &temp1.at(0), 1);

      // call blas level-2 function 
      formic::dgemv('N', Ns, Nind, 1.0, &_le_der_appro.at(0, 0), Ns, &x.at(0), 1, 0.0, &temp2.at(0), 1);

      // combine these two temp vector together 
      temp1 -= temp2;

      // left multiply by _der_rat_^T
      formic::dgemv('T', Ns, Nind, _appro_factor, &_der_rat_appro.at(0, 0), Ns, &temp1.at(0), 1, 0.0, &y.at(0), 1);
      
      return;
    }

    // if we want H^T * x
    else {
      
      // temp vector that store _der_rat * x
      formic::ColVec<double> temp1(Ns);

      // call blas level-2 function 
      formic::dgemv('N', Ns, Nind, _appro_factor, &_der_rat_appro.at(0, 0), Ns, &x.at(0), 1, 0.0, &temp1.at(0), 1);

      // temp vector that store _le_der^T * _der_rat * x
      formic::ColVec<double> temp2(Nind);

      // call blas level-2 function 
      formic::dgemv('T', Ns, Nind, 1.0, &_le_der_appro.at(0, 0), Ns, &temp1.at(0), 1, 0.0, &temp2.at(0), 1);

      // call bals level-2 function 
      formic::dgemv('T', Ns, Nind, _hd_shift, &_der_rat_appro.at(0, 0), Ns, &temp1.at(0), 1, 0.0, &y.at(0), 1);

      // get the resulting vector 
      y -= temp2;

      return;
    }
  }
}

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

void cqmc::engine::SpamLMHD::HMatMatOp(const formic::Matrix<double> & x, formic::Matrix<double> & y, const bool transpose, const bool approximate)
{
  
  int my_rank = formic::mpi::rank(); 
  int num_rank = formic::mpi::size();
  //int my_rank;
  //int num_rank;
  //MPI_Comm_rank(MPI_COMM_WORLD, & my_rank);
  //MPI_Comm_size(MPI_COMM_WORLD, & num_rank);

  // size the resulting matrix correctly
  y.reset(x.rows(), x.cols());

  // the number of samples on each process
  int Ns = _le_der.rows();

  // the number of independent variables 
  int Nind = _le_der.cols();

  // the number of new krylov vectors
  int Nnew = x.cols();

  // if the approximate flag is set to be true, throw out an error
  if ( approximate ) 
    throw formic::Exception("Matrix-Matrix multiplication doesn't support appriximate matrix");

  // check to see whether derivative vector matrices have the same size
  if ( _le_der.rows() != _der_rat.rows() && _le_der.cols() != _der_rat.cols())
    throw formic::Exception("the input derivative vector %d by %d and %d by %d matrices are of different size!") % _le_der.rows() % _le_der.cols() % _der_rat.rows() % _der_rat.cols();

  // if we are doing ground state calculation 
  if ( _ground ) {
    
    // if we do H*x
    if ( !transpose ) {
      
      // temp matrix to store le_der * x
      formic::Matrix<double> temp(Ns, Nnew);
      
      // call blas level-3 function
      formic::dgemm('N', 'N', Ns, Nnew, Nind, 1.0, &_le_der.at(0, 0), Ns, &x.at(0, 0), Nind, 0.0, &temp.at(0, 0), Ns);

      // call blas level-3 function 
      formic::dgemm('T', 'N', Nind, Nnew, Ns, 1.0,  &_der_rat.at(0, 0), Ns, &temp.at(0, 0), Ns, 0.0, &y.at(0, 0), Nind);

      return;
    }

    // if we do H^T * x
    else {
     
      // temp mattrix that stores der_rat * x
      formic::Matrix<double> temp(Ns, Nnew);

      // call blas level-3 function 
      formic::dgemm('N', 'N', Ns, Nnew, Nind, 1.0, &_der_rat.at(0, 0), Ns, &x.at(0, 0), Nind, 0.0, &temp.at(0, 0), Ns);

      // call blas level-3 function 
      formic::dgemm('T', 'N', Nind, Nnew, Ns, 1.0, &_le_der.at(0, 0), Ns, &temp.at(0, 0), Ns, 0.0, &y.at(0, 0), Nind);

      return;
    }
  }

  // if we are doing excited state calculation 
  else if ( !_ground ) {
    
    // if we want H*x
    if ( !transpose ) {
      
      // temp matrix that stores omega * der_rat * x
      formic::Matrix<double> temp1(Ns, Nnew);

      // temp matrix that stores le_der * x
      formic::Matrix<double> temp2(Ns, Nnew);

      // call blas level-3 function 
      formic::dgemm('N', 'N', Ns, Nnew, Nind, _hd_shift, &_der_rat.at(0, 0), Ns, &x.at(0, 0), Nind, 0.0, &temp1.at(0, 0), Ns);

      // call blas level-3 function 
      formic::dgemm('N', 'N', Ns, Nnew, Nind, 1.0, &_le_der.at(0, 0), Ns, &x.at(0, 0), Nind, 0.0, &temp2.at(0, 0), Ns);

      // combine these two temp matrices together
      temp1 -= temp2;

      // left multiply by _der_rat^T
      formic::dgemm('T', 'N', Nind, Nnew, Ns, 1.0, &_der_rat.at(0, 0), Ns, &temp1.at(0, 0), Ns, 0.0, &y.at(0, 0), Nind);

      return;
    }

    // if we want H^T*x
    else {
      
      // temp vector that stored _der_rat * x
      formic::Matrix<double> temp1(Ns, Nnew);

      // call blas level-3 function 
      formic::dgemm('N', 'N', Ns, Nnew, Nind, 1.0, &_der_rat.at(0, 0), Ns, &x.at(0, 0), Nind, 0.0, &temp1.at(0, 0), Ns);

      // temp matrix that stores _le_der^T * _der_rat * x
      formic::Matrix<double> temp2(Nind, Nnew);

      // call blas level-3 function
      formic::dgemm('T', 'N', Nind, Nnew, Ns, 1.0, &_le_der.at(0, 0), Ns, &temp1.at(0, 0), Ns, 0.0, &temp2.at(0, 0), Nind);

      // call blas level-3 function 
      formic::dgemm('T', 'N', Nind, Nnew, Ns, _hd_shift, &_der_rat.at(0, 0), Ns, &temp1.at(0, 0), Ns, 0.0, &y.at(0, 0), Nind);

      // get the resulting vector 
      y -= temp2;

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

void cqmc::engine::SpamLMHD::SMatVecOp(const formic::ColVec<double> & x, formic::ColVec<double> & y, const bool approximate)
{
  
  // size the resulting vector correctly
  y.reset(x.size());

  // since we do not have the matrix, then we need to do der_rat * le_der * x on each process
  // the matrix free implementation is a little bit complicated, I will explain it here 
  // we have two derivative vector matrox, L(i, j) = <i|H|[psi_j>/<i|psi>, D(i, j) = <i|psi_j>/<i|psi> 
  // i denotes configuration(electron position in real space and number vector in Hilbert space) 
  // psi_j denotes wfn derivative w.r.t. jth variable, and j=0 means undifferentiated wfn
  // note that |value/guiding|^2 and weights should be absorbed in L and D matrix
  // in ground state calculation, H = D^T * D, Hx = D^T * Dx
  // in excited state calculation, H = (omega * D - L)^T * (omega * D - L), temp1 = omage * D * x
  // temp2 = Lx

  // number of samples
  int Ns = _le_der.rows();

  // the number of independent variables + 1
  const int Nind = _le_der.cols();

  // check to see whether derivative vector matrices have the same size
  if ( _le_der.rows() != _der_rat.rows() && _le_der.cols() != _der_rat.cols() )
    throw formic::Exception("input derivative vectors are of different size");

  // modify number of samples based on approximate degree
  if ( approximate )
    Ns /= _appro_degree;

  // if we are doing ground state calculation
  if ( _ground && !approximate ) {
    
    // temp vector that store _der_rat * x
    formic::ColVec<double> temp(Ns);

    // call blas level-2 function 
    formic::dgemv('N', Ns, Nind, 1.0, &_der_rat.at(0, 0), Ns, &x.at(0), 1, 0.0, &temp.at(0), 1);

    // call blas levec-2 function 
    formic::dgemv('T', Ns, Nind, 1.0, &_der_rat.at(0, 0), Ns, &temp.at(0), 1, 0.0, &y.at(0), 1);

    return;
  }

  // if we multiply by approximate matrix 
  else if ( _ground && approximate ) {
    
    // temp vector that stores _der_rat * x
    formic::ColVec<double> temp(Ns);

    // call blas level-2 function 
    formic::dgemv('N', Ns, Nind, 1.0, &_der_rat_appro.at(0, 0), Ns, &x.at(0), 1, 0.0, &temp.at(0), 1);

    // call blas level-2 function 
    formic::dgemv('T', Ns, Nind, _appro_factor, &_der_rat_appro.at(0, 0), Ns, &temp.at(0), 1, 0.0, &y.at(0), 1);

    return;
  }

  else if ( !_ground && !approximate ) {
    
    // temp vectpr to store omega * der_rat * x
    formic::ColVec<double> temp1(Ns);

    // temp vector to store le_der * x
    formic::ColVec<double> temp2(Ns);

    // temp vector that store omega * der_rat^T * (omega * der_rat - le_der) * x
    formic::ColVec<double> temp3(x.size());

    // call blas level-2 function 
    formic::dgemv('N', Ns, Nind, _hd_shift, &_der_rat.at(0, 0), Ns, &x.at(0), 1, 0.0, &temp1.at(0), 1);

    // call blas level-2 function 
    formic::dgemv('N', Ns, Nind, 1.0, &_le_der.at(0, 0), Ns, &x.at(0), 1, 0.0, &temp2.at(0), 1);

    // combine these two temp vector together 
    temp1 -= temp2;

    // omega * D^T * (omega * D - L) * x  
    formic::dgemv('T', Ns, Nind, _hd_shift, &_der_rat.at(0, 0), Ns, &temp1.at(0), 1, 0.0, &y.at(0), 1);

    // L^T * (omega * D - L) * x
    formic::dgemv('T', Ns, Nind, 1.0, &_le_der.at(0, 0), Ns, &temp1.at(0), 1, 0.0, &temp3.at(0), 1);

    // (omega * D^T - L^T) * (omega * D - L) * x
    y -= temp3;

    return;
  }

  // if we multiply by approximate matrix 
  else if ( !_ground && approximate ) {
    
    // temp vectpr to store omega * der_rat * x
    formic::ColVec<double> temp1(Ns);

    // temp vector to store le_der * x
    formic::ColVec<double> temp2(Ns);

    // temp vector that store omega * der_rat^T * (omega * der_rat - le_der) * x
    formic::ColVec<double> temp3(x.size());

    // call blas level-2 function 
    formic::dgemv('N', Ns, Nind, _hd_shift, &_der_rat_appro.at(0, 0), Ns, &x.at(0), 1, 0.0, &temp1.at(0), 1);

    // call blas level-2 function 
    formic::dgemv('N', Ns, Nind, 1.0, &_le_der_appro.at(0, 0), Ns, &x.at(0), 1, 0.0, &temp2.at(0), 1);

    // combine these two temp vector together 
    temp1 -= temp2;

    // omega * D^T * (omega * D - L) * x  
    formic::dgemv('T', Ns, Nind, _hd_shift, &_der_rat_appro.at(0, 0), Ns, &temp1.at(0), 1, 0.0, &y.at(0), 1);

    // L^T * (omega * D - L) * x
    formic::dgemv('T', Ns, Nind, 1.0, &_le_der_appro.at(0, 0), Ns, &temp1.at(0), 1, 0.0, &temp3.at(0), 1);

    // (omega * D^T - L^T) * (omega * D - L) * x
    y -= temp3;

    // account for approximation prefactor
    y *= _appro_factor;

    return;
  }
}

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

void cqmc::engine::SpamLMHD::SMatMatOp(const formic::Matrix<double> & x, formic::Matrix<double> & y, const bool approximate)
{
  
  int my_rank = formic::mpi::rank(); 
  int num_rank = formic::mpi::size();
  //int my_rank;
  //int num_rank;
  //MPI_Comm_rank(MPI_COMM_WORLD, & my_rank);
  //MPI_Comm_size(MPI_COMM_WORLD, & num_rank);

  // size the resulting matrix correctly
  y.reset(x.rows(), x.cols());

  // the number of samples on each process
  int Ns = _le_der.rows();

  // the number of independent variables 
  int Nind = _le_der.cols();

  // the number of new krylov vectors
  int Nnew = x.cols();

  // if the approximate flag is set to be true, throw out an error
  if ( approximate ) 
    throw formic::Exception("Matrix-Matrix multiplication doesn't support approximate matrix");

  // check to see whether derivative vector matrices have the same size
  if ( _le_der.rows() != _der_rat.rows() && _le_der.cols() != _der_rat.cols())
    throw formic::Exception("the input derivative vector matrices are of different size!");

  // if we are doing ground state calculation
  if ( _ground ) {
    
    // temp matrix that stores _der_rat * x
    formic::Matrix<double> temp(Ns, Nnew);

    // call blas level-3 function 
    formic::dgemm('N', 'N', Ns, Nnew, Nind, 1.0, &_der_rat.at(0, 0), Ns, &x.at(0, 0), Nind, 0.0, &temp.at(0, 0), Ns);

    // call blas level-3 function 
    formic::dgemm('T', 'N', Nind, Nnew, Ns, 1.0, &_der_rat.at(0, 0), Ns, &temp.at(0, 0), Ns, 0.0, &y.at(0, 0), Nind);

    return;
  }

  // if we are doing excited state calculation
  else {
    
    // temp matrix that stores omega * der_rat * x
    formic::Matrix<double> temp1(Ns, Nnew);

    // temp matrix that stores le_der * x
    formic::Matrix<double> temp2(Ns, Nnew);

    // temp matrix that stores omega * der_rat^T * (omega * der_rat - le_der) * x
    formic::Matrix<double> temp3(Nind, Nnew);
    
    // call blas level-3 function 
    formic::dgemm('N', 'N', Ns, Nnew, Nind, _hd_shift, &_der_rat.at(0, 0), Ns, &x.at(0, 0), Nind, 0.0, &temp1.at(0, 0), Ns);

    // call blas level-3 function 
    formic::dgemm('N', 'N', Ns, Nnew, Nind, 1.0, &_le_der.at(0, 0), Ns, &x.at(0, 0), Nind, 0.0, &temp2.at(0, 0), Ns);

    // combine these two temp vectors together 
    temp1 -= temp2;

    // omega * D^T * (omega * D - L) * x
    formic::dgemm('T', 'N', Nind, Nnew, Ns, _hd_shift, &_der_rat.at(0, 0), Ns, &temp1.at(0, 0), Ns, 0.0, &y.at(0, 0), Nind);

    // L^T * (omega * D - L) * x
    formic::dgemm('T', 'N', Nind, Nnew, Ns, 1.0, &_le_der.at(0, 0), Ns, &temp1.at(0, 0), Ns, 0.0, &temp3.at(0, 0), Nind);

    // (omega * D^T - L^T) * (omega * D - L) * x
    y -= temp3;

    return;
  }
}


/////////////////////////////////////////////////////////////////////////////////
// \brief constructor
//
//
//
/////////////////////////////////////////////////////////////////////////////////
	
cqmc::engine::SpamLMHD::SpamLMHD(const formic::VarDeps * dep_ptr,
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
                                 formic::Matrix<double> & le_der_appro)
:EigenSolver(dep_ptr,
             nfds,
             lm_eigen_thresh,
             var_deps_use,
             chase_lowest,
             chase_closest,
             ground,
             false,
             vf,
             init_energy,
             0.0,
             hd_shift,
             0.0,
             lm_max_e_change,
             total_weight,
             vgsa,
             der_rat,
             le_der),
_nkry(0),
_nkry_full(0),
_appro_degree(appro_degree),
_n_max_iter(lm_krylov_iter),
_inner_maxIter(inner_maxIter),
_inner_print(inner_print),
_singular_value_threshold(lm_min_S_eval),
_appro_factor(appro_factor),
_init_energy(init_energy),
_energy_outer(init_energy),
_energy_inner(init_energy),
_der_rat_appro(der_rat_appro),
_le_der_appro(le_der_appro),
_smallest_sin_value_inner(0.0),
_smallest_sin_value_outer(0.0)
{
  // get rank number and number of ranks 
  int my_rank = formic::mpi::rank(); 
  int num_rank = formic::mpi::size();
  //int my_rank; 
  //int num_rank;
  //MPI_Comm_rank(MPI_COMM_WORLD, & my_rank);
  //MPI_Comm_size(MPI_COMM_WORLD, & num_rank);
		
}

/////////////////////////////////////////////////////////////////////////////////////
// \brief solves the eigenvalue problem via the normal davidson method
//
//
//
/////////////////////////////////////////////////////////////////////////////////////
  	
bool cqmc::engine::SpamLMHD::iterative_solve(double & eval, std::ostream & output)
{

  // get rank number and number of rank 
  int my_rank = formic::mpi::rank(); 
  int num_rank = formic::mpi::size();
  //int my_rank;
  //int num_rank;
  //MPI_Comm_rank(MPI_COMM_WORLD, & my_rank);
  //MPI_Comm_size(MPI_COMM_WORLD, & num_rank);

  // initialize the solution vector to the unit vector along the first direction 
  _evecs.reset( ( _var_deps_use ? 1 + _dep_ptr->n_tot() : _nfds ), 0.0 );
  _evecs.at(0) = 1.0;

  // ensure that at least one vector is in the outer krylov subspace 
  if ( my_rank == 0 && _nkry == 0) 
    throw formic::Exception("Empty krylov subspace upon entry to iterative_solve. Did you forget to add the initial vector?");

  // return value, whether the solver is successful or not 
  bool retval = true;

  // converge flag(outer)
  bool converged_outer = false;

  // converge flag(inner)
  bool converged_inner = false;

  // best outer residual
  double _best_residual_outer = 1.0e100;

  // best inner residual 
  double _best_residual_inner = 1.0e100;

  // times of outer iteration 
  int iter_outer = 0;

  // times of inner iteration
  int iter_inner = 0;

  // print out that we have started the iteration 
  if (my_rank == 0)
    output << boost::format("iteration solving starts here(engine) \n") << std::endl << std::endl;

  while(true) {
    
    // smallest singular value 
    double smallest_sin_value_outer = 0.0;

    // solve subspace eigenvalue problem on root process 
    if ( my_rank == 0 ) { 
      this -> solve_subspace_nonsymmetric(true);
    }

    // send resulting eigenvalues to all processes and record it as the new best estimate
    formic::mpi::bcast(&_energy_outer, 1);
    eval = _energy_outer;

    // check if the energy has an imaginary component and stop iteration when it does
    formic::mpi::bcast(&_eval_was_complex, 1);
    if ( _eval_was_complex ) {
      if ( my_rank == 0 ) 
        output << boost::format("spam iteration %4i: stopping due to imaginary component in energy") % iter_outer << std::endl;
      break;
    }

    // if energy change is unreasonable, stop iterating and set bad solve flag
    if (std::abs(_energy_outer - _init_energy) > _max_energy_change) {
      retval = false;
      if ( my_rank == 0 )
        output << boost::format("spam iteration %4i stopping due to too large eigenvalue change") % iter_outer << std::endl;
      break;
    }

    // if the overlap matrix becomes singular, stop iterating
    formic::mpi::bcast(&_smallest_sin_value_outer, 1);
    if (std::abs(_smallest_sin_value_outer) < _singular_value_threshold) {
      if (my_rank == 0) 
        output << boost::format("spam iteration %4i stopping due to small subspace S singular value of %.2e") % iter_outer % _smallest_sin_value_outer << std::endl;
      break;
    }

    // construct new Krylov vector from subspace eigenvector
    _wv1.reset(_nfds);
    if (my_rank == 0) {
      _wv6.reset(_nfds);
      _wv6 = _kvecs * _sub_evec_outer;
      // get normalization factor
      const double temp_norm = std::sqrt(_wv6.norm2());

      // normalize this new vector
      _wv6 /= temp_norm;

      // add up linear combination of Hamiltonian and overlap products to make the new residual vector
      _wv1.reset(_nfds);
      _wv1 = _hvecs * _sub_evec_outer;
      _wv1 = _wv1 - _sub_eval_outer * _svecs * _sub_evec_outer;
    }

    // send this new vector to all processes
    formic::mpi::bcast(&_wv1.at(0), _wv1.size());

    // compute the residual norm and send it to all processes
    double current_outer_residual;
    current_outer_residual = _wv1.norm2();
    formic::mpi::bcast(&current_outer_residual, 1);

    // if this is the best residual, save it and save the new eigenvector estimate
    if (my_rank == 0 && current_outer_residual < _best_residual_outer) {
      _best_residual_outer = current_outer_residual;

      // get current eigenvector estimate, which corresponds to the set of independent variables
      formic::ColVec<double> ind_evecs;
      ind_evecs = _wv6.clone();
    
      // if our actual variables are dependent on the set of independent variables worded with here, expand the eigenvector into the full set of variables
      if ( _var_deps_use ) {
      
        // size the eigenvector correctly
        _evecs.reset( (1 + _dep_ptr -> n_tot()), 0.0);
        _evecs.at(0) = ind_evecs.at(0);

        // get some temporary vectors
        formic::ColVec<double> _evec_temp(_dep_ptr -> n_tot());
        formic::ColVec<double> _ind_temp(_dep_ptr -> n_ind());
        for (int i = 0; i < _ind_temp.size(); i++) {
          _ind_temp.at(i) = ind_evecs.at(i+1);
        }
        _dep_ptr -> expand_ind_to_all(&_ind_temp.at(0), &_evec_temp.at(0));
      
        for ( int i = 0; i < _evec_temp.size(); i++) {
          _evecs.at(i+1) = _evec_temp.at(i);
        }
      }

      // otherwise just copy the eigenvector into output since the independent and total variable sets are the same 
      else {
        _evecs = ind_evecs.clone();
      }
    }

    // print iteration results
    if (my_rank == 0) {

      // if we are doing ground state calculation, then print out energy
      if ( _ground ) 
        output << boost::format("spam outer iteration %4i:   krylov dim = %3i   energy = %20.12f       residual = %.2e           smallest_sin_value = %.2e")
        % iter_outer
        % _nkry
        % _energy_outer
        % current_outer_residual
        % _smallest_sin_value_outer
        << std::endl;

      // if we are doing excited state calculation, then print out target function value 
      else 
        output << boost::format("spam outer iteration %4i:   krylov dim = %3i   tar_fn = %20.12f       residual = %.2e           smallest_sin_value = %.2e")
        % iter_outer
        % _nkry
        % _sub_eval_outer
        % current_outer_residual
        % _smallest_sin_value_outer
        << std::endl;
    }

    // check for convergence
    converged_outer = current_outer_residual < _residual_threshold;

    // if iteration has already converged, we exit iteration 
    if ( converged_outer ) 
      break;

    // if iteration hasn't converged, we increment the iteration count by 1 and stop if maximum number of iterations has been reached
    if ( iter_outer ++ >= _n_max_iter)
      break;

    // now this is important, we add the new krylov basis vector to inner loop 
    this -> add_krylov_vector_inner(_wv1);

    // now enter inner loop 
    while (true) {
     
      // average of smallest singular value
      double smallest_sin_val_avg_inner = 0.0;

      // solve subspace eigenvalue problem on root process 
      if ( my_rank == 0 ) {
        this -> solve_subspace_nonsymmetric(false);
      }

      // send resulting eigenvalues to all processes and record it as the new best estimate
      formic::mpi::bcast(&_energy_inner, 1);
      
      // check if the energy(or target function) has an imaginary component and stop iteration when it does
      formic::mpi::bcast(&_eval_was_complex, 1);
      if ( _eval_was_complex ) {
        if ( my_rank == 0 ) 
          output << boost::format("spam outer iteration %4i inner iteration %4i: stopping due to imaginary component in energy") % iter_outer % iter_inner << std::endl;
        break;
      }

      // if energy(or target function) change is unreasonable, stop iterating but don't set bag solve flag
      if (std::abs(_energy_inner - _init_energy) > _max_energy_change) { 
        //retval = false;
        if ( my_rank == 0 ) 
          output << boost::format("spam outer iteration %4i inner iteration %4i: stopping due to too large eigenvalue change") % iter_outer % iter_inner << std::endl;
        break;
      }

      // if the overlap matrix becomes singular, stop iterating
      formic::mpi::bcast(&_smallest_sin_value_inner, 1);
      if (std::abs(_smallest_sin_value_inner) < _singular_value_threshold) {
        if (my_rank == 0) 
          output << boost::format("spam outer iteration %4i inner iteration %4i: stopping due to too small S singular value of %.2e") % iter_outer % iter_inner % _smallest_sin_value_inner << std::endl;
        break;
      }

      // construct new krylov vector from subspace eigenvector
      if (my_rank == 0) {
        _wv6.reset(_nfds);
        _wv6 = _kvecs * _sub_evec_inner;
        // get normalization factor 
        const double temp_norm = std::sqrt(_wv6.norm2());

        // normalize this new vector
        _wv6 /= temp_norm;
        
        // add up linear combination of Hamiltonian and overlap products to make the new residual vector
        _wv1.reset(_nfds);
        _wv1 = _hhvecs * _sub_evec_inner;
        _wv1 = _wv1 - _sub_eval_inner * _hsvecs * _sub_evec_inner;
      }

      // send this new vector to all processes
      formic::mpi::bcast(&_wv1.at(0), _wv1.size());
      
      // compute the residual norm and send it to all processes
      double current_inner_residual;
      current_inner_residual = _wv1.norm2();
      formic::mpi::bcast(&current_inner_residual, 1);

      // if this is the best residual, save it and save the new eigenvector estimate
      if (my_rank == 0 && current_inner_residual < _best_residual_inner) 
        _best_residual_inner = current_inner_residual;
      
      // print iteration results if request
      if (my_rank == 0 && _inner_print) {
        
        // if we are doing ground state calculation, then print out energy
        if ( _ground ) 
          output << boost::format("spam outer iteration %4i inner iteration %4i: krylov dim = %3i  energy = %20.12f      residual = %.2e        smallest_sin_value = %.2e")
          % iter_outer
          % iter_inner
          % _kvecs.cols()
          % _energy_inner
          % current_inner_residual
          % _smallest_sin_value_inner
          << std::endl;

        // if we are doing excited state calculation, then print out target function value
        else 
          output << boost::format("spam outer iteration %4i inner iteration %4i: krylov dim = %3i  energy = %20.12f      residual = %.2e        smallest_sin_value = %.2e") 
          % iter_outer
          % iter_inner
          % _kvecs.cols()
          % _sub_eval_inner
          % current_inner_residual
          % _smallest_sin_value_inner
          << std::endl;
      }

      // check for convergence
      converged_inner = current_inner_residual < _residual_threshold;

      // if iteration has already converged, we exit iteration 
      if ( converged_inner ) 
        break;

      // if iteration hasn't converged, we increment the iteration count by 1 and stop if maximum number of iterations has been reached
      if ( iter_inner ++ >= _inner_maxIter ) 
        break;

      // add this new krylov basis vector to inner loop
      this -> add_krylov_vector_inner(_wv1);

    // end of spam inner loop 
    }
    
    // size the intermediate vectors correctly for non-root process
    if ( my_rank != 0 ) 
      _kvecs_about_to_add.reset(_nfds, _nkry - _nkry_full);
      
    // broadcast intermediate vectors that will be added to full krylov space 
    formic::mpi::bcast(&_kvecs_about_to_add.at(0, 0), _kvecs_about_to_add.size());

    // add these new krylov vectors to outer loop 
    this -> add_krylov_vectors_outer(_kvecs_about_to_add);

    // clear these intermediate vectors 
    _kvecs_about_to_add.reset(0, 0);
  
    // reset the number of inner iterations
    iter_inner = 0;

  // end of spam outer loop
  }

  // print iteration information 
  if (converged_outer && my_rank == 0) 
    output << boost::format("spam solver converged in %10i iterations") % iter_outer << std::endl << std::endl;

  else if (my_rank == 0) 
    output << boost::format("spam solver did not converge after %10i iterations") % iter_outer << std::endl << std::endl;

  return retval;

}

///////////////////////////////////////////////////////////////////////////////////////////////
// \brief solves the eigenvalue problem
//
//
//
///////////////////////////////////////////////////////////////////////////////////////////////

bool cqmc::engine::SpamLMHD::solve(double & eval, std::ostream & output)
{
  return this -> iterative_solve(eval, output);
}

////////////////////////////////////////////////////////////////////////////////////
// \brief updates hamiltonian * krylov vector and hamiltonian projection based on 
//        new shift 
//
//
////////////////////////////////////////////////////////////////////////////////////

void cqmc::engine::SpamLMHD::update_hvecs_sub(const double new_i_shift, const double new_s_shift)
{
  
  // get rank number and number of ranks
  int my_rank = formic::mpi::rank();
  int num_rank = formic::mpi::size();
  //MPI_Comm_rank(MPI_COMM_WORLD, & my_rank);
  //MPI_Comm_size(MPI_COMM_WORLD, & num_rank);

  // get the different between new shift and old shift
  const double diff_shift_i = new_i_shift - _hshift_i;
  const double diff_shift_s = new_s_shift - _hshift_s;

  if (my_rank == 0) {
    // update "identity shift" for the hamiltonian product
    for (int j = 0; j < _nkry; j ++) {
      for (int i = 1; i < _nfds; i ++) {
        _hvecs.at(i, j) += diff_shift_i * _kvecs.at(i, j);
        _thvecs.at(i, j) += diff_shift_i * _kvecs.at(i, j);
      }
    }

    // update "overlap shift" for the hamiltonian product 
    for (int j = 1; j < _nkry; j++) {
      for (int i = 0; i < _nfds; i++) {
        _hvecs.at(i, j) += diff_shift_s * _svecs.at(i, j);
        _thvecs.at(i, j) += diff_shift_s * _svecs.at(i, j);
      }
    }

    _hhvecs = _hvecs;

    // update projection of hamiltonian matrix 
    _subH = _kvecs.t() * _hvecs;
    _hy_subH = _subH.clone();
  }
}

////////////////////////////////////////////////////////////////////////////////////
// \brief reset the eigen solver 
// 
// \brief clear subspace projection of Hamiltonian and overlap matrix, clear Krylov
//        subspace and action of Hamiltonian and overlap matrix
////////////////////////////////////////////////////////////////////////////////////

void cqmc::engine::SpamLMHD::child_reset()
{
  // clear subspace projection of Hamiltonian and overlap matrix
  _subH.reset(0, 0);
  _hy_subH.reset(0, 0);
  _subS.reset(0, 0);
  _hy_subS.reset(0, 0);

  // clear Krylov subspace
  _nkry = 0;
  _nkry_full = 0;
  _kvecs.reset(0, 0);
  _kvecs_about_to_add.reset(0, 0);

  // clear Hamiltonian and overlap matrix's action on krylov subspace
  _hvecs.reset(0, 0);
  _thvecs.reset(0, 0);
  _ahvecs.reset(0, 0);
  _athvecs.reset(0, 0);
  _hhvecs.reset(0, 0);

  _svecs.reset(0, 0);
  _asvecs.reset(0, 0);
  _hsvecs.reset(0, 0);

  // clear eigenvector and wavefunction coefficients
  _sub_evec_outer.reset(0);
  _sub_evec_inner.reset(0);
  _evecs_inner.reset(0);

  // clear all values calculated from last solve
  _sub_eval_outer = 0.0;
  _sub_eval_inner = 0.0;
  _smallest_sin_value_outer = 1e10;
  _smallest_sin_value_inner = 1e10;

  // set the inner and outer energy to be initial energy
  _energy_inner = _init_energy;
  _energy_outer = _init_energy;
}

///////////////////////////////////////////////////////////////////////////////////////////////
// \brief converts eigenvectors into wave function coefficients
//        solving this question Hc = Sd for d, S is the ground overlap matrix 
//
//
///////////////////////////////////////////////////////////////////////////////////////////////

void cqmc::engine::SpamLMHD::child_convert_to_wf_coeff()
{
  
  return;

}
