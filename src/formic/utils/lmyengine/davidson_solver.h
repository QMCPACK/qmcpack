//////////////////////////////////////////////////////////////////////////////////////////////////
// \file davidson_solver.h       a class to perform generalized eigenvalue problem
//
//
//
//  solve with normal davidson method(either build or not build the matrix)
//////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef EIGINE_DAVIDSON_SOLVER_HEADER
#define EIGINE_DAVIDSON_SOLVER_HEADER

#include <vector>
#include <string>
#include <numeric>
#include <cassert>
#include <algorithm>
#include <cmath>
#include <complex>
#include <iostream>

#include <boost/format.hpp>
#include <boost/shared_ptr.hpp>

//#include<mpi.h>

#include "formic/utils/zero_one.h"
#include "formic/utils/numeric.h"
#include "formic/utils/lmyengine/var_dependencies.h"
#include "formic/utils/lmyengine/eigen_solver.h"
#include "formic/utils/matrix.h"


namespace cqmc
{
namespace engine
{
template<typename S>
class DavidsonLMHD : public EigenSolver<S>
{
private:
  /// \brief number of krylov vectors
  int _nkry;

  /// \brief maximum number of iteration allowed in the solver
  int _n_max_iter;

  /// \brief the smallest subspace overlap matrix singular value
  double _smallest_sin_value;

  /// \brief the subspace eigenvalue(in harmonic Davidson this is not the energy)
#ifndef QMC_COMPLEX
  double _sub_eval;
#else
  std::complex<double> _sub_eval_complex;
#endif

  /// \brief the imaginary part of subspace eigenvalue
  double _sub_eval_imag;

  /// \brief subspace overlap matrix singular value threshold
  double _singular_value_threshold;

  // \brief best residual
  double _best_residual;

  /// \brief difference between new and old tau
  double _tau_diff;

  // \brief flag to tell whether we have built the matrix or not
  bool _build_lm_matrix;

  /// \brief the Hamiltonian matrix
  formic::Matrix<S>& _hmat;

  /// \brief the overlap matrix
  formic::Matrix<S>& _smat;

  /// \brief the overlap matrix of normal linear method
  formic::Matrix<S>& _lmsmat;

  /// \brief the preconditioning matrix
  //formic::Matrix<double> _mmat;

  /// \brief inverse-square-root of the diagonal of the overlap matrix
  //std::vector<double> _ovl_diag_neg_half;

  /// \brief normalized krylov space basis vectors
  formic::Matrix<S> _kvecs;

  /// \brief vectors resulting from Hamiltonian action on the krylov basis
  formic::Matrix<S> _hvecs;

  /// \brief vectors resulting from Hamiltonian transpose action on the krylov basis
  formic::Matrix<S> _htvecs;

  /// \brief vectors resulting from overlap matrix action on the krylov basis
  formic::Matrix<S> _svecs;

  /// \brief subspace Hamiltonian matrix
  formic::Matrix<S> _subH;

  /// \brief subspace overlap matrix
  formic::Matrix<S> _subS;

  /// \brief work vector
  formic::ColVec<S> _wv1;

  /// \brief work vector
  formic::ColVec<S> _wv2;

  /// \brief work vector
  formic::RowVec<S> _wv3;

  /// \brief work vector
  formic::ColVec<std::complex<double>> _wv4;

  /// \brief work vector
  formic::ColVec<std::complex<double>> _wv5;

  /// \brief work vector
  formic::ColVec<S> _wv6;

  /// \brief work vector
  formic::ColVec<S> _wvX;

  /// \brief work vector
  formic::ColVec<S> _wvY;

  /// \brief subspace eigenvector
#ifndef QMC_COMPLEX
  formic::ColVec<double> _sub_evec;
#else
  formic::ColVec<std::complex<double>> _sub_evec_complex;
#endif


public:
  ///////////////////////////////////////////////////////////////////
  // \brief solve the subspace generalized eigenvalue problem
  //        with nonsymmetric H and S
  //
  //
  ///////////////////////////////////////////////////////////////////

  void solve_subspace_nonsymmetric(const bool outer) override
  {
    // get rank number and number of ranks
    int my_rank = formic::mpi::rank();

    //one and zero in complex form
    const std::complex<double> complex_one(1.0, 0.0);
    const std::complex<double> complex_zero(0.0, 0.0);

    //get subspace dimension
    const int m = _nkry;

    //create vectors and matrices used in svd routine
    formic::Matrix<S> u, v, vt;
    formic::ColVec<S> sin_vals;
    int truncate_index = 0;

    // make sure the subspace matrix is not empty
    if (_subS.rows() == 0 || _subH.rows() == 0)
      throw formic::Exception("subspace matrix is empty upon solving subspace eigenvalue problem");

    //perform svd to subspace overlap matrix
    _subS.svd(u, sin_vals, vt);
    v = vt.clone();
    v.cip();

    //record the smallest singular value of subspace overlap matrix
    _smallest_sin_value = std::abs(formic::real(sin_vals.at(0)));
    for (int i = 0; i < m; i++)
    {
      _smallest_sin_value = std::min(_smallest_sin_value, std::abs(formic::real(sin_vals(i))));
      truncate_index      = i;

      //check if singular value is smaller that the singular value threshold
      if (_smallest_sin_value < _singular_value_threshold)
      {
        break;
      }
    }

    //get the number of colums of new U and V matrix by add 1 to truncate index
    truncate_index++;

    //throw away those columns in U and V matrix which corresponds to singular values below the threshold
    u.conservativeResize(m, truncate_index);
    v.conservativeResize(m, truncate_index);

    //throw away those small singular values
    sin_vals.conservativeResize(truncate_index);

    //convert the truncated singular value vector to a truncated_index * truncated_index diagonal matrix
    formic::Matrix<S> trun_sin_val_matrix(truncate_index, truncate_index, formic::zero(S()));
    for (int i = 0; i < truncate_index; i++)
      trun_sin_val_matrix.at(i, i) = sin_vals.at(i);

    //calculate the inverse of this matrix
    for (int i = 0; i < truncate_index; i++)
    {
      for (int j = 0; j < truncate_index; j++)
      {
        trun_sin_val_matrix.at(i, j) =
            (trun_sin_val_matrix.at(i, j) == formic::zero(S()) ? trun_sin_val_matrix.at(i, j)
                                                               : formic::unity(S()) / trun_sin_val_matrix.at(i, j));
      }
    }

    //calculate matrix S_trun^-1 * U^-1 * H * V
    formic::Matrix<S> new_sub_H(truncate_index, truncate_index, formic::zero(S()));
    new_sub_H = trun_sin_val_matrix * u.c() * _subH * v;
    new_sub_H = _subH;

    //solve this standard eigenvalue problem ( new_H * y = lambda * y, y = V^T * x, x is the original eigenvector)
    formic::ColVec<std::complex<double>> e_evals;
    formic::Matrix<std::complex<double>> e_evecs;
    new_sub_H.nonsym_eig(e_evals, e_evecs);

    // set an eigen array to hold energies
    formic::ColVec<std::complex<double>> energy_list = e_evals.clone();

    int selected = 0;
    //if we want to chase the closest, selected the eigenvalue that is most similar to the previous eigenvalue( this is essestially an attempt to stay in the same solution)
    if (this->_chase_closest)
    {
      std::complex<double> closest_cost = e_evals.at(0);
      //for (int j = 1; j < truncate_index; j++){
      //  if(std::abs( complex_one * this->_cost - e_evals.at(j) ) < std::abs(complex_one * this->_cost - closest_cost )){
      //    selected = j;
      //    closest_cost = e_evals.at(j);
      //  }
      //}
      //formic::Matrix<S> new_e_evecs = v * e_evecs;
      double small_diff = std::abs(1.0 - std::abs(e_evecs.col_as_vec(selected).at(0)));
      for (int j = 1; j < truncate_index; j++)
      {
        double diff = std::abs(1.0 - std::abs(e_evecs.col_as_vec(j).at(0)));
        if (diff < small_diff)
        {
          selected     = j;
          closest_cost = e_evals.at(j);
          small_diff   = diff;
        }
      }

      // if the eigenvalue has an imaginary component, we abort
      this->_eval_was_complex = false;
      if (std::abs(closest_cost.imag()) > 1.0)
      {
        this->_eval_was_complex = true;
        return;
      }

      // if the eigenvalue is real, we record it and the corresponding eigenvector
      // record energy
      this->_cost = closest_cost.real();

// record the eigenvalue
#ifndef QMC_COMPLEX
      _sub_eval = this->_cost;
#else
      _sub_eval_complex = closest_cost;
#endif
      _sub_eval_imag = closest_cost.imag();

      // record the eigenvector y
      _wv4 = e_evecs.col_as_vec(selected);

#ifndef QMC_COMPLEX
      //formic::ColVec<double> _wv4_dup(_wv4.size(), 0.0);
      //for (int i = 0; i < _wv4.size(); i++)
      //  _wv4_dup.at(i) = formic::real(_wv4.at(i));

      // convert y to x
      //_wv5.reset(m);
      //formic::Matrix<std::complex<double> > v_complex(v.rows(), v.cols(), complex_zero);
      //for (int i = 0; i < v_complex.rows(); i++) {
      //  for (int j = 0; j < v_complex.cols(); j++) {
      //    v_complex.at(i,j) = std::complex<double>(v.at(i,j), 0.0);
      //  }
      //}

      //formic::xgemm('N', 'N', m, 1, truncate_index, complex_one, &v_complex.at(0,0), m, &_wv4.at(0), truncate_index, complex_zero, &_wv5.at(0), m);

      _sub_evec.reset(m);
      for (int i = 0; i < m; i++)
      {
        _sub_evec.at(i) = _wv4.at(i).real();
      }
#else
      _sub_evec_complex.reset(m);
      _sub_evec_complex = _wv4;
//std::cout << _sub_evec_complex.print("%12.6f", "sub_evec_complex");
//formic::xgemm('N', 'N', m, 1, truncate_index, complex_one, &v.at(0,0), m, &_wv4.at(0), truncate_index, complex_zero, &_sub_evec_complex.at(0), m);
#endif
    }

    // if we want to chase the lowest, selected the lowest eigenvalue
    if (this->_chase_lowest)
    {
      std::complex<double> lowest_eval = e_evals.at(0);
      for (int j = 1; j < truncate_index; j++)
      {
        if (e_evals.at(j).real() < lowest_eval.real())
        {
          selected    = j;
          lowest_eval = e_evals.at(j);
        }
      }

      // if the eigenvalue has an imaginary component, we abort
      EigenSolver<S>::_eval_was_complex = false;
      if (std::abs(lowest_eval.imag()) > 1.0)
      {
        this->_eval_was_complex = true;
        return;
      }

      // record the cost function
      this->_cost = lowest_eval.real();

// record the eigenvalue
#ifndef QMC_COMPLEX
      _sub_eval = this->_cost;
#else
      _sub_eval_complex = lowest_eval;
#endif
      _sub_eval_imag = lowest_eval.imag();

      // record the eigenvector y
      _wv4 = e_evecs.col_as_vec(selected);
//std::cout << _wv4.print("%12.6f", "wv4");
#ifndef QMC_COMPLEX
      //formic::ColVec<double> _wv4_dup(_wv4.size(), 0.0);
      //for (int i = 0; i < _wv4.size(); i++)
      //  _wv4_dup.at(i) = formic::real(_wv4.at(i));

      // convert y to x
      //_wv5.reset(m);

      //// convert v to complex form to make sure that all quantities in dgemm call are of the same type
      //formic::Matrix<std::complex<double> > v_complex(v.rows(), v.cols(), complex_zero);
      //for (int i = 0; i < v_complex.rows(); i++) {
      //  for (int j = 0; j < v_complex.cols(); j++) {
      //    v_complex.at(i,j) = std::complex<double>(v.at(i,j), 0.0);
      //  }
      //}

      //formic::xgemm('N', 'N', m, 1, truncate_index, complex_one, &v_complex.at(0,0), m, &_wv4.at(0), truncate_index, complex_zero, &_wv5.at(0), m);
      _sub_evec.reset(m);
      for (int i = 0; i < m; i++)
        _sub_evec.at(i) = _wv4.at(i).real();
#else
      _sub_evec_complex.reset(m);
      _sub_evec_complex = _wv4;
//formic::xgemm('N', 'N', m, 1, truncate_index, complex_one, &v.at(0,0), m, &_wv4.at(0), truncate_index, complex_zero, &_sub_evec_complex.at(0), m);
#endif
    }
  }

  ///////////////////////////////////////////////////////////////////////////////
  // \brief solves the subspace eigenvalue problem
  //
  //
  //
  ///////////////////////////////////////////////////////////////////////////////

  void solve_subspace(const bool outer) override
  {
    this->solve_subspace_nonsymmetric(outer);
    return;
  }

public:
  ////////////////////////////////////////////////////////////////////////////////
  // \brief constructor with given parameters
  //
  //
  //
  /////////////////////////////////////////////////////////////////////////////////

  DavidsonLMHD(const formic::VarDeps* dep_ptr,
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
               const std::vector<S>& vf,
               const double init_cost,
               const double init_variance,
               const double hd_shift,
               const double var_weight,
               const double lm_max_e_change,
               const double total_weight,
               const double vgsa,
               formic::Matrix<S>& der_rat,
               formic::Matrix<S>& le_der,
               formic::Matrix<S>& hmat,
               formic::Matrix<S>& smat,
               formic::Matrix<S>& lmsmat)
      : EigenSolver<S>(dep_ptr,
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
        _nkry(0),
        _n_max_iter(lm_krylov_iter),
        _smallest_sin_value(0.0),
        _singular_value_threshold(lm_min_S_eval),
        _build_lm_matrix(build_lm_matrix),
        _hmat(hmat),
        _smat(smat),
        _lmsmat(lmsmat)
  {
    // get rank number and number of ranks
    int my_rank = formic::mpi::rank();

    // initialize the eigenvector
    this->_evecs.reset((this->_var_deps_use ? 1 + this->_dep_ptr->n_tot() : this->_nfds), formic::zero(S()));
    this->_evecs.at(0) = formic::unity(S());

    // initialize the independent eigenvector
    this->_ind_evecs.reset((this->_var_deps_use ? 1 + this->_dep_ptr->n_ind() : this->_nfds), formic::zero(S()));
    this->_ind_evecs.at(0) = formic::unity(S());
  }

  ////////////////////////////////////////////////////////////////////////////////////
  // \brief evaluate tau and compute the "variance corrected" Hamiltonian
  //
  //
  ////////////////////////////////////////////////////////////////////////////////////

  void tau_and_correct_ham() override
  {
    // get rank number and number of ranks
    int my_rank = formic::mpi::rank();

    if (my_rank == 0)
    {
      // first remember the old tau
      double old_tau = this->_tau;

      // evaluate new tau
      this->_tau = formic::real(formic::dotc(this->_ind_evecs, _smat * this->_ind_evecs) /
                                formic::dotc(this->_ind_evecs, _lmsmat * this->_ind_evecs));

      // evaluate the difference between old and new tau
      _tau_diff = std::abs(old_tau - this->_tau);

      // modify hamiltonian matrix based on the difference between new and old tau
      double square_diff = this->_tau * this->_tau - old_tau * old_tau;

      // update variance modified hamiltonian matrix
      _hmat += -1.0 * this->_var_weight * 10.0 * square_diff * _lmsmat;
    }
  }

  ////////////////////////////////////////////////////////////////////////////////////
  // \brief adds a new Krylov basis vector for normal davidson method
  //
  //
  //
  ////////////////////////////////////////////////////////////////////////////////////

  void add_krylov_vector(const formic::ColVec<S>& v) override
  {
    // get rank number and number of ranks
    int my_rank = formic::mpi::rank();

    // check vector length
    if (my_rank == 0 && !_build_lm_matrix && v.size() != this->_der_rat.cols())
      throw formic::Exception("bad vector length of %d in DavidsonLMHD::add_krylov_vector: expected length of %d") %
          v.size() % this->_der_rat.cols();

    else if (my_rank == 0 && _build_lm_matrix && v.size() != this->_nfds)
      throw formic::Exception("bad vector length of %d in DavidsonLMHD::add_krylov_vector: expected length of %d") %
          v.size() % this->_nfds;

    // increment krylov subspace size and remember the old size
    const int nold = _nkry++;

    // copy the new vector to work vector
    _wv1.reset(v.size());
    _wv1 = v.clone();

    // perform graham schmidt orthogonolization against the existing krylov vectors
    if (my_rank == 0)
    {
      for (int i = 0; i < nold; i++)
      {
        _wv1 -= (formic::unity(S()) * formic::dotc(_svecs.col_as_vec(i), _wv1) * _kvecs.col_as_vec(i));
      }
    }

    // if we don't build the matrix, we need to broadcast this vector to all processes
    if (!_build_lm_matrix)
      formic::mpi::bcast(&_wv1.at(0), _wv1.size());

    // temp vector holding normal LM overlap matrix times this vector
    formic::ColVec<S> _wv_temp;
    // compute the product of the overlap matrix times the new krylov vector on the root process
    if (_build_lm_matrix)
    {
      if (my_rank == 0)
      {
        this->SMatVecOp(_wv1, _wv2);
        if (!this->_ground)
          this->LSMatVecOp(_wv1, _wv_temp);
      }
    }

    else
    {
      this->SMatVecOp(_wv1, _wv2);
      // reduce the result vector
      formic::ColVec<S> _wv2_avg(_wv2.size());
      formic::mpi::reduce(&_wv2.at(0), &_wv2_avg.at(0), _wv2.size(), MPI_SUM);
      _wv2 = _wv2_avg.clone();
    }

    // compute the product of the hamiltonian matrix times the new krylov vector
    formic::ColVec<S> hs(this->_nfds);
    if (_build_lm_matrix)
    {
      if (my_rank == 0)
      {
        this->HMatVecOp(_wv1, hs);
      }
    }
    else
    {
      this->HMatVecOp(_wv1, hs);
      // reduce the result vector
      formic::ColVec<S> hs_avg(hs.size());
      formic::mpi::reduce(&hs.at(0), &hs_avg.at(0), hs.size(), MPI_SUM);
      hs = hs_avg.clone();
    }

    // modify hamiltonian product to account for "identity" shift
    if (my_rank == 0)
    {
      for (int i = 1; i < hs.size(); i++)
      {
        hs.at(i) += (formic::unity(S()) * this->_hshift_i) * _wv1.at(i);
      }
    }

    // modify hamiltonian product to account for "overlap" shift
    if (my_rank == 0 && nold > 0)
    {
      if (this->_ground)
        hs += (formic::unity(S()) * this->_hshift_s) * _wv2;
      else
#ifdef QMC_COMPLEX
        hs += (formic::unity(S()) * this->_hshift_s) * _wv_temp;
#else
        hs += (formic::unity(S()) * this->_hshift_s) * _wv2;
#endif
    }

    // normalize the new krylov vector and save the vector and its operation on overlap matrix
    if (my_rank == 0)
    {
      const double norm = std::sqrt(formic::real(formic::dotc(_wv1, _wv2)));
      _wv1 /= norm;
      _wv2 /= norm;
      hs /= norm;

      _kvecs.conservativeResize(this->_nfds, _nkry);
      std::copy(_wv1.begin(), _wv1.end(), _kvecs.col_begin(_nkry - 1));
      _svecs.conservativeResize(this->_nfds, _nkry);
      std::copy(_wv2.begin(), _wv2.end(), _svecs.col_begin(_nkry - 1));

      // multiply new vector by the Hamiltonian
      _hvecs.conservativeResize(this->_nfds, _nkry);
      std::copy(hs.begin(), hs.end(), _hvecs.col_begin(_nkry - 1));

      // update subspace projection of the Hamiltonian
      _subH.conservativeResize(_nkry, _nkry);
      _wv3 = _wv1.c() * _hvecs;
      _wv6 = _kvecs.c() * hs;
      std::copy(_wv3.begin(), _wv3.end(), _subH.row_begin(_nkry - 1));
      std::copy(_wv6.begin(), _wv6.end(), _subH.col_begin(_nkry - 1));

      // update subspace projection of the overlap
      _subS.conservativeResize(_nkry, _nkry);
      _wv3 = _wv1.c() * _svecs;
      _wv6 = _kvecs.c() * _wv2;
      std::copy(_wv3.begin(), _wv3.end(), _subS.row_begin(_nkry - 1));
      std::copy(_wv6.begin(), _wv6.end(), _subS.row_begin(_nkry - 1));
    }
  }

  ////////////////////////////////////////////////////////////////////////////////
  // \brief adds a new Krylov basis vector for spam inner loop
  //
  //
  //
  ////////////////////////////////////////////////////////////////////////////////

  void add_krylov_vector_inner(const formic::ColVec<S>& v) override {}

  ////////////////////////////////////////////////////////////////////////////////
  // \brief adds a bunch of new Krylov vectors(which is a matrix) for spam outer
  //        loop
  //
  //
  ////////////////////////////////////////////////////////////////////////////////

  void add_krylov_vectors_outer(const formic::Matrix<S>& m) override {}

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

  void HMatVecOp(const formic::ColVec<S>& x,
                 formic::ColVec<S>& y,
                 const bool transpose   = false,
                 const bool approximate = false) override
  {
    int my_rank = formic::mpi::rank();

    // size the resulting vector correctly
    y.reset(x.size());

    // if we have already built the matrix, we just multiply it
    if (_build_lm_matrix)
    {
      // if the approximate flag is set to be true, then throw out a exception
      if (approximate)
        throw formic::Exception(
            "Davidson solver's matrix-vector multiplication routine doesn't support approximate matrix!");

      // call level-2 blas function
      formic::xgemv('N', _hmat.rows(), _hmat.cols(), formic::unity(S()), &_hmat.at(0, 0), _hmat.rows(), &x.at(0), 1,
                    formic::zero(S()), &y.at(0), 1);

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

    else
    {
      // the number of samples on each process
      const int Ns = this->_le_der.rows();

      // the number of independent variables
      const int Nind = this->_le_der.cols();

      // check whether derivative vector matrices have the same size
      if (this->_le_der.rows() != this->_der_rat.rows() && this->_le_der.cols() != this->_der_rat.cols())
        throw formic::Exception("the input derivative vector matrices are of different size!");

      // if the approximate flag is set to be true, then throw out an exception
      if (approximate)
        throw formic::Exception(
            "Davidson solver's matrix-vector multiplication routine doesn't support approximate matrix!");

      // if we are doing ground state calculation
      if (this->_ground)
      {
        // temp vector to store le_der * x
        formic::ColVec<S> temp(this->_le_der.rows());

        // call blas level-2 function
        formic::xgemv('N', Ns, Nind, formic::unity(S()), &this->_le_der.at(0, 0), Ns, &x.at(0), 1, formic::zero(S()),
                      &temp.at(0), 1);

        // call blas level-2 function
        formic::xgemv('T', Ns, Nind, formic::unity(S()), &this->_der_rat.at(0, 0), Ns, &temp.at(0), 1,
                      formic::zero(S()), &y.at(0), 1);

        return;
      }

      else if (!this->_ground)
      {
        // temp vectpr to store omega * der_rat * x
        formic::ColVec<S> temp1(Ns);

        // temp vector to store le_der * x
        formic::ColVec<S> temp2(Ns);

        // call blas level-2 function
        formic::xgemv('N', Ns, Nind, formic::unity(S()) * this->_hd_shift, &this->_der_rat.at(0, 0), Ns, &x.at(0), 1,
                      formic::zero(S()), &temp1.at(0), 1);

        // call blas level-2 function
        formic::xgemv('N', Ns, Nind, formic::unity(S()), &this->_le_der.at(0, 0), Ns, &x.at(0), 1, formic::zero(S()),
                      &temp2.at(0), 1);

        // combine these two temp vector together
        temp1 -= temp2;

        // left multiply by _der_rat^T
        formic::xgemv('T', Ns, Nind, formic::unity(S()), &this->_der_rat.at(0, 0), Ns, &temp1.at(0), 1,
                      formic::zero(S()), &y.at(0), 1);

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

  void HMatMatOp(const formic::Matrix<double>& x,
                 formic::Matrix<double>& y,
                 const bool transpose,
                 const bool approximate = false) override
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

  void SMatVecOp(const formic::ColVec<S>& x, formic::ColVec<S>& y, const bool approximate = false) override
  {
    // size the resulting vector correctly
    y.reset(x.size());

    // if we have already built the matrix, we just multiply it
    if (_build_lm_matrix)
    {
      // if the approximate flag is set to be true, then throw out a exception
      if (approximate)
        throw formic::Exception(
            "Davidson solver's matrix-vector multiplication routine doesn't support approximate matrix!");

      // call level-2 blas function
      formic::xgemv('N', _smat.rows(), _smat.cols(), formic::unity(S()), &_smat.at(0, 0), _smat.rows(), &x.at(0), 1,
                    formic::zero(S()), &y.at(0), 1);

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

    else
    {
      // the number of samples on each process
      const int Ns = this->_le_der.rows();

      // the number of independent variables + 1
      const int Nind = this->_le_der.cols();

      // check to see whether derivative vector matrices have same size
      if (this->_le_der.rows() != this->_der_rat.rows() && this->_le_der.cols() != this->_der_rat.cols())
        throw formic::Exception("input derivative vectors are of different size");

      // if the approximate flag is set to be true, then throw out an exception
      if (approximate)
        throw formic::Exception(
            "Davidson solver's matrix-vector multiplication routine doesn't support approximate matrix!");

      // if we are doing ground state calculation
      if (this->_ground)
      {
        // temp vector to store der_rat * x
        formic::ColVec<S> temp(Ns);

        // call blas level-2 function
        formic::xgemv('N', Ns, Nind, formic::unity(S()), &this->_der_rat.at(0, 0), Ns, &x.at(0), 1, formic::zero(S()),
                      &temp.at(0), 1);

        // call blas level-2 function
        formic::xgemv('T', Ns, Nind, formic::unity(S()), &this->_der_rat.at(0, 0), Ns, &temp.at(0), 1,
                      formic::zero(S()), &y.at(0), 1);

        return;
      }

      else if (!this->_ground)
      {
        // temp vectpr to store omega * der_rat * x
        formic::ColVec<S> temp1(Ns);

        // temp vector to store le_der * x
        formic::ColVec<S> temp2(Ns);

        // temp vector that store omega * der_rat^T * (omega * der_rat - le_der) * x
        formic::ColVec<S> temp3(x.size());

        // call blas level-2 function
        formic::xgemv('N', Ns, Nind, formic::unity(S()) * this->_hd_shift, &this->_der_rat.at(0, 0), Ns, &x.at(0), 1,
                      formic::zero(S()), &temp1.at(0), 1);

        // call blas level-2 function
        formic::xgemv('N', Ns, Nind, formic::unity(S()), &this->_le_der.at(0, 0), Ns, &x.at(0), 1, formic::zero(S()),
                      &temp2.at(0), 1);

        // combine these two temp vector together
        temp1 -= temp2;

        // omega * D^T * (omega * D - L) * x
        formic::xgemv('T', Ns, Nind, formic::unity(S()) * this->_hd_shift, &this->_der_rat.at(0, 0), Ns, &temp1.at(0),
                      1, formic::zero(S()), &y.at(0), 1);

        // L^T * (omega * D - L) * x
        formic::xgemv('T', Ns, Nind, formic::unity(S()), &this->_le_der.at(0, 0), Ns, &temp1.at(0), 1,
                      formic::zero(S()), &temp3.at(0), 1);

        // (omega * D^T - L^T) * (omega * D - L) * x
        y -= temp3;

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

  void LSMatVecOp(const formic::ColVec<S>& x, formic::ColVec<S>& y, const bool approximate = false) override
  {
    // size the resulting vector correctly
    y.reset(x.size());

    // if we have already built the matrix, we just multiply it
    if (_build_lm_matrix)
    {
      // if the approximate flag is set to be true, then throw out a exception
      if (approximate)
        throw formic::Exception(
            "Davidson solver's matrix-vector multiplication routine doesn't support approximate matrix!");

      // call level-2 blas function
      formic::xgemv('N', _lmsmat.rows(), _lmsmat.cols(), formic::unity(S()), &_lmsmat.at(0, 0), _lmsmat.rows(),
                    &x.at(0), 1, formic::zero(S()), &y.at(0), 1);

      // return product vector
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

  void SMatMatOp(const formic::Matrix<double>& x, formic::Matrix<double>& y, const bool approximate = false) override {}

  //void MMatVecOp(const formic::ColVec<double> & x, formic::ColVec<double> & y, const int mmat_first);

  ////////////////////////////////////////////////////////////////////////////////////
  // \brief updates hamiltonian * krylov vector and hamiltonian projection based on
  //        new shift
  //
  //
  ////////////////////////////////////////////////////////////////////////////////////

  void update_hvecs_sub(const double new_i_shift, const double new_s_shift) override
  {
    // get rank number and number of ranks
    int my_rank = formic::mpi::rank();

    // get the different between new shift and old shift
    const double diff_shift_i = new_i_shift - this->_hshift_i;
    const double diff_shift_s = new_s_shift - this->_hshift_s;

    if (my_rank == 0)
    {
      // update "identity shift" for the hamiltonian product
      for (int j = 0; j < _nkry; j++)
      {
        for (int i = 1; i < this->_nfds; i++)
        {
          _hvecs.at(i, j) += diff_shift_i * _kvecs.at(i, j);
        }
      }

      // update "overlap shift" for the hamiltonian product
      for (int j = 1; j < _nkry; j++)
      {
        for (int i = 0; i < this->_nfds; i++)
        {
          _hvecs.at(i, j) += diff_shift_s * _svecs.at(i, j);
        }
      }

      // update projection of hamiltonian matrix
      _subH = _kvecs.c() * _hvecs;
    }
  }

  ////////////////////////////////////////////////////////////////////////////////////
  // \brief update Hamiltonian and overlap matrix for the solver
  //
  //
  //
  ////////////////////////////////////////////////////////////////////////////////////

  void update_hamovlp(formic::Matrix<S>& hmat, formic::Matrix<S>& smat)
  {
    // get rank number and number of ranks
    int my_rank = formic::mpi::rank();

    // check if the input matrices have the right size
    bool right_size =
        (hmat.rows() == this->_nfds && hmat.cols() == this->_nfds && smat.rows() == this->_nfds && smat.cols());

    if (my_rank == 0 && !right_size)
      throw formic::Exception("bad matrix size found in update hamiltonian and overlap function");

    else
    {
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

  void child_reset() override
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
#ifndef QMC_COMPLEX
    _sub_evec.reset(0);
#else
    _sub_evec_complex.reset(0);
#endif

// clear all values calculated from last solve
#ifndef QMC_COMPLEX
    _sub_eval = 0.0;
#else
    _sub_eval_complex = 0.0;
#endif
    _smallest_sin_value = 1e10;

    // set the eigenvector to initial guess
    this->_evecs.reset((this->_var_deps_use ? 1 + this->_dep_ptr->n_tot() : this->_nfds), formic::zero(S()));
    this->_evecs.at(0) = formic::unity(S());

    // set the independent eigenvector to initial guess
    this->_ind_evecs.reset((this->_var_deps_use ? 1 + this->_dep_ptr->n_ind() : this->_nfds), formic::zero(S()));
    this->_ind_evecs.at(0) = formic::unity(S());
  }

  /////////////////////////////////////////////////////////////////////////////////////
  // \brief solves the eigenvalue problem via the normal davidson method
  //
  //
  //
  /////////////////////////////////////////////////////////////////////////////////////

  bool iterative_solve(double& eval, std::ostream& output) override
  {
    // get rank number and number of ranks
    int my_rank = formic::mpi::rank();

    // initialize the solution vector to the unit vector along the first direction
    this->_evecs.reset((this->_var_deps_use ? 1 + this->_dep_ptr->n_tot() : this->_nfds), formic::zero(S()));
    this->_evecs.at(0) = formic::unity(S());

    // ensure that at least one vector is in the krylov subspace
    if (my_rank == 0 && _nkry == 0)
      throw formic::Exception(
          "Empty krylov subspace upon entry to iterative_solve. Did you forget to add the initial vector?");

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

    while (true)
    {
      // average of smallest singular value
      double smallest_sin_val_avg = 0.0;

      // solve subspace eigenvalue problem on root process
      if (my_rank == 0)
      {
        this->solve_subspace(true);
      }

      // send resulting eigenvalues to all processes and record it as the new best estimate(if we build the matrix)
      formic::mpi::bcast(&this->_cost, 1);
      eval = this->_cost;

      // check if the cost function has imaginary component and stop iterating when it does
      formic::mpi::bcast(&this->_eval_was_complex, 1);
      if (this->_eval_was_complex)
      {
        if (my_rank == 0)
          output << boost::format("davidson iteration %4i: stopping due to imaginary component in cost function") % iter
                 << std::endl;
        break;
      }

      // if cost function change is unreasonable, stop iterating and set bad solve flag
      if (std::abs(this->_cost - this->_init_cost) > this->_max_energy_change)
      {
        retval = false;
        if (my_rank == 0)
          output << boost::format("davidson iteration %4i stopping due to too large eigenvalue change") % iter
                 << std::endl;
        break;
      }

      // if the overlap matrix becomes singular, stop iterating
      formic::mpi::bcast(&_smallest_sin_value, 1);
      if (std::abs(_smallest_sin_value) < _singular_value_threshold)
      {
        if (my_rank == 0)
          output << boost::format("davidson iteration %4i stopping due to small subspace S singular value of %.2e") %
                  iter % _smallest_sin_value
                 << std::endl;
        break;
      }

      // construct new krylov vectors from subspace eigenvector
      if (my_rank == 0)
      {
        _wv6.reset(this->_nfds);
#ifndef QMC_COMPLEX
        _wv6 = _kvecs * _sub_evec;
#else
        _wv6 = _kvecs * _sub_evec_complex;
#endif

        // get the normalization factor
        const double temp_norm = std::sqrt(formic::real(_wv6.norm2()));

        // normalize this new vector
        _wv6 /= temp_norm;
        // add up linear combination of Hamiltonian and overlap products to make the new residual vector
        _wv1.reset(this->_nfds);

#ifndef QMC_COMPLEX
        _wv1 = _hvecs * _sub_evec;
        _wv1 = _wv1 - _sub_eval * _svecs * _sub_evec;
#else
        _wv1 = _hvecs * _sub_evec_complex;
        _wv1 = _wv1 - _sub_eval_complex * _svecs * _sub_evec_complex;
#endif
      }
      // if we don't build matrix, then send this new vector to all processes
      if (!_build_lm_matrix)
        formic::mpi::bcast(&_wv1.at(0), _wv1.size());

      // compute the residual norm and send it to all processes
      double current_residual;
      double current_residual_avg;

#ifdef QMC_COMPLEX
      current_residual = std::sqrt(formic::real(_wv1.norm2()));
#else
      current_residual = std::sqrt(_wv1.norm2());
#endif
      formic::mpi::bcast(&current_residual, 1);

      // if this is the best residual, save it and save the new eigenvector estimate
      if (my_rank == 0 && current_residual < _best_residual)
      {
        _best_residual = current_residual;

        // get current eigenvector estimate, which corresponds to the set of independent variables
        this->_ind_evecs = _wv6;

        // if our actual variables are dependent on the set of independent variables worked with here, expand the eigenvector into the full set of variables
        if (this->_var_deps_use)
        {
          // size the eigenvector correctly
          this->_evecs.reset(1 + this->_dep_ptr->n_tot(), formic::zero(S()));
          this->_evecs.at(0) = this->_ind_evecs.at(0);

          // get some temporary vectors
          formic::ColVec<S> _evec_temp(this->_dep_ptr->n_tot());
          formic::ColVec<S> _ind_temp(this->_dep_ptr->n_ind());
          for (int i = 0; i < _ind_temp.size(); i++)
          {
            _ind_temp.at(i) = this->_ind_evecs.at(i + 1);
          }
          this->_dep_ptr->expand_ind_to_all(&_ind_temp.at(0), &_evec_temp.at(0));

          for (int i = 0; i < _evec_temp.size(); i++)
          {
            this->_evecs.at(i + 1) = _evec_temp.at(i);
          }
        }

        // otherwise just copy the eigenvector into output since the independent and total variable sets are the same
        else
        {
          this->_evecs = this->_ind_evecs.clone();
        }
      }

      // print itertion result
      if (my_rank == 0)
      {
        // if we are doing ground state calculation, then print out energy
        if (this->_ground)
          output << boost::format("davidson iteration %4i:   krylov dim = %3i   energy = %20.12f       residual = %.2e "
                                  "          smallest_sin_value = %.2e") %
                  iter % _nkry % this->_cost % current_residual % _smallest_sin_value
                 << std::endl;

        // if we are doing excited state calculation, then print out target function value
        else
          output << boost::format("davidson iteration %4i:   krylov dim = %3i   tar_fn = %20.12f       residual = %.2e "
                                  "          smallest_sin_value = %.2e") %
                  iter % _nkry % this->_cost % current_residual % _smallest_sin_value
                 << std::endl;
      }

      // check for convergence
      converged = current_residual < this->_residual_threshold;

      // if we build matrix and iteration has already converged, we exit iteration
      if (converged)
        break;

      // if iteration hasn't converged, we increment the iteration count and stop if maximum number of iterations has been reached
      if (iter++ >= _n_max_iter)
        break;

      // add the new krylov basis vector
      this->add_krylov_vector(_wv1);
    }

    // print iteration information
    if (converged && my_rank == 0)
      output << boost::format("davidson solver converged in %10i iterations") % iter << std::endl << std::endl;

    else if (my_rank == 0)
      output << boost::format("davidson solver did not converge after %10i iterations") % iter << std::endl
             << std::endl;

    return retval;
  }

  ///////////////////////////////////////////////////////////////////////////////////////////////
  // \brief solves the eigenvalue problem
  //
  //
  //
  ///////////////////////////////////////////////////////////////////////////////////////////////

  bool solve(double& eval, std::ostream& output) override
  {
    // get rank number and number of ranks
    int my_rank = formic::mpi::rank();

    if (!this->_variance_correct)
      return this->iterative_solve(eval, output);

    else
    {
      // first evaluate tau and correct hamiltonian matrix
      this->tau_and_correct_ham();

      // bcast tau
      formic::mpi::bcast(&this->_tau, 1);

      // flag to tell whether tau is converged
      bool tau_converged = false;

      // number of iteration of tau
      int iter_tau = 0;

      // return value
      bool retval = false;

      // reset eigensolver
      this->child_reset();

      // set the current energy as initial wavefunction's energy
      this->_cost = this->_init_cost;

      // add initial guess
      // get the dimension of matrix
      const int N = _hmat.cols();
      formic::ColVec<S> temp(N);
      for (int j = 0; j < temp.size(); j++)
        temp.at(j) = (j == 0 ? formic::unity(S()) : formic::zero(S()));
      this->add_krylov_vector(temp);

      // iteratively solve the generalized eigenvalue problem
      retval = this->iterative_solve(eval, output);

      // return solve result
      return retval;
    }
  }

  ///////////////////////////////////////////////////////////////////////////////////////////////
  // \brief converts eigenvectors into wave function coefficients
  //        solving this question Hc = Sd for d, S is the ground overlap matrix
  //
  //
  ///////////////////////////////////////////////////////////////////////////////////////////////

  void child_convert_to_wf_coeff() override { return; }
};

} // namespace engine

} // namespace cqmc

#endif
