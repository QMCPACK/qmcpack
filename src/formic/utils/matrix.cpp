///////////////////////////////////////////////////////////////////////////////////////////////////
/// \file formic/utils/matrix.cpp
///
/// \brief   implementation file for the formic::Matrix class
///
///////////////////////////////////////////////////////////////////////////////////////////////////

#include<cmath>
#include<complex>
#include<stack>
#include<numeric>
//#include<iostream>

#include<formic/utils/matrix.h>
#include<formic/utils/numeric.h>
#include<formic/utils/random.h>

///////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief Computes and returns the exponential of a matrix by evaluating its tailor series until
///        all elements of the next term in the series are smaller than the provided threshold.
///
/// \param[in]      m        the matrix to exponentiate
/// \param[in]      tol      the tolerance below which further terms in the series are ignored
///
///////////////////////////////////////////////////////////////////////////////////////////////////
template<class S>
formic::Matrix<S> formic::matrix_exponent(const formic::ConstMatrix<S> & m, const double tol) {

  // check tolerance sanity
  if ( tol <= 0.0 )
    throw formic::Exception("tolerance must be positive in formic::matrix_exponent");

  // check for square matrix
  if ( m.rows() != m.cols() )
    throw formic::Exception("matrix must be square in formic::matrix_exponent");

  // initialize the return value to the identity matrix
  formic::Matrix<S> r = formic::identity_matrix<S>(m.rows());

  // get matrix to hold each tailor series term
  formic::Matrix<S> p = r.clone();

  // compute the tailor series of the matrix exponential until
  // the element-wise tolerance is reached
  double x = 1.0;
  for (double max_diff = 2.0 * tol; max_diff > tol; x += 1.0) {

    // get the next term in the exponent tailor series
    p = m * p;
    p *= formic::unity(S()) / x;

    // add the term to the total
    r += p;

    // get the maximum magnitude of the elements of p
    max_diff = -1.0;
    for (int j = 0; j < p.cols(); j++)
    for (int i = 0; i < p.rows(); i++)
      max_diff = std::max(std::abs(p(i,j)), max_diff);

    //std::cout << "matrix exponent iter " << x << "    max_diff = " << max_diff << std::endl;

  }
  //std::cout << std::endl;

  // return the exponentiated matrix
  return r;

}

template formic::Matrix<double> formic::matrix_exponent(const formic::ConstMatrix<double> &, const double);
template formic::Matrix<std::complex<double> > formic::matrix_exponent(const formic::ConstMatrix<std::complex<double> > &, const double);

///////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief Computes and returns the derivative adjoint of a matrix for the matrix
///        exponential function given the original matrix and the derivative adjoints 
///        of the exponentiated matrix.
///
/// \param[in]      m        the original (not-exponentiated) matrix
/// \param[in]      a        the derivative adjoints of the exponentiated matrix
/// \param[in]      tol      the tolerance below which further terms in the series are ignored
///
///////////////////////////////////////////////////////////////////////////////////////////////////
template<class S>
formic::Matrix<S> formic::matrix_exponent_der_adj(const formic::ConstMatrix<S> & m, const formic::ConstMatrix<S> & a, const double tol) {

  //return a.clone();
  //return a + 0.5 * formic::unity(S()) * a * m.t() + 0.5 * formic::unity(S()) * m.t() * a;

  // check tolerance sanity
  if ( tol <= 0.0 )
    throw formic::Exception("tolerance must be positive in formic::matrix_exponent_der_adj");

  // check for square matrices
  if ( m.rows() != m.cols() )
    throw formic::Exception("original matrix must be square in formic::matrix_exponent_der_adj");
  if ( a.rows() != a.cols() )
    throw formic::Exception("exponentiated matrix's derivative adjoint must be square in formic::matrix_exponent_der_adj");

  // check that matrices are the same size
  if ( m.rows() != a.rows() )
    throw formic::Exception("original matrix must be the same size as the exponentiated matrix's derivative adjoint in formic::matrix_exponent_der_adj");

  // get matrix to hold each tailor series term
  formic::Matrix<S> p = formic::identity_matrix<S>(m.rows());

  // get container to hold tailor series history
  std::stack<formic::ConstMatrix<S> > pstack;

  // compute the tailor series terms until
  // the element-wise tolerance is reached
  double x = 1.0;
  for (double max_diff = 2.0 * tol; max_diff > tol; x += 1.0) {

    // save the tailor series term
    pstack.push(p);

    // get the next term in the exponent tailor series
    p = m * p;
    p *= formic::unity(S()) / x;

    // get the maximum magnitude of the elements of p
    max_diff = -1.0;
    for (int j = 0; j < p.cols(); j++)
    for (int i = 0; i < p.rows(); i++)
      max_diff = std::max(std::abs(p(i,j)), max_diff);

  }

  // the matrix p will now hold the adjoint of each tailor series term
  p <<= a;

  // initialize the return value to zero
  formic::Matrix<S> r(m.rows(), m.cols(), formic::zero(S()));

  // get the adjoint contribution from each term in the tailor series
  x -= 1.0;
  for ( ; x > 0.5; x -= 1.0) {

    // make sure we didn't mess up
    if ( pstack.empty() )
      throw formic::Exception("unexpectedly empty pstack in formic::matrix_exponent_der_adj");

    // add this tailor series term's contribution to the overall derivative
    r += ( p * pstack.top().t() ) * ( formic::unity(S()) / x );

    // discard the intermediate we no longer need
    pstack.pop();

    // update the accumulation of the intermediates' adjoints
    p = a + ( formic::unity(S()) / x ) * m.t() * p;

  }

  // make sure we didn't mess up
  if ( ! pstack.empty() )
    throw formic::Exception("unexpectedly non-empty pstack in formic::matrix_exponent_der_adj");

  // return the derivative adjoint
  return r;

}

template formic::Matrix<double> formic::matrix_exponent_der_adj(const formic::ConstMatrix<double> &,
                                                                const formic::ConstMatrix<double> &,
                                                                const double);
template formic::Matrix<std::complex<double> > formic::matrix_exponent_der_adj(const formic::ConstMatrix<std::complex<double> > &,
                                                                               const formic::ConstMatrix<std::complex<double> > &,
                                                                               const double);

///////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief Computes the singular value decomposition of the matrix.
///        A = U s VT
///
/// \param[out]     u        on exit, the left singular vectors
/// \param[out]     s        on exit, the singular values in descending order
/// \param[out]     vt       on exit, the conjugate transpose of the right singular vectors
///
///////////////////////////////////////////////////////////////////////////////////////////////////
template<class S> void formic::ConstMatrix<S>::svd(formic::Matrix<S> & u, formic::ColVec<S> & s, formic::Matrix<S> & vt) const {

  // check that the matrix is not empty
  if ( this->size() == 0 )
    throw formic::Exception("cannot take the svd of a %i by %i matrix") % m_n % m_m;

  // get large and small dimensions
  const size_t sd = std::min(m_n, m_m);
  const size_t ld = std::max(m_n, m_m);

  // size the output vectors
  u.reset(m_n, m_n, formic::zero(S()));
  s.reset(sd, formic::zero(S()));
  vt.reset(m_m, m_m, formic::zero(S()));

  // get a copy of the matrix
  std::vector<S> mat(this->size());
  *this >>= &mat.at(0);

  // get space to hold singular values in real form
  std::vector<double> t(sd, 0.0);

  // get workspace
  const size_t lwork = 30 * ld;
  std::vector<S> work(lwork, formic::zero(S()));
  std::vector<double> rwork(5 * sd, 0.0);

  // perform SVD
  int info = 0;
  formic::xgesvd('A', 'A', m_n, m_m, &mat.at(0), m_n, &t.at(0), u.begin(), m_n, vt.begin(), m_m, &work.at(0), lwork, &rwork.at(0), info);
  if ( info != 0 )
    throw formic::Exception("xgesvd failed with info = %i in formic::ConstMatrix::svd") % info;

  // convert the real singular values to type S
  for (int i = 0; i < sd; i++)
    s.at(i) = formic::unity(S()) * t.at(i);

}

template void formic::ConstMatrix<double>::svd(formic::Matrix<double> &, formic::ColVec<double> &, formic::Matrix<double> &) const;
template void formic::ConstMatrix<std::complex<double> >::svd(formic::Matrix<std::complex<double> > &, formic::ColVec<std::complex<double> > &, formic::Matrix<std::complex<double> > &) const;

///////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief Returns a new matrix that is the inverse or pseudo-inverse of this matrix
///
/// \param[in]      thresh   Threshold below which singular values are treated as zeros.
///                          If any zeros are encountered, the pseudo-inverse is computed instead.
///
///////////////////////////////////////////////////////////////////////////////////////////////////
template<class S> formic::Matrix<S> formic::ConstMatrix<S>::inv(const double thresh) const {

  // ensure matrix is square and not empty
  if ( m_n != m_m || this->size() == 0 )
    throw formic::Exception("cannot invert a %i by %i matrix") % m_n % m_m;

  // compute SVD
  formic::Matrix<S> u, vt;
  formic::ColVec<S> s;
  this->svd(u, s, vt);

  // compute and return inverse
  vt.cip();
  for (size_t i = 0; i < vt.cols(); i++)
    vt.scale_col_by(i, formic::unity(S()) * ( formic::real(s.at(i)) > thresh ? 1.0 / formic::real(s.at(i)) : 0.0 ));
  return vt * u.cip();

}

template formic::Matrix<double> formic::ConstMatrix<double>::inv(const double) const;
template formic::Matrix<std::complex<double> > formic::ConstMatrix<std::complex<double> >::inv(const double) const;

///////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief Computes the eigenvalues and right eigenvectors of the matrix, which is not assumed
///        to be symmetric.
///
/// \param[out]     w        on exit, the eigenvalues
/// \param[out]     v        on exit, the right eigenvectors are stored in the columns of v
///
///////////////////////////////////////////////////////////////////////////////////////////////////
template<class S> void formic::ConstMatrix<S>::nonsym_eig(formic::ColVec<std::complex<double> > & w, formic::Matrix<std::complex<double> > & v) const {

  // ensure matrix is square and not empty
  if ( m_n != m_m || this->size() == 0 )
    throw formic::Exception("formic::ConstMatrix::nonsym_eig cannot get eigenvalues and eigenvectors of a %i by %i matrix") % m_n % m_m;

  // size the output vectors
  w.reset(m_n, formic::zero(S()));
  v.reset(m_n, m_m, formic::zero(S()));

  // one and zero in complex form
  const std::complex<double> one = formic::unity(std::complex<double>());
  const std::complex<double> zero = formic::zero(std::complex<double>());

  // create arrays and other variables used in diagonalization routines
  const size_t lwork = 30 * m_n;
  std::vector<std::complex<double> > work(lwork, zero);
  std::vector<double> rwork(2 * m_n, 0.0);

  // get a copy of the matrix in complex form
  std::vector<std::complex<double> > mat(this->size(), zero);
  auto itthis = this->begin();
  std::for_each(mat.begin(), mat.end(), [&] (std::complex<double> & a) { a = one * (*itthis++); });

  // compute the eigenvalues and right eigenvectors
  int info = 0;
  formic::zgeev('N', 'V', m_n, &mat.at(0), m_n, w.begin(), v.begin(), m_n, v.begin(), m_n, &work.at(0), lwork, &rwork.at(0), info);
  if ( info != 0 )
    throw formic::Exception("formic::zgeev failed formic::ConstMatrix::nonsym_eig") % info;

}

template void formic::ConstMatrix<double>::nonsym_eig(formic::ColVec<std::complex<double> > & w, formic::Matrix<std::complex<double> > & v) const;
template void formic::ConstMatrix<std::complex<double> >::nonsym_eig(formic::ColVec<std::complex<double> > & w, formic::Matrix<std::complex<double> > & v) const;

namespace formic {
  
  ////////////////////////////////////////////////////////////////////////////////////////////////
  /// \brief Computes the eigenvalues and eigenvectors of the real symmetric matrix
  ///
  /// \param[out]   w       on exit, the eigenvalues in ascending order
  /// \param[ut]    v       on exit, the eigenvectors are stored in the columns of v
  ///
  ////////////////////////////////////////////////////////////////////////////////////////////////
  template<> void ConstMatrix<double>::sym_eig(formic::ColVec<double> & w, formic::Matrix<double> & v) const { 
    
    // ensure matrix is square and not empty
    if ( m_n != m_m || this->size() == 0 ) 
      throw formic::Exception("formic::ConstMatrix::sym_eig cannot get eigenvalues and eigenvectors of a %i by %i matrix") % m_n % m_m;

    // size the output eigenvalue vector
    w.reset(m_n, 0.0);

    // copy the matrix into the space in which it will be diagonalized into eigenvectors
    v <<= *this;

    // create arrays and other variables used in diagonalization routines
    const size_t lwork = 10 * m_n;
    std::vector<double> work(lwork,0.0);

    // compute the eigenvalues and eigenvectors
    int info = 0;
    formic::dsyev('V', 'U', m_n, v.begin(), m_n, w.begin(), &work.at(0), lwork, info);
    if ( info != 0 ) 
      throw formic::Exception("formic::dsyev failed formic::ConstMatrix::sym_eig") % info;

  }

  ////////////////////////////////////////////////////////////////////////////////////////////////
  /// \brief Computes the eigenvalues and eigenvectors of the Hermitian matrix
  ///
  /// \param[out]   w       on exit, the eigenvalues in ascending order
  /// \param[ut]    v       on exit, the eigenvectors are stored in the columns of v
  ///
  ////////////////////////////////////////////////////////////////////////////////////////////////
  template<> void ConstMatrix<std::complex<double> >::sym_eig(formic::ColVec<double> & w, formic::Matrix<std::complex<double> > & v) const { 
    
    // ensure matrix is square and not empty
    if ( m_n != m_m || this->size() == 0 ) 
      throw formic::Exception("formic::ConstMatrix::sym_eig cannot get eigenvalues and eigenvectors of a %i by %i matrix") % m_n % m_m;

    // size the output eigenvalue vector
    w.reset(m_n, 0.0);

    // copy the matrix into the space in which it will be diagonalized into eigenvectors
    v <<= *this;

    // create arrays and other variables used in diagonalization routines
    const size_t lwork = 10 * m_n;
    std::vector<std::complex<double> > work(lwork, formic::zero(std::complex<double>()));
    std::vector<double> rwork(3*m_n, 0.0);

    // compute the eigenvalues and eigenvectors
    int info = 0;
    formic::zheev('V', 'U', m_n, v.begin(), m_n, w.begin(), &work.at(0), lwork, &rwork.at(0), info);
    if ( info != 0 ) 
      throw formic::Exception("formic::zheev failed formic::ConstMatrix::sym_eig") % info;

  }
}
///////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief Returns an iterator to the element with the largest absolute value
///
///////////////////////////////////////////////////////////////////////////////////////////////////
template<class S> const S * formic::ConstMatrix<S>::max_abs_elem_iter() const {

  // invalid if we are empty
  if ( this->size() == 0 )
    throw formic::Exception("formic::ConstMatrix::max_abs_elem_iter encountered an illegal matrix dimension of %i by %i") % m_n % m_m;

  // return an iterator to the element with the largest absolute value
  return std::accumulate(this->begin(),
                         this->end(),
                         this->begin(),
                         [] (const S * p, const S & a) { return ( std::abs(*p) > std::abs(a) ? p : &a ); });

}

template const int * formic::ConstMatrix<int>::max_abs_elem_iter() const;
template const double * formic::ConstMatrix<double>::max_abs_elem_iter() const;
template const std::complex<double> * formic::ConstMatrix<std::complex<double> >::max_abs_elem_iter() const;

///////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief Returns an iterator to the element with the smallest absolute value
///
///////////////////////////////////////////////////////////////////////////////////////////////////
template<class S> const S * formic::ConstMatrix<S>::min_abs_elem_iter() const {

  // invalid if we are empty
  if ( this->size() == 0 )
    throw formic::Exception("formic::ConstMatrix::min_abs_elem_iter encountered an illegal matrix dimension of %i by %i") % m_n % m_m;

  // return an iterator to the element with the smallest absolute value
  return std::accumulate(this->begin(),
                         this->end(),
                         this->begin(),
                         [] (const S * p, const S & a) { return ( std::abs(*p) < std::abs(a) ? p : &a ); });

}

template const int * formic::ConstMatrix<int>::min_abs_elem_iter() const;
template const double * formic::ConstMatrix<double>::min_abs_elem_iter() const;
template const std::complex<double> * formic::ConstMatrix<std::complex<double> >::min_abs_elem_iter() const;

///////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief  returns a new matrix in which each element is the complex<double> equivalent of this
///         matrix's corresponding element
///
///////////////////////////////////////////////////////////////////////////////////////////////////
template<class S> formic::Matrix<std::complex<double> > formic::ConstMatrix<S>::get_complex_form() const {
  formic::Matrix<std::complex<double> > retval(this->rows(), this->cols());
  auto y = retval.begin();
  std::for_each(this->begin(),
                this->end(),
                [&] (const S & a) {
                  *y++ =   double(formic::real(a)) * std::complex<double>(1.0, 0.0)
                         + double(formic::imag(a)) * std::complex<double>(0.0, 1.0);
                }
               );
  return retval;
}

template formic::Matrix<std::complex<double> > formic::ConstMatrix<int>::get_complex_form() const;
template formic::Matrix<std::complex<double> > formic::ConstMatrix<double>::get_complex_form() const;
template formic::Matrix<std::complex<double> > formic::ConstMatrix<std::complex<double> >::get_complex_form() const;

///////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief  returns a new column vector whose elements are copied from the specified column
///         of this matrix
///
/// \param[in]      j        index of the desired column
///
///////////////////////////////////////////////////////////////////////////////////////////////////
template<class S> formic::ColVec<S> formic::ConstMatrix<S>::col_as_vec(const size_t j) const {
  if ( j >= m_m )
    throw formic::Exception("formic::ConstMatrix::col_as_vec received out-of-bound index %i for a matrix of dimension %i by %i") % j % m_n % m_m;
  formic::ColVec<S> retval(m_n);
  std::copy(this->col_begin(j), this->col_end(j), retval.begin());
  return retval;
}

template formic::ColVec<int> formic::ConstMatrix<int>::col_as_vec(const size_t j) const;
template formic::ColVec<double> formic::ConstMatrix<double>::col_as_vec(const size_t j) const;
template formic::ColVec<std::complex<double> > formic::ConstMatrix<std::complex<double> >::col_as_vec(const size_t j) const;

///////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief  returns a new row vector whose elements are copied from the specified row
///         of this matrix
///
/// \param[in]      i        index of the desired row
///
///////////////////////////////////////////////////////////////////////////////////////////////////
template<class S> formic::RowVec<S> formic::ConstMatrix<S>::row_as_vec(const size_t i) const {
  if ( i >= m_n )
    throw formic::Exception("formic::ConstMatrix::row_as_vec received out-of-bound index %i for a matrix of dimension %i by %i") % i % m_n % m_m;
  formic::RowVec<S> retval(m_m);
  std::copy(this->row_begin(i), this->row_end(i), retval.begin());
  return retval;
}

template formic::RowVec<int> formic::ConstMatrix<int>::row_as_vec(const size_t i) const;
template formic::RowVec<double> formic::ConstMatrix<double>::row_as_vec(const size_t i) const;
template formic::RowVec<std::complex<double> > formic::ConstMatrix<std::complex<double> >::row_as_vec(const size_t i) const;

///////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief  places the elements of the supplied vector in the specified column
///
/// \param[in]      j        index of the desired column
/// \param[in]      v        the vector
///
///////////////////////////////////////////////////////////////////////////////////////////////////
template<class S> void formic::Matrix<S>::put_vec_in_col(const size_t j, const formic::ConstMatrix<S> & v) {
  if ( j >= m_m )
    throw formic::Exception("formic::Matrix::put_vec_in_col received out-of-bound index %i for a matrix of dimension %i by %i") % j % m_n % m_m;
  if ( v.size() != m_n || ( v.rows() != 1 && v.cols() != 1 ) )
    throw formic::Exception("formic::Matrix::put_vec_in_col cannot put a vector of dimension %i by %i in the column of a matrix of dimension %i by %i")
          % v.rows() % v.cols() % m_n % m_m;
  std::copy(v.begin(), v.end(), this->col_begin(j));
}

template void formic::Matrix<int>::put_vec_in_col(const size_t j, const formic::ConstMatrix<int> & v);
template void formic::Matrix<double>::put_vec_in_col(const size_t j, const formic::ConstMatrix<double> & v);
template void formic::Matrix<std::complex<double> >::put_vec_in_col(const size_t j, const formic::ConstMatrix<std::complex<double> > & v);

///////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief  places the elements of the supplied vector in the specified row
///
/// \param[in]      i        index of the desired row
/// \param[in]      v        the vector
///
///////////////////////////////////////////////////////////////////////////////////////////////////
template<class S> void formic::Matrix<S>::put_vec_in_row(const size_t i, const formic::ConstMatrix<S> & v) {
  if ( i >= m_n )
    throw formic::Exception("formic::Matrix::put_vec_in_row received out-of-bound index %i for a matrix of dimension %i by %i") % i % m_n % m_m;
  if ( v.size() != m_m || ( v.rows() != 1 && v.cols() != 1 ) )
    throw formic::Exception("formic::Matrix::put_vec_in_row cannot put a vector of dimension %i by %i in the row of a matrix of dimension %i by %i")
          % v.rows() % v.cols() % m_n % m_m;
  std::copy(v.begin(), v.end(), this->row_begin(i));
}

template void formic::Matrix<int>::put_vec_in_row(const size_t i, const formic::ConstMatrix<int> & v);
template void formic::Matrix<double>::put_vec_in_row(const size_t i, const formic::ConstMatrix<double> & v);
template void formic::Matrix<std::complex<double> >::put_vec_in_row(const size_t i, const formic::ConstMatrix<std::complex<double> > & v);

///////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief  return the element by element dot product of two matrices in which the
///         complex conjugates of the first matrix's elements are used
///
/// \param[in]      mat1     the 1st matrix
/// \param[in]      mat2     the 2nd matrix
///
///////////////////////////////////////////////////////////////////////////////////////////////////
template<class S> S formic::dotc(const formic::ConstMatrix<S> & mat1, const formic::ConstMatrix<S> & mat2) {
  if ( mat1.dimensions_differ_from(mat2) )
    throw formic::Exception("cannot take dot product of a %i by %i matrix and a %i by %i matrix in formic::dotc");
  const S * v = mat1.begin();
  return std::accumulate(mat2.begin(), mat2.end(), formic::zero(S()), [&] (S init, S b) { return init + formic::conj(*v++) * b; } );
}

template int formic::dotc(const formic::ConstMatrix<int> & mat1, const formic::ConstMatrix<int> & mat2);
template double formic::dotc(const formic::ConstMatrix<double> & mat1, const formic::ConstMatrix<double> & mat2);
template std::complex<double> formic::dotc(const formic::ConstMatrix<std::complex<double> > & mat1, const formic::ConstMatrix<std::complex<double> > & mat2);

///////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief  Fills the matrix with uniformly distributed random numbers
/// 
/// \param[in]    lb      lower bound for random number range (default is -1.0)
/// \param[in]    ub      upper bound for random number range (default is  1.0)
///
///////////////////////////////////////////////////////////////////////////////////////////////////
template<class S> void formic::Matrix<S>::random_fill(const double lb, const double ub) {
  if ( ub < lb )
    throw formic::Exception("ub was less than lb in formic::Matrix::random_fill");
  for (auto it = this->begin(); it != this->end(); it++) 
    *it = formic::unity(S()) * ( lb + formic::random_number<double>() * ( ub - lb) ) + formic::imaginary_unity(S()) * ( lb + formic::random_number<double>() * ( ub - lb ) );
}
template void formic::Matrix<double>::random_fill(const double lb, const double ub);
template void formic::Matrix<std::complex<double> >::random_fill(const double lb, const double ub);
