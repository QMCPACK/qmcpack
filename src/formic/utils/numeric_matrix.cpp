///////////////////////////////////////////////////////////////////////////////////////////////////
/// \file formic/utils/numeric_matrix.cpp
///
/// \brief   implementation for some matrix mathematics functions
///
///////////////////////////////////////////////////////////////////////////////////////////////////

#include<algorithm>
#include<complex>
#include<cstdlib>
#include<cmath>
#include<vector>

#include<formic/utils/numeric.h>
#include<formic/utils/lapack_interface.h>
#include<formic/utils/exception.h>

///////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief   Computes a matrix's inverse and determinant using an LU decomposition.
///
/// \param[in]     n                       the dimension of the matrix A
/// \param[out]    det                     on exit, the matrix's determinant
/// \param[in,out] A          size n*n     On entry, the matrix to be inverted.
///                                        On exit, the inverse.
/// \param[out]    work       size n*n     work space for computing the inverse
/// \param[out]    iwork      size 2*n     integer work space
///
///////////////////////////////////////////////////////////////////////////////////////////////////
template <class S> void formic::matrix_inverse_lu(const int n,
                                                  S & det,
                                                  S * const A,
                                                  S * const work,
                                                  int * const iwork) {

  // use some of the integer work space for pivoting during the LU decomposition
  int * const ipiv = iwork;

  // use some of the integer work space for inverting the LU pivoting
  int * const inv_perm = iwork + n;

  // copy the matrix into the workspace
  formic::xcopy(n*n, A, 1, work, 1);

  // compute the LU decomposition of the matrix
  int info;
  formic::xgetrf(n, n, work, n, ipiv, info);
  if (info < 0)
    throw formic::Exception("formic::xgetrf failed with error code %i in formic::matrix_inverse_lu") % info;
  if (info > 0)
    throw formic::Exception("U(%i,%i) is exactly zero in formic::matrix_inverse_lu, thus the matrix is not invertable") % info % info;

  // since FORTRAN counts from 1, subtract 1 from each of the pivot array values
  for (int i = 0; i < n; i++)
    ipiv[i] -= 1;

  // compute the determinant
  det = formic::unity(S());
  for (int i = 0; i < n; i++)
    if ( i != ipiv[i] )
      det *= - work[i*n+i];
    else
      det *=   work[i*n+i];
  //if (std::abs(det) < 1.0e-12)
  //  throw pcps::SingularMatrixException( (boost::format("matrix determinant is zero in pcps::matrix_inverse_lu, thus the matrix is not invertable")).str() );

  // determine the inverse of the pivot array's permutation
  for (int i = 0; i < n; i++)
    inv_perm[i] = i;
  for (int i = n-1; i >= 0; i--)
    std::swap(inv_perm[i], inv_perm[ipiv[i]]);

  // compute the inverse of the upper triangle matrix
  for (int i = 0; i < n; i++)
  for (int j = i; j < n; j++) {
    A[i*n+j] = ( i == j ? formic::unity(S()) : formic::zero(S()) );
    for (int k = i; k < j; k++)
      A[i*n+j] -= A[i*n+k] * work[j*n+k];  // work indices transposed due to FORTRAN
    A[i*n+j] /= work[j*n+j];
  }

  // compute the inverse of the lower triangle matrix (whose diagonal elements are one and are not stored)
  for (int i = 1; i < n; i++)
  for (int j = i-1; j >= 0; j--) {
    A[i*n+j] = -work[j*n+i];  // work indices transposed due to FORTRAN
    for (int k = j+1; k < i; k++)
      A[i*n+j] -= A[k*n+j] * work[k*n+i];  // work indices transposed due to FORTRAN
  }

  // compute the product of the inverse triangular matrices
  {
    const S * const U_inv = A;
    const S * const L_inv = A;
    for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++) {
      // term with diagonal of L_inv (which equals one)
      work[i*n+j] = ( i <= j ? U_inv[i*n+j] : formic::zero(S()) );
      // other terms
      for (int k = std::max(i, j+1); k < n; k++)
        work[i*n+j] += U_inv[i*n+k] * L_inv[k*n+j];
    }
  }

  // Multiply by the inverse of the permutation matrix (column reordering)
  // We also transpose the final matrix to account for FORTRAN indexing
  for (int i = 0; i < n; i++)
  for (int j = 0; j < n; j++)
    A[j*n+i] = work[i*n+inv_perm[j]];

}

// explicitly instantiate the function
template void formic::matrix_inverse_lu(const int n,
                                        double & det,
                                        double * const A,
                                        double * const work,
                                        int * const iwork);
template void formic::matrix_inverse_lu(const int n,
                                        std::complex<double> & det,
                                        std::complex<double> * const A,
                                        std::complex<double> * const work,
                                        int * const iwork);

///////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief   computes the determinant of a matrix using LU decomposition
///
/// \param[in]       n       The matrix dimension.
/// \param[in,out]   A       On entry, the matrix whose determinant we want.  Overwritten on exit. Size n*n.
/// \param[out]      iwork   Work space for pivoting during the LU decomposition, size n.
///
/// \return the determinant in mantissa/exponent format
///
///////////////////////////////////////////////////////////////////////////////////////////////////
template <class S> formic::MantExp<S> formic::matrix_determinant_lu(const int n, S * const A, int * const iwork) {

  // compute the LU decomposition of the matrix
  int info;
  formic::xgetrf(n, n, A, n, iwork, info);
  if (info < 0)
    throw formic::Exception("formic::xgetrf failed with error code %i in formic::matrix_determinant_lu") % info;

  // since FORTRAN counts from 1, subtract 1 from each of the pivot array values
  for (int i = 0; i < n; i++)
    iwork[i] -= 1;

  // compute the determinant
  formic::MantExp<S> det(formic::unity(S()));
  for (int i = 0; i < n; i++)
    if ( i != iwork[i] )
      det *= - A[i*n+i];
    else
      det *=   A[i*n+i];

  // return the determinant value
  return det;

}

// explicitly instantiate the function
template formic::MantExp<double               > formic::matrix_determinant_lu(const int, double               * const, int * const);
template formic::MantExp<std::complex<double> > formic::matrix_determinant_lu(const int, std::complex<double> * const, int * const);


///////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief   Computes the hermitian matrix A^(-1/2) for a hermitian, positive definite matrix A
///
/// \param[in]       n                     the dimension of the matrix A
/// \param[in,out]   A        size n*n     On entry, the matrix to be negative halved.
///                                        On exit, A^(-1/2).
///
///////////////////////////////////////////////////////////////////////////////////////////////////
template <class S> void formic::matrix_neg_half(const int n, S * const A) {

  int info;
  const int lwork = std::max(5,n) * n;
  std::vector<S> work(lwork, formic::zero(S()));
  std::vector<S> work2(lwork, formic::zero(S()));
  std::vector<double> rwork(3*n, 0.0);
  std::vector<double> eval(n, 0.0);

  // diagonalize the matrix, after which the eigenvectors will be in the rows the matrix
  formic::xsyev('V', 'U', n, A, n, &eval.at(0), &work.at(0), lwork, &rwork.at(0), info);
  if (info != 0)
    throw formic::Exception("formic::xsyev failed with error code %i in formic::matrix_neg_half") % info;

  // compute the A^(-1/2) matrix
  for (int i = 0; i < n; i++) {
    const double x = 1.0 / std::sqrt(eval[i]);
    for (int j = 0; j < n; j++) {
      work[i*n+j] = A[j*n+i];
      work2[i*n+j] = A[i*n+j] * x;
    }
  }
  formic::xgemm('N', 'N', n, n, n, formic::unity(S()), &work2.at(0), n, &work.at(0), n, formic::zero(S()), A, n);

}

template void formic::matrix_neg_half(const int n, double * const A);
template void formic::matrix_neg_half(const int n, std::complex<double> * const A);

///////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief   Computes the eigenvalues and right eigenvectors for the eigensystem H x = lambda S x.
///          H and S can be any matrices with complex elements: they need not be symmetric,
///          positive definite, or even positive semi-definite.
///          When S is singular, the first m vectors returned are eigenvectors that lie outside
///          the kernal of S.  In this case, the remaining vectors that are returned form an
///          orthonormal basis spanning the kernal of S and their corresponding eigenvalues are
///          set to zero.  The function returns m (the rank of S) as its return value.
///
/// \param[in]       n                     The dimensions of the H and S matrices.
/// \param[in]       hmat                  A length n*n array containing H in column-major order.
/// \param[in]       smat                  A length n*n array containing S in column-major order.
/// \param[out]      evals                 On input, a length n array.
///                                        On exit, the eigenvalues.
/// \param[out]      evecs                 On input, a length n*n array.
///                                        On exit, the normalized eigenvectors and/or kernal basis
///                                        vectors, one per column, in column-major order.
/// \param[in]       thresh                Threshold below which singular values of S are treated
///                                        as zeros.
/// \return  The rank of the matrix S.
///
///////////////////////////////////////////////////////////////////////////////////////////////////
int formic::solve_general_nonsymmetric_eigensystem(const int n,
                                                   const std::complex<double> * const hmat,
                                                   const std::complex<double> * const smat,
                                                   std::complex<double> * const evals,
                                                   std::complex<double> * const evecs,
                                                   const double thresh) {

  // ensure threshold is not negative
  if ( thresh < 0.0 )
    throw formic::Exception("threshold argument must not be negative in formic::solve_general_nonsymmetric_eigensystem");

  // get convenient typedef for complex
  typedef std::complex<double> Cplx;

  // get zero and one
  const Cplx c0 = std::complex<double>(0.0, 0.0);
  const Cplx c1 = std::complex<double>(1.0, 0.0);

  // initialize output to zero
  std::fill(evals, evals + n, c0);
  std::fill(evecs, evecs + n*n, c0);

  // get SVD of S
  std::vector<double> bigs(n, 0.0);
  std::vector<Cplx> bigu(n*n, c0);
  std::vector<Cplx> bigvt(n*n, c0);
  {
    int info = 0;
    int lwork = 20*n;
    std::vector<Cplx> scopy(smat, smat + n*n);
    std::vector<Cplx> work(lwork, c0);
    std::vector<double> rwork(5*n, 0.0);
    formic::zgesvd('A', 'A', n, n, &scopy.at(0), n, &bigs.at(0), &bigu.at(0), n, &bigvt.at(0), n, &work.at(0), lwork, &rwork.at(0), info);
    if ( info != 0 )
      throw formic::Exception("zgesvd failed with info = %i in formic::solve_general_nonsymmetric_eigensystem") % info;
  }

  //// print singular values
  //formic::of << std::endl;
  //for (int i = 0; i < n; i++)
  //  formic::of << boost::format("singular value %4i  = %16.4e\n") % i % bigs.at(i);
  //formic::of << std::endl;

  //// print singular vectors and values
  //formic::of << std::endl;
  //formic::of << boost::format("printing singular values and vectors") << std::endl;
  //for (int j = 0; j < n; j++)
  //  formic::of << boost::format("  %12.2e") % bigs.at(j);
  //formic::of << std::endl;
  //for (int j = 0; j < n; j++)
  //  formic::of << boost::format("--------------");
  //formic::of << std::endl;
  //for (int i = 0; i < n; i++) {
  //  for (int j = 0; j < n; j++)
  //    formic::of << boost::format("  %12.4f") % bigvt.at(j+i*n).real();
  //  formic::of << std::endl;
  //}
  //formic::of << std::endl;

  // count zero singular values
  int nz = 0;
  for (int i = 0; i < n; i++)
    nz += ( bigs.at(i) < thresh ? 1 : 0 );

  // set and print reduced dimension
  const int m = n - nz;

  // get U' and V'
  std::vector<Cplx> uprm(n*m, c0);
  std::vector<Cplx> vprm(n*m, c0);
  for (int i = 0; i < n; i++)
  for (int j = 0; j < m; j++) {
    uprm.at(i+j*n) =            bigu.at(i+j*n);
    vprm.at(i+j*n) = std::conj(bigvt.at(j+i*n));
  }

  // get inverted nonzero singular values
  std::vector<Cplx> inverted_sv(m, c0);
  for (int j = 0; j < m; j++)
    inverted_sv.at(j) = c1 / bigs.at(j);

  //formic::of << boost::format("printing vprm") << std::endl;
  //for (int i = 0; i < n; i++) {
  //  for (int j = 0; j < m; j++)
  //    formic::of << boost::format("  %12.4f") % vprm.at(i+j*m);
  //  formic::of << std::endl;
  //}
  //formic::of << std::endl;

  // Project H and multiply it by inverse singular values.
  // i.e. set gmat = diagonal_matrix_of(inverted_sv) * conjugate_transpose(uprm) * hmat * vprm
  std::vector<Cplx> ct_uprm_times_hmat_times_vprm(m*m, c0);
  {
    std::vector<Cplx> hmat_times_vprm(n*m, c0);
    formic::xgemm('N', 'N', n, m, n, c1,        hmat, n,            &vprm.at(0), n, c0,               &hmat_times_vprm.at(0), n);
    formic::xgemm('C', 'N', m, m, n, c1, &uprm.at(0), n, &hmat_times_vprm.at(0), n, c0, &ct_uprm_times_hmat_times_vprm.at(0), m);
  }
  for (int j = 0; j < m; j++)
  for (int i = 0; i < m; i++)
    ct_uprm_times_hmat_times_vprm.at(i+j*m) *= inverted_sv.at(i);
  Cplx * const gmat = &ct_uprm_times_hmat_times_vprm.at(0);

  //std::vector<Cplx> gmat_copy(gmat, gmat+m*m);

  // compute eigenvalues and eigenvectors of this matrix (evecs placed in columns of evecs_m)
  std::vector<Cplx> evecs_m(m*m, c0);
  const bool use_eigen_diagonalizer = false;
  {
    int info = 0;
    int lwork = 20*m;
    std::vector<Cplx> work(lwork, c0);
    std::vector<double> rwork(2*m, 0.0);
    formic::zgeev('N', 'V', m, gmat, m, evals, &evecs_m.at(0), m, &evecs_m.at(0), m, &work.at(0), lwork, &rwork.at(0), info);
    if ( info != 0 )
      throw formic::Exception("zgeev failed with info = %i in formic::solve_general_nonsymmetric_eigensystem") % info;
  }

  //std::vector<Cplx> resid(m, c0);
  //for (int i = 0; i < m; i++) {
  //  formic::xgemm('N', 'N', m, 1, m, c1, &gmat_copy.at(0), m, &evecs_m.at(i*m), m, c0, &resid.at(0), m);
  //  formic::xaxpy(m, -evals[i], &evecs_m.at(i*m), 1, &resid.at(0), 1);
  //  formic::of << boost::format("error %i = %16.4e") % i % std::sqrt(std::abs(formic::xdotc(m, &resid.at(0), 1, &resid.at(0), 1))) << std::endl;
  //}

  //formic::of << boost::format("printing vprm") << std::endl;
  //for (int i = 0; i < n; i++) {
  //  for (int j = 0; j < m; j++)
  //    formic::of << boost::format("  %12.4f") % vprm.at(i+j*m);
  //  formic::of << std::endl;
  //}
  //formic::of << std::endl;

  // compute eigenvectors of original eigensystem corresponding to non-zero singular values
  formic::xgemm('N', 'N', n, m, m, c1, &vprm.at(0), n, &evecs_m.at(0), m, c0, evecs, n);

  //formic::of << std::endl;
  //std::vector<Cplx> resid0(n, c0);
  //std::vector<Cplx> resid1(n, c0);
  //std::vector<Cplx> resid2(n, c0);
  //for (int i = 0; i < m; i++) {
  //  formic::xgemm('N', 'N', n, 1, m, c1, &vprm.at(0), n, &evecs_m.at(i*m), m, c0, &resid0.at(0), n);
  //  formic::xgemm('N', 'N', n, 1, n, c1,        hmat, n,    &resid0.at(0), n, c0, &resid1.at(0), n);
  //  formic::xgemm('C', 'N', m, 1, n, c1, &uprm.at(0), n,    &resid1.at(0), n, c0, &resid2.at(0), m);
  //  for (int j = 0; j < m; j++)
  //    resid2.at(j) *= inverted_sv.at(j);
  //  formic::xaxpy(m, -evals[i], &evecs_m.at(i*m), 1, &resid2.at(0), 1);
  //  formic::of << boost::format("after xgemm, error %4i  = %16.4e") % i % std::sqrt(std::abs(formic::xdotc(m, &resid2.at(0), 1, &resid2.at(0), 1))) << std::endl;
  //  //formic::xgemm('N', 'N', n, 1, n, c1, hmat, n, evecs+i*n, n, c0, &resid.at(0), n);
  //  //formic::xgemm('N', 'N', n, 1, n, -evals[i], smat, n, evecs+i*n, n, c1, &resid.at(0), n);
  //  //formic::xgemm('C', 'N', m, 1, n, c1, &uprm.at(0), n, &resid.at(0), n, c0, &resid2.at(0), m);
  //  //formic::of << boost::format("error %i = %16.4e") % i % std::sqrt(std::abs(formic::xdotc(m, &resid2.at(0), 1, &resid2.at(0), 1))) << std::endl;
  //}
  //formic::of << std::endl;

  // get eigenvectors of original eigensystem corresponding to zero singular values
  for (int j = m; j < n; j++)
  for (int i = 0; i < n; i++)
    evecs[i+j*n] = std::conj(bigvt.at(j+i*n));

  // return the number of eigenvectors corresponding to non-zero singular values
  return m;

}
