///////////////////////////////////////////////////////////////////////////////////////////////////
/// \file formic/utils/numeric.h
///
/// \brief   header file for miscellaneous functions related to numbers
///
///////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef FORMIC_UTILS_NUMERIC_HEADER
#define FORMIC_UTILS_NUMERIC_HEADER

#include<string>
#include<complex>
#include<cassert>
#include<cmath>
#include<vector>
#include<set>

#include<formic/utils/exception.h>
#include<formic/utils/zero_one.h>
#include<formic/utils/mantissa_exponent.h>

namespace formic {

  void unbiased_ratio_of_means(const int n, const double * const p, const double * const f, const double * const g, double & r, double & v);
  void mpi_unbiased_ratio_of_means(const int n, const double * const p, const double * const f, const double * const g, double & r, double & v);

  void get_pair_ioff(const int n, std::vector<int> & ioff);

  int binom_coeff(int n, int m);

  int n_integer_solutions(const int n, const int r, const int k, int * work = 0);

  std::string format_number(const std::string & f, const double value);

  std::string format_number(const std::string & f, const std::complex<double> value);

  ///////////////////////////////////////////////////////////////////////////////////////////////////
  /// \brief   computes integer powers
  ///
  /// \param[in]     a        number to be raised to a power
  /// \param[in]     b        the integer power to raise the number to, assumed >= 0
  ///
  /// \return a^b
  ///
  ///////////////////////////////////////////////////////////////////////////////////////////////////
  template <class T> inline T pow(const T a, int b) {
    assert(b >= 0);
    T retval = 1;
    for ( ; b > 0; b--) retval *= a;
    return retval;
  };

  ///////////////////////////////////////////////////////////////////////////////////////////////////
  /// \brief   returns the complex conjugate of a real number
  ///
  /// \param[in]     x    the number to be conjugated
  ///
  /// \return the complex conjugate of x
  ///
  ///////////////////////////////////////////////////////////////////////////////////////////////////
  template <class S> inline S conj(S x) { return x; }

  ///////////////////////////////////////////////////////////////////////////////////////////////////
  /// \brief   returns the complex conjugate of a complex number
  ///
  /// \param[in]     x    the number to be conjugated
  ///
  /// \return the complex conjugate of x
  ///
  ///////////////////////////////////////////////////////////////////////////////////////////////////
  template <class S> inline std::complex<S> conj(std::complex<S> x) { return std::conj(x); }

  ///////////////////////////////////////////////////////////////////////////////////////////////////
  /// \brief   returns the real part of a real number
  ///
  /// \param[in]     x    the number whose real part we want
  ///
  /// \return the real part of x
  ///
  ///////////////////////////////////////////////////////////////////////////////////////////////////
  template <class S> inline S real(S x) { return x; }

  ///////////////////////////////////////////////////////////////////////////////////////////////////
  /// \brief   returns the real part of a complex number
  ///
  /// \param[in]     x    the number whose real part we want
  ///
  /// \return the real part of x
  ///
  ///////////////////////////////////////////////////////////////////////////////////////////////////
  template <class S> inline S real(std::complex<S> x) { return x.real(); }

  ///////////////////////////////////////////////////////////////////////////////////////////////////
  /// \brief   returns the imaginary part of a real number
  ///
  /// \param[in]     x    the number whose imaginary part we want
  ///
  /// \return the imaginary part of x
  ///
  ///////////////////////////////////////////////////////////////////////////////////////////////////
  template <class S> inline S imag(S x) { return 0; }

  ///////////////////////////////////////////////////////////////////////////////////////////////////
  /// \brief   returns the imaginary part of a complex number
  ///
  /// \param[in]     x    the number whose imaginary part we want
  ///
  /// \return the imaginary part of x
  ///
  ///////////////////////////////////////////////////////////////////////////////////////////////////
  template <class S> inline S imag(std::complex<S> x) { return x.imag(); }

  ///////////////////////////////////////////////////////////////////////////////////////////////////
  /// \brief   returns the square of the norm of a real number
  ///
  /// \param[in]     x        the number whose square norm we want
  ///
  /// \return the square norm of x
  ///
  ///////////////////////////////////////////////////////////////////////////////////////////////////
  template <class S> inline S square_norm(S x) { return x*x; }

  ///////////////////////////////////////////////////////////////////////////////////////////////////
  /// \brief   returns the square of the norm of a complex number
  ///
  /// \param[in]     x        the number whose square norm we want
  ///
  /// \return the square norm of x
  ///
  ///////////////////////////////////////////////////////////////////////////////////////////////////
  template <class S> inline S square_norm(std::complex<S> x) { return x.real() * x.real() + x.imag() * x.imag(); }

  ///////////////////////////////////////////////////////////////////////////////////////////////////
  /// \brief   rounds a double to the nearest integer
  ///
  /// \param[in]     d        the number to round
  ///
  /// \return the integer nearest to d
  ///
  ///////////////////////////////////////////////////////////////////////////////////////////////////
  inline int round(const double d) { return ( d >= 0 ? int(d+0.5) : int(d-0.5) ); }

  ///////////////////////////////////////////////////////////////////////////////////////////////////
  /// \brief   rounds a double to a given digit precision
  ///
  /// \param[in]     x        the number to round
  /// \param[in]     digits   the number of digits to round to
  ///
  /// \return the rounded number
  ///
  ///////////////////////////////////////////////////////////////////////////////////////////////////
  inline double round(const double x, const int digits) {
    const double factor = std::pow(double(10.0), digits);
    const double y = factor * x;
    return double( y >= 0.0 ? (long int)(y+0.5) : (long int)(y-0.5) ) / factor;
  }

  ///////////////////////////////////////////////////////////////////////////////////////////////////
  /// \brief   rounds a complex double to a given digit precision
  ///
  /// \param[in]     x        the number to round
  /// \param[in]     digits   the number of digits to round to
  ///
  /// \return the rounded number
  ///
  ///////////////////////////////////////////////////////////////////////////////////////////////////
  inline std::complex<double> round(const std::complex<double> x, const int digits) {
    return std::complex<double>( formic::round(x.real(), digits), formic::round(x.imag(), digits) );
  }

  ///////////////////////////////////////////////////////////////////////////////////////////////////
  /// \brief   checks whether a double is equal to some integer divided by two
  ///
  /// \param[in]     x        the number to check
  ///
  /// \return whether the number is half of an integer
  ///
  ///////////////////////////////////////////////////////////////////////////////////////////////////
  inline bool is_half_integer(const double x) {
    return ( std::fabs( formic::round(2.0 * x) - 2.0 * x ) < 1.0e-9 );
  }

  ///////////////////////////////////////////////////////////////////////////////////////////////////
  /// \brief   computes a compound index for a pair of distinct numbers i < j in the
  ///          range 0, 1, 2, ..., n-1
  ///
  /// \param[in]     i             the 1st element of the pair
  /// \param[in]     j             the 2nd element of the pair
  /// \param[in]     n             the number of elements in the range
  ///
  /// \return the compound index
  ///
  ///////////////////////////////////////////////////////////////////////////////////////////////////
  inline int cmpd_pair_index(const int i, const int j, const int n) {
    assert( i < j );
    assert( i >= 0 && i < n );
    assert( j >= 0 && j < n );
    return ( 2*n - i - 3 ) * i / 2 + j - 1;
  }

  ///////////////////////////////////////////////////////////////////////////////////////////////////
  /// \brief   Computes a compound index for a pair of distinct numbers i,j in the
  ///          range 0, 1, 2, ..., n-1.  Works for both i < j and i > j.
  ///
  /// \param[in]     i             the 1st element of the pair
  /// \param[in]     j             the 2nd element of the pair
  /// \param[in]     n             the number of elements in the range
  ///
  /// \return the compound index
  ///
  ///////////////////////////////////////////////////////////////////////////////////////////////////
  inline int cmpd_pair_index_unordered(const int i, const int j, const int n) {
    assert( i != j );
    assert( i >= 0 && i < n );
    assert( j >= 0 && j < n );
    return formic::cmpd_pair_index(std::min(i,j), std::max(i,j), n);
  }

  ///////////////////////////////////////////////////////////////////////////////////////////////////
  /// \brief   computes a compound index for a pair of numbers i <= j in the
  ///          range 0, 1, 2, ..., n-1
  ///
  /// \param[in]     i             the 1st element of the pair
  /// \param[in]     j             the 2nd element of the pair
  /// \param[in]     n             the number of elements in the range
  ///
  /// \return the compound index
  ///
  ///////////////////////////////////////////////////////////////////////////////////////////////////
  inline int cmpd_pair_index_incl(const int i, const int j, const int n) {
    assert( i <= j );
    assert( i >= 0 && i < n );
    assert( j >= 0 && j < n );
    return ( 2*n - i - 1 ) * i / 2 + j;
  }

  ///////////////////////////////////////////////////////////////////////////////////////////////////
  /// \brief   computes a compound index for any pair of numbers i,j in the
  ///          range 0, 1, 2, ..., n-1.
  ///
  ///          Note that both (i,j) and (j,i) are allowed as inputs but
  ///          that they will both return the same compound index.
  ///
  /// \param[in]     i             the 1st element of the pair
  /// \param[in]     j             the 2nd element of the pair
  /// \param[in]     n             the number of elements in the range
  ///
  /// \return the compound index
  ///
  ///////////////////////////////////////////////////////////////////////////////////////////////////
  inline int cmpd_pair_index_incl_unordered(const int i, const int j, const int n) {
    return cmpd_pair_index_incl(std::min(i,j), std::max(i,j), n);
  }

  double porabola_min_max(const double * const xvec, const double * const yvec);

  ///////////////////////////////////////////////////////////////////////////////////////////////////
  /// \brief   checks whether all elements in the range [first,last) are distinct
  ///          by performing a double loop over all elements, checking each pairwise equality
  ///
  /// \param[in]     first         iterator to the beginning of the range
  /// \param[in]     last          iterator to the end of the range
  ///
  /// \return true if all elements are different, false otherwise
  ///
  ///////////////////////////////////////////////////////////////////////////////////////////////////
  template <class Iterator>
  inline bool all_different_by_loops(Iterator first, Iterator last) {
    for ( ; first != last; first++) {
      Iterator i = first;
      i++;
      for ( ; i != last; i++)
        if (*first == *i)
          return false;
    }
    return true;
  }

  ///////////////////////////////////////////////////////////////////////////////////////////////////
  /// \brief   given the number of indices and a container of indices, returns a vector of
  ///          the indices that are not in the container.
  ///
  /// \param[in]     n       the number of indices ( indices count from 0 to n-1 )
  /// \param[in]     c       a container holding indices
  ///
  /// \return  a vector of the indices that were not in the container
  ///
  ///////////////////////////////////////////////////////////////////////////////////////////////////
  template <class Container>
  inline std::vector<int> other_indices(const int n, const Container & c) {

    // check sanity
    if ( n < 0 )
      throw formic::Exception("n must be non-negative in formic::other_indices");

    // get a set of the indices in the container
    std::set<int> s(c.begin(), c.end());

    // place the indices that are not in the container in a vector
    std::vector<int> others;
    for (int i = 0; i < n; i++)
      if ( s.count(i) == 0 )
        others.push_back(i);

    // return the vector of indices not in the container
    return others;

  }

  ///////////////////////////////////////////////////////////////////////////////////////////////////
  /// \brief   Binary operator for use with std::sort that will sort elements in ascending order
  ///          by their absolute values.
  ///
  ///          Usage:  std::sort(x.begin(), x.end(), formic::AscendingNormSorter<S>());
  ///
  ///////////////////////////////////////////////////////////////////////////////////////////////////
  template<class S> struct AscendingNormSorter {
    bool operator()(const S a, const S b) { return ( std::abs(a) > std::abs(b) ? true : false ); }
  };

  ///////////////////////////////////////////////////////////////////////////////////////////////////
  /// \brief   Binary operator for use with std::accumulate that picks out the element with the
  ///          largest absolute value.
  ///
  ///          Usage:  std::accumulate(x.begin(), x.end(), 0.0, formic::LargestAbsAccumulator<S>());
  ///
  ///////////////////////////////////////////////////////////////////////////////////////////////////
  template<class S> struct LargestAbsAccumulator {
    S operator()(const S a, const S b) { return ( std::abs(a) > std::abs(b) ? a : b ); }
  };

  ///////////////////////////////////////////////////////////////////////////////////////////////////
  /// \brief   returns the number of real parameters needed to hold the array's elements
  ///
  /// \param[in]     n    the array length
  /// \param[in]     x    the array (used only to know the variable type)
  ///
  ///////////////////////////////////////////////////////////////////////////////////////////////////
  inline int n_real_params(const int n, const double * const x) {
    return n;
  }

  ///////////////////////////////////////////////////////////////////////////////////////////////////
  /// \brief   returns the number of real parameters needed to hold the array's elements
  ///
  /// \param[in]     n    the array length
  /// \param[in]     x    the array (used only to know the variable type)
  ///
  ///////////////////////////////////////////////////////////////////////////////////////////////////
  inline int n_real_params(const int n, const std::complex<double> * const x) {
    return 2*n;
  }

  ///////////////////////////////////////////////////////////////////////////////////////////////////
  /// \brief   Converts a real array of length 2n to a complex array of length n.
  ///
  ///////////////////////////////////////////////////////////////////////////////////////////////////
  inline void compress_to_complex(const int n, const double * const r, std::complex<double> * const c) {
    for (int i = 0; i < n; i++)
      c[i] = std::complex<double>(r[i], r[i+n]);
  }

  ///////////////////////////////////////////////////////////////////////////////////////////////////
  /// \brief   Copies one real array to another
  ///
  ///////////////////////////////////////////////////////////////////////////////////////////////////
  inline void compress_to_complex(const int n, const double * const r, double * const c) {
    for (int i = 0; i < n; i++)
      c[i] = r[i];
  }

  ///////////////////////////////////////////////////////////////////////////////////////////////////
  /// \brief   Copies one real array to another
  ///
  ///////////////////////////////////////////////////////////////////////////////////////////////////
  inline void expand_to_real(const int n, const double * const c, double * const r) {
    for (int i = 0; i < n; i++)
      r[i] = c[i];
  }

  ///////////////////////////////////////////////////////////////////////////////////////////////////
  /// \brief   Converts a complex array of length n to a real array of length 2n.
  ///
  ///////////////////////////////////////////////////////////////////////////////////////////////////
  inline void expand_to_real(const int n, const std::complex<double> * const c, double * const r) {
    for (int i = 0; i < n; i++)
      r[i] = c[i].real();
    for (int i = 0; i < n; i++)
      r[i+n] = c[i].imag();
  }

  ///////////////////////////////////////////////////////////////////////////////////////////////////
  /// \brief   finds the integer closest to the square root of a given integer
  ///
  /// \param[in]     a        a non-negative integer
  ///
  /// \return  the integer closest to the square root of a
  ///
  ///////////////////////////////////////////////////////////////////////////////////////////////////
  template <class T> inline T int_sqrt(const T a) {

    // ensure argument is non negative
    if ( a < 0 )
      throw formic::Exception("formic::int_sqrt requires a non negative argument.");

    // treat zero as a special case
    if ( a == 0 ) return 0;

    // get a non-negative integer near the square root
    T low = T(std::abs(std::sqrt(double(a))));

    // if guess is exact then we're done
    if ( low*low == a ) return low;

    // make sure this integer is not larger than the square root
    for ( ; low*low > a; low--) {}

    // if value is now exact then we're done
    if ( low*low == a ) return low;

    // The integer is now strictly less than the square root.
    // Now make it the one directly below the square root.
    for ( ; (low+1)*(low+1) < a; low++) {}

    // now get the integer equal to or directly above the square root
    const T high = low + 1;

    // return the integer closest to the square root
    const double ldist = std::abs( double(low ) - std::abs(std::sqrt(double(a))) );
    const double hdist = std::abs( double(high) - std::abs(std::sqrt(double(a))) );
    return ( ldist < hdist ? low : high );

  };

  ///////////////////////////////////////////////////////////////////////////////////////////////////
  /// \brief   determines if the supplied integer is a perfect square and if so returns its sqrt,
  ///          returning -1 in the case that the integer is not a perfect square
  ///
  /// \param[in]     a        a non-negative integer
  ///
  /// \return  sqrt(a) if a is a perfect square, -1 otherwise
  ///
  ///////////////////////////////////////////////////////////////////////////////////////////////////
  inline int perfect_sqrt(const int a) {

    // ensure argument is non negative
    if ( a < 0 )
      throw formic::Exception("formic::perfect_sqrt requires a non negative argument.");

    // if a is a perfect square, then return the sqrt
    const int s = formic::int_sqrt(a);
    if ( s * s == a )
      return s;

    // otherwise return -1 to signal that a is not a perfect square
    return -1;

  }

  template <class S> void matrix_inverse_lu(const int n, S & det, S * const A, S * const work, int * const iwork);

  template <class S> formic::MantExp<S> matrix_determinant_lu(const int n, S * const A, int * const iwork);

  template <class S> void matrix_neg_half(const int n, S * const A);

  int solve_general_nonsymmetric_eigensystem(const int n,
                                             const std::complex<double> * const hmat,
                                             const std::complex<double> * const smat,
                                             std::complex<double> * const evals,
                                             std::complex<double> * const evecs,
                                             const double thresh = 1.0e-9);

}

#endif
