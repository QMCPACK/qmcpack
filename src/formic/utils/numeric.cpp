///////////////////////////////////////////////////////////////////////////////////////////////////
/// \file formic/utils/numeric.cpp
///
/// \brief   implementation file for miscellaneous functions related to numbers
///
///////////////////////////////////////////////////////////////////////////////////////////////////

#include<formic/utils/numeric.h>
#include<formic/utils/mpi_interface.h>

#include<boost/scoped_array.hpp>
#include<boost/format.hpp>

///////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief   get an offset array used in compounding pairs of distinct indices
///
///   For two indices i,j with i < j, we have:    compound(i,j) = i + ioff[j];
///
/// \param[in]       n        desired length of the array
/// \param[in,out]   ioff     on exit, the offset array
///
///////////////////////////////////////////////////////////////////////////////////////////////////
void formic::get_pair_ioff(int n, std::vector<int> & ioff) {

  if (ioff.size() != n)
    ioff.resize(n);

  if ( n <= 0 )
    return;

  ioff.at(0) = 0;
  for (int i = 1; i < n; i++)
    ioff.at(i) = ioff.at(i-1) + i - 1;

}

///////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief   computes the binomial coefficient
///
/// \param[in]     n        number of things
/// \param[in]     m        how many things to take at a time
///
/// \return the number of ways n things can be taken m at a time
///
///////////////////////////////////////////////////////////////////////////////////////////////////
int formic::binom_coeff(int n, int m) {
  if (n < 0 || m < 0 || m > n) return 0;
  double retval = 1.0;
  while (m > 0) retval = ( retval * (n--) ) / (m--);
  return int(retval+0.5);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief   Returns the number of solutions to the equation
///          x(1) + x(2) + ... + x(n) = r
///          when the variables x(i) are constrained to be
///          integers in the range (0, 1, 2, ..., k)
///
/// \param[in]     n        number of variables
/// \param[in]     r        sum of varaibles
/// \param[in]     k        range of each variable
/// \param[out]    work     integer workspace, either null or size >= k+1
///
///////////////////////////////////////////////////////////////////////////////////////////////////
int formic::n_integer_solutions(const int n, const int r, const int k, int * work) {

  assert( n >= 0 );
  assert( r >= 0 );
  assert( k >= 0 );

  // if requested, dynamically allocate the work space
  boost::scoped_array<int> dynamic_work;
  if (work == 0) {
    dynamic_work.reset( new int[k+1] );
    work = dynamic_work.get();
  }

  // initialize an array to hold the number of variables having each allowed value
  int * const n_with_value = work;
  work += (k+1);
  for (int i = 0; i <= k; i++)
    n_with_value[i] = 0;

  // initialize the return value
  int retval = 0;

  // Loop over all possible distributions of variables among the values.
  // Note that we do not directly track of how many variables are equal to zero,
  // as this is known by how many variables take on other values.
  while (true) {

    // compute the number of nonzero variables
    int n_nonzero = 0;
    for (int i = 1; i <= k; i++)
      n_nonzero += n_with_value[i];

    // compute the sum of the variables
    int sum = 0;
    for (int i = 1; i <= k; i++)
      sum += i * n_with_value[i];

    // if this distribution solves the equation, count how many ways it can occur
    if (sum == r && n_nonzero <= n) {

      // determine how many variables are nonzero
      int t = 0;
      for (int i = 1; i <= k; i++)
        t += n_with_value[i];

      // count how many ways the variables can satisfy this distribution
      int occurances = formic::binom_coeff(n, t);
      for (int i = 1; i < k; i++) {
        occurances *= formic::binom_coeff(t, n_with_value[i]);
        t -= n_with_value[i]; // t is now equal to the number of variables greater than i
      }

      // record how many ways the variables satisfy this distribution
      retval += occurances;

    }

    // increment to the next distribution of variables
    int p;
    for (p = k; p > 0; p--)
      if (++n_with_value[p] > n)
        n_with_value[p] = 0;
      else
        break;

    // stop iterating if all distributions have been processed
    if (p == 0) break;

  }

  // return the result
  return retval;

}

///////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief   formats a real number into a string
///
/// \param[in]     f        the formatting string used by boost::format
/// \param[in]     value    the number to be formatted
///
/// \return the string containing the formatted number
///
///////////////////////////////////////////////////////////////////////////////////////////////////
std::string formic::format_number(const std::string & f, const double value) {
  return (boost::format(f) % value).str();
}

///////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief   formats a complex number into a string
///
/// \param[in]     f        the formatting string used by boost::format
/// \param[in]     value    the number to be formatted
///
/// \return the string containing the formatted number
///
///////////////////////////////////////////////////////////////////////////////////////////////////
std::string formic::format_number(const std::string & f, const std::complex<double> value) {
  std::string retval;
  retval.append("( ");
  retval.append( (boost::format(f) % value.real()).str() );
  retval.append(", ");
  retval.append( (boost::format(f) % value.imag()).str() );
  retval.append(" )");
  return retval;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief   computes the unbiased estimate of a ratio of means:  <f>_p / <g>_p  in which the
///          numerator and denominator values are sampled from the same probability distribution p
///
/// \param[in]     n        the number of samples
/// \param[in]     p        the probability weight for each sample
/// \param[in]     f        the numerator samples
/// \param[in]     g        the denominator samples
/// \param[out]    r        on exit, the estimate of the ratio <f>_p / <g>_p
/// \param[out]    v        on exit, the estimate of the variance in the ratio
///
///////////////////////////////////////////////////////////////////////////////////////////////////
void formic::unbiased_ratio_of_means(const int n, const double * const p, const double * const f, const double * const g, double & r, double & v) {

  // compute the normalization, the numerator and denominator means, the means of the squares, and the mean of the products
  double nm = 0.0; // normalization constant
  double mf = 0.0; // mean of numerator
  double mg = 0.0; // mean of denominator
  double sf = 0.0; // mean of the square of the numerator terms
  double sg = 0.0; // mean of the square of the denominator terms
  double mp = 0.0; // mean of the product of numerator times denominator
  for (int i = 0; i < n; i++) {
    nm += p[i];
    double x = p[i] * f[i];
    mf += x;
    sf += x * f[i];
    mp += x * g[i];
    x = p[i] * g[i];
    mg += x;
    sg += x * g[i];
  }
  mf /= nm;
  mg /= nm;
  sf /= nm;
  sg /= nm;
  mp /= nm;

  // compute the numerator and denominator variances and the covariance
  const double vf = ( sf - mf * mf ) * double(n) / double(n-1);
  const double vg = ( sg - mg * mg ) * double(n) / double(n-1);
  const double cv = ( mp - mf * mg ) * double(n) / double(n-1);

  // compute the unbiased estimate of the ratio of means
  r = ( mf / mg ) / ( 1.0 + ( vg / mg / mg - cv / mf / mg ) / double(n) );

  // compute the unbiased estimate of the variance of the ratio of means
  v = ( mf * mf / mg / mg / double(n) ) * ( vf / mf / mf + vg / mg / mg - 2.0 * cv / mf / mg );

}

///////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief   computes the unbiased estimate of a ratio of means:  <f>_p / <g>_p  in which the
///          numerator and denominator values are sampled from the same probability distribution p
///          and samples are combined across all processors
///
/// \param[in]     n        the number of samples on this process
/// \param[in]     p        the probability weight for each sample
/// \param[in]     f        the numerator samples
/// \param[in]     g        the denominator samples
/// \param[out]    r        on exit, the estimate of the ratio <f>_p / <g>_p
/// \param[out]    v        on exit, the estimate of the variance in the ratio
///
///////////////////////////////////////////////////////////////////////////////////////////////////
void formic::mpi_unbiased_ratio_of_means(const int n, const double * const p, const double * const f, const double * const g, double & r, double & v) {

  // compute the normalization, the numerator and denominator means, the means of the squares, and the mean of the products
  double y[7];
  y[0] = 0.0; // normalization constant
  y[1] = 0.0; // mean of numerator
  y[2] = 0.0; // mean of denominator
  y[3] = 0.0; // mean of the square of the numerator terms
  y[4] = 0.0; // mean of the square of the denominator terms
  y[5] = 0.0; // mean of the product of numerator times denominator
  y[6] = double(n); // number of samples
  for (int i = 0; i < n; i++) {
    y[0] += p[i];
    double x = p[i] * f[i];
    y[1] += x;
    y[3] += x * f[i];
    y[5] += x * g[i];
    x = p[i] * g[i];
    y[2] += x;
    y[4] += x * g[i];
  }
  double z[7];
  formic::mpi::allreduce(&y[0], &z[0], 7, MPI_SUM);
  const double mf = z[1] / z[0]; // mean of numerator
  const double mg = z[2] / z[0]; // mean of denominator
  const double sf = z[3] / z[0]; // mean of the square of the numerator terms
  const double sg = z[4] / z[0]; // mean of the square of the denominator terms
  const double mp = z[5] / z[0]; // mean of the product of numerator times denominator
  const double ns = z[6];        // number of samples

  // compute the numerator and denominator variances and the covariance
  const double vf = ( sf - mf * mf ) * ns / ( ns - 1.0 );
  const double vg = ( sg - mg * mg ) * ns / ( ns - 1.0 );
  const double cv = ( mp - mf * mg ) * ns / ( ns - 1.0 );

  // compute the unbiased estimate of the ratio of means
  r = ( mf / mg ) / ( 1.0 + ( vg / mg / mg - cv / mf / mg ) / ns );

  // compute the unbiased estimate of the variance of the ratio of means
  v = ( mf * mf / mg / mg ) * ( vf / mf / mf + vg / mg / mg - 2.0 * cv / mf / mg );

}
