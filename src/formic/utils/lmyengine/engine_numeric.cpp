///////////////////////////////////////////////////////////////////////////////////////////////////
/// \file cqmc/numeric/numeric.cpp
///
/// \brief   implementation file for miscellaneous functions related to numbers
///
///////////////////////////////////////////////////////////////////////////////////////////////////

//#include <mpi.h>
#include "formic/utils/mpi_interface.h"

#include "formic/utils/zero_one.h"
#include "formic/utils/numeric.h"
#include "engine_numeric.h"


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
template<typename S, typename T> void cqmc::unbiased_ratio_of_means(const int n, const double * const p, const S * const f, const T * const g, const bool correct, double & r, double & v) {

  // compute the normalization, the numerator and denominator means, the means of the squares, and the mean of the products
  double nm = 0.0; // normalization constant
  double mf = 0.0; // mean of numerator
  double mg = 0.0; // mean of denominator
  double sf = 0.0; // mean of the square of the numerator terms
  double sfnc = 0.0; // mean of the square of the numerator used in exact sampling 
  double sg = 0.0; // mean of the square of the denominator terms
  double mp = 0.0; // mean of the product of numerator times denominator
  for (int i = 0; i < n; i++) {
    nm += formic::real(p[i]);
    S x = p[i] * f[i];
    mf += formic::real(x);
    if ( correct ) {
      sf += formic::real(x * f[i]);
      mp += formic::real(x * g[i]);
    }
    else 
      sfnc += formic::real(p[i] * (f[i] / g[i]) * (f[i] / g[i]) * g[i]);
    x = p[i] * g[i];
    mg += formic::real(x);
    if ( correct )
      sg += formic::real(x * g[i]);
  }
  mf /= nm;
  mg /= nm;
  sf /= nm;
  sfnc /= nm;
  sg /= nm;
  mp /= nm;

  // compute the numerator and denominator variances and the covariance
  const double vf = ( sf - mf * mf ) * static_cast<double>(n) / static_cast<double>(n-1);
  const double vg = ( sg - mg * mg ) * static_cast<double>(n) / static_cast<double>(n-1);
  const double cv = ( mp - mf * mg ) * static_cast<double>(n) / static_cast<double>(n-1);

  double z = ( correct ? 1.0 : 0.0);

  // compute the unbiased estimate of the ratio of means
  r = ( mf / mg ) / ( 1.0 + z * ( vg / mg / mg - cv / mf / mg ) / static_cast<double>(n) );

  // compute the unbiased estimate of the variance of the ratio of means
  if ( correct ) 
    v = ( mf * mf / mg / mg / static_cast<double>(n) ) * ( vf / mf / mf + vg / mg / mg - 2.0 * cv / mf / mg );
  else 
    v = sfnc - mf * mf; 

}

template void cqmc::unbiased_ratio_of_means(const int n, const double * const p, const double * const f, const double * const g, const bool correct, double & r, double & v);
template void cqmc::unbiased_ratio_of_means(const int n, 
                                            const double * const p, 
                                            const std::complex<double> * const f, 
                                            const std::complex<double> * const g, 
                                            const bool correct, 
                                            double & r, 
                                            double & v);
template void cqmc::unbiased_ratio_of_means(const int n, 
                                            const double * const p, 
                                            const std::complex<double> * const f, 
                                            const double * const g, 
                                            const bool correct, 
                                            double & r, 
                                            double & v);


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
template<typename S, typename T> void cqmc::mpi_unbiased_ratio_of_means(const int n, const double * const p, const S * const f, const T * const g, const bool correct, S & r, double & v) {

  // compute the normalization, the numerator and denominator means, the means of the squares, and the mean of the products
  S y[8];
  y[0] = 0.0; // normalization constant
  y[1] = 0.0; // mean of numerator
  y[2] = 0.0; // mean of denominator
  y[3] = 0.0; // mean of the square of the numerator terms
  y[4] = 0.0; // mean of the square of the denominator terms
  y[5] = 0.0; // mean of the product of numerator times denominator
  y[6] = static_cast<S>(n); // number of samples
  y[7] = 0.0; // mean of square of numerator terms used in exact sampling(when correct is false) 

  // normalization constant 
  double norm = 0.0;

  // accumulate data
  for (int i = 0; i < n; i++) {
    norm += p[i];
    S x = p[i] * f[i];
    y[1] += x;
    if ( correct ) {
      y[3] += x * formic::conj(f[i]);
      y[5] += x * g[i];
    }
    else 
      y[7] += p[i] * f[i] * formic::conj(f[i]);
    x = p[i] * g[i];
    y[2] += x;
    if ( correct )
      y[4] += x * formic::conj(g[i]);
  }
  S z[8];
  double total_norm = 0.0;
  formic::mpi::allreduce(&y[0], &z[0], 8, MPI_SUM);
  formic::mpi::allreduce(&norm, &total_norm, 1, MPI_SUM);

  const S mf = z[1] / total_norm; // mean of numerator
  const S mg = z[2] / total_norm; // mean of denominator
  const S sf = z[3] / total_norm; // mean of the square of the numerator terms
  const S sg = z[4] / total_norm; // mean of the square of the denominator terms
  const S mp = z[5] / total_norm; // mean of the product of numerator times denominator
  const S sfnc = z[7] / total_norm; // mean of square of the numerator used in exact sampling  
  const S ns = z[6];        // number of samples

  // compute the numerator and denominator variances and the covariance
  const S vf = ( sf - mf * mf ) * ns / ( ns - 1.0 );
  const S vg = ( sg - mg * mg ) * ns / ( ns - 1.0 );
  const S cv = ( mp - mf * mg ) * ns / ( ns - 1.0 );

  double c = (correct ? 1.0 : 0.0);

  // compute the unbiased estimate of the ratio of means
  r = ( mf / mg ) / ( 1.0 + c * ( vg / mg / mg - cv / mf / mg ) / ns );
  //std::cout << vg / mg / mg << "   " << cv / mf / mg << std::endl;

  // compute the unbiased estimate of the variance of the ratio of means
  if ( correct ) 
    v = formic::real(( mf * formic::conj(mf) / mg / mg ) * ( vf / mf / formic::conj(mf) + vg / mg / mg - 2.0 * cv / mf / mg ));
  else 
    v = formic::real(sfnc - mf * formic::conj(mf));

}

template void cqmc::mpi_unbiased_ratio_of_means(const int n, const double * const p, const double * const f, const double * const g, const bool correct, double & r, double & v);
template void cqmc::mpi_unbiased_ratio_of_means(const int n, 
                                                const double * const p, 
                                                const std::complex<double> * const f, 
                                                const std::complex<double> * const g, 
                                                const bool correct, 
                                                std::complex<double> & r, 
                                                double & v);
template void cqmc::mpi_unbiased_ratio_of_means(const int n, 
                                                const double * const p, 
                                                const std::complex<double> * const f, 
                                                const double * const g, 
                                                const bool correct, 
                                                std::complex<double> & r, 
                                                double & v);


///////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief   rounds a double to the nearest integer
///
/// \param[in]    d          the number to round
///
/// \return the integer nearst to d
///////////////////////////////////////////////////////////////////////////////////////////////////
int cqmc::my_round(const double d) { return ( d >= 0 ? int(d + 0.5) : int(d - 0.5) ); }
