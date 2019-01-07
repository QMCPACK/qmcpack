//////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source
// License.  See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by:
//    Lawrence Livermore National Laboratory 
//
// File created by:
// Miguel A. Morales, moralessilva2@llnl.gov 
//    Lawrence Livermore National Laboratory 
////////////////////////////////////////////////////////////////////////////////

#include<cassert>
#include <complex>
#include <type_traits>
#include <thrust/complex.h>
#include <thrust/device_ptr.h>
#include <thrust/fill.h>
#define QMC_CUDA 1
#include "Numerics/detail/cuda_utilities.hpp"
//#include "AFQMC/Kernels/strided_range.hpp"

namespace kernels 
{

void fill_n(int * first, int N, int incx, int const value)
{ thrust::fill_n(thrust::device_ptr<int>(first),N,value); }
/*
{ 
 thrust::device_ptr<int> x(first);
 strided_range<thrust::device_ptr<int> > strided(x, x+N, incx);
 thrust::fill(strided.begin(),strided.end(),value); 
}
*/

void fill_n(float * first, int N, int incx, float const value)
{ thrust::fill_n(thrust::device_ptr<float>(first),N,value); }
/*
{ 
 thrust::device_ptr<float> x(first);
 strided_range<thrust::device_ptr<float> > strided(x, x+N, incx);
 thrust::fill(strided.begin(),strided.end(),value); 
}
*/

void fill_n(double * first, int N, int incx, double const value)
{ thrust::fill_n(thrust::device_ptr<double>(first),N,value); }
/*
{ 
 thrust::device_ptr<double> x(first);
 strided_range<thrust::device_ptr<double> > strided(x, x+N, incx);
 thrust::fill(strided.begin(),strided.end(),value); 
}
*/

void fill_n(std::complex<float> * first, int N, int incx, std::complex<float> const value)
/*
{ 
 thrust::device_ptr<thrust::complex<float> > x(reinterpret_cast<thrust::complex<float> *>(first));
 strided_range<thrust::device_ptr<thrust::complex<float> > > strided(x, x+N, incx);
 thrust::fill(strided.begin(),strided.end(),static_cast<thrust::complex<float> const>(value));
}
*/
{ thrust::fill_n(thrust::device_ptr<thrust::complex<float> >(
                    reinterpret_cast<thrust::complex<float> *>(first)),N,
                    static_cast<thrust::complex<float> const >(value)); }


void fill_n(std::complex<double> * first, int N, int stride, std::complex<double> const value)
/*
{ 
 thrust::device_ptr<thrust::complex<double> > x(reinterpret_cast<thrust::complex<double> *>(first));
 strided_range<thrust::device_ptr<thrust::complex<double> > > strided(x, x+N, incx);
 thrust::fill(strided.begin(),strided.end(),static_cast<thrust::complex<double> const>(value));
}
*/
{ thrust::fill_n(thrust::device_ptr<thrust::complex<double> >(
                    reinterpret_cast<thrust::complex<double> *>(first)),N,
                    static_cast<thrust::complex<double> const >(value)); }


void fill_n(int * first, int N, int const value)
{ thrust::fill_n(thrust::device_ptr<int>(first),N,value); }
void fill_n(float * first, int N, float const value)
{ thrust::fill_n(thrust::device_ptr<float>(first),N,value); }
void fill_n(double * first, int N, double const value)
{ thrust::fill_n(thrust::device_ptr<double>(first),N,value); }
void fill_n(std::complex<float> * first, int N, std::complex<float> const value)
{ thrust::fill_n(thrust::device_ptr<thrust::complex<float> >(
                    reinterpret_cast<thrust::complex<float> *>(first)),N,
                    static_cast<thrust::complex<float> const >(value)); }
void fill_n(std::complex<double> * first, int N, std::complex<double> const value)
{ thrust::fill_n(thrust::device_ptr<thrust::complex<double> >(
                    reinterpret_cast<thrust::complex<double> *>(first)),N,
                    static_cast<thrust::complex<double> const >(value)); }


}

