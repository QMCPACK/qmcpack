//////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source
// License.  See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by:
// Alfredo Correa, correaa@llnl.gov 
//    Lawrence Livermore National Laboratory 
// Miguel A. Morales, moralessilva2@llnl.gov 
//    Lawrence Livermore National Laboratory 
//
// File created by:
// Alfredo Correa, correaa@llnl.gov 
//    Lawrence Livermore National Laboratory 
////////////////////////////////////////////////////////////////////////////////

#ifndef MA_BLAS_HPP
#define MA_BLAS_HPP

#include "AFQMC/Numerics/detail/blas.hpp"
#include<utility> //std::enable_if
#include<cassert>
#include<iostream>

#include "AFQMC/Numerics/ma_blas_extensions.hpp"

namespace ma{

template<class MultiArray1DX,
         class MultiArray1DY,
         typename = typename std::enable_if_t<std::decay<MultiArray1DX>::type::dimensionality == 1>,
         typename = typename std::enable_if_t<std::decay<MultiArray1DY>::type::dimensionality == 1>
        >
MultiArray1DY copy(MultiArray1DX&& x, MultiArray1DY&& y){
        copy(x.size(), x.origin(), x.strides()[0], y.origin(), y.strides()[0]);
        return std::forward<MultiArray1DY>(y);
}

template<class MultiArray2DX,
         class MultiArray2DY,
         typename = typename std::enable_if_t<std::decay<MultiArray2DX>::type::dimensionality == 2>,
         typename = typename std::enable_if_t<std::decay<MultiArray2DY>::type::dimensionality == 2>,
         typename = void
        >
MultiArray2DY copy(MultiArray2DX&& x, MultiArray2DY&& y){
        assert( x.strides()[0] == x.shape()[1] ); // only on contiguous arrays 
        assert( x.strides()[1] == 1 );            // only on contiguous arrays 
        copy(x.num_elements(), x.origin(), 1, y.origin(), 1);
        return std::forward<MultiArray2DY>(y);
}

template<class MultiArray1Dx, 
         class MultiArray1Dy,
         typename = typename std::enable_if<std::decay<MultiArray1Dx>::type::dimensionality == 1>::type,
         typename = typename std::enable_if<std::decay<MultiArray1Dy>::type::dimensionality == 1>::type 
>
//auto 
typename std::decay<MultiArray1Dx>::type::element
dot(MultiArray1Dx&& x, MultiArray1Dy&& y){
        assert(x.size() == y.size());
// FIX FIX FIX
        return dot(x.size(), std::addressof(*x.origin()), x.strides()[0], std::addressof(*y.origin()), y.strides()[0]);
}

template<class T, class MultiArray1D, typename = typename std::enable_if<std::decay<MultiArray1D>::type::dimensionality == 1>::type >
MultiArray1D scal(T a, MultiArray1D&& x){
	scal(x.size(), a, x.origin(), x.strides()[0]);
	return std::forward<MultiArray1D>(x);
}

template<class T, class MultiArray1D>
auto operator*=(MultiArray1D&& x, T a) -> decltype(scal(a, std::forward<MultiArray1D>(x))){
	return scal(a, std::forward<MultiArray1D>(x));
}

template<class T, class MultiArray1DA, class MultiArray1DB, 
	typename = typename std::enable_if<std::decay<MultiArray1DA>::type::dimensionality == 1 and std::decay<MultiArray1DB>::type::dimensionality == 1>::type
>
MultiArray1DB axpy(T x, MultiArray1DA const& a, MultiArray1DB&& b){
	assert( a.shape()[0] == b.shape()[0] );
	axpy(a.shape()[0], x, a.origin(), a.strides()[0], b.origin(), b.strides()[0]);
	return std::forward<MultiArray1DB>(b);
}

template<class T, class MultiArray2DA, class MultiArray2DB,
        typename = typename std::enable_if<std::decay<MultiArray2DA>::type::dimensionality == 2 and std::decay<MultiArray2DB>::type::dimensionality == 2>::type,
        typename = void // TODO change to use dispatch 
>
MultiArray2DB axpy(T x, MultiArray2DA const& a, MultiArray2DB&& b){
        assert( a.num_elements() == b.num_elements() );
        assert( a.strides()[0] == a.shape()[1] ); // only on contiguous arrays 
        assert( a.strides()[1] == 1 );            // only on contiguous arrays 
        assert( b.strides()[0] == b.shape()[1] ); // only on contiguous arrays 
        assert( b.strides()[1] == 1 );            // only on contiguous arrays 
        axpy(a.num_elements(), x, a.origin(), 1, b.origin(), 1);
        return std::forward<MultiArray2DB>(b);
}

template<char IN, class T, class MultiArray2DA, class MultiArray1DX, class MultiArray1DY,
	typename = typename std::enable_if< MultiArray2DA::dimensionality == 2 and MultiArray1DX::dimensionality == 1 and std::decay<MultiArray1DY>::type::dimensionality == 1>::type
>
MultiArray1DY gemv(T alpha, MultiArray2DA const& A, MultiArray1DX const& x, T beta, MultiArray1DY&& y){
        assert( (IN == 'N') || (IN == 'T') || (IN == 'C')  );
	if(IN == 'T' or IN == 'C') assert( x.shape()[0] == A.shape()[1] and y.shape()[0] == A.shape()[0]);
	else if(IN == 'N') assert( x.shape()[0] == A.shape()[0] and y.shape()[0] == A.shape()[1]);
	assert( A.strides()[1] == 1 ); // gemv is not implemented for arrays with non-leading stride != 1
	int M = A.shape()[1];
	int N = A.shape()[0];
	gemv(IN, M, N, alpha, A.origin(), A.strides()[0], 
                   x.origin(), x.strides()[0], beta, 
                   y.origin(), y.strides()[0]);
	return std::forward<MultiArray1DY>(y);
} //y := alpha*A*x + beta*y,

template<char IN, class MultiArray2DA, class MultiArray1DX, class MultiArray1DY>
MultiArray1DY gemv(MultiArray2DA const& A, MultiArray1DX const& x, MultiArray1DY&& y){
	return gemv<IN>(1., A, x, 0., std::forward<MultiArray1DY>(y));
} //y := alpha*A*x

//	gemm<'T', 'T'>(1., A, B, 0., C); // C = T(A*B) = T(B)*T(A) or T(C) = A*B
//	gemm<'N', 'N'>(1., A, B, 0., C); // C = B*A = T(T(A)*T(B)) or T(C) = T(A)*T(B)
//	gemm<'T', 'N'>(1., A, B, 0., C); // C = T(A*T(B)) = B*T(A) or T(C) = A*T(B)
//	gemm<'N', 'T'>(1., A, B, 0., C); // C =  T(T(A)*B) = T(B)*A or T(C) = T(A)*B

template<char TA, char TB, class T, class MultiArray2DA, class MultiArray2DB, class MultiArray2DC, 
	typename = typename std::enable_if< MultiArray2DA::dimensionality == 2 and MultiArray2DB::dimensionality == 2 and std::decay<MultiArray2DC>::type::dimensionality == 2>::type
>
MultiArray2DC gemm(T alpha, MultiArray2DA const& a, MultiArray2DB const& b, T beta, MultiArray2DC&& c){
	assert( a.strides()[1] == 1 );
	assert( b.strides()[1] == 1 );
	assert( c.strides()[1] == 1 );
	assert( (TA == 'N') || (TA == 'T') || (TA == 'C')  );
	assert( (TB == 'N') || (TB == 'T') || (TB == 'C')  );
	int M = -1;
	int N = -1;
	int K = -1;
	if(TA == 'N' and TB == 'N'){
		M = a.shape()[1];
		N = b.shape()[0];
		K = a.shape()[0];
		assert(a.shape()[0] == b.shape()[1] and c.shape()[0] == b.shape()[0] and c.shape()[1] == a.shape()[1]);
	}
	if((TA == 'T' or TA == 'C') and (TB == 'T' or TB == 'C')){
		M = a.shape()[0];
		N = b.shape()[1];
		K = a.shape()[1];
		assert(a.shape()[1] == b.shape()[0] and c.shape()[0] == b.shape()[1] and c.shape()[1] == a.shape()[0]);
	}
	if((TA == 'T' or TA == 'C') and TB == 'N'){
		M = a.shape()[0];
		N = b.shape()[0];
		K = a.shape()[1];
		assert(a.shape()[1] == b.shape()[1] and c.shape()[0] == b.shape()[0] and c.shape()[1] == a.shape()[0]);
	}
	if(TA == 'N' and (TB == 'T' or TB == 'C')){
		M = a.shape()[1];
		N = b.shape()[1];
		K = a.shape()[0];
		assert(a.shape()[0] == b.shape()[0] and c.shape()[0] == b.shape()[1] and c.shape()[1] == a.shape()[1]);
	}
	gemm(
		TA, TB, 
		M, N, K, alpha, 
		a.origin(), a.strides()[0], 
		b.origin(), b.strides()[0],
		beta, 
		c.origin(), c.strides()[0]
	);
	return std::forward<MultiArray2DC>(c);
}

// Expect: A[nbatch][nrow][ncol]
template<char TA, char TB, class T, class MultiArray3DA, class MultiArray3DB, class MultiArray3DC,
        typename = typename std::enable_if< MultiArray3DA::dimensionality == 3 and
                                    MultiArray3DB::dimensionality == 3 and
                                    std::decay<MultiArray3DC>::type::dimensionality == 3>::type
        >
MultiArray3DC gemmStridedBatched(T alpha, MultiArray3DA const& a, MultiArray3DB const& b,
                                 T beta, MultiArray3DC&& c){
        assert( a.strides()[2] == 1 );
        assert( b.strides()[2] == 1 );
        assert( c.strides()[2] == 1 );
        assert( a.shape()[0] == b.shape()[0] );
        assert( a.shape()[0] == c.shape()[0] );
        assert( (TA == 'N') || (TA == 'T') || (TA == 'C')  );
        assert( (TB == 'N') || (TB == 'T') || (TB == 'C')  );
        int M = -1;
        int N = -1;
        int K = -1;
        if(TA == 'N' and TB == 'N'){
                M = a.shape()[2];
                N = b.shape()[1];
                K = a.shape()[1];
                assert(a.shape()[1] == b.shape()[2] and c.shape()[1] == b.shape()[1] and c.shape()[2] == a.shape()[2]);
        }
        if((TA == 'T' or TA == 'C') and (TB == 'T' or TB == 'C')){
                M = a.shape()[1];
                N = b.shape()[2];
                K = a.shape()[2];
                assert(a.shape()[2] == b.shape()[1] and c.shape()[1] == b.shape()[2] and c.shape()[2] == a.shape()[1]);
        }
        if((TA == 'T' or TA == 'C') and TB == 'N'){
                M = a.shape()[1];
                N = b.shape()[1];
                K = a.shape()[2];
                assert(a.shape()[2] == b.shape()[2] and c.shape()[1] == b.shape()[1] and c.shape()[2] == a.shape()[1]);
        }
        if(TA == 'N' and (TB == 'T' or TB == 'C')){
                M = a.shape()[2];
                N = b.shape()[2];
                K = a.shape()[1];
                assert(a.shape()[1] == b.shape()[1] and c.shape()[1] == b.shape()[2] and c.shape()[2] == a.shape()[2]);
        }
        gemmStridedBatched(
                TA, TB,
                M, N, K,
                alpha,
                a.origin(), a.strides()[1],a.strides()[0],
                b.origin(), b.strides()[1],b.strides()[0],
                beta,
                c.origin(), c.strides()[1],c.strides()[0],
                a.shape()[0]
        );
        return std::forward<MultiArray3DC>(c);
}

template<char TA, char TB, class T, class MultiArray2DA, class MultiArray2DB, class MultiArray2DC>
MultiArray2DC gemm(MultiArray2DA const& a, MultiArray2DB const& b, MultiArray2DC&& c){
	return gemm(1., a, b, 0., std::forward<MultiArray2DC>(c));
}

template<char TA, char TB, class T, class MultiArray2DA, class MultiArray2DB, class MultiArray2DC,
        typename = typename std::enable_if< MultiArray2DA::dimensionality == 2 and
                                            MultiArray2DB::dimensionality == 2 and
                                            std::decay<MultiArray2DC>::type::dimensionality == 2>::type
>
MultiArray2DC geam(T alpha, MultiArray2DA const& a, T beta, MultiArray2DB const& b, MultiArray2DC&& c){
        assert( a.strides()[1] == 1 );
        assert( b.strides()[1] == 1 );
        assert( c.strides()[1] == 1 );
        assert( (TA == 'N') || (TA == 'T') || (TA == 'C')  );
        assert( (TB == 'N') || (TB == 'T') || (TB == 'C')  );
        if(TA == 'N' and TB == 'N'){
                assert(a.shape()[0] == c.shape()[0] and a.shape()[1] == c.shape()[1]);
                assert(b.shape()[0] == c.shape()[0] and b.shape()[1] == c.shape()[1]);
        }
        if((TA == 'T' or TA == 'C') and (TB == 'T' or TB == 'C')){
                assert(a.shape()[1] == c.shape()[0] and a.shape()[0] == c.shape()[1]);
                assert(b.shape()[1] == c.shape()[0] and b.shape()[0] == c.shape()[1]);
        }
        if((TA == 'T' or TA == 'C') and TB == 'N'){
                assert(a.shape()[1] == c.shape()[0] and a.shape()[0] == c.shape()[1]);
                assert(b.shape()[0] == c.shape()[0] and b.shape()[1] == c.shape()[1]);
        }
        if(TA == 'N' and (TB == 'T' or TB == 'C')){
                assert(a.shape()[0] == c.shape()[0] and a.shape()[1] == c.shape()[1]);
                assert(b.shape()[1] == c.shape()[0] and b.shape()[0] == c.shape()[1]);
        }
        geam(   TA, TB, c.shape()[1], c.shape()[0],
                alpha, a.origin(), a.strides()[0],
                beta, b.origin(), b.strides()[0],
                c.origin(), c.strides()[0]
        );
        return std::forward<MultiArray2DC>(c);
}

template<char TA, class T, class MultiArray2DA, class MultiArray2DC,
        typename = typename std::enable_if< MultiArray2DA::dimensionality == 2 and
                                            std::decay<MultiArray2DC>::type::dimensionality == 2>::type
>
MultiArray2DC geam(T alpha, MultiArray2DA const& a, MultiArray2DC&& c){
        assert( a.strides()[1] == 1 );
        assert( c.strides()[1] == 1 );
        assert( (TA == 'N') || (TA == 'T') || (TA == 'C')  );
        if(TA == 'N'){
                assert(a.shape()[0] == c.shape()[0] and a.shape()[1] == c.shape()[1]);
        }
        if((TA == 'T' or TA == 'C')) {
                assert(a.shape()[1] == c.shape()[0] and a.shape()[0] == c.shape()[1]);
        }
        auto aorg(a.origin());
        geam(   TA, TA, c.shape()[1], c.shape()[0],
                alpha, a.origin(), a.strides()[0],
                T(0), a.origin(), a.strides()[0],
                c.origin(), c.strides()[0]
        );
        return std::forward<MultiArray2DC>(c);
}

}

#endif

