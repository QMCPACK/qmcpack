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
        copy(x.size(), pointer_dispatch(x.origin()), x.stride(0), pointer_dispatch(y.origin()), y.stride(0));
        return std::forward<MultiArray1DY>(y);
}

template<class MultiArray2DX,
         class MultiArray2DY,
         typename = typename std::enable_if_t<std::decay<MultiArray2DX>::type::dimensionality == 2>,
         typename = typename std::enable_if_t<std::decay<MultiArray2DY>::type::dimensionality == 2>,
         typename = void
        >
MultiArray2DY copy(MultiArray2DX&& x, MultiArray2DY&& y){
        assert( x.stride(0) == x.size(1) ); // only on contiguous arrays 
        assert( x.stride(1) == 1 );            // only on contiguous arrays 
        copy(x.num_elements(), pointer_dispatch(x.origin()), 1, pointer_dispatch(y.origin()), 1);
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
        return dot(x.size(), pointer_dispatch(x.origin()), x.stride(0), pointer_dispatch(y.origin()), y.stride(0));
}

template<class MultiArray2Dx,
         class MultiArray2Dy,
         typename = typename std::enable_if<std::decay<MultiArray2Dx>::type::dimensionality == 2>::type,
         typename = typename std::enable_if<std::decay<MultiArray2Dy>::type::dimensionality == 2>::type,
         typename = void
>
typename std::decay<MultiArray2Dx>::type::element
dot(MultiArray2Dx&& x, MultiArray2Dy&& y){
        assert( x.stride(0) == x.size(1) ); // only on contiguous arrays 
        assert( x.stride(1) == 1 );            // only on contiguous arrays 
        assert( y.stride(0) == y.size(1) ); // only on contiguous arrays 
        assert( y.stride(1) == 1 );            // only on contiguous arrays 
        assert( x.num_elements() == y.num_elements());
        return dot(x.num_elements(), pointer_dispatch(x.origin()), 1, pointer_dispatch(y.origin()), 1);
}

template<class T, class MultiArray1D, typename = typename std::enable_if<std::decay<MultiArray1D>::type::dimensionality == 1>::type >
MultiArray1D scal(T a, MultiArray1D&& x){
	scal(x.size(), a, pointer_dispatch(x.origin()), x.stride(0));
	return std::forward<MultiArray1D>(x);
}

template<class T,
        class MultiArrayND,
        typename = typename std::enable_if<(std::decay<MultiArrayND>::type::dimensionality > 1)>::type,
        typename = void // TODO change to use dispatch 
    >
MultiArrayND scal(T a, MultiArrayND&& x){
#ifdef NDEBUG
        long sz(x.size(0));
        for(int i=1; i<int(std::decay<MultiArrayND>::type::dimensionality); ++i) 
          sz *= x.size(i);
        assert( x.num_elements() == sz ); 
        assert( x.stride(std::decay<MultiArrayND>::type::dimensionality-1) == 1 );            // only on contiguous arrays 
#endif
        scal(x.num_elements(), a, pointer_dispatch(x.origin()), 1);
        return std::forward<MultiArrayND>(x);
}

template<class T,
        class MultiArray3D,
        typename = typename std::enable_if<std::decay<MultiArray3D>::type::dimensionality == 3>::type,
        typename = void, // TODO change to use dispatch 
        typename = void // TODO change to use dispatch 
    >
MultiArray3D scal(T a, MultiArray3D&& x){
        assert( x.stride(0) == x.size(1)*x.size(2) ); // only on contiguous arrays 
        assert( x.stride(1) == x.size(2) ); // only on contiguous arrays 
        assert( x.stride(2) == 1 );            // only on contiguous arrays 
        scal(x.num_elements(), a, pointer_dispatch(x.origin()), 1);
        return std::forward<MultiArray3D>(x);
}

template<class T, class MultiArray1D>
auto operator*=(MultiArray1D&& x, T a) -> decltype(scal(a, std::forward<MultiArray1D>(x))){
	return scal(a, std::forward<MultiArray1D>(x));
}

template<class T, class MultiArray1DA, class MultiArray1DB, 
	typename = typename std::enable_if<std::decay<MultiArray1DA>::type::dimensionality == 1 and std::decay<MultiArray1DB>::type::dimensionality == 1>::type
>
MultiArray1DB axpy(T x, MultiArray1DA const& a, MultiArray1DB&& b){
	assert( a.size(0) == b.size(0) );
	axpy(a.size(0), x, pointer_dispatch(a.origin()), a.stride(0), pointer_dispatch(b.origin()), b.stride(0));
	return std::forward<MultiArray1DB>(b);
}

template<class T, class MultiArray2DA, class MultiArray2DB,
        typename = typename std::enable_if<std::decay<MultiArray2DA>::type::dimensionality == 2 and std::decay<MultiArray2DB>::type::dimensionality == 2>::type,
        typename = void // TODO change to use dispatch 
>
MultiArray2DB axpy(T x, MultiArray2DA const& a, MultiArray2DB&& b){
        assert( a.num_elements() == b.num_elements() );
        assert( a.stride(0) == a.size(1) ); // only on contiguous arrays 
        assert( a.stride(1) == 1 );            // only on contiguous arrays 
        assert( b.stride(0) == b.size(1) ); // only on contiguous arrays 
        assert( b.stride(1) == 1 );            // only on contiguous arrays 
        axpy(a.num_elements(), x, pointer_dispatch(a.origin()), 1, pointer_dispatch(b.origin()), 1);
        return std::forward<MultiArray2DB>(b);
}

template<char IN, class T, class MultiArray2DA, class MultiArray1DX, class MultiArray1DY,
	typename = typename std::enable_if< MultiArray2DA::dimensionality == 2 and MultiArray1DX::dimensionality == 1 and std::decay<MultiArray1DY>::type::dimensionality == 1>::type
>
MultiArray1DY gemv(T alpha, MultiArray2DA const& A, MultiArray1DX const& x, T beta, MultiArray1DY&& y){
        assert( (IN == 'N') || (IN == 'T') || (IN == 'C')  );
	if(IN == 'T' or IN == 'C') assert( x.size(0) == A.size(1) and y.size(0) == A.size(0));
	else if(IN == 'N') assert( x.size(0) == A.size(0) and y.size(0) == A.size(1));
	assert( A.stride(1) == 1 ); // gemv is not implemented for arrays with non-leading stride != 1
	int M = A.size(1);
	int N = A.size(0);
	gemv(IN, M, N, alpha, pointer_dispatch(A.origin()), A.stride(0), 
                   pointer_dispatch(x.origin()), x.stride(0), beta, 
                   pointer_dispatch(y.origin()), y.stride(0));
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
	assert( a.stride(1) == 1 );
	assert( b.stride(1) == 1 );
	assert( c.stride(1) == 1 );
	assert( (TA == 'N') || (TA == 'T') || (TA == 'C')  );
	assert( (TB == 'N') || (TB == 'T') || (TB == 'C')  );
	int M = -1;
	int N = -1;
	int K = -1;
	if(TA == 'N' and TB == 'N'){
		M = a.size(1);
		N = b.size(0);
		K = a.size(0);
		assert(a.size(0) == b.size(1) and c.size(0) == b.size(0) and c.size(1) == a.size(1));
	}
	if((TA == 'T' or TA == 'C') and (TB == 'T' or TB == 'C')){
		M = a.size(0);
		N = b.size(1);
		K = a.size(1);
		assert(a.size(1) == b.size(0) and c.size(0) == b.size(1) and c.size(1) == a.size(0));
	}
	if((TA == 'T' or TA == 'C') and TB == 'N'){
		M = a.size(0);
		N = b.size(0);
		K = a.size(1);
		assert(a.size(1) == b.size(1) and c.size(0) == b.size(0) and c.size(1) == a.size(0));
	}
	if(TA == 'N' and (TB == 'T' or TB == 'C')){
		M = a.size(1);
		N = b.size(1);
		K = a.size(0);
		assert(a.size(0) == b.size(0) and c.size(0) == b.size(1) and c.size(1) == a.size(1));
	}
	gemm(
		TA, TB, 
		M, N, K, alpha, 
		pointer_dispatch(a.origin()), a.stride(0), 
		pointer_dispatch(b.origin()), b.stride(0),
		beta, 
		pointer_dispatch(c.origin()), c.stride(0)
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
        assert( a.stride(2) == 1 );
        assert( b.stride(2) == 1 );
        assert( c.stride(2) == 1 );
        assert( a.size(0) == b.size(0) );
        assert( a.size(0) == c.size(0) );
        assert( (TA == 'N') || (TA == 'T') || (TA == 'C')  );
        assert( (TB == 'N') || (TB == 'T') || (TB == 'C')  );
        int M = -1;
        int N = -1;
        int K = -1;
        if(TA == 'N' and TB == 'N'){
                M = a.size(2);
                N = b.size(1);
                K = a.size(1);
                assert(a.size(1) == b.size(2) and c.size(1) == b.size(1) and c.size(2) == a.size(2));
        }
        if((TA == 'T' or TA == 'C') and (TB == 'T' or TB == 'C')){
                M = a.size(1);
                N = b.size(2);
                K = a.size(2);
                assert(a.size(2) == b.size(1) and c.size(1) == b.size(2) and c.size(2) == a.size(1));
        }
        if((TA == 'T' or TA == 'C') and TB == 'N'){
                M = a.size(1);
                N = b.size(1);
                K = a.size(2);
                assert(a.size(2) == b.size(2) and c.size(1) == b.size(1) and c.size(2) == a.size(1));
        }
        if(TA == 'N' and (TB == 'T' or TB == 'C')){
                M = a.size(2);
                N = b.size(2);
                K = a.size(1);
                assert(a.size(1) == b.size(1) and c.size(1) == b.size(2) and c.size(2) == a.size(2));
        }
        gemmStridedBatched(
                TA, TB,
                M, N, K,
                alpha,
                pointer_dispatch(a.origin()), a.stride(1),a.stride(0),
                pointer_dispatch(b.origin()), b.stride(1),b.stride(0),
                beta,
                pointer_dispatch(c.origin()), c.stride(1),c.stride(0),
                a.size(0)
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
        assert( a.stride(1) == 1 );
        assert( b.stride(1) == 1 );
        assert( c.stride(1) == 1 );
        assert( (TA == 'N') || (TA == 'T') || (TA == 'C')  );
        assert( (TB == 'N') || (TB == 'T') || (TB == 'C')  );
        if(TA == 'N' and TB == 'N'){
                assert(a.size(0) == c.size(0) and a.size(1) == c.size(1));
                assert(b.size(0) == c.size(0) and b.size(1) == c.size(1));
        }
        if((TA == 'T' or TA == 'C') and (TB == 'T' or TB == 'C')){
                assert(a.size(1) == c.size(0) and a.size(0) == c.size(1));
                assert(b.size(1) == c.size(0) and b.size(0) == c.size(1));
        }
        if((TA == 'T' or TA == 'C') and TB == 'N'){
                assert(a.size(1) == c.size(0) and a.size(0) == c.size(1));
                assert(b.size(0) == c.size(0) and b.size(1) == c.size(1));
        }
        if(TA == 'N' and (TB == 'T' or TB == 'C')){
                assert(a.size(0) == c.size(0) and a.size(1) == c.size(1));
                assert(b.size(1) == c.size(0) and b.size(0) == c.size(1));
        }
        geam(   TA, TB, c.size(1), c.size(0),
                alpha, pointer_dispatch(a.origin()), a.stride(0),
                beta, pointer_dispatch(b.origin()), b.stride(0),
                pointer_dispatch(c.origin()), c.stride(0)
        );
        return std::forward<MultiArray2DC>(c);
}

template<char TA, class T, class MultiArray2DA, class MultiArray2DC,
        typename = typename std::enable_if< MultiArray2DA::dimensionality == 2 and
                                            std::decay<MultiArray2DC>::type::dimensionality == 2>::type
>
MultiArray2DC geam(T alpha, MultiArray2DA const& a, MultiArray2DC&& c){
        assert( a.stride(1) == 1 );
        assert( c.stride(1) == 1 );
        assert( (TA == 'N') || (TA == 'T') || (TA == 'C')  );
        if(TA == 'N'){
                assert(a.size(0) == c.size(0) and a.size(1) == c.size(1));
        }
        if((TA == 'T' or TA == 'C')) {
                assert(a.size(1) == c.size(0) and a.size(0) == c.size(1));
        }
        geam(   TA, TA, c.size(1), c.size(0),
                alpha, pointer_dispatch(a.origin()), a.stride(0),
                T(0), pointer_dispatch(a.origin()), a.stride(0),
                pointer_dispatch(c.origin()), c.stride(0)
        );
        return std::forward<MultiArray2DC>(c);
}

}

#endif

