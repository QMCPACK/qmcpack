#ifdef COMPILATION_INSTRUCTIONS
(echo '#include"'$0'"'>$0.cpp)&&c++ -Wall -D_TEST_MULTI_ADAPTORS_BLAS_FILLING $0.cpp -o $0x `pkg-config --cflags --libs blas` -lboost_unit_test_framework&&$0x&&rm $0x $0.cpp;exit
#endif
// © Alfredo A. Correa 2019

#ifndef MULTI_ADAPTORS_BLAS_FILLING_HPP
#define MULTI_ADAPTORS_BLAS_FILLING_HPP

#include    "../blas/core.hpp"
#include    "../blas/operations.hpp"
#include "../../array_ref.hpp"

namespace boost{
namespace multi{
namespace blas{

enum class uplo  : char{L='L', U='U'};

enum class filling : char{
	lower = static_cast<char>(uplo::U),
	upper = static_cast<char>(uplo::L)
};

filling flip(filling side){
	switch(side){
		case filling::lower: return filling::upper;
		case filling::upper: return filling::lower;
	} __builtin_unreachable();
}

filling operator-(filling side){return flip(side);}
filling operator+(filling side){return side;}

template<class A2D>
filling detect_triangular_aux(A2D const& A, std::false_type){
	{
		for(auto i = size(A); i != 0; --i){
			auto const asum_up = blas::asum(begin(A[i-1])+i, end(A[i-1]));
			if(asum_up!=asum_up) return filling::lower;
			else if(asum_up!=0.) return filling::upper;

			auto const asum_lo = blas::asum(begin(rotated(A)[i-1])+i, end(rotated(A)[i-1]));
			if(asum_lo!=asum_lo) return filling::upper;
			else if(asum_lo!=0.) return filling::lower;
		}
	}
	return filling::lower;
}

template<class A2D>
filling detect_triangular_aux(A2D const& A, std::true_type);

template<class A2D>
filling detect_triangular(A2D const& A){
#if __cpp_if_constexpr>=201606
	if constexpr(not is_hermitized<A2D>()){
		using blas::asum;
		for(auto i = size(A); i != 0; --i){
			auto const asum_up = asum(A[i-1]({i, A[i-1].size()}));
			if(asum_up!=asum_up) return filling::lower;
			else if(asum_up!=0.) return filling::upper;

			auto const asum_lo = asum(rotated(A)[i-1]({i, rotated(A)[i-1].size()}));
			if(asum_lo!=asum_lo) return filling::upper;
			else if(asum_lo!=0.) return filling::lower;
		}
	}else{
		return flip(detect_triangular(hermitized(A)));
	}
	return filling::lower;
#else
	return detect_triangular_aux(A, is_conjugated_t<A2D>{});//std::integral_constant<bool, not is_hermitized<A2D>()>{});
#endif
}

template<class A2D>
filling detect_triangular_aux(A2D const& A, std::true_type){
	return flip(detect_triangular(hermitized(A)));
}

}}

}

#if _TEST_MULTI_ADAPTORS_BLAS_FILLING

#define BOOST_TEST_MODULE "C++ Unit Tests for Multi adaptors side"
#define BOOST_TEST_DYN_LINK
#include<boost/test/unit_test.hpp>

#include "../../array.hpp"
#include "../../utility.hpp"
#include "../blas/nrm2.hpp"

#include<complex>
#include<cassert>
#include<iostream>
#include<numeric>
#include<algorithm>

using std::cout;

template<class M> 
decltype(auto) print(M const& C){
	using boost::multi::size;
	for(int i = 0; i != size(C); ++i){
		for(int j = 0; j != size(C[i]); ++j) cout<< C[i][j] <<' ';
		cout<<std::endl;
	}
	return cout<<"---"<<std::endl;
}

namespace multi = boost::multi;
using complex = std::complex<double>;
auto const I = complex(0., 1.);

template<class T> void what();

BOOST_AUTO_TEST_CASE(multi_adaptors_blas_side){
	return;
}

#endif
#endif

