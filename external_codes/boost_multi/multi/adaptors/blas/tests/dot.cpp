#ifdef COMPILATION_INSTRUCTIONS
c++ -std=c++14 -Wall -Wextra -Wpedantic -lblas -DADD_ $0 -o $0x.x  && time $0x.x $@ && rm -f $0x.x; exit
#endif

#include "../../blas.hpp"

#include "../../../array.hpp"

#include<complex>
#include<cassert>

using std::cout;
namespace multi = boost::multi;

int main(){
	{
		multi::array<double, 2> CA = {
			{1.,  2.,  3.,  4.},
			{5.,  6.,  7.,  8.},
			{9., 10., 11., 12.}
		};
		using multi::blas::dot;
		auto d = dot(CA[1], CA[2]); // ddot
		assert(d == std::inner_product(begin(CA[1]), begin(CA[2]), end(CA[1]), 0.));
	}
	{
		multi::array<float, 1> const A = {1.,2.,3.};
		multi::array<float, 1> const B = {1.,2.,3.};
		using multi::blas::dot;
		auto f = dot(A, B); // sdot
		assert( f == std::inner_product(begin(A), end(A), begin(B), float{0}) );
	}
	{
		multi::array<float, 1> const A = {1.,2.,3.};
		multi::array<float, 1> const B = {1.,2.,3.};
		using multi::blas::dot;
		auto d = dot<double>(A, B); // dsdot
		assert( d == std::inner_product(begin(A), end(A), begin(B), double{0}) );
	}
	{
		multi::array<float, 1> const A = {1.,2.,3.};
		multi::array<float, 1> const B = {1.,2.,3.};
		using multi::blas::dot;
		auto d = dot<double>(begin(A) + 1, end(A), begin(B) + 1); // dsdot, mixed precision can be called explicitly
		static_assert( std::is_same<decltype(d), double>{}, "!");
		assert( d == std::inner_product(begin(A) + 1, end(A), begin(B) + 1, double{0}) );
	}
	{
		multi::array<float, 1> const A = {1.,2.,3.};
		multi::array<float, 1> const B = {1.,2.,3.};
		using multi::blas::dot;
		auto f = dot(begin(A) + 1, end(A), begin(B) + 1); // sdot
		static_assert( std::is_same<decltype(f), float>{}, "!" ); 
		assert( f == std::inner_product(begin(A) + 1, end(A), begin(B) + 1, float{0}) );
	}
	{
		using Z = std::complex<double>; Z const I(0.,1.);
		multi::array<Z, 1> const A = {I,2.*I,3.*I};
		using multi::blas::dotu; // zdotc
		assert( dotu(A, A) == std::inner_product(begin(A), end(A), begin(A), std::complex<double>(0)) );
	}
	{
		using Z = std::complex<double>; Z const I(0.,1.);
		multi::array<Z, 1> const A = {I,2.*I,3.*I};
		using multi::blas::dotc;
		std::cout << dotc(A, A) << std::endl; // zdotc
		assert( dotc(A, A) == std::inner_product(begin(A), end(A), begin(A), std::complex<double>(0), std::plus<>{}, [](auto&& a, auto&& b){return a*conj(b);}) );
	}
	{
		using Z = std::complex<double>; Z const I(0.,1.);
		multi::array<Z, 1> const A = {I,2.*I,3.*I};
		using multi::blas::dot; 
		std::cout << dot(A, A) << std::endl; // zdotc (dot defaults to dotc for complex)
		assert( dot(A, A) == std::inner_product(begin(A), end(A), begin(A), std::complex<double>(0), std::plus<>{}, [](auto&& a, auto&& b){return a*conj(b);}) );
	}
	{
		multi::array<float, 1> const A = {1.,2.,3.};
		multi::array<float, 1> const B = {1.,2.,3.};
		using multi::blas::dot;
		auto f = dot(1.2, A, B); // sdsdot, 1.2 is demoted to 1.2f
		assert( f == 1.2f + std::inner_product(begin(A), end(A), begin(B), float{0}) );
	}
	{
	//	multi::array<double, 1> const A = {1.,2.,3.};
	//	multi::array<double, 1> const B = {1.,2.,3.};
	//	using multi::blas::dot;
	//	auto f = dot(1.2, A, B); // this strange function is only implement for floats
	}

}

