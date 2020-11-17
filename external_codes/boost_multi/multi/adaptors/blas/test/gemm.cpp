#ifdef COMPILATION// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;-*-
$CXX -D_MULTI_CUBLAS_ALWAYS_SYNC $0 -o $0x `pkg-config --libs blas` -lcudart -lcublas -lboost_unit_test_framework&&$0x&&rm $0x; exit
#endif
// © Alfredo A. Correa 2019-2020

#define BOOST_TEST_MODULE "C++ Unit Tests for Multi cuBLAS gemm"
#define BOOST_TEST_DYN_LINK
#include<boost/test/unit_test.hpp>

#include "../../../memory/adaptors/cuda/managed/ptr.hpp"

#include "../../../adaptors/blas.hpp"
#include "../../../adaptors/blas/cuda.hpp"

#include "../../../adaptors/cuda.hpp"
#include "../../../array.hpp"

namespace multi = boost::multi;

template<class M> decltype(auto) print(M const& C){
	using multi::size;
	for(int i = 0; i != size(C); ++i){
		for(int j = 0; j != size(C[i]); ++j) std::cout<< C[i][j] <<' ';
		std::cout<<std::endl;
	}
	return std::cout<<std::endl;
}

BOOST_AUTO_TEST_CASE(multi_blas_gemm_square_real){
	namespace blas = multi::blas;
	multi::array<double, 2> const a = {
		{1, 3, 4},
		{9, 7, 1},
		{1, 2, 3}
	};
	multi::array<double, 2> const b = {	
		{11, 12, 4},
		{ 7, 19, 1},
		{11, 12, 4}
	};
	{
		multi::array<double, 2> c({size(a), size(rotated(b))}, 9999);
		blas::gemm(1., a, b, 0., c);
		BOOST_REQUIRE( c[2][1] == 86 );
	}
	{
		multi::array<double, 2> c({size(a), size(rotated(b))}, 9999);
		blas::gemm(1., a, blas::T(b), 0., c);
		BOOST_REQUIRE( c[2][1] == 48 );
	}
	{
		multi::array<double, 2> c({size(a), size(rotated(b))}, 9999);
		blas::gemm(1., blas::T(a), b, 0., c);
		BOOST_REQUIRE( c[2][1] == 103 );
	}
	{
		multi::array<double, 2> c({size(a), size(rotated(b))}, 9999);
		using multi::blas::transposed; using multi::blas::gemm;
		blas::gemm(1., blas::T(a), blas::T(b), 0., c);
		BOOST_REQUIRE( c[2][1] == 50 );		
	}
	{
		multi::array<double, 2> c({size(a), size(rotated(b))}, 9999);
		using multi::blas::gemm;
		blas::gemm(1., a, b, 0., c);
		BOOST_REQUIRE( c[2][1] == 86 );
	}
	{
		multi::array<double, 2> c({size(a), size(rotated(b))}, 9999);
		blas::gemm(1., a, blas::T(b), 0., c);
		BOOST_REQUIRE( c[2][1] == 48 );
	}
	{
		multi::array<double, 2> c({size(a), size(rotated(b))}, 9999);
		blas::gemm(1., blas::T(a), b, 0., c);
		BOOST_REQUIRE( c[2][1] == 103 );		
	}
	{
		multi::array<double, 2> c({size(a), size(rotated(b))}, 9999);
		blas::gemm(1., blas::T(a), blas::T(b), 0., c);
		BOOST_REQUIRE( c[2][1] == 50 );		
	}
	{
		multi::array<double, 2> c({size(a), size(rotated(b))}, 9999);
		blas::gemm(2., blas::H(a), blas::H(b), 0., c);
		BOOST_REQUIRE( c[2][1] == 100 );

		multi::array<double, 2> const c_copy = blas::gemm(2., blas::H(a), blas::H(b));
		BOOST_REQUIRE( c == c_copy );
		multi::array<double, 2> const c_copy2 = blas::gemm(blas::H(a), blas::H(b));
		BOOST_REQUIRE( c_copy2[2][1] == 50 );
	}
}

namespace utf = boost::unit_test;

BOOST_AUTO_TEST_CASE(multi_adaptors_blas_gemm_real_nonsquare_automatic, *utf::tolerance(0.00001)){
	namespace blas = multi::blas;
	multi::array<double, 2> const a = {
		{ 1., 3., 1.},
		{ 9., 7., 1.},
	};
	multi::array<double, 2> const b = {	
		{ 11., 12., 4., 8.},
		{  7., 19., 2., 7.},
		{  5.,  3., 3., 1.}
	};
	{
		multi::array<double, 2> c({2, 4});
		blas::gemm(1., a, b, 0., c); // c=ab, c⸆=b⸆a⸆
		BOOST_REQUIRE( c[1][2] == 53 );
	}
	{
		multi::array<double, 2> c({2, 4});
		blas::gemm(0.1, a, b, 0., c); // c=ab, c⸆=b⸆a⸆
		BOOST_TEST_REQUIRE( c[1][2] == 5.3 );
	}
	{
		multi::array<double, 2> c({2, 4});
		blas::gemm(0.1, a, b, 0., c); // c=ab, c⸆=b⸆a⸆
		BOOST_TEST_REQUIRE( c[1][2] == 5.3 );
	}
	{
		auto c = blas::gemm(0.1, a, b); // c=ab, c⸆=b⸆a⸆
		BOOST_TEST_REQUIRE( c[1][2] == 5.3 );
	}
}

BOOST_AUTO_TEST_CASE(multi_adaptors_blas_gemm_real_nonsquare_automatic_cuda, *utf::tolerance(0.00001)){
	namespace cuda = multi::cuda;
	namespace blas = multi::blas;
	cuda::array<double, 2> const a = {
		{ 1., 3., 1.},
		{ 9., 7., 1.},
	};
	cuda::array<double, 2> const b = {	
		{ 11., 12., 4., 8.},
		{  7., 19., 2., 7.},
		{  5.,  3., 3., 1.}
	};
	{
		cuda::array<double, 2> c({2, 4});
		blas::gemm(1., a, b, 0., c); // c=ab, c⸆=b⸆a⸆
		BOOST_REQUIRE( c[1][2] == 53 );
	}
	{
		cuda::array<double, 2> c({2, 4});
		blas::gemm(0.1, a, b, 0., c); // c=ab, c⸆=b⸆a⸆
		BOOST_REQUIRE( c[1][2] == 5.3 );
	}
	{
		cuda::array<double, 2> c({2, 4});
		blas::gemm(0.1, a, b, 0., c); // c=ab, c⸆=b⸆a⸆
		BOOST_REQUIRE( c[1][2] == 5.3 );
	}
	{
		auto c = blas::gemm(0.1, a, b); // c=ab, c⸆=b⸆a⸆
		BOOST_REQUIRE( c[1][2] == 5.3 );
	}
}

BOOST_AUTO_TEST_CASE(multi_adaptors_blas_gemm_real_nonsquare_hermitized_second, *utf::tolerance(0.00001)){
	namespace blas = multi::blas;
	multi::array<double, 2> const a = {
		{1, 3, 1},
		{9, 7, 1},
	};
	multi::array<double, 2> const b = {	
		{11,  7, 5},
		{12, 19, 3},
		{ 4,  2, 3},
		{ 8,  7, 1}
	};
	{
		multi::array<double, 2> c({2, 4});
		blas::gemm(1., a, blas::H(b), 0., c); // c=ab, c⸆=b⸆a⸆
		BOOST_REQUIRE( c[1][2] == 53 );
	}
	{
		multi::array<double, 2> c({2, 4});
		blas::gemm(0.1, a, blas::H(b), 0., c); // c=ab, c⸆=b⸆a⸆
		BOOST_TEST( c[1][2] == 5.3 );
	}
	{
		multi::array<double, 2> c({2, 4});
		blas::gemm(0.1, a, blas::H(b), 0., c); // c=ab, c⸆=b⸆a⸆
		BOOST_TEST( c[1][2] == 5.3 );
	}
	{
		auto c = blas::gemm(0.1, a, blas::H(b)); // c=ab, c⸆=b⸆a⸆
		BOOST_TEST( c[1][2] == 5.3 );
	}
}

BOOST_AUTO_TEST_CASE(multi_adaptors_blas_gemm_real_nonsquare_hermitized_second_gpu, *utf::tolerance(0.00001)){
	namespace cuda = multi::cuda;
	namespace blas = multi::blas;
	cuda::array<double, 2> const a = {
		{1, 3, 1},
		{9, 7, 1},
	};
	cuda::array<double, 2> const b = {	
		{11,  7, 5},
		{12, 19, 3},
		{ 4,  2, 3},
		{ 8,  7, 1}
	};
	using multi::blas::gemm;using multi::blas::hermitized;
	{
		cuda::array<double, 2> c({2, 4});
		blas::gemm(1., a, blas::H(b), 0., c); // c=ab, c⸆=b⸆a⸆
		BOOST_REQUIRE( c[1][2] == 53 );
	}
	{
		cuda::array<double, 2> c({2, 4});
		blas::gemm(0.1, a, blas::H(b), 0., c); // c=ab, c⸆=b⸆a⸆
		BOOST_REQUIRE( c[1][2] == 5.3 );
	}
	{
		cuda::array<double, 2> c({2, 4});
		blas::gemm(0.1, a, blas::H(b), 0., c); // c=ab, c⸆=b⸆a⸆
		BOOST_REQUIRE( c[1][2] == 5.3 );
	}
	{
		cuda::array<double, 2> c({2, 4});
		auto c_copy = blas::gemm(0.1, a, blas::H(b), 0., c); // c=ab, c⸆=b⸆a⸆
		BOOST_REQUIRE( c_copy[1][2] == 5.3 );
	}
	{
		multi::cuda::array<double, 2> c({2, 4});
		auto c_copy = blas::gemm(0.1, a, blas::H(b), 0., c); // c=ab, c⸆=b⸆a⸆
		BOOST_REQUIRE( c_copy[1][2] == 5.3 );
	}
	{
	//	auto f = [](auto&& a, auto&& b){return blas::gemm(0.1, a, blas::H(b));};
	//	auto c = f(a, b);
	//	BOOST_REQUIRE( c[1][2] == 5.3 );
	}
	{
	//	auto f = [](auto&& a, auto&& b){return blas::gemm(0.1, a, blas::H(b));};
	//	cuda::array<double, 2> c;
	//	c = f(a, b);
	//	BOOST_REQUIRE( c[1][2] == 5.3 );
	}
	{
	//	auto c = blas::gemm(0.1, a, blas::H(b)); // c=ab, c⸆=b⸆a⸆
	//	BOOST_REQUIRE( c[1][2] == 5.3 );
	}
}

BOOST_AUTO_TEST_CASE(multi_adaptors_blas_gemm_real_nonsquare_hermitized_second_managed, *utf::tolerance(0.00001)){
	namespace cuda = multi::cuda;
	namespace blas = multi::blas;
	cuda::managed::array<double, 2> const a = {
		{1, 3, 1},
		{9, 7, 1},
	};
	cuda::managed::array<double, 2> const b = {	
		{11,  7, 5},
		{12, 19, 3},
		{ 4,  2, 3},
		{ 8,  7, 1}
	};
	{
		cuda::managed::array<double, 2> c({2, 4});
		blas::gemm(1., a, blas::H(b), 0., c); // c=ab, c⸆=b⸆a⸆
		BOOST_REQUIRE( c[1][2] == 53 );
	}
	{
		cuda::managed::array<double, 2> c({2, 4});
		blas::gemm(0.1, a, blas::H(b), 0., c); // c=ab, c⸆=b⸆a⸆
		BOOST_REQUIRE( c[1][2] == 5.3 );
	}
	{
		cuda::managed::array<double, 2> c({2, 4});
		blas::gemm(0.1, a, blas::H(b), 0., c); // c=ab, c⸆=b⸆a⸆
		BOOST_REQUIRE( c[1][2] == 5.3 );
	}
	{
		cuda::managed::array<double, 2> c({2, 4});
		auto c_copy = blas::gemm(0.1, a, blas::H(b), 0., c); // c=ab, c⸆=b⸆a⸆
		BOOST_REQUIRE( c_copy[1][2] == 5.3 );
	}
	{
		multi::cuda::managed::array<double, 2> c({2, 4});
		auto c_copy = blas::gemm(0.1, a, blas::H(b), 0., c); // c=ab, c⸆=b⸆a⸆
		BOOST_REQUIRE( c_copy[1][2] == 5.3 );
	}
	{
	//	auto f = [](auto&& a, auto&& b){return blas::gemm(0.1, a, blas::H(b));};
	//	auto c = f(a, b);
	//	BOOST_REQUIRE( c[1][2] == 5.3 );
	}
	{
	//	auto f = [](auto&& a, auto&& b){return blas::gemm(0.1, a, blas::H(b));};
	//	cuda::managed::array<double, 2> c;
	//	c = f(a, b);
	//	BOOST_REQUIRE( c[1][2] == 5.3 );
	}
	{
		auto f = [](auto&& a, auto&& b){return blas::gemm(0.1, a, blas::H(b));};
		multi::cuda::managed::array<double, 2> c = a; BOOST_REQUIRE(size(c) == 2 and size(rotated(c)) == 3);
		c = f(a, b);
		BOOST_REQUIRE( c[1][2] == 5.3 ); BOOST_REQUIRE(size(c) == 2 and size(rotated(c)) == 4);
	}
	{
		auto c = blas::gemm(0.1, a, blas::H(b)); // c=ab, c⸆=b⸆a⸆
		BOOST_REQUIRE( c[1][2] == 5.3 );
	}
}

using complex = std::complex<double>; complex const I{0, 1};

BOOST_AUTO_TEST_CASE(multi_adaptors_blas_gemm_complex_nonsquare_automatic){
	namespace cuda = multi::cuda;
	namespace blas = multi::blas;
	multi::array<complex, 2> const a = {
		{ 1. + 2.*I, 3. - 3.*I, 1.-9.*I},
		{ 9. + 1.*I, 7. + 4.*I, 1.-8.*I},
	};
	multi::array<complex, 2> const b = {	
		{ 11.+1.*I, 12.+1.*I, 4.+1.*I, 8.-2.*I},
		{  7.+8.*I, 19.-2.*I, 2.+1.*I, 7.+1.*I},
		{  5.+1.*I,  3.-1.*I, 3.+8.*I, 1.+1.*I}
	};
	{
		multi::array<complex, 2> c({2, 4});
		blas::gemm(1., a, b, 0., c); // c=ab, c⸆=b⸆a⸆
		BOOST_REQUIRE( c[1][2] == complex(112, 12) );
	}
	{
		cuda::array<complex, 2> const acu = a;
		cuda::array<complex, 2> const bcu = b;
		cuda::array<complex, 2> ccu({2, 4});
		blas::gemm(complex(1.), acu, bcu, complex(0.), ccu);
		BOOST_REQUIRE( ccu[1][2] == complex(112, 12) );
	}
	{
		cuda::managed::array<complex, 2> const amcu = a;
		cuda::managed::array<complex, 2> const bmcu = b;
		cuda::managed::array<complex, 2> cmcu({2, 4});
		blas::gemm(1., amcu, bmcu, 0., cmcu);
		BOOST_REQUIRE( cmcu[1][2] == complex(112, 12) );
	}
}

struct multiplies_bind1st{
	multiplies_bind1st(multi::cuda::managed::array<complex, 2>&& m) : m_(std::move(m)){}
	template<class A>
	auto operator()(A const& a) const{
		using multi::blas::gemm;
		return gemm(m_, a);
	}
private:
	multi::cuda::managed::array<complex, 2> m_;
};

BOOST_AUTO_TEST_CASE(multi_constructors_inqnvcc_bug){
	namespace cuda = multi::cuda;
	namespace blas = multi::blas;
	cuda::managed::array<complex, 2> m = {
		{ 1. + 2.*I, 3. - 3.*I, 1.-9.*I},
		{ 9. + 1.*I, 7. + 4.*I, 1.-8.*I},
	};
	cuda::managed::array<complex, 2> const b = {	
		{ 11.+1.*I, 12.+1.*I, 4.+1.*I, 8.-2.*I},
		{  7.+8.*I, 19.-2.*I, 2.+1.*I, 7.+1.*I},
		{  5.+1.*I,  3.-1.*I, 3.+8.*I, 1.+1.*I}
	};
	auto c = blas::gemm(m, b);
	BOOST_REQUIRE( c[1][2] == complex(112, 12) );
	BOOST_REQUIRE( b[1][2] == 2.+1.*I );

	auto m_as_operator2 = [&](auto const& B){return blas::gemm(m, B);};
	auto c2 = m_as_operator2(b);
	BOOST_REQUIRE( c == c2 );

	auto m_as_operator3 = [=](auto const& B){return blas::gemm(m, B);};
	auto c3 = m_as_operator3(b);
	BOOST_REQUIRE( c == c3 );

	multiplies_bind1st m_as_operator4(std::move(m));
	auto c4 = m_as_operator4(b);
	BOOST_REQUIRE( c == c4 );
	BOOST_REQUIRE( is_empty(m) );
}

BOOST_AUTO_TEST_CASE(multi_blas_gemm_elongated){
	namespace blas = multi::blas;
	multi::array<complex, 2> const a = {
		{1.-2.*I, 9.-1.*I}
	};
	BOOST_REQUIRE( size(a) == 1 and size(a[0]) == 2 );
	multi::array<complex, 2> const b = {
		{2. + 3.*I}, 
		{19.+11.*I}
	};
	BOOST_REQUIRE( size(b) == 2 and size(b[0]) == 1 );
	{
		multi::array<complex, 2> c({1, 1});
		blas::gemm(1., a, b, 0., c); // c=ab, c⸆=b⸆a⸆
		BOOST_REQUIRE( c[0][0] == a[0][0]*b[0][0] + a[0][1]*b[1][0] );
	}
	{
		multi::array<complex, 2> c({2, 2});
		blas::gemm(1., b, a, 0., c); // c=ba, c⸆=a⸆b⸆
		BOOST_REQUIRE( c[0][0] == b[0][0]*a[0][0] );
		BOOST_TEST( c[1][1] == b[1][0]*a[0][1] );
		auto const c_copy = blas::gemm(1., b, a);
		BOOST_REQUIRE( c_copy == c );
	}
	{
		multi::array<complex, 2> c({1, 1});
		blas::gemm(1., a, blas::H(a), 0., c); // c=ab, c⸆=b⸆a⸆
		BOOST_REQUIRE( c[0][0] == a[0][0]*conj(a[0][0]) + a[0][1]*conj(a[0][1]) );
		auto const c_copy = blas::gemm(1., a, blas::H(a));
		BOOST_REQUIRE( c_copy == c );
	}
	{
		multi::array<complex, 2> c({2, 2});
		blas::gemm(1., blas::H(a), a, 0., c); // c=ab, c⸆=b⸆a⸆
		BOOST_REQUIRE( c[0][0] == a[0][0]*conj(a[0][0]) );
		auto const c_copy = blas::gemm(blas::H(a), a);
		BOOST_REQUIRE( c_copy == c );
	}
	{
		multi::array<complex, 2> const a = {
			{1.-2.*I, 9.-1.*I}
		};
		multi::array<complex, 2> const b = {
			{2.+3.*I}, 
			{19.+11.*I}
		};
		multi::array<complex, 2> c({1, 1});
		blas::gemm(1., a, b, 0., c);
	}
	{
		multi::array<double, 2> const a = {{2., 3.}};
		multi::array<double, 2> const b = {{4., 5.}};
		multi::array<double, 2> c({1, 1});
		blas::gemm(1., a, blas::T(b), 0., c); // blas error
	}
	{
		multi::array<double, 2> const a = {
			{2.},
			{3.},
			{5.}
		};
		multi::array<double, 2> const b = {
			{4.},
			{5.},
			{6.}
		};
		multi::array<double, 2> c1({1, 1}), c2({1, 1});
		auto ra = rotated(a).decay();
		blas::gemm(1., ra, b, 0., c1); // ok
		BOOST_REQUIRE( c1[0][0] == a[0][0]*b[0][0] + a[1][0]*b[1][0] + a[2][0]*b[2][0] );

	//	gemm(1., rotated(a), b, 0., c2); // was blas error
	//	BOOST_REQUIRE(c1 == c2); // not reached
	}
	{
		multi::array<double, 2> const a = {
			{2.},
			{3.},
			{5.}
		};
		multi::array<double, 2> const b = {
			{4., 2.},
			{5., 1.},
			{6., 2.}
		};
		multi::array<double, 2> c1({1, 2}), c2({1, 2});
		auto ra = rotated(a).decay();

		blas::gemm(1., ra, b, 0., c1); // ok
		BOOST_REQUIRE( c1[0][0] == a[0][0]*b[0][0] + a[1][0]*b[1][0] + a[2][0]*b[2][0] );

		blas::gemm(1., blas::T(a), b, 0., c2);
		BOOST_REQUIRE(c1 == c2);
	}
	if(0){
		multi::array<double, 2> const a = {
			{2.},
			{3.},
			{5.}
		};
		multi::array<double, 2> const b = {
			{4.},
			{5.},
			{6.}
		};
		multi::array<double, 2> c1({1, 1}), c2({1, 1});
		auto ra = rotated(a).decay();
		blas::gemm(1., ra, b, 0., c1); // ok
		BOOST_REQUIRE( c1[0][0] == a[0][0]*b[0][0] + a[1][0]*b[1][0] + a[2][0]*b[2][0] );

		blas::gemm(1., blas::T(a), b, 0., c2);
		BOOST_REQUIRE(c1 == c2);
	}
	if(0){
		multi::array<complex, 2> const a = {
			{2. + 1.*I},
			{3. + 2.*I}
		};
		multi::array<complex, 2> const b = {
			{4. + 3.*I}, 
			{5. + 4.*I}
		};
		multi::array<complex, 2> c1({1, 1}), c2({1, 1});
		auto ha = blas::hermitized(a).decay(); 
		blas::gemm(1., ha, b, 0., c1); // ok
		blas::gemm(1., blas::H(a), b, 0., c2); // was blas error
		print(c1);
		print(c2);
	//	BOOST_REQUIRE(c1 == c2);
	}
	{
		multi::array<complex, 2> const a = {
			{1.-2.*I}, 
			{9.-1.*I}
		};
		multi::array<complex, 2> const b = {
			{2.+3.*I, 2. + 999.*I}, 
			{19.+11.*I, 1. + 999.*I}
		};
		multi::array<complex, 2> c1({1, 2}), c2({2, 1}, 999.);
		auto ha = blas::hermitized(a).decay();
		blas::gemm(1., ha, b, 0., c1);
		blas::gemm(1., blas::H(b), a, 0., c2);
//		print(c1);
//		print(c2);
	//	std::cout << std::endl;
		BOOST_REQUIRE( c1 == blas::hermitized(c2) );
	}
	{
	//	multi::array<complex, 2> c({1, 1});
	//	using multi::blas::gemm;
	//	using multi::blas::hermitized;
//		gemm(1., hermitized(b), b, 0., c);
	//	BOOST_REQUIRE( c[0][0] == b[0][0]*conj(b[0][0]) + b[1][0]*conj(b[1][0]) );
	}
}

template<class A> void what(A&&) = delete;
template<class... A> void what() = delete;

BOOST_AUTO_TEST_CASE(multi_adaptors_blas_gemm_complex_nonsquare_automatic2){
	namespace cuda = multi::cuda;
	namespace blas = multi::blas;
	multi::array<complex, 2> const a = {
		{1.-2.*I, 9.-1.*I},
		{3.+3.*I, 7.-4.*I},
		{1.+9.*I, 1.+8.*I}
	};
	multi::array<complex, 2> const b = {	
		{ 11.+1.*I, 12.+1.*I, 4.+1.*I, 8.-2.*I},
		{  7.+8.*I, 19.-2.*I, 2.+1.*I, 7.+1.*I},
		{  5.+1.*I,  3.-1.*I, 3.+8.*I, 1.+1.*I}
	};
	{
		multi::array<complex, 2> c({2, 4});
		blas::gemm(1., blas::H(a), b, 0., c); // c=ab, c⸆=b⸆a⸆
		BOOST_REQUIRE( c[1][2] == complex(112, 12) );

		multi::array<complex, 2> const c_copy  = blas::gemm(1., blas::H(a), b);
		multi::array<complex, 2> const c_copy2 = blas::gemm(blas::H(a), b);
		BOOST_REQUIRE(( c == c_copy and c == c_copy2 ));
	}
	{
		cuda::array<complex, 2> const acu = a;
		cuda::array<complex, 2> const bcu = b;
		cuda::array<complex, 2> ccu({2, 4}, acu.get_allocator());
		blas::gemm(1., blas::H(acu), bcu, 0., ccu);
		BOOST_REQUIRE( ccu[1][2] == complex(112, 12) );

	//	what(base(acu));
	//	what(blas::H(acu).get_allocator());
	//	what<decltype(blas::H(acu))::decay_type, decltype(acu)::decay_type>();
	//	cuda::array<complex, 2> const ccu_copy  = blas::gemm(1., blas::H(acu), bcu);
	//	cuda::array<complex, 2> const ccu_copy2 = blas::gemm(blas::H(acu), bcu);
	//	BOOST_REQUIRE(( ccu_copy == ccu and ccu_copy2 == ccu ));
	}
#if 0
	{
		cuda::managed::array<complex, 2> const amcu = a;
		cuda::managed::array<complex, 2> const bmcu = b;
		cuda::managed::array<complex, 2> cmcu({2, 4});
		blas::gemm(1., blas::H(amcu), bmcu, 0., cmcu);
		BOOST_REQUIRE( cmcu[1][2] == complex(112, 12) );

	//	[](void*){}();

		cuda::managed::array<complex, 2> const cmcu_copy  = blas::gemm(1., blas::H(amcu), bmcu);
		cuda::managed::array<complex, 2> const cmcu_copy2 = blas::gemm(blas::H(amcu), bmcu);
		BOOST_REQUIRE(( cmcu_copy == cmcu and cmcu_copy2 == cmcu ));
	}
#endif
}

BOOST_AUTO_TEST_CASE(multi_adaptors_blas_gemm_complex_nonsquare_automatic3){
	namespace cuda = multi::cuda;
	namespace blas = multi::blas;
	multi::array<complex, 2> const a = {
		{1.-2.*I, 9.-1.*I},
		{3.+3.*I, 7.-4.*I},
		{1.+9.*I, 1.+8.*I}
	};
	multi::array<complex, 2> const bH = {	
		{ 11.+1.*I, 12.+1.*I, 4.+1.*I, 8.-2.*I},
		{  7.+8.*I, 19.-2.*I, 2.+1.*I, 7.+1.*I},
		{  5.+1.*I,  3.-1.*I, 3.+8.*I, 1.+1.*I}
	};
	multi::array<complex, 2> const b = multi::blas::hermitized(bH);
	{
		multi::array<complex, 2> c({2, 4});
		blas::gemm(1., blas::H(a), blas::H(b), 0., c); // c=ab, c⸆=b⸆a⸆
		BOOST_REQUIRE( c[1][2] == complex(112, 12) );
	}
	{
		cuda::array<complex, 2> const acu = a;
		cuda::array<complex, 2> const bcu = b;
		cuda::array<complex, 2> ccu({2, 4});
		blas::gemm(1., blas::H(acu), blas::H(bcu), 0., ccu);
		BOOST_REQUIRE( ccu[1][2] == complex(112, 12) );
	}
	{
		cuda::managed::array<complex, 2> const amcu = a;
		cuda::managed::array<complex, 2> const bmcu = b;
		cuda::managed::array<complex, 2> cmcu({2, 4});
		blas::gemm(1., blas::H(amcu), blas::H(bmcu), 0., cmcu);
		BOOST_REQUIRE( cmcu[1][2] == complex(112, 12) );
	}
}

BOOST_AUTO_TEST_CASE(multi_adaptors_blas_gemm_complex_nonsquare_automatic4){
	namespace blas = multi::blas;

	multi::array<complex, 2> c({12, 12});
	{
		multi::array<complex, 2> const a({12, 100}, 1.+2.*I);
		multi::array<complex, 2> const b({12, 100}, 1.+2.*I);
		using multi::blas::hermitized;
		using multi::blas::gemm;
		blas::gemm(1., a, blas::H(b), 0., c);
		BOOST_REQUIRE( real(c[0][0]) > 0);

		auto c_copy = blas::gemm(1., a, blas::H(b));
		BOOST_REQUIRE( c_copy == c );
	}
	{
		multi::array<complex, 2> const a_block({24, 100}, 1.+2.*I);
		multi::array<complex, 2> const b({12, 100}, 1.+2.*I);
		multi::array<complex, 2> c2({12, 12});

		blas::gemm(1., a_block.strided(2), blas::H(b), 0., c2);

		BOOST_REQUIRE( real(c[0][0]) > 0);
		BOOST_REQUIRE( c == c2 );

		auto c2_copy = blas::gemm(1., a_block.strided(2), blas::H(b));
		BOOST_REQUIRE( c2_copy == c2 );
	}
}

BOOST_AUTO_TEST_CASE(multi_blas_gemm_complex_issue68){
	namespace cuda = multi::cuda;
	namespace blas = multi::blas;
	cuda::managed::array<complex, 2> const a = {
		{1.-2.*I},
		{3.+3.*I},
		{1.+9.*I}
	};
	{
		cuda::managed::array<complex, 2> c({1, 1});
		blas::gemm(1., blas::H(a), a, 0., c);
		BOOST_REQUIRE( c[0][0] == 105. + 0.*I );
	}
#if 0
	{
		auto c = blas::gemm(2., blas::H(a), a);
		BOOST_REQUIRE( c[0][0] == 210. + 0.*I );
	}
	{
		auto c = blas::gemm(blas::H(a), a);

		BOOST_REQUIRE( c[0][0] == 105. + 0.*I );
	}
#endif
}

