#ifdef COMPILATION// -*-indent-tabs-mode: t; c-basic-offset: 4; tab-width: 4;-*-
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
		using multi::blas::gemm;
		gemm(1., a, b, 0., c);
		BOOST_REQUIRE( c[2][1] == 86 );
	}
	{
		multi::array<double, 2> c({size(a), size(rotated(b))}, 9999);
		using multi::blas::gemm;
		using multi::blas::transposed;
		gemm(1., a, transposed(b), 0., c);
		BOOST_REQUIRE( c[2][1] == 48 );
	}
	{
		multi::array<double, 2> c({size(a), size(rotated(b))}, 9999);
		using multi::blas::transposed;
		using multi::blas::gemm;
		gemm(1., transposed(a), b, 0., c);
		BOOST_REQUIRE( c[2][1] == 103 );
	}
	{
		multi::array<double, 2> c({size(a), size(rotated(b))}, 9999);
		using multi::blas::transposed; using multi::blas::gemm;
		gemm(1., transposed(a), transposed(b), 0., c);
		BOOST_REQUIRE( c[2][1] == 50 );		
	}
	{
		multi::array<double, 2> c({size(a), size(rotated(b))}, 9999);
		using multi::blas::gemm;
		gemm(1., a, b, 0., c);
		BOOST_REQUIRE( c[2][1] == 86 );
	}
	{
		multi::array<double, 2> c({size(a), size(rotated(b))}, 9999);
		using multi::blas::gemm;
		gemm(1., a, rotated(b), 0., c);
		BOOST_REQUIRE( c[2][1] == 48 );
	}
	{
		multi::array<double, 2> c({size(a), size(rotated(b))}, 9999);
		using multi::blas::gemm;
		gemm(1., rotated(a), b, 0., c);
		BOOST_REQUIRE( c[2][1] == 103 );		
	}
	{
		multi::array<double, 2> c({size(a), size(rotated(b))}, 9999);
		using multi::blas::gemm;
		gemm(1., rotated(a), rotated(b), 0., c);
		BOOST_REQUIRE( c[2][1] == 50 );		
	}
	{
		multi::array<double, 2> c({size(a), size(rotated(b))}, 9999);
		using multi::blas::gemm;
		using multi::blas::hermitized;
		gemm(2., hermitized(a), hermitized(b), 0., c);
		BOOST_REQUIRE( c[2][1] == 100 );

		multi::array<double, 2> const c_copy = gemm(2., hermitized(a), hermitized(b));
		BOOST_REQUIRE( c == c_copy );
		multi::array<double, 2> const c_copy2 = gemm(hermitized(a), hermitized(b));
		BOOST_REQUIRE( c_copy2[2][1] == 50 );
	}
}

namespace utf = boost::unit_test;

BOOST_AUTO_TEST_CASE(multi_adaptors_blas_gemm_real_nonsquare_automatic, *utf::tolerance(0.00001)){
	multi::array<double, 2> const a = {
		{ 1., 3., 1.},
		{ 9., 7., 1.},
	};
	multi::array<double, 2> const b = {	
		{ 11., 12., 4., 8.},
		{  7., 19., 2., 7.},
		{  5.,  3., 3., 1.}
	};
	using multi::blas::gemm;
	{
		multi::array<double, 2> c({2, 4});
		gemm(1., a, b, 0., c); // c=ab, c⸆=b⸆a⸆
		BOOST_REQUIRE( c[1][2] == 53 );
	}
	{
		multi::array<double, 2> c({2, 4});
		gemm(0.1, a, b, 0., c); // c=ab, c⸆=b⸆a⸆
		BOOST_TEST( c[1][2] == 5.3 );
	}
	{
		multi::array<double, 2> c({2, 4});
		gemm(0.1, a, b, 0., c); // c=ab, c⸆=b⸆a⸆
		BOOST_TEST( c[1][2] == 5.3 );
	}
	{
		auto c = gemm(0.1, a, b); // c=ab, c⸆=b⸆a⸆
		BOOST_TEST( c[1][2] == 5.3 );
	}
}

BOOST_AUTO_TEST_CASE(multi_adaptors_blas_gemm_real_nonsquare_automatic_cuda, *utf::tolerance(0.00001)){
	multi::cuda::array<double, 2> const a = {
		{ 1., 3., 1.},
		{ 9., 7., 1.},
	};
	multi::cuda::array<double, 2> const b = {	
		{ 11., 12., 4., 8.},
		{  7., 19., 2., 7.},
		{  5.,  3., 3., 1.}
	};
	using multi::blas::gemm;
	{
		multi::cuda::array<double, 2> c({2, 4});
		gemm(1., a, b, 0., c); // c=ab, c⸆=b⸆a⸆
		BOOST_TEST( c[1][2] == 53 );
	}
	{
		multi::cuda::array<double, 2> c({2, 4});
		gemm(0.1, a, b, 0., c); // c=ab, c⸆=b⸆a⸆
		BOOST_TEST( c[1][2] == 5.3 );
	}
	{
		multi::cuda::array<double, 2> c({2, 4});
		gemm(0.1, a, b, 0., c); // c=ab, c⸆=b⸆a⸆
		BOOST_TEST( c[1][2] == 5.3 );
	}
	{
		auto c = gemm(0.1, a, b); // c=ab, c⸆=b⸆a⸆
		BOOST_TEST( c[1][2] == 5.3 );
	}
}

BOOST_AUTO_TEST_CASE(multi_adaptors_blas_gemm_real_nonsquare_hermitized_second, *utf::tolerance(0.00001)){
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
	using multi::blas::gemm;using multi::blas::hermitized;
	{
		multi::array<double, 2> c({2, 4});
		gemm(1., a, hermitized(b), 0., c); // c=ab, c⸆=b⸆a⸆
		BOOST_REQUIRE( c[1][2] == 53 );
	}
	{
		multi::array<double, 2> c({2, 4});
		gemm(0.1, a, hermitized(b), 0., c); // c=ab, c⸆=b⸆a⸆
		BOOST_TEST( c[1][2] == 5.3 );
	}
	{
		multi::array<double, 2> c({2, 4});
		gemm(0.1, a, hermitized(b), 0., c); // c=ab, c⸆=b⸆a⸆
		BOOST_TEST( c[1][2] == 5.3 );
	}
	{
		auto c = gemm(0.1, a, hermitized(b)); // c=ab, c⸆=b⸆a⸆
		BOOST_TEST( c[1][2] == 5.3 );
	}
}

BOOST_AUTO_TEST_CASE(multi_adaptors_blas_gemm_real_nonsquare_hermitized_second_gpu, *utf::tolerance(0.00001)){
	multi::cuda::array<double, 2> const a = {
		{1, 3, 1},
		{9, 7, 1},
	};
	multi::cuda::array<double, 2> const b = {	
		{11,  7, 5},
		{12, 19, 3},
		{ 4,  2, 3},
		{ 8,  7, 1}
	};
	using multi::blas::gemm;using multi::blas::hermitized;
	{
		multi::cuda::array<double, 2> c({2, 4});
		gemm(1., a, hermitized(b), 0., c); // c=ab, c⸆=b⸆a⸆
		BOOST_REQUIRE( c[1][2] == 53 );
	}
	{
		multi::cuda::array<double, 2> c({2, 4});
		gemm(0.1, a, hermitized(b), 0., c); // c=ab, c⸆=b⸆a⸆
		BOOST_TEST( c[1][2] == 5.3 );
	}
	{
		multi::cuda::array<double, 2> c({2, 4});
		gemm(0.1, a, hermitized(b), 0., c); // c=ab, c⸆=b⸆a⸆
		BOOST_TEST( c[1][2] == 5.3 );
	}
	{
		multi::cuda::array<double, 2> c({2, 4});
		auto c_copy = gemm(0.1, a, hermitized(b), 0., c); // c=ab, c⸆=b⸆a⸆
		BOOST_TEST( c_copy[1][2] == 5.3 );
	}
	{
		multi::cuda::array<double, 2> c({2, 4});
		auto c_copy = gemm(0.1, a, hermitized(b), 0., c); // c=ab, c⸆=b⸆a⸆
		BOOST_TEST( c_copy[1][2] == 5.3 );
	}
	{
		auto f = [](auto&& a, auto&& b){return gemm(0.1, a, hermitized(b));};
		auto c = f(a, b);
		BOOST_TEST( c[1][2] == 5.3 );
	}
	{
		auto f = [](auto&& a, auto&& b){return gemm(0.1, a, hermitized(b));};
		multi::cuda::array<double, 2> c;
		c = f(a, b);
		BOOST_TEST( c[1][2] == 5.3 );
	}
	{
		auto c = gemm(0.1, a, hermitized(b)); // c=ab, c⸆=b⸆a⸆
		BOOST_TEST( c[1][2] == 5.3 );
	}
}

BOOST_AUTO_TEST_CASE(multi_adaptors_blas_gemm_real_nonsquare_hermitized_second_managed, *utf::tolerance(0.00001)){
	multi::cuda::managed::array<double, 2> const a = {
		{1, 3, 1},
		{9, 7, 1},
	};
	multi::cuda::managed::array<double, 2> const b = {	
		{11,  7, 5},
		{12, 19, 3},
		{ 4,  2, 3},
		{ 8,  7, 1}
	};
	using multi::blas::gemm; using multi::blas::hermitized;
	{
		multi::cuda::managed::array<double, 2> c({2, 4});
		gemm(1., a, hermitized(b), 0., c); // c=ab, c⸆=b⸆a⸆
		BOOST_REQUIRE( c[1][2] == 53 );
	}
	{
		multi::cuda::managed::array<double, 2> c({2, 4});
		gemm(0.1, a, hermitized(b), 0., c); // c=ab, c⸆=b⸆a⸆
		BOOST_TEST( c[1][2] == 5.3 );
	}
	{
		multi::cuda::managed::array<double, 2> c({2, 4});
		gemm(0.1, a, hermitized(b), 0., c); // c=ab, c⸆=b⸆a⸆
		BOOST_TEST( c[1][2] == 5.3 );
	}
	{
		multi::cuda::managed::array<double, 2> c({2, 4});
		auto c_copy = gemm(0.1, a, hermitized(b), 0., c); // c=ab, c⸆=b⸆a⸆
		BOOST_TEST( c_copy[1][2] == 5.3 );
	}
	{
		multi::cuda::managed::array<double, 2> c({2, 4});
		auto c_copy = gemm(0.1, a, hermitized(b), 0., c); // c=ab, c⸆=b⸆a⸆
		BOOST_TEST( c_copy[1][2] == 5.3 );
	}
	{
		auto f = [](auto&& a, auto&& b){return gemm(0.1, a, hermitized(b));};
		auto c = f(a, b);
		BOOST_TEST( c[1][2] == 5.3 );
	}
	{
		auto f = [](auto&& a, auto&& b){return gemm(0.1, a, hermitized(b));};
		multi::cuda::managed::array<double, 2> c;
		c = f(a, b);
		BOOST_TEST( c[1][2] == 5.3 );
	}
	{
		auto f = [](auto&& a, auto&& b){return gemm(0.1, a, hermitized(b));};
		multi::cuda::managed::array<double, 2> c = a; BOOST_REQUIRE(size(c) == 2 and size(rotated(c)) == 3);
		c = f(a, b);
		BOOST_TEST( c[1][2] == 5.3 ); BOOST_REQUIRE(size(c) == 2 and size(rotated(c)) == 4);
	}
	{
		auto c = gemm(0.1, a, hermitized(b)); // c=ab, c⸆=b⸆a⸆
		BOOST_TEST( c[1][2] == 5.3 );
	}
}



using complex = std::complex<double>; complex const I{0, 1};

BOOST_AUTO_TEST_CASE(multi_adaptors_blas_gemm_complex_nonsquare_automatic){
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
		using multi::blas::gemm;
		gemm(1., a, b, 0., c); // c=ab, c⸆=b⸆a⸆
		BOOST_REQUIRE( c[1][2] == complex(112, 12) );
	}
	{
		namespace cuda = multi::cuda;
		cuda::array<complex, 2> const acu = a;
		cuda::array<complex, 2> const bcu = b;
		cuda::array<complex, 2> ccu({2, 4});
		using multi::blas::gemm;
		gemm(complex(1.), acu, bcu, complex(0.), ccu);
		BOOST_REQUIRE( ccu[1][2] == complex(112, 12) );
	}
	{
		namespace cuda = multi::cuda;
		cuda::managed::array<complex, 2> const amcu = a;
		cuda::managed::array<complex, 2> const bmcu = b;
		cuda::managed::array<complex, 2> cmcu({2, 4});
		using multi::blas::gemm;
		gemm(1., amcu, bmcu, 0., cmcu);
		BOOST_REQUIRE( cmcu[1][2] == complex(112, 12) );
	}
}




using complex = std::complex<double>;

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
	multi::cuda::managed::array<complex, 2> m = {
		{ 1. + 2.*I, 3. - 3.*I, 1.-9.*I},
		{ 9. + 1.*I, 7. + 4.*I, 1.-8.*I},
	};
	multi::cuda::managed::array<complex, 2> const b = {	
		{ 11.+1.*I, 12.+1.*I, 4.+1.*I, 8.-2.*I},
		{  7.+8.*I, 19.-2.*I, 2.+1.*I, 7.+1.*I},
		{  5.+1.*I,  3.-1.*I, 3.+8.*I, 1.+1.*I}
	};
	auto c = multi::blas::gemm(m, b);
	BOOST_REQUIRE( c[1][2] == complex(112, 12) );
	BOOST_REQUIRE( b[1][2] == 2.+1.*I );

	auto m_as_operator2 = [&](auto const& B){
		using multi::blas::gemm; return gemm(m, B);
	};
	auto c2 = m_as_operator2(b);
	BOOST_REQUIRE( c == c2 );

	auto m_as_operator3 = [=](auto const& B){
		using multi::blas::gemm; 
		return gemm(m, B);
	};
	auto c3 = m_as_operator3(b);
	BOOST_REQUIRE( c == c3 );	

	multiplies_bind1st m_as_operator4(std::move(m));
	auto c4 = m_as_operator4(b);
	BOOST_REQUIRE( c == c4 );
	BOOST_REQUIRE( is_empty(m) );
}

BOOST_AUTO_TEST_CASE(multi_blas_gemm_elongated){
	using complex = std::complex<double>; complex const I{0,1};
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
		using multi::blas::gemm;
		using multi::blas::hermitized;
		gemm(1., a, b, 0., c); // c=ab, c⸆=b⸆a⸆
		BOOST_REQUIRE( c[0][0] == a[0][0]*b[0][0] + a[0][1]*b[1][0] );
	}
	{
		multi::array<complex, 2> c({2, 2});
		using multi::blas::gemm;
		using multi::blas::hermitized;
		gemm(1., b, a, 0., c); // c=ba, c⸆=a⸆b⸆
		BOOST_REQUIRE( c[0][0] == b[0][0]*a[0][0] );
		BOOST_TEST( c[1][1] == b[1][0]*a[0][1] );
		auto const c_copy = gemm(1., b, a);
		BOOST_REQUIRE( c_copy == c );
	}
	{
		multi::array<complex, 2> c({1, 1});
		using multi::blas::gemm;
		using multi::blas::hermitized;
		gemm(1., a, hermitized(a), 0., c); // c=ab, c⸆=b⸆a⸆
		BOOST_REQUIRE( c[0][0] == a[0][0]*conj(a[0][0]) + a[0][1]*conj(a[0][1]) );
		auto const c_copy = gemm(1., a, hermitized(a));
		BOOST_REQUIRE( c_copy == c );
	}
	{
		multi::array<complex, 2> c({2, 2});
		using multi::blas::gemm;
		using multi::blas::hermitized;
		gemm(1., hermitized(a), a, 0., c); // c=ab, c⸆=b⸆a⸆
		BOOST_REQUIRE( c[0][0] == a[0][0]*conj(a[0][0]) );
		auto const c_copy = gemm(hermitized(a), a);
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
		using multi::blas::gemm;
		using multi::blas::hermitized;
		gemm(1., a, b, 0., c);
	}
	{
		multi::array<double, 2> const a = {{2., 3.}};
		multi::array<double, 2> const b = {{4., 5.}};
		multi::array<double, 2> c({1, 1});
		using multi::blas::gemm; using multi::blas::hermitized;
		gemm(1., a, rotated(b), 0., c); // blas error
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
		using multi::blas::gemm; using multi::blas::hermitized;
		auto ra = rotated(a).decay();
		gemm(1., ra, b, 0., c1); // ok
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
		using multi::blas::gemm; using multi::blas::hermitized;
		auto ra = rotated(a).decay();
		gemm(1., ra, b, 0., c1); // ok
		BOOST_REQUIRE( c1[0][0] == a[0][0]*b[0][0] + a[1][0]*b[1][0] + a[2][0]*b[2][0] );
		gemm(1., rotated(a), b, 0., c2);
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
		using multi::blas::gemm; using multi::blas::hermitized;
		auto ra = rotated(a).decay();
		gemm(1., ra, b, 0., c1); // ok
		BOOST_REQUIRE( c1[0][0] == a[0][0]*b[0][0] + a[1][0]*b[1][0] + a[2][0]*b[2][0] );
		gemm(1., rotated(a), b, 0., c2);
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
		using multi::blas::gemm; using multi::blas::hermitized;
		auto ha = hermitized(a).decay(); 
		gemm(1., ha, b, 0., c1); // ok
		gemm(1., hermitized(a), b, 0., c2); // was blas error
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
		using multi::blas::gemm;
		using multi::blas::hermitized;
		auto ha = hermitized(a).decay();
		gemm(1., ha, b, 0., c1);
		gemm(1., hermitized(b), a, 0., c2);
//		print(c1);
//		print(c2);
	//	std::cout << std::endl;
		BOOST_REQUIRE( c1 == hermitized(c2) );
	}
	{
	//	multi::array<complex, 2> c({1, 1});
	//	using multi::blas::gemm;
	//	using multi::blas::hermitized;
//		gemm(1., hermitized(b), b, 0., c);
	//	BOOST_REQUIRE( c[0][0] == b[0][0]*conj(b[0][0]) + b[1][0]*conj(b[1][0]) );
	}
}

BOOST_AUTO_TEST_CASE(multi_adaptors_blas_gemm_complex_nonsquare_automatic2){
	using complex = std::complex<double>; complex const I{0,1};
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
		using multi::blas::gemm;
		using multi::blas::hermitized;
		gemm(1., hermitized(a), b, 0., c); // c=ab, c⸆=b⸆a⸆
		BOOST_REQUIRE( c[1][2] == complex(112, 12) );

		multi::array<complex, 2> const c_copy = gemm(1., hermitized(a), b);
		multi::array<complex, 2> const c_copy2 = gemm(hermitized(a), b);
		BOOST_REQUIRE(( c == c_copy and c == c_copy2 ));
	}
	{
		namespace cuda = multi::cuda;
		cuda::array<complex, 2> const acu = a;
		cuda::array<complex, 2> const bcu = b;
		cuda::array<complex, 2> ccu({2, 4});
		using multi::blas::gemm;
		using multi::blas::hermitized;
		gemm(1., hermitized(acu), bcu, 0., ccu);
		BOOST_REQUIRE( ccu[1][2] == complex(112, 12) );

		cuda::array<complex, 2> const ccu_copy = gemm(1., hermitized(acu), bcu);
		cuda::array<complex, 2> const ccu_copy2 = gemm(hermitized(acu), bcu);
		BOOST_REQUIRE(( ccu_copy == ccu and ccu_copy2 == ccu ));
	}
	{
		namespace cuda = multi::cuda;
		cuda::managed::array<complex, 2> const amcu = a;
		cuda::managed::array<complex, 2> const bmcu = b;
		cuda::managed::array<complex, 2> cmcu({2, 4});
		using multi::blas::gemm;
		using multi::blas::hermitized;
		gemm(1., hermitized(amcu), bmcu, 0., cmcu);
		BOOST_REQUIRE( cmcu[1][2] == complex(112, 12) );

	//	[](void*){}();

		cuda::managed::array<complex, 2> const cmcu_copy = gemm(1., hermitized(amcu), bmcu);
		cuda::managed::array<complex, 2> const cmcu_copy2 = gemm(hermitized(amcu), bmcu);
		BOOST_REQUIRE(( cmcu_copy == cmcu and cmcu_copy2 == cmcu ));
	}
}

BOOST_AUTO_TEST_CASE(multi_adaptors_blas_gemm_complex_nonsquare_automatic3){
	using complex = std::complex<double>; complex const I{0,1};
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
		using multi::blas::gemm;
		using multi::blas::hermitized;
		gemm(1., hermitized(a), hermitized(b), 0., c); // c=ab, c⸆=b⸆a⸆
		BOOST_REQUIRE( c[1][2] == complex(112, 12) );
	}
	{
		namespace cuda = multi::cuda;
		cuda::array<complex, 2> const acu = a;
		cuda::array<complex, 2> const bcu = b;
		cuda::array<complex, 2> ccu({2, 4});
		using multi::blas::gemm;
		using multi::blas::hermitized;
		gemm(1., hermitized(acu), hermitized(bcu), 0., ccu);
		BOOST_REQUIRE( ccu[1][2] == complex(112, 12) );
	}
	{
		namespace cuda = multi::cuda;
		cuda::managed::array<complex, 2> const amcu = a;
		cuda::managed::array<complex, 2> const bmcu = b;
		cuda::managed::array<complex, 2> cmcu({2, 4});
		using multi::blas::gemm;
		using multi::blas::hermitized;
		gemm(1., hermitized(amcu), hermitized(bmcu), 0., cmcu);
		BOOST_REQUIRE( cmcu[1][2] == complex(112, 12) );
	}
}

BOOST_AUTO_TEST_CASE(multi_adaptors_blas_gemm_complex_nonsquare_automatic4){
	using complex = std::complex<double>; complex const I{0,1};
	multi::array<complex, 2> c({12, 12});
	{
		multi::array<complex, 2> const a({12, 100}, 1.+2.*I);
		multi::array<complex, 2> const b({12, 100}, 1.+2.*I);
		using multi::blas::hermitized;
		using multi::blas::gemm;
		gemm(1., a, hermitized(b), 0., c);
		BOOST_REQUIRE( real(c[0][0]) > 0);

		auto c_copy = gemm(1., a, hermitized(b));
		BOOST_REQUIRE( c_copy == c );
	}
	{
		multi::array<complex, 2> const a_block({24, 100}, 1.+2.*I);
		multi::array<complex, 2> const b({12, 100}, 1.+2.*I);
		multi::array<complex, 2> c2({12, 12});

		using multi::blas::hermitized;
		using multi::blas::gemm;
		gemm(1., a_block.strided(2), hermitized(b), 0., c2);

		BOOST_REQUIRE( real(c[0][0]) > 0);
		BOOST_REQUIRE( c == c2 );

		auto c2_copy = gemm(1., a_block.strided(2), hermitized(b));
		BOOST_REQUIRE( c2_copy == c2 );
	}
}

template<class... T> void what(T&&...) = delete;

BOOST_AUTO_TEST_CASE(multi_blas_gemm_complex_issue68){
	using complex = std::complex<double>; complex const I{0,1};
	multi::cuda::managed::array<complex, 2> const a = {
		{1.-2.*I},
		{3.+3.*I},
		{1.+9.*I}
	};
	{
		multi::cuda::managed::array<complex, 2> c({1, 1});
		using multi::blas::gemm;
		using multi::blas::hermitized;
		gemm(1., hermitized(a), a, 0., c);
		BOOST_REQUIRE( c[0][0] == 105. + 0.*I );
	}
	{
		using multi::blas::gemm;
		using multi::blas::hermitized;
		auto c = gemm(2., hermitized(a), a);
		BOOST_REQUIRE( c[0][0] == 210. + 0.*I );
	}
	{
		using multi::blas::gemm;
		using multi::blas::hermitized;
		auto c = gemm(hermitized(a), a);
		
	//	what<decltype(hermitized(a)),  decltype(base(hermitized(a))), decltype(underlying(base(hermitized(a)))), multi::pointer_traits<decltype(base(a))>::default_allocator_type, decltype(hermitized(a))::decay_type>();

//[with T = {boost::multi::basic_array<std::complex<double>, 2, boost::multi::blas::involuter<boost::multi::memory::cuda::managed::ptr<std::complex<double>, std::complex<double>*>, boost::multi::blas::conjugate>, boost::multi::layout_t<2, long int> >,
//boost::multi::blas::involuter<boost::multi::memory::cuda::managed::ptr<std::complex<double>, std::complex<double>*>, boost::multi::blas::conjugate>, 
//boost::multi::memory::cuda::managed::ptr<std::complex<double>, std::complex<double>*>, 
//boost::multi::memory::cuda::managed::allocator<std::complex<double> >, 
//boost::multi::array<std::complex<double>, 2, std::allocator<std::complex<double> > >}]’

	//	boost::multi::basic_array<std::complex<double>, 2, boost::multi::blas::involuter<boost::multi::memory::cuda::managed::ptr<std::complex<double>, std::complex<double>*>, boost::multi::blas::conjugate>, boost::multi::layout_t<2, long int> >
	//	boost::multi::array<std::complex<double>, 2, std::allocator<std::complex<double> > >
		BOOST_REQUIRE( c[0][0] == 105. + 0.*I );
	}
}

