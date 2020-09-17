#ifdef COMPILATION_INSTRUCTIONS
nvcc -std=c++14 -x cu --expt-relaxed-constexpr $0 -o $0x -Wno-deprecated-declarations -lboost_unit_test_framework&&$0x&&rm $0x; exit
#endif

#define BOOST_TEST_MODULE "C++ Unit Tests for Multi CUDA allocators"
#define BOOST_TEST_DYN_LINK
#include<boost/test/unit_test.hpp>

#include "../../../../../multi/array.hpp"

template<class ItA, class ItB, class ItC, class Size, class ItD>
__global__ void transform3_plus_CU(ItA first_a, ItB first_b, ItC first_c, Size max_n, ItD first_d){
	auto const ii = blockIdx.x*blockDim.x + threadIdx.x; (void)ii;
	if(ii < max_n) [[likely]] first_d[ii] = first_a[ii] + first_b[ii] + first_c[ii];
}

template<class ItA, class ItB, class ItC, typename Size, class ItD>
ItD transform3_n_plus_cuda(ItA first_a, ItB first_b, ItC first_c, Size n, ItD first_d){
	static constexpr Size blockSize = 256;
	transform3_plus_CU<<<(n - 1 + blockSize)/blockSize, blockSize>>>(first_a, first_b, first_c, n, first_d);
	cudaDeviceSynchronize();
	return first_d + n;
}

template<class ArrA, class ArrB, class ArrC, class ArrD>
ArrD&& plus3_cuda(ArrA const& a, ArrB const& b, ArrC const& c, ArrD&& d){
	assert( size(a) == size(b) and size(b) == size(c) and size(c)==size(d) );
	auto const last = transform3_n_plus_cuda(begin(a), begin(b), begin(c), size(a), begin(d));
	assert( end(d) == last );
	return std::forward<ArrD>(d);
}

template<class A, class B, class C, 
	class Ret = boost::multi::array<
		decltype(typename A::element_type{}+typename B::element_type{}+typename C::element_type{}), 
		A::dimensionality, typename A::allocator_type
	>
>
[[nodiscard]] 
Ret plus3_cuda(A const& a, B const& b, C const& c){
	Ret ret(size(a), 0, get_allocator(a));
	plus3_cuda(a, b, c, ret);
	return ret;
}

#include "../../cuda/allocator.hpp"
#include "../../cuda/managed/allocator.hpp"

#include<thrust/complex.h>

namespace multi = boost::multi;
namespace cuda = multi::memory::cuda;

#if 0
BOOST_AUTO_TEST_CASE(cuda_transform_plus_complex){
	using complex = thrust::complex<double>;
	using T = complex;
	using allocator = cuda::allocator<T>;
	using array1d = multi::array<T, 1, allocator>;
	array1d A1(200, T{0.});
	array1d B1(200, T{0.});
	array1d C1(200, T{0.});

	A1[100] = T{1., 1.};
	B1[100] = T{2., 1.};
	C1[100] = T{3., 1.};

	{
		multi::array<T, 1, allocator> D1(200, T{0.});

		transform3_plus_CU<<<size(A1), 1>>>(cbegin(A1), cbegin(B1), cbegin(C1), size(A1), begin(D1));

		BOOST_REQUIRE(( T{6., 3.} == D1[100] )); //complex{6., 6.} )); 
	}
}
#endif

template<class T> T what() = delete;
template<class T> T what(T&&) = delete;

BOOST_AUTO_TEST_CASE(cuda_transform_plus_double){

	using allocator = cuda::allocator<double>;
	multi::array<double, 1, allocator> A1(200, 0.);
	multi::array<double, 1, allocator> B1(200, 0.);
	multi::array<double, 1, allocator> C1(200, 0.);
	BOOST_REQUIRE( size(A1) == 200 );

	multi::array<double, 1, allocator> const Z1(200, 0.);
	auto p = Z1.data();
	p = nullptr;
//	decltype(Z1)::element_ptr p = nullptr;
//	what( Z1.cdata() );
//	assert( Z1[10] == 0. );

	A1[100] = 1.;
	B1[100] = 2.;
	C1[100] = 3.;

	{
		multi::array<double, 1, allocator> D1(200, 0.);

		transform3_plus_CU<<<size(A1), 1>>>(begin(A1), begin(B1), begin(C1), size(A1), begin(D1));

		BOOST_REQUIRE( D1[100] == 6. ); 
	}
	{
		multi::array<double, 1, cuda::allocator<double>> D1(200, 0.);

		auto e = transform3_n_plus_cuda(begin(A1), begin(B1), begin(C1), size(A1), begin(D1));	

		BOOST_REQUIRE( D1[100] == 6. );
		BOOST_REQUIRE( e == end(D1) );
	}
	{
		multi::array<double, 1, cuda::allocator<double>> D1(200, 0.);

		plus3_cuda(A1, B1, C1, D1);	

		BOOST_REQUIRE( D1[100] == 6. );
	}
	{
		auto DD1 = plus3_cuda(A1, B1, C1);	

		BOOST_REQUIRE( DD1[100] == 6. );
	}
	{
		multi::array<double, 2, cuda::allocator<double>> A2({100, 100}, 0.);
		multi::array<double, 2, cuda::allocator<double>> B2({100, 100}, 0.);
		multi::array<double, 2, cuda::allocator<double>> C2({100, 100}, 0.);

		A2[50][50] = 11.;
		B2[50][50] = 22.;
		C2[50][50] = 33.;

		multi::array<double, 2, cuda::allocator<double>> D2({100, 100}, 0.);

		auto e = transform3_n_plus_cuda(data_elements(A2), data_elements(B2), data_elements(C2), num_elements(A2), data_elements(D2));

		BOOST_REQUIRE( D2[50][50] == 66. );
	}
	{
		multi::array<double, 1, cuda::allocator<double>> D1(200, 0.);

		auto e = transform3_n_plus_cuda(begin(A1), begin(B1), begin(C1), 0, begin(D1));	

		BOOST_REQUIRE( e == begin(D1) );
	}
}

