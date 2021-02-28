#ifdef COMPILATION// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4-*-
$CXX $0 -o $0x -lcudart -lboost_unit_test_framework -lboost_timer -ldl&&$0x&&rm $0x;exit
#endif
// © Alfredo A. Correa 2019-2020

#ifndef MULTI_ADAPTORS_CUDA_HPP
#define MULTI_ADAPTORS_CUDA_HPP

#include "../memory/adaptors/cuda/allocator.hpp"
#include "../memory/adaptors/cuda/managed/allocator.hpp"
#include "../adaptors/cuda/algorithms/copy.hpp"

#include "../array.hpp"

namespace boost{
namespace multi{
namespace cuda{

	template<class T>
	using allocator = multi::memory::cuda::allocator<T>;

	template<class T> using ptr = multi::memory::cuda::ptr<T>;

	template<class T, multi::dimensionality_type D>
	using array = multi::array<T, D, cuda::allocator<T>>;

	template<class T, multi::dimensionality_type D>
	using array_ref = multi::array_ref<T, D, cuda::ptr<T>>;

	template<class T, multi::dimensionality_type D>
	using static_array = multi::static_array<T, D, cuda::allocator<T>>;

//	template<class A> auto raw_array_cast(A&& a)
//	->decltype(static_array_cast<typename A::element_type, decltype(raw_pointer_cast(base(std::forward<A>(a))))>(std::forward<A>(a))){
//		return static_array_cast<typename A::element_type, decltype(raw_pointer_cast(base(std::forward<A>(a))))>(std::forward<A>(a));}

	template<class A> auto raw_array_cast(A&& a)
	->decltype(std::forward<A>(a).template static_array_cast<typename A::element_type, decltype(raw_pointer_cast(base(std::forward<A>(a))))>()){
		return std::forward<A>(a).template static_array_cast<typename A::element_type, decltype(raw_pointer_cast(base(std::forward<A>(a))))>();}

	namespace managed{
		template<class T>
		using allocator = multi::memory::cuda::managed::allocator<T>;

		template<class T> using ptr = multi::memory::cuda::managed::ptr<T>;

		template<class T, multi::dimensionality_type D>
		using array = multi::array<T, D, cuda::managed::allocator<T>>;

		template<class T, multi::dimensionality_type D>
		using array_ref = multi::array<T, D, multi::memory::cuda::managed::ptr<T>>;

		template<class T, multi::dimensionality_type D>
		using static_array = multi::array<T, D, multi::memory::cuda::managed::ptr<T>>;
	}

}

/*
auto copy(const double* first, const double* last, boost::multi::array_iterator<double, 1, boost::multi::memory::cuda::managed::ptr<double, double*>, double&> d_first){
	return copy(
		boost::multi::array_iterator<double, 1, double const*, double const&>(first), 
		boost::multi::array_iterator<double, 1, double const*, double const&>(last), 
		d_first
	);
}*/

}}

#if defined(__INCLUDE_LEVEL__) and not __INCLUDE_LEVEL__
#define BOOST_TEST_MODULE "C++ Unit Tests for Multi CUDA adaptor"
#define BOOST_TEST_DYN_LINK
#include<boost/test/unit_test.hpp>
#include<boost/timer/timer.hpp>

#include<chrono>
#include<numeric>

template <class T>
__attribute__((always_inline)) inline void DoNotOptimize(const T &value) {
  asm volatile("" : "+m"(const_cast<T &>(value)));
}

struct watch : private std::chrono::high_resolution_clock{
	std::string label_; time_point  start_;
	watch(std::string label ="") : label_{label}, start_{now()}{}
	~watch(){
		std::cerr<< label_<<": "<< std::chrono::duration<double>(now() - start_).count() <<" sec"<<std::endl;
	}
};


namespace multi = boost::multi;
namespace cuda = multi::cuda;
namespace utf = boost::unit_test::framework;

BOOST_AUTO_TEST_CASE(multi_adaptors_cuda_construct_1d){
	multi::array<double, 1> A(4, 99.);
	cuda::array<double, 1> Agpu{A};
	BOOST_REQUIRE( extensions(A) == extensions(Agpu) );
	BOOST_REQUIRE( Agpu == A );
}

BOOST_AUTO_TEST_CASE(multi_adaptors_cuda_copy_1d){
	multi::array<double, 1> A(4, 99.);
	cuda::array<double, 1> Agpu(4);
	BOOST_REQUIRE( extensions(A) == extensions(Agpu) );
	Agpu({0, 4}) = A({0, 4});
	BOOST_REQUIRE( Agpu == A );
}

BOOST_AUTO_TEST_CASE(multi_adaptors_cuda_construct_2d){
	multi::array<double, 2> A({4, 6}, 99.);
	cuda::array<double, 2> Agpu{A};
	BOOST_REQUIRE( extensions(A) == extensions(Agpu) );
	BOOST_REQUIRE( Agpu == A );
	
	A[1][1] = Agpu[1][1];
}

BOOST_AUTO_TEST_CASE(multi_adaptors_cuda_copy_2d){
	multi::array<double, 2> A({4, 6}); std::iota(A.data_elements(), A.data_elements() + A.num_elements(), 1.);
	cuda::array<double, 2> Agpu({4, 6}, 99.);
	BOOST_REQUIRE( extensions(A) == extensions(Agpu) );
	Agpu({0, 4}, {1, 6}) = A({0, 4}, {1, 6});
	BOOST_REQUIRE( Agpu != A );
	Agpu = A;
	BOOST_REQUIRE( Agpu == A );
}

BOOST_AUTO_TEST_CASE(multi_adaptors_cuda_1d_initializer_list){
	cuda::array<double, 1> Bgpu = {1., 2., 3., 4.};
	BOOST_REQUIRE( Bgpu[1] == 2. );
}

BOOST_AUTO_TEST_CASE(multi_adaptors_cuda_2d_initializer_list){
	cuda::array<double, 2> Bgpu = {
		{1., 2., 3., 4.},
		{5., 6., 7., 8.},
		{9., 10., 11., 12.},
	};
	BOOST_REQUIRE( size(Bgpu) == 3 );
	BOOST_REQUIRE( Bgpu[1][1] == 6. );
}

BOOST_AUTO_TEST_CASE(multi_adaptors_cuda_2d_initializer_list_bis){
	multi::array<double, 2> A({3, 4}); std::iota(A.data_elements(), A.data_elements() + A.num_elements(), 1.);

	cuda::array<double, 2> Agpu({3, 4}, 99.);
	BOOST_REQUIRE( extensions(A) == extensions(Agpu) );
	Agpu({0, 3}, {1, 4}) = A({0, 3}, {1, 4});
	BOOST_REQUIRE( Agpu != A );
	Agpu = A;
	BOOST_REQUIRE( Agpu == A );
}

BOOST_AUTO_TEST_CASE(multi_adaptors_cuda_copy_vs_move){
	cuda::array<double, 4> Agpu({30, 100, 100, 100}, 99.);
	[&, _ = watch{utf::current_test_case().full_name()+" COPY"}]{
		cuda::array<double, 4> Agpu_cpy = Agpu;
		BOOST_REQUIRE( &Agpu_cpy[1][2][3][4] != &Agpu[1][2][3][4] );
		BOOST_REQUIRE( Agpu_cpy[1][2][3][4] == Agpu[1][2][3][4] );
	}();
	[&, _ = watch{utf::current_test_case().full_name()+" MOVE"}]{
		cuda::array<double, 4> Agpu_mov = std::move(Agpu);
		BOOST_REQUIRE( Agpu.empty() );
		BOOST_REQUIRE( Agpu_mov.size() == 30 );
	}();
}

BOOST_AUTO_TEST_CASE(multi_adaptors_cuda_copy_vs_move_complex){
	using complex = std::complex<double>;
	cuda::array<complex, 4> Agpu({30, 100, 100, 100}, 99.);
	[&, _ = watch{utf::current_test_case().full_name()+" COPY"}]{
		cuda::array<complex, 4> Agpu_cpy = Agpu;
		BOOST_REQUIRE( &Agpu_cpy[1][2][3][4] != &Agpu[1][2][3][4] );
		BOOST_REQUIRE( Agpu_cpy[1][2][3][4] == Agpu[1][2][3][4] );
	}();
	[&, _ = watch{utf::current_test_case().full_name()+" MOVE"}]{
		cuda::array<complex, 4> Agpu_mov = std::move(Agpu);
		BOOST_REQUIRE( Agpu.empty() );
		BOOST_REQUIRE( Agpu_mov.size() == 30 );
	}();
}

BOOST_AUTO_TEST_CASE(multi_adaptors_cuda_managed_double){
	cuda::managed::array<double, 4> A({2,3,4,5});
	cuda::managed::array<double, 4> B({2,3,4,5}, 0.);
	cuda::managed::array<double, 4> C({2,3,4,5}, 5.);
}

BOOST_AUTO_TEST_CASE(multi_adaptors_cuda_managed_ai3){
	using ai3 = std::array<int, 3>;
	cuda::managed::array<ai3, 4> A({2,3,4,5}); // default initialize elements
	cuda::managed::array<ai3, 4> B({2,3,4,5}, ai3{} ); // value initialize elements
	cuda::managed::array<ai3, 4> C({2,3,4,5}, ai3{11, 22, 33} ); // value initialize elements
}

BOOST_AUTO_TEST_CASE(multi_adaptor_cuda_decay){
	cuda::array<double, 2> A = {
		{1., 2., 3., 4.},
		{5., 6., 7., 8.},
		{1., 2., 3., 4.}
	};
	cuda::array<double, 1> A1 = A[1];

//	cuda::array<complex, 2> A = {
//		{1. + 2.*I, 2. + 3.*I, 3. + 4.*I, 4. + 5.*I},
//		{5. + 2.*I, 6. + 3.*I, 7. + 4.*I, 8. + 5.*I},
//		{1. + 1.*I, 2. + 2.*I, 3. + 3.*I, 4. + 4.*I}
//	};
//	cuda::array<complex
}

#if 0
BOOST_AUTO_TEST_CASE(multi_adaptors_cuda_copy){

	multi::array<double, 2> A({4, 5}, 99.);
	cuda::array<double, 2> Agpu = A;

}

BOOST_AUTO_TEST_CASE(multi_adaptors_cuda){

	multi::array<double, 2> A({4, 5}, 99.);
	cuda::array<double, 2> Agpu = A;
	assert( Agpu == A );

	cuda::managed::array<double, 2> Amng = A;
	assert( Amng == Agpu );

	cuda::array_ref<double, 2> Rgpu(data_elements(Agpu), extensions(Agpu));

	{std::allocator<double> a = get_allocator(A);}

	{
		cuda::ptr<double> p;
		using multi::get_allocator;
		cuda::allocator<double> a = get_allocator(p); (void)a;
	}
	{
		cuda::managed::ptr<double> p;
		using multi::get_allocator;
		cuda::managed::allocator<double> a = get_allocator(p); (void)a;
	}
	{
		double* p = nullptr;
		using multi::get_allocator;
		std::allocator<double> a = get_allocator(p); (void)a;
	}
	{
		multi::array<double, 2> arr;
		std::allocator<double> a = get_allocator(arr);
	}
	{
		cuda::array<double, 2> arr;
		cuda::allocator<double> a = get_allocator(arr); (void)a;
	}
	{
//		cuda::array<double, 0> arr = 45.;
//		BOOST_REQUIRE( arr() == 45. );
	}
	{
//		cuda::managed::array<double, 0> arr = 45.;
//		BOOST_REQUIRE( arr() == 45. );
	}
	{
		cuda::managed::array<double, 1> arr = {1.2, 3.4, 4.5};
	}
	{
		using complex = std::complex<double>;
		cuda::managed::array<complex, 2> a({1000, 1000}, 99.);
		BOOST_REQUIRE( size(a) == 1000 );
		cuda::managed::array<complex, 2> b;
		b = std::move(a);
		BOOST_REQUIRE( size(b) == 1000 );
		BOOST_REQUIRE( size(a) == 0 );
	}
}
#endif
#endif
#endif

