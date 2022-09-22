#define BOOST_TEST_MODULE "C++ Unit Tests for Multi CUDA thrust"
#define BOOST_TEST_DYN_LINK
#include<boost/test/unit_test.hpp>

#include "../../cuda/managed.hpp"

namespace multi = boost::multi;

void set_one(double* p){
	*p = 1.;
}

void set_two_gpu(thrust::cuda::pointer<double> p){
	*p = 2.;
}

void set_three_ref(double& p){
	p = 3.;
}

template<class Pointer, class V = typename std::iterator_traits<Pointer>::value_type, class = std::enable_if_t<std::is_same<V, double>{} and std::is_convertible<Pointer, thrust::cuda::pointer<V>>{}> >
void some_fun(Pointer p){}

template<class Pointer, class V = typename std::iterator_traits<Pointer>::value_type, class = std::enable_if_t<std::is_same<V, double>{} and std::is_convertible<Pointer, V*>{}> >
void some_other_fun(Pointer p){}

template<int N> class prio : std::conditional_t<N!=0, prio<N-1>, std::false_type>{};

template<class Pointer, class V = typename std::iterator_traits<Pointer>::value_type, std::enable_if_t<std::is_same<V, double>{} and std::is_convertible<Pointer, thrust::cuda::pointer<V>>{}, int> =0>
int overload_aux(Pointer p, prio<0>){return 0;}

template<class Pointer, class V = typename std::iterator_traits<Pointer>::value_type, std::enable_if_t<std::is_same<V, double>{} and std::is_convertible<Pointer, V*>{}, int> =0>
int overload_aux(Pointer p, prio<1>){return 1;}

template<class Pointer> int overload(Pointer p){return overload_aux(p, prio<1>{});}

BOOST_AUTO_TEST_CASE(vector){

	multi::thrust::cuda::managed::allocator<double> alloc;
	multi::thrust::cuda::managed::pointer<double> p = alloc.allocate(100);

	p[17] = 3.;
	BOOST_TEST_REQUIRE( p[17] == 3. );

	set_one(p);
	BOOST_TEST_REQUIRE( p[0] == 1. );

	set_two_gpu(p);
	BOOST_TEST_REQUIRE( p[0] == 2. );

	set_three_ref( p[1] );
	BOOST_TEST_REQUIRE( p[1] == 3. );

	some_fun(p);

	BOOST_TEST_REQUIRE(overload(p) == 1);

	alloc.deallocate(p, 100);

}
