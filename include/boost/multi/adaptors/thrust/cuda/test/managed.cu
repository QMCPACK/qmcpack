
auto main() -> int {
}
// #include<boost/test/unit_test.hpp>

// #include "../../cuda/managed.hpp"

// namespace multi = boost::multi;

// void set_one(double* p){
// 	*p = 1.;
// }

// void set_two_gpu(thrust::cuda::pointer<double> p){
// 	*p = 2.;
// }

// void set_three_ref(double& p){
// 	p = 3.;
// }

// template<class Pointer, class V = typename std::iterator_traits<Pointer>::value_type, class = std::enable_if_t<std::is_same<V, double>{} and std::is_convertible<Pointer, thrust::cuda::pointer<V>>{}> >
// void some_fun(Pointer p){}

// template<class Pointer, class V = typename std::iterator_traits<Pointer>::value_type, class = std::enable_if_t<std::is_same<V, double>{} and std::is_convertible<Pointer, V*>{}> >
// void some_other_fun(Pointer p){}

// template<int N> class prio : std::conditional_t<N!=0, prio<N-1>, std::false_type>{};

// template<class Pointer, class V = typename std::iterator_traits<Pointer>::value_type, std::enable_if_t<std::is_same<V, double>{} and std::is_convertible<Pointer, thrust::cuda::pointer<V>>{}, int> =0>
// int overload_aux(Pointer p, prio<0>){return 0;}

// template<class Pointer, class V = typename std::iterator_traits<Pointer>::value_type, std::enable_if_t<std::is_same<V, double>{} and std::is_convertible<Pointer, V*>{}, int> =0>
// int overload_aux(Pointer p, prio<1>){return 1;}

// template<class Pointer> int overload(Pointer p){return overload_aux(p, prio<1>{});}

// BOOST_AUTO_TEST_CASE(vector){

// 	multi::thrust::cuda::managed::allocator<double> alloc;
// 	multi::thrust::cuda::managed::pointer<double> p = alloc.allocate(100);

// 	p[17] = 3.;
// 	BOOST_TEST_REQUIRE( p[17] == 3. );

// 	set_one(p);
// 	BOOST_TEST_REQUIRE( p[0] == 1. );

// 	set_two_gpu(p);
// 	BOOST_TEST_REQUIRE( p[0] == 2. );

// 	set_three_ref( p[1] );
// 	BOOST_TEST_REQUIRE( p[1] == 3. );

// 	some_fun(p);

// 	BOOST_TEST_REQUIRE(overload(p) == 1);

// 	alloc.deallocate(p, 100);

// }

// BOOST_AUTO_TEST_CASE(vector)
// {
// 	static_assert(std::is_same<typename std::pointer_traits<cuda::ptr<double>>::element_type, double>{}, "!");
// 	cuda::allocator<double> calloc;
// 	cuda::ptr<double> p = calloc.allocate(100);
// 	cuda::ptr<void> v = p;
// 	cuda::ptr<void const> vc{v};
// 	v = const_pointer_cast(vc);
// 	assert( vc == v );
// 	std::pointer_traits<decltype(p)>::rebind<double const> pc = p; // cuda::ptr<double const> pc = p;
// 	assert( pc == p );
// 	using cuda::const_pointer_cast;
// 	auto end = p + 100;
// 	auto rbegin = std::make_reverse_iterator(end);
// 	auto rend = std::make_reverse_iterator(p);
// 	std::transform(rbegin, rend, rbegin, [](auto&& e){return std::forward<decltype(e)>(e) + 99.;});
// 	assert( p[11] == 99. );
// 	p[33] = 123.;
// 	p[99] = 321.;
// //  p[33] += 1;
// 	add_one(p[33]);
// 	double p33 = p[33];
// 	assert( p33 == 124. );
// 	assert( p[33] == 124. );
// 	assert( p[33] == p[33] );
// 	swap(p[33], p[99]);
// 	assert( p[99] == 124. );
// 	assert( p[33] == 321. );
// 	std::cout << p[33] << std::endl;
// 	calloc.deallocate(p, 100);

// 	multi::array<double, 1, cuda::allocator<double>> arr2(multi::array<double, 1>::extensions_type{100l}, 999.);

// 	assert(size(arr2) == 100);
// }

// #ifdef COMPILATION_INSTRUCTIONS
// nvcc -ccbin cuda-c++ -std=c++14 $0 -o $0x && $0x && rm -f $0x; exit
// #endif

// #include "../../../../multi/array.hpp"
// #include "../../../../multi/detail/stack_allocator.hpp"
// #include "../../../../multi/detail/cuda/allocator.hpp"

// #include<iostream>

// namespace multi = boost::multi;
// namespace cuda = multi::detail::memory::cuda;

// using std::cout;

// int main(){
// 	{
// 		std::size_t stack_size = 4000;
// 		multi::stack_buffer<cuda::allocator<>> buf{stack_size};
// 		for(int i = 0; i != 3; ++i){
// 			cout<<"pass "<< i << std::endl;
// 			{
// 				multi::array<double, 2, multi::stack_allocator<double, cuda::allocator<>>> A({2, 10}, &buf);
// 				multi::array<int,    2, multi::stack_allocator<double, cuda::allocator<>>> B({3, 10}, &buf);
// 				multi::array<double, 2, multi::stack_allocator<double, cuda::allocator<>>> C({4, 10}, &buf);
// 				for(int j = 0; j != 100; ++j)
// 					multi::array<double, 2, multi::stack_allocator<double, cuda::allocator<>>> D({4, 10}, &buf);
// 				B[1][1] = 33.;
// 				B[2][2] = 33.;
// 				assert( B[1][1] == B[2][2] );
// 			}
// 			cout
// 				<<"  size: "<< buf.size() 
// 				<<"\n  hits: "<< buf.hits() 
// 				<<"\n  misses "<< buf.misses() 
// 				<<"\n  allocated(bytes) "<< buf.allocated_bytes() 
// 				<<"\n  deallocated(bytes) "<< buf.deallocated_bytes()
// 				<<"\n  max_needed(bytes) "<< buf.max_needed()
// 				<<"\n  stack recovered(bytes) " << buf.stack_recovered()
// 				<< std::endl
// 			;
// 			assert( buf.allocated_bytes() == buf.deallocated_bytes() );
// 			if(buf.max_needed() > buf.size()) buf.reset(buf.max_needed());
// 		}
// 	}
// 	assert( cuda::allocation_counter::n_allocations == 1 );
// }
