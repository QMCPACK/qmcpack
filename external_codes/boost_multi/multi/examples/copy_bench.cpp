#ifdef COMPILATION_INSTRUCTIONS //sudo cpupower frequency-set --governor performance
clang++ -O3 $0 -o$0x -DBOOST_TEST_DYN_LINK -lboost_unit_test_framework -lboost_timer -lbenchmark -lpthread&&$0x&&rm $0x;exit
#endif

//#define BOOST_TEST_MODULE "C++ Unit Tests for Multi move"
//#include<boost/test/unit_test.hpp>

#include "../array.hpp"

#include <boost/timer/timer.hpp>

#include<vector>
#include<complex>
#include<numeric>

#include <benchmark/benchmark.h>

namespace multi = boost::multi;

#if 0
BOOST_AUTO_TEST_CASE(multi_array_copy){
	std::vector<multi::array<double, 2> > A(10, multi::array<double, 2>({4, 5}, 99.));
	std::vector<multi::array<double, 2> > B;
	B.reserve(10);

	for(auto& v: A)	B.emplace_back(v);

	BOOST_REQUIRE( size(B) == size(A) );
	BOOST_REQUIRE( size(A[4]) == 4 );
	BOOST_REQUIRE( size(B[5]) == 4 );
	BOOST_REQUIRE( B[5][1][2] == 99. );
}
#endif

//BOOST_AUTO_TEST_CASE(multi_copy_vs_vector_copy){
namespace boost{
namespace multi{
	template<class T1, class P1, class T2, class P2>
	auto copy(array_iterator<T1, 1, P1> f, array_iterator<T1, 1, P1> l, array_iterator<T2, 1, P2> d)
	->decltype(array_iterator<T2, 1, P2>{(stride(f)==1 and stride(d)==1)?std::copy(base(f), base(l), base(d)):std::copy(f, l, d)}){assert(stride(f)==stride(l));
		return array_iterator<T2, 1, P2>{(stride(f)==1 and stride(d)==1)?std::copy(base(f), base(l), base(d)):std::copy(f, l, d)};}
}}

class A{
	double i;
public:
	A(double const& i) : i{i}{}
	A(A const& o) = default;
	A& operator=(A const& a) = default;
};

static void CopyMulti(benchmark::State& state) {
	multi::array<A, 1> a(1 << 16, {1.1}); std::iota(begin(a), end(a), 1.1111);
	multi::array<A, 1> b(1 << 16, {99.});
	for(auto _ : state){
		using std::copy;
		copy(begin(a), end(a), begin(b));
		benchmark::DoNotOptimize(b);
	}
}
BENCHMARK(CopyMulti);

static void CopyMulti2(benchmark::State& state) {
	multi::array<A, 1> a(1 << 16, {1.1}); std::iota(begin(a), end(a), 1.1111);
	multi::array<A, 1> b(1 << 16, {99.});
	for(auto _ : state){
		std::copy(begin(a), end(a), begin(b));
		benchmark::DoNotOptimize(b);
	}
}
BENCHMARK(CopyMulti2);


static void Copy(benchmark::State& state) {
	std::vector<A> v(1 << 16, {1.1}); std::iota(begin(v), end(v), 1.1111);
	std::vector<A> w(1 << 16, {99.});
	for(auto _ : state){
		std::copy(begin(v), end(v), begin(w));
		benchmark::DoNotOptimize(w);
	}
}
BENCHMARK(Copy);

BENCHMARK_MAIN();


