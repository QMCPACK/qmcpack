#define ENABLE_GPU 1

#include <boost/core/lightweight_test.hpp>

#include <boost/multi/array.hpp>

#include <thrust/universal_allocator.h>

#include <boost/multi/adaptors/thrust.hpp>
#include <boost/multi/array.hpp>

namespace multi = boost::multi;

__global__ void add(int a, int b, int* c) {
	*c = a + b;
}

auto main() -> int {
	multi::thrust::cuda::array<int, 0> dev_c(multi::extents_t<0>{}, 0);

	add<<<1, 1>>(2, 7, dev_c.base());

	multi::array<int, 0> c(dev_c);

	BOOST_TEST( *c.base() == 9 );

	return boost::report_errors();
}
