#include <boost/multi/array.hpp>
#include <cassert>

auto main() -> int {
    namespace multi = boost::multi;
    multi::array<double, 2> A({3, 3});
    A[1][1] = 42.0;
    assert(A[1][1] == 42.0);
}
