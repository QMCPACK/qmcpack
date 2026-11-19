// Copyright 2018-2025 Alfredo A. Correa

#include <mpi3/communicator.hpp>
#include <mpi3/environment.hpp>

#include <boost/core/lightweight_test.hpp>
#include <cstddef>
#include <iterator>
#include <tuple>
#include <utility>
#include <vector>

namespace mpi3 = boost::mpi3;

namespace {
void part1(mpi3::communicator& world) {
	using T = std::tuple<double, double>;
	std::vector<T> v_local(10, T{world.rank(), world.rank()});
	std::vector<T> v(v_local.size() * static_cast<std::size_t>(world.size()));
	auto           end = world.all_gather_n(v_local.begin(), v_local.size(), v.begin());
	BOOST_TEST(end == v.end());
	BOOST_TEST((v[ 0] == T{0.0, 0.0}));
	BOOST_TEST((v[10] == T{1.0, 1.0}));
	BOOST_TEST((v[20] == T{2.0, 2.0}));
}

void part2(mpi3::communicator& world) {
	using T = std::tuple<double, double>;
	std::vector<T> v_local(10, T{world.rank(), world.rank()});
	std::vector<T> v(v_local.size() * static_cast<std::size_t>(world.size()));
	auto           d_last = world.all_gather(begin(v_local), end(v_local), begin(v));
	BOOST_TEST(d_last == end(v));
	BOOST_TEST((v[0] == T{0.0, 0.0}));
	BOOST_TEST((v[10] == T{1.0, 1.0}));
	BOOST_TEST((v[20] == T{2.0, 2.0}));
}

void part3(mpi3::communicator& world) {
	using T = std::pair<double, int>;
	std::vector<T> v_local(10, T{world.rank(), world.rank()});
	std::vector<T> v(v_local.size() * static_cast<std::size_t>(world.size()));
	auto           end = world.all_gather_n(v_local.begin(), v_local.size(), v.begin());
	BOOST_TEST(end == v.end());
	BOOST_TEST((v[0] == T{0.0, 0}));
	BOOST_TEST((v[10] == T{1.0, 1}));
	BOOST_TEST((v[20] == T{2.0, 2}));
}

void part4(mpi3::communicator& world) {
	using T = std::pair<double, int>;
	std::vector<T> v_local(10, T{world.rank(), world.rank()});
	std::vector<T> v(v_local.size() * static_cast<std::size_t>(world.size()));
	auto           d_last = world.all_gather(begin(v_local), end(v_local), begin(v));
	BOOST_TEST(d_last == end(v));
	BOOST_TEST((v[0] == T{0.0, 0}));
	BOOST_TEST((v[10] == T{1.0, 1}));
	BOOST_TEST((v[20] == T{2.0, 2}));
}

void part5(mpi3::communicator& world) {
	using T = std::pair<double, int>;
	std::vector<T> v_local(static_cast<std::size_t>(world.rank()) + 10, T{world.rank(), world.rank()});
	std::vector<T> v(v_local.size() * static_cast<std::size_t>(world.size()));
	auto           d_last = world.all_gather(begin(v_local), begin(v_local) + 4, begin(v));
	BOOST_TEST(std::distance(begin(v), d_last) == 4L * world.size());
	// BOOST_TEST(e == end(v));
	BOOST_TEST((v[0] == T{0.0, 0}));
	BOOST_TEST((v[4] == T{1.0, 1}));
	// BOOST_TEST(( v[10] == T{1.0, 1} ));
	// BOOST_TEST(( v[20] == T{2.0, 2} ));
}

void part6(mpi3::communicator& world) {
	auto cs = world.all_gather_as<std::vector<int>>(world.rank());
	BOOST_TEST(cs[0] == 0);
	BOOST_TEST(cs[1] == 1);
	BOOST_TEST(cs[2] == 2);
}

void part7(mpi3::communicator& world) {
	using T = double;
	std::vector<T> v_local(static_cast<std::size_t>(world.rank() + 1), world.rank());
	std::vector<T> v(1 + 2 + 3);
	auto           end = world.all_gatherv_n(v_local.begin(), v_local.size(), v.begin());
	BOOST_TEST(end == v.end());
	BOOST_TEST((v[0] == 0.0));
	BOOST_TEST((v[1] == 1.0));
	BOOST_TEST((v[2] == 1.0));
	BOOST_TEST((v[3] == 2.0));
	BOOST_TEST((v[4] == 2.0));
	BOOST_TEST((v[5] == 2.0));
}
}  // end namespace

auto main(int argc, char** argv) -> int try {
	mpi3::environment env(argc, argv);

	auto world = env.world();

	part1(world);
	part2(world);
	part3(world);
	part4(world);
	part5(world);
	part6(world);
	part7(world);

	return boost::report_errors();
} catch(...) {
	return 1;
}
