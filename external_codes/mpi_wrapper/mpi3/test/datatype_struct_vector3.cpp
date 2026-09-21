// Copyright 2023-2024 Alfredo A. Correa

#include <mpi3/communicator.hpp>
#include <mpi3/detail/datatype.hpp>
#include <mpi3/environment.hpp>
#include <mpi3/type.hpp>

#include <boost/core/lightweight_test.hpp>
#include <vector>

namespace mpi3 = boost::mpi3;

struct vec3 {
	double x;  // NOLINT(misc-non-private-member-variables-in-classes)
	double y;  // NOLINT(misc-non-private-member-variables-in-classes)
	double z;  // NOLINT(misc-non-private-member-variables-in-classes)

	constexpr bool operator==(vec3 const& other) const { return x == other.x and y == other.y and z == other.z; }
	constexpr bool operator!=(vec3 const& other) const { return x != other.x or y != other.y or z != other.z; }

	constexpr auto operator+(vec3 const& other) const { return vec3{x + other.x, y + other.y, z + other.z}; }
};

template<> struct mpi3::datatype<vec3> : mpi3::struct_<double, double, double> {};

namespace {
mpi3::environment mpienv;  // NOLINT(fuchsia-statically-constructed-objects,cert-err58-cpp,cppcoreguidelines-avoid-non-const-global-variables) experiment with global environment
}  // end namespace

auto main(int /*argc*/, char** /*argv*/) -> int try {

	mpi3::communicator world{mpienv.world()};

	switch(world.rank()) {
	case 0: {
		std::vector<vec3> v(5);
		v[2] = vec3{1.0, 2.0, 3.0};
		world.send_n(begin(v), 5, 1);
		break;
	};
	case 1: {
		std::vector<vec3> v(5);
		world.receive_n(begin(v), 5, 0);
		BOOST_TEST((v[2] == vec3{1.0, 2.0, 3.0}));
		break;
	};
	default:;
	}

	{
		std::vector<vec3> w = {
			vec3{1.0, 2.0, 3.0},
			vec3{4.0, 5.0, 6.0}
		};
		std::vector<vec3> sum(2);

		world.all_reduce_n(w.begin(), w.size(), sum.begin());

		BOOST_TEST(sum[0].x == w[0].x * world.size());
		BOOST_TEST(sum[1].y == w[1].y * world.size());
	}
	{
		std::vector<vec3> w = {
			vec3{1.0, 2.0, 3.0},
			vec3{4.0, 5.0, 6.0}
		};
		std::vector<vec3> sum(2);

		world.all_reduce(w.begin(), w.end(), sum.begin());

		BOOST_TEST(sum[0].x == w[0].x * world.size());
		BOOST_TEST(sum[1].y == w[1].y * world.size());
	}
	{
		std::vector<vec3> w = {
			vec3{1.0, 2.0, 3.0},
			vec3{4.0, 5.0, 6.0}
		};
		std::vector<vec3> sum(2);

		world.all_reduce_n(w.begin(), w.size(), sum.begin());

		BOOST_TEST(sum[0].x == w[0].x * world.size());
		BOOST_TEST(sum[1].y == w[1].y * world.size());
	}
	{
		std::vector<vec3> w = {
			vec3{1.0, 2.0, 3.0},
			vec3{4.0, 5.0, 6.0}
		};
		std::vector<vec3> sum(2);

		world.all_reduce(w.begin(), w.end(), sum.begin());

		BOOST_TEST(sum[0].x == w[0].x * world.size());
		BOOST_TEST(sum[1].y == w[1].y * world.size());
	}
	// {
	// 	std::vector<vec3> w = {
	// 		vec3{1.0, 2.0, 3.0},
	// 		vec3{4.0, 5.0, 6.0}
	// 	};
	// 	auto w_copy = w;

	// 	world.all_reduce_n(w.begin(), w.size());

	// 	BOOST_TEST(w[0].x == w_copy[0].x * world.size());
	// 	BOOST_TEST(w[1].y == w_copy[1].y * world.size());
	// }
	{
		std::vector<vec3> w = {
			vec3{1.0, 2.0, 3.0},
			vec3{4.0, 5.0, 6.0}
		};
		auto w_copy = w;

		world.all_reduce(w.begin(), w.end());

		BOOST_TEST(w[0].x == w_copy[0].x * world.size());
		BOOST_TEST(w[1].y == w_copy[1].y * world.size());
	}
	static_assert(boost::mpi3::has_datatype<vec3>{});
	return boost::report_errors();
} catch(...) {
	return 1;
}
