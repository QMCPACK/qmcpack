// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;autowrap:nil;-*-
// Copyright 2023 Alfredo A. Correa

#include <mpi3/communicator.hpp>
#include <mpi3/environment.hpp>

// #include <mpi3/main.hpp>

namespace mpi3 = boost::mpi3;

struct vec3 {
	double x;  // NOLINT(misc-non-private-member-variables-in-classes)
	double y;  // NOLINT(misc-non-private-member-variables-in-classes)
	double z;  // NOLINT(misc-non-private-member-variables-in-classes)

	constexpr bool operator==(vec3 const& other) const {return x == other.x and y == other.y and z == other.z;}
	constexpr bool operator!=(vec3 const& other) const {return x != other.x or  y != other.y or  z != other.z;}

	constexpr auto operator+(vec3 const& other) const {return vec3{x + other.x, y + other.y, z + other.z};}
};

template<> struct mpi3::datatype<vec3> : mpi3::struct_<
	double,
	double,
	double
> {};

mpi3::environment mpienv;  // NOLINT(fuchsia-statically-constructed-objects,cert-err58-cpp,cppcoreguidelines-avoid-non-const-global-variables) experiment with global environment

auto main(int /*argc*/, char** /*argv*/) -> int try { // NOLINT(bugprone-exception-escape)

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
		assert(( v[2] == vec3{1.0, 2.0, 3.0} ));
		break;
	};
	}

	{
		std::vector<vec3> w = { vec3{1.0, 2.0, 3.0}, vec3{4.0, 5.0, 6.0} };
		std::vector<vec3> sum(2);

		world.all_reduce_n(w.begin(), w.size(), sum.begin());

		assert(sum[0].x == w[0].x*world.size() );
		assert(sum[1].y == w[1].y*world.size() );
	}
	{
		std::vector<vec3> w = { vec3{1.0, 2.0, 3.0}, vec3{4.0, 5.0, 6.0} };
		std::vector<vec3> sum(2);

		world.all_reduce(w.begin(), w.end(), sum.begin());

		assert(sum[0].x == w[0].x*world.size() );
		assert(sum[1].y == w[1].y*world.size() );
	}
	{
		std::vector<vec3> w = { vec3{1.0, 2.0, 3.0}, vec3{4.0, 5.0, 6.0} };
		std::vector<vec3> sum(2);

		world.all_reduce_n(w.begin(), w.size(), sum.begin());

		assert(sum[0].x == w[0].x*world.size() );
		assert(sum[1].y == w[1].y*world.size() );
	}
	{
		std::vector<vec3> w = { vec3{1.0, 2.0, 3.0}, vec3{4.0, 5.0, 6.0} };
		std::vector<vec3> sum(2);

		world.all_reduce(w.begin(), w.end(), sum.begin());

		assert(sum[0].x == w[0].x*world.size() );
		assert(sum[1].y == w[1].y*world.size() );
	}
	{
		std::vector<vec3> w = { vec3{1.0, 2.0, 3.0}, vec3{4.0, 5.0, 6.0} };
		auto w_copy = w;

		world.all_reduce_n(w.begin(), w.size());

		assert(w[0].x == w_copy[0].x*world.size() );
		assert(w[1].y == w_copy[1].y*world.size() );
	}
	{
		std::vector<vec3> w = { vec3{1.0, 2.0, 3.0}, vec3{4.0, 5.0, 6.0} };
		auto w_copy = w;

		world.all_reduce(w.begin(), w.end());

		assert(w[0].x == w_copy[0].x*world.size() );
		assert(w[1].y == w_copy[1].y*world.size() );
	}
	static_assert(boost::mpi3::has_datatype<vec3>{});
} catch(...) {
	throw;
	return 1;
}
