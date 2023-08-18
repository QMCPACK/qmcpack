// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;autowrap:nil;-*-
// Copyright 2023 Alfredo A. Correa

// #include <mpi3/communicator.hpp>

#include <mpi3/main_environment.hpp>
#include <mpi3/operation.hpp>

namespace mpi3 = boost::mpi3;

struct singleton {
	int value;  // NOLINT(misc-non-private-member-variables-in-classes)
	bool operator==(singleton const& other) const {return value == other.value;}
	bool operator!=(singleton const& other) const {return value != other.value;}
	constexpr auto operator+(singleton const& other) const {return singleton{value + other.value};}
};

template<> struct mpi3::datatype<singleton> : mpi3::struct_<int> {};

int mpi3::main(int /*argc*/, char** /*argv*/, mpi3::environment& mpie) {  // NOLINT(bugprone-exception-escape)
	auto world{mpie.world()};

	{
		auto data = singleton{world.rank()};

		singleton result{};
		world.reduce_n(&data, 1, &result);
		world.broadcast_n(&result, 1);

		singleton correct_result{0};
		for(int i = 0; i != world.size(); ++i) {correct_result = correct_result + singleton{i};}  // NOLINT(altera-unroll-loops,altera-id-dependent-backward-branch)
		assert(result == correct_result);
	}
	{
		auto data = std::vector<singleton>{{ {world.rank()}, {world.rank()}, {world.rank()} }};

		std::vector<singleton> result(3);
		world.reduce_n(data.begin(), 3, result.begin());
		world.broadcast_n(result.begin(), 3);

		singleton correct_result{0};
		for(int i = 0; i != world.size(); ++i) {correct_result = correct_result + singleton{i};}  // NOLINT(altera-unroll-loops,altera-id-dependent-backward-branch)

		assert(result[1].value == correct_result.value);
	}
	{
		auto data = std::vector<singleton>{{ {world.rank()}, {world.rank()}, {world.rank()} }};

		std::vector<singleton> result(3);
		world.reduce(data.begin(), data.end(), result.begin());
		world.broadcast(result.begin(), result.end());

		singleton correct_result{0};
		for(int i = 0; i != world.size(); ++i) {correct_result = correct_result + singleton{i};}  // NOLINT(altera-unroll-loops,altera-id-dependent-backward-branch)

		assert(result[1].value == correct_result.value);
	}

	return 0;
}