// Copyright 2018-2025 Alfredo A. Correa

#include <mpi3/communicator.hpp>
#include <mpi3/environment.hpp>
#include <mpi3/operation.hpp>

#include <boost/core/lightweight_test.hpp>

#include <cstddef>
#include <functional>
#include <iterator>
#include <numeric>
#include <vector>

namespace mpi3 = boost::mpi3;

int main(int argc, char** argv) {  // NOLINT(bugprone-exception-escape)
	mpi3::environment env(argc, argv);

	auto world = env.world();

	BOOST_TEST(world.size() > 1);

	// {
	// 	std::vector<std::size_t> global(120);

	// 	std::vector<std::size_t> local(global.size());
	// 	iota(begin(local), end(local), 0);  // NOLINT(boost-use-ranges) for C++20, use std::ranges::iota

	// 	auto last = world.all_reduce_n(local.begin(), local.size(), global.begin());
	// 	BOOST_TEST(last == global.end());

	// 	BOOST_TEST(
	// 		std::inner_product(
	// 			global.begin(), global.end(), local.begin(),
	// 			true, std::logical_and<bool>{},
	// 			[sz = static_cast<std::size_t>(world.size())](auto const& e1, auto const& e2) { return e1 == e2 * sz; }
	// 		)
	// 	);

	// 	auto const sum_of_ranks = (world += world.rank());
	// 	BOOST_TEST(sum_of_ranks == world.size() * (world.size() - 1) / 2);

	// 	auto rank     = world.rank();
	// 	auto sum_rank = 0;
	// 	world.all_reduce_n(&rank, 1, &sum_rank);
	// 	//  sum_rank = (world += rank);
	// 	BOOST_TEST(sum_rank == world.size() * (world.size() - 1) / 2);

	// 	auto max_rank = -1;
	// 	world.all_reduce_n(&rank, 1, &max_rank, mpi3::max<>{});
	// 	BOOST_TEST(max_rank == world.size() - 1);

	// 	auto min_rank = -1;
	// 	world.all_reduce_n(&rank, 1, &min_rank, mpi3::min<>{});
	// 	BOOST_TEST(min_rank == 0);
	// }
	{
		std::vector<int> local(20, 1);
		if(world.rank() == 2) {
			local[1] = 0;
		}

		std::vector<int> global(local.size());
		world.all_reduce_n(local.begin(), local.size(), global.begin() /*, std::plus<>{}*/);

		BOOST_TEST(global[0] == world.size());
		BOOST_TEST(global[1] == world.size() - 1);
	}
	{
		std::vector<int> local(20, 1);
		if(world.rank() == 2) {
			local[1] = 9;
		}

		std::vector<int> global(local.size());
		world.all_reduce_n(local.begin(), local.size(), global.begin(), std::logical_and<>{});

		BOOST_TEST(global[0] != 0);
		BOOST_TEST(global[1] == 1);
	}
	{
		int b = 1;
		if(world.rank() == 2) {
			b = 0;
		}
		int const all = (world += b);
		BOOST_TEST(all == world.size() - 1);
	}
	{
		int const b   = world.rank() == 2 ? 0 : 1;
		int const all = (world &= b);
		BOOST_TEST(all == false);
	}
	{
		bool const b      = not(world.rank() == 1);
		bool const all_of = (world &= b);
		BOOST_TEST(all_of == false);
	}
	{
		bool const b      = (world.rank() == 1);
		bool       any_of = false;
		world.all_reduce_n(&b, 1, &any_of, std::logical_or<>{});
		BOOST_TEST(any_of);

		bool all_of = true;
		world.all_reduce_n(&b, 1, &all_of, std::logical_and<>{});
		BOOST_TEST(not all_of);
	}
	{
		auto const local_rnk = world.rank();
		auto const sum = world.operator+(local_rnk);
		BOOST_TEST( sum == world.size()*(world.size() - 1)/2 );
	}
	{
		bool const b      = (world.rank() == 1);
		bool const any_of = (world |= b);
		BOOST_TEST(any_of);

		bool const all_of = (world &= b);  // cppcheck-suppress unreadVariable
		BOOST_TEST(not all_of);
	}

	return boost::report_errors();
}