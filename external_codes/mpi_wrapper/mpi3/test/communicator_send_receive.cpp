// Copyright 2017-2025 Alfredo A. Correa

#include <mpi3/cartesian_communicator.hpp>
#include <mpi3/communicator.hpp>
#include <mpi3/environment.hpp>

#include <algorithm>
#include <array>
#include <boost/core/lightweight_test.hpp>
#include <complex>
#include <exception>
#include <iostream>
#include <list>
#include <string>
#include <vector>

namespace mpi3 = boost::mpi3;


auto main(int argc, char** argv) -> int try {  // NOLINT(readability-function-cognitive-complexity)
	mpi3::environment env(argc, argv);

	auto world = env.world();

	{
		BOOST_TEST( not world.is_empty() );
		auto const s = world.size();
		BOOST_TEST( s != 0 );

		auto const right = (world.rank() + 1 + s) % s;
		auto const left  = (world.rank() - 1 + s) % s;

		{
			using Container = std::vector<double>;
			Container c(10, world.rank());
			world.send_receive_n(c.begin(), c.size(), left, right);
			BOOST_TEST( c.front() == right );
		}
		{
			using Container = std::vector<double>;
			Container c(10, world.rank());
			world.send_receive(c.begin(), c.end(), left, right);
			BOOST_TEST( c.front() == right );
		}
		{
			using Container = std::list<double>;
			Container c(10, world.rank());
			world.send_receive(c.begin(), c.end(), left, right);
			BOOST_TEST( c.front() == right );
		}
		{
			using Container = std::list<std::string>;
			Container c(10, std::to_string(world.rank()));
			world.send_receive_n(c.begin(), c.size(), left, right);
			BOOST_TEST( c.front() == std::to_string(right) );
		}
		{
			std::array<int, 10> buffer{};
			buffer[5] = world.rank();
			std::array<int, 10> buffer2{};
			buffer2[5] = -1;
			world.send_receive_n(
				buffer.data(), 10, left,
				buffer2.data(), right
			);
			BOOST_TEST(buffer2[5] == right);
		}
		{
			std::vector<int> buffer(10);
			buffer[5] = world.rank();
			std::vector<int> buffer2(10);
			buffer2[5] = -1;
			world.send_receive(
				buffer.begin(), buffer.end(), left,
				buffer2.begin(), buffer2.end(), right
			);
			BOOST_TEST(buffer2[5] == right);
		}
		{
			std::list<std::complex<double>> b(10, std::complex<double>{});  // std::to_string(1.*world.rank()));
			world.send_receive_n(b.begin(), b.size(), left, right);
			//  BOOST_TEST( *b.begin() == std::to_string(1.*right) );
		}
		{
			mpi3::circular_communicator circle{world};
			std::vector<int>            buffer(10);
			buffer[5] = circle.coordinate();
			circle.send_receive(buffer.begin(), buffer.end(), circle.rank(circle.coordinate() - 1), circle.rank(circle.coordinate() + 1));
			BOOST_TEST( buffer[5] == right );
		}
		try {
			mpi3::ring       circle{world};
			std::vector<int> buffer(10);
			buffer[5] = circle.coordinate();
			circle.rotate(buffer.begin(), buffer.end());
			//  BOOST_TEST( buffer[5] == circle.rank(circle.coordinate() - 1) );
		} catch(std::exception& e) {
			std::cout << e.what() << '\n';
		}
	}
	{
		std::vector<int> buffer(10);
		buffer[5]        = world.rank();
		auto right = world.rank() + 1;
		if(right >= world.size()) {
			right = 0;
		}
		auto left = world.rank() - 1;
		if(left < 0) {
			left = world.size() - 1;
		}
		world.send_receive_replace(buffer.begin(), buffer.end(), left, right);
		//  MPI_Sendrecv_replace(buffer, 10, MPI_INT, left, 123, right, 123, MPI_COMM_WORLD, &status);
		BOOST_TEST( buffer[5] == right );
	}
	{
		std::vector<int> buffer(10);
		buffer[5] = world.rank();
		world.send_receive(buffer.begin(), buffer.end(), world.rank(), world.rank());
		BOOST_TEST( buffer[5] == world.rank() );
	}
	{
		std::vector<int> buffer(10);
		buffer[5]  = world.rank();
		auto right = world.rank() + 1;
		if(right >= world.size()) {
			right = 0;
		}
		auto left = world.rank() - 1;
		if(left < 0) {
			left = world.size() - 1;
		}
		world.send_receive(buffer.begin(), buffer.end(), left, right);
		BOOST_TEST( buffer[5] == right );
	}
	{
		std::vector<int> const vs = {1, 2, 3};
		switch(world.rank()) {
			case 0: world.send(vs.begin(), vs.end(), 1); break;
			case 1: {
				std::vector<int> vr(vs.size());
				auto end = world.receive(vr.begin());  // NOLINT(clang-diagnostic-deprecated-declarations)
				BOOST_TEST( end == vr.end() );
				BOOST_TEST( std::equal(vs.begin(), vs.end(), vr.begin(), vr.end()) );  // NOLINT(boost-use-ranges)
				break;
			}
			default:;
		}
	}
	{
		std::vector<int> const vs = {1, 2, 3};
		switch(world.rank()) {
			case 0: world.send(vs.begin(), vs.end(), 1); break;
			case 1: {
				std::vector<int> vr(vs.size());
				world.receive(vr.begin(), vr.end());
				BOOST_TEST( std::equal(vs.begin(), vs.end(), vr.begin(), vr.end()) );  // NOLINT(boost-use-ranges)
				break;
			}
			default:;
		}
	}
	{
		std::vector<int> const vs = {1, 2, 3};
		switch(world.rank()) {
			case 0: world.send(vs.begin(), vs.end(), 1); break;
			case 1: {
				std::list<int> vr(vs.size());
				world.receive(vr.begin());  // NOLINT(clang-diagnostic-deprecated-declarations)
				BOOST_TEST( std::equal(vs.begin(), vs.end(), vr.begin(), vr.end()) );  // NOLINT(boost-use-ranges)
				break;
			}
			default:;
		}
	}
	{
		std::vector<int> const vs = {1, 2, 3};
		switch(world.rank()) {
			case 0: world.send(vs.begin(), vs.end(), 1); break;
			case 1: {
				std::list<int> vr(vs.size());
				world.receive(vr.begin(), vr.end());
				BOOST_TEST( std::equal(vs.begin(), vs.end(), vr.begin(), vr.end()) );  // NOLINT(boost-use-ranges)
				break;
			}
			default:;
		}
	}

	return boost::report_errors();
} catch(...) {
	return 1;
}
