#if COMPILATION_INSTRUCTIONS
mpic++ -O3 -std=c++14 -O3 -Wall -Wextra `#-Wfatal-errors` $0 -o $0x.x -lboost_serialization && time mpirun -n 1 $0x.x $@ && rm -f $0x.x; exit
#endif

#include "../../mpi3/main.hpp"
#include "../../mpi3/cartesian_communicator.hpp"

namespace mpi3 = boost::mpi3;

auto mpi3::main(int/*argc*/, char**/*argv*/, mpi3::communicator world) -> int try{

	assert( not world.is_empty() );
	auto const s = world.size();
	assert( s != 0 );

	auto right = (world.rank() + 1 + s ) % s;
	auto left  = (world.rank() - 1 + s ) % s;

	{
		using Container = std::vector<double>;
		Container c(10, world.rank());
		world.send_receive_n(c.begin(), c.size(), left, right);
		assert( c.front() == right );
	}
	{
		using Container = std::vector<double>;
		Container c(10, world.rank());
		world.send_receive(c.begin(), c.end(), left, right);
		assert( c.front() == right );
	}
	{
		using Container = std::list<double>;
		Container c(10, world.rank());
		world.send_receive(c.begin(), c.end(), left, right);
		assert( c.front() == right );
	}
	{
		using Container = std::list<std::string>;
		Container c(10, std::to_string(world.rank()));
		world.send_receive_n(c.begin(), c.size(), left, right);
		assert( c.front() == std::to_string(right) );
	}
	{
		std::array<int, 10> buffer {}; buffer [5] = world.rank();
		std::array<int, 10> buffer2{}; buffer2[5] = -1;
		world.send_receive_n(
			buffer .data(), 10, left ,
			buffer2.data()    , right
		);
		assert(buffer2[5] == right);
	}
	{
		std::vector<int> buffer (10); buffer [5] = world.rank();
		std::vector<int> buffer2(10); buffer2[5] = -1;
		world.send_receive(
			buffer .begin(), buffer .end(), left , 
			buffer2.begin(), buffer2.end(), right
		);
		assert(buffer2[5] == right);
	}
	{
		std::list<std::complex<double>> b(10, std::complex<double>{});//std::to_string(1.*world.rank()));
		world.send_receive_n(b.begin(), b.size(), left, right);
	//	assert( *b.begin() == std::to_string(1.*right) );
	}
	{
		std::vector<int> buffer (10); buffer [5] = world.rank();
	    auto right =  world.rank() + 1; if(right >= world.size()) {right = 0               ;}
    	auto left  =  world.rank() - 1; if(left  <             0) {left  = world.size() - 1;}
		world.send_receive_replace(buffer.begin(), buffer.end(), left, right);
	//  MPI_Sendrecv_replace(buffer, 10, MPI_INT, left, 123, right, 123, MPI_COMM_WORLD, &status);
		assert( buffer[5] == right );
	}
	{
		std::vector<int> buffer (10); buffer [5] = world.rank();
		world.send_receive(buffer.begin(), buffer.end(), world.rank(), world.rank());
		assert( buffer[5] == world.rank() );
	}
	{
		std::vector<int> buffer (10); buffer [5] = world.rank();
	    auto right =  world.rank() + 1; if(right >= world.size()) {right = 0               ;}
    	auto left  =  world.rank() - 1; if(left  <             0) {left  = world.size() - 1;}
		world.send_receive(buffer.begin(), buffer.end(), left, right);
		assert( buffer[5] == right );
	}
	{
		mpi3::circular_communicator circle{world};
		std::vector<int> buffer(10); buffer [5] = circle.coordinate();
		circle.send_receive(buffer.begin(), buffer.end(), circle.rank(circle.coordinate() - 1), circle.rank(circle.coordinate() + 1));
		assert( buffer[5] == right );
	}
	try {
		mpi3::ring circle{world};
		std::vector<int> buffer(10); buffer [5] = circle.coordinate();
		circle.rotate(buffer.begin(), buffer.end());
	//  assert( buffer[5] == circle.rank(circle.coordinate() - 1) );
	} catch(std::exception& e) {std::cout << e.what() <<std::endl;}
	return 0;
}catch(...) {
	return 1;
}

