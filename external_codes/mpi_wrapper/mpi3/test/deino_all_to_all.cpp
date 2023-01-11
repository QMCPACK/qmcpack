// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;autowrap:nil;-*-
// Copyright 2022 Alfredo A. Correa

// based on http://mpi.deino.net/mpi_functions/MPI_Alltoall.html

#include <mpi3/main.hpp>
#include <mpi3/ostream.hpp>

#include "mpi.h"

namespace mpi3 = boost::mpi3;

//template<class InputIt, class Size, class OutputIt>
//void do_all_to_all_n(mpi3::communicator& comm, InputIt first, Size count, OutputIt result) {
// #ifndef MPICH_VERSION
// 	comm.all_to_all_n(first, count, result);
// #else
//     if(first != result) {
//         comm.all_to_all_n(first, count, result);
//     } else {
// 		std::vector<MPI_Request> reqs(comm.size(), MPI_REQUEST_NULL);
		
//         assert(count % comm.size() == 0);
// 		for(int iproc = 0; iproc < comm.size(); iproc++) {
// 			MPI_Isendrecv_replace(first, count/comm.size(), mpi3::detail::basic_datatype<typename std::iterator_traits<InputIt>::value_type>(), iproc, 0, MPI_ANY_SOURCE, MPI_ANY_TAG, &comm, &reqs[iproc]);
// 		}

// 		MPI_Waitall(static_cast<int>(reqs.size()), reqs.data(), MPI_STATUSES_IGNORE);
//     }
// #endif
//}

int mpi3::main(int /*argc*/, char** /*argv*/, mpi3::communicator world) try {
	std::size_t chunk = 5;

	auto sb = std::vector<int>(static_cast<std::size_t>(world.size()) * chunk);
    std::iota(sb.begin(), sb.end(), 40000 + world.rank()*100);

	auto rb = std::vector<int>(static_cast<std::size_t>(world.size()) * chunk);

	auto sz = static_cast<std::size_t>(world.size()); assert( sz != 0 );
	assert( sb.size() % sz == 0);

	world.all_to_all_n(sb.data(), sb.size()/sz, rb.data());
	// do_all_to_all_n(world, sb.data(), sb.size(), rb.data());
	// do_all_to_all_n(world, sb.data(), sb.size(), rb.data());

    mpi3::ostream wout(world);
    std::copy(sb.begin(), sb.end(), std::ostream_iterator<int>(wout<<"sb = ", ", ")); wout<<std::endl;
    std::copy(rb.begin(), rb.end(), std::ostream_iterator<int>(wout<<"rb = ", ", ")); wout<<std::endl;

	world.all_to_all_inplace_n(sb.data(), sb.size()/sz);
	// do_all_to_all_n(world, sb.data(), sb.size(), sb.data());
	// world.all_to_all_n(sb.data(), sb.size()); //  , sb.data());
	std::copy(sb.begin(), sb.end(), std::ostream_iterator<int>(wout<<"sb (inplace) = ", ", ")); wout<<std::endl;

	assert(sb == rb);

	return 0;
} catch(...) {return 0;}

#if 0
#include "../../mpi3/main.hpp"
#include "../../mpi3/process.hpp"

namespace mpi3 = boost::mpi3;
using std::cout;


int mpi3::main(int, char*[], mpi3::communicator world) {
	assert(world.size() == 4);

	std::vector<int> send_buff = {
		world.rank() * 10 + 0,
		world.rank() * 10 + 1,
		world.rank() * 10 + 2,
		world.rank() * 10 + 3,
	};

	assert((int)send_buff.size() == world.size());

	std::vector<int> recv_buff(world.size(), 99.);
	{

		assert((int)recv_buff.size() == world.size());
		world.all_to_all(send_buff.begin(), recv_buff.begin());

		if(world.rank() == 0)
			assert(recv_buff[0] == 00 and recv_buff[1] == 10 and recv_buff[2] == 20 and recv_buff[3] == 30);

		if(world.rank() == 1)
			assert(recv_buff[0] == 01 and recv_buff[1] == 11 and recv_buff[2] == 21 and recv_buff[3] == 31);
	}
	std::vector<int> recv_buff2(world.size(), 99.);
	for(std::size_t i = 0; i != send_buff.size(); ++i)
		world[i] << send_buff[i] >> recv_buff2[i];

	assert(recv_buff == recv_buff2);

	{
		auto const& csend_buff = send_buff;
		auto const  r          = world & csend_buff;
		assert(recv_buff == r);
	}
	{
		auto const& csend_buff = send_buff;
		auto const  r          = world(csend_buff);
		assert(recv_buff == r);
	}
	{
		world& send_buff;
		assert(recv_buff == send_buff);
	}

	return 0;
}
#endif