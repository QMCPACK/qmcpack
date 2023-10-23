#include "../../../mpi3/nccl/communicator.hpp"
#include "../../../mpi3/nccl/universal_communicator.hpp"

#include "../../../mpi3/main.hpp"

#include <thrust/system/cuda/memory.h>
#include <thrust/device_vector.h>
#include <thrust/host_vector.h>

#include <thrust/complex.h>

namespace mpi3 = boost::mpi3;

int mpi3::main(int /*argc*/, char** /*argv*/, mpi3::communicator WORLD) {
	assert(WORLD.size() == 4);

	cudaSetDevice(WORLD.rank());

	auto HEMI = WORLD / 2;

	using Communicator = mpi3::nccl::universal_communicator;

	Communicator magnesium{HEMI};
	assert(magnesium.rank() == HEMI.rank());

	using T = thrust::complex<double>;  // int64_t;
//  thust::device_vector<T, thrust::cuda::universal_allocator<T>> A(1000, world.rank());
	thrust::device_vector<T, thrust::cuda::allocator<T>> A(1000, T{1.*WORLD.rank()});

	magnesium.all_reduce_n(A.data(), A.size(), A.data());
	thrust::host_vector<T> H = A;

	std::cout<<"[WORLD rank"<< WORLD.rank() <<" HEMI rank"<< HEMI.rank() <<"] result:"<< H[0] <<std::endl;

	assert( magnesium.size() == 2 );

	std::vector<double> V(100, 1.);
	magnesium.all_reduce_n(V.data(), V.size(), V.data());
	assert( V[10] == magnesium.size() );

	switch(magnesium.rank()) {
	case 0: {
		magnesium.send_n(A.data(), A.size(), 1);
	}
	case 1: {
		thrust::device_vector<T, thrust::cuda::allocator<T>> B(1000, T{});
		magnesium.receive_n(B.data(), B.size(), 0);
		assert( A == B );
	}
	}

//	thrust::device_vector<int, thrust::cuda::allocator<int>> singleton(1);
	std::vector<int> singleton(1);
	if(magnesium.rank() == 0) { singleton[0] = 99; }
	magnesium.broadcast_n(singleton.data(), 1);
	cudaStreamSynchronize(NULL);
	assert( singleton[0] == 99 );

	auto magnesium2{magnesium};
	assert( magnesium2.size() == magnesium.size() );

	return 0;
}
