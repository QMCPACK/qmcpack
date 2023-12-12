#ifdef COMPILATION_INSTRUCTIONS
nvcc -ccbin cuda-c++ -std=c++14 $0 -o $0x && $0x && rm -f $0x; exit
#endif

#include "../../../../multi/array.hpp"
#include "../../../../multi/detail/stack_allocator.hpp"
#include "../../../../multi/detail/cuda/allocator.hpp"

#include<iostream>

namespace multi = boost::multi;
namespace cuda = multi::detail::memory::cuda;

using std::cout;

int main(){
	{
		std::size_t stack_size = 4000;
		multi::stack_buffer<cuda::allocator<>> buf{stack_size};
		for(int i = 0; i != 3; ++i){
			cout<<"pass "<< i << std::endl;
			{
				multi::array<double, 2, multi::stack_allocator<double, cuda::allocator<>>> A({2, 10}, &buf);
				multi::array<int,    2, multi::stack_allocator<double, cuda::allocator<>>> B({3, 10}, &buf);
				multi::array<double, 2, multi::stack_allocator<double, cuda::allocator<>>> C({4, 10}, &buf);
				for(int j = 0; j != 100; ++j)
					multi::array<double, 2, multi::stack_allocator<double, cuda::allocator<>>> D({4, 10}, &buf);
				B[1][1] = 33.;
				B[2][2] = 33.;
				assert( B[1][1] == B[2][2] );
			}
			cout
				<<"  size: "<< buf.size() 
				<<"\n  hits: "<< buf.hits() 
				<<"\n  misses "<< buf.misses() 
				<<"\n  allocated(bytes) "<< buf.allocated_bytes() 
				<<"\n  deallocated(bytes) "<< buf.deallocated_bytes()
				<<"\n  max_needed(bytes) "<< buf.max_needed()
				<<"\n  stack recovered(bytes) " << buf.stack_recovered()
				<< std::endl
			;
			assert( buf.allocated_bytes() == buf.deallocated_bytes() );
			if(buf.max_needed() > buf.size()) buf.reset(buf.max_needed());
		}
	}
	assert( cuda::allocation_counter::n_allocations == 1 );
}

