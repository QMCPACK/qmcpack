#include <thrust/device_vector.h>

int main(){
//	thrust::device_vector<int> D(5);
//	assert( D.size() == 5 );

//	cudaDeviceSynchronize();
	std::allocator<int> alloc;
	int* p = alloc.allocate(10);
	p[0] = 2;
	return p[0] + 1;
}
