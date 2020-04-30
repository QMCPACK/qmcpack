#ifdef COMPILATION_INSTRUCTIONS
(echo "#include\""$0"\"" > $0.cpp) && c++ -std=c++17 -Wall -Wextra -Wfatal-errors -D_TEST_BOOST_MULTI_MEMORY_INSTRUMENTED $0.cpp -o $0x && valgrind $0x && rm $0x $0.cpp; exit
#endif
#ifndef BOOST_MULTI_MEMORY_INSTRUMENTED_HPP
#define BOOST_MULTI_MEMORY_INSTRUMENTED_HPP

#include<cstddef> // max_align_t
#include<memory>
#include<stack>
#include<cassert>
#include<stdexcept>
#include<map>

#include<memory_resource>

//#include "../memory/block.hpp"
//#include "../memory/allocator.hpp"

#include<numeric> // accumulate
#include<iostream>
namespace boost{
namespace multi{
namespace memory{

template<
	class MemoryResource = std::pmr::memory_resource, 
	typename SizeType = std::size_t, 
	typename VoidPointer = void*
>
class instrumented{
	MemoryResource* back_ = std::pmr::get_default_resource();
	using void_pointer = VoidPointer;
	using size_type = SizeType;
	std::map<void_pointer, std::ptrdiff_t> blocks_;
public:
	instrumented() = default;
	instrumented(instrumented const&) = delete;
	std::map<void_pointer, std::ptrdiff_t> const& blocks() const{return blocks_;}
	auto leak(){
		return std::accumulate(
			blocks_.begin(), blocks_.end(), size_type{}, [](auto a, auto&& e){return a+e.second;}
		);
	}
	typename instrumented::void_pointer
	allocate(size_type required_bytes, size_type align = alignof(std::max_align_t)){
		std::cout << "allocating " << required_bytes << std::endl;
		auto ret = back_->allocate(required_bytes, align);
		blocks_[ret] += required_bytes;
		return ret;
	}
	void deallocate(typename instrumented::void_pointer p, typename instrumented::size_type discarded_bytes, size_type align = alignof(std::max_align_t)){
		std::cout << "deallocating " << discarded_bytes << std::endl;
		back_->deallocate(p, discarded_bytes, align);
		blocks_[p] -= discarded_bytes;
	}
};
}}}

#include "../../multi/memory/allocator.hpp"

namespace boost{
namespace multi{
namespace memory{

template<class T>
using instrumented_allocator = multi::memory::allocator<T, multi::memory::instrumented<>>;

}}}

#if _TEST_BOOST_MULTI_MEMORY_INSTRUMENTED

#include "../../multi/array.hpp"
#include "../memory/monotonic.hpp"

#include<iostream>
#include<vector>
#include<cmath>

namespace multi = boost::multi;
using std::cout;

int main(){
{
	multi::memory::instrumented<> im;
	auto p1 = im.allocate(1*sizeof(double), alignof(double));
	auto p2 = im.allocate(255*sizeof(double), alignof(double));
	im.deallocate(p2, 255*sizeof(double));
	im.deallocate(p1, 1*sizeof(double));
	assert( im.blocks().size() == 2 );
	assert( not im.leak() );
	{
		multi::memory::instrumented<> im;
		multi::memory::allocator<double, multi::memory::instrumented<> > A(&im);
		double* p = A.allocate(1);
		A.construct(p, 8.);
		assert( *p == 8. );
		double* arr = A.allocate(255);
		A.construct(arr, 81.);
		assert( *arr == 81. );
		A.deallocate(arr, 255);
		A.deallocate(p, 1);
		assert( not im.leak() );
	}
	{
		multi::memory::instrumented<> im;
		multi::memory::allocator<double, multi::memory::instrumented<> > A(&im);
		{
			std::vector<double, decltype(A)> v(A);
			v.push_back(99);
			v.push_back(10);
			v.push_back(12);
			v.resize(1);
			v.push_back(10);
		}
		assert( not im.leak() );
	}
	{
		multi::memory::instrumented<> im;
		using alloc = multi::memory::instrumented_allocator<double>;
		alloc A(&im);
		alloc B(A);
		{
			multi::static_array<double, 1, alloc> arr1({10}, 99., A);
			multi::static_array<double, 2, alloc> arr2({10, 20}, 99., A);
			multi::static_array<double, 3, alloc> arr3({10, 20, 30}, 99., B);
			multi::static_array<double, 3, alloc> brr3(arr3);
		}
		assert( not im.leak() );
	}
	{
		multi::memory::instrumented<> im;
		multi::memory::allocator<double, multi::memory::instrumented<> > A(&im);
		{
			multi::array<double, 1, decltype(A)> arr1({10}, 99., A);
			multi::array<double, 2, decltype(A)> arr2({10, 20}, 99., A);
			multi::array<double, 3, decltype(A)> arr3({10, 20, 30}, 99., A);
		}
		assert( not im.leak() );
	}
	{
		multi::memory::instrumented<> im;
		{
			using alloc = multi::memory::instrumented_allocator<double>;
			multi::array<double, 3, alloc> A({10, 20, 30}, 99., &im);
			multi::array<double, 3, alloc> B({10, 20, 30}, 11., &im);
			B = std::move(A); assert( empty(A) );
		}
		assert( not im.leak() );
	}
	std::cout << "-------------------" << std::endl;
	{
		multi::memory::instrumented<> im;
		{
			using alloc = multi::memory::instrumented_allocator<double>;
		//	using alloc = std::allocator<double>;
			multi::array<double, 1, alloc> arr1({10}, 99., &im);
			multi::array<double, 2, alloc> arr2({10, 20}, 99., &im);
			multi::array<double, 3, alloc> arr3({10, 20, 30}, 99., &im);
			arr1.reextent({20});
			arr2.reextent({200, 10});
			arr3.reextent({201, 10, 100});
			arr2.clear();
			multi::array<double, 1, alloc> brr1 = arr1;
			multi::array<double, 2, alloc> brr2(arr2);
			assert( arr3.num_elements() == 201*10*100 );
			multi::array<double, 3, alloc> brr3(std::move(arr3));
			assert( arr3.num_elements() == 0 );
			assert( brr3.num_elements() == 201*10*100 );
		}
		assert( not im.leak() );
	}
	return 0;
//	return 0;
}
}
#endif 
#endif

