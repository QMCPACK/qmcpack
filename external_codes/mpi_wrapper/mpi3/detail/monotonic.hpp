#if COMPILATION_INSTRUCTIONS
(echo "#include\""$0"\"" > $0x.cpp) && clang++ -O3 -std=c++14 -Wall -Wfatal-errors -D_TEST_MPI3_DETAIL_MONOTONIC $0x.cpp -o $0x.x && $0x.x $@ && rm -f $0x.x $0x.cpp; exit
#endif
#ifndef MPI3_DETAIL_MONOTONIC_ALLOCATOR_HPP
#define MPI3_DETAIL_MONOTONIC_ALLOCATOR_HPP

#include<memory>

namespace boost{
namespace mpi3{
namespace detail{

template<class VoidPtr = void*, typename Size = std::size_t>
struct block{
	VoidPtr ptr;
	Size sz;
};

// not thread safe
template<class A>
struct monotonic : std::allocator_traits<A>{
	using allocator_type = monotonic;
	using typename std::allocator_traits<A>::value_type;
	using typename std::allocator_traits<A>::pointer;
	using typename std::allocator_traits<A>::void_pointer;
	using typename std::allocator_traits<A>::size_type;
	template<class U> struct rebind{
		using other = monotonic<typename A::template rebind<U>::other>;
	};
private:
	std::shared_ptr<block<void_pointer, size_type>> bP_; //	std::shared_ptr<void_pointer> markP_; //	std::shared_ptr<std::size_t> szP_;
	using char_ptr = typename A::template rebind<char>::other::pointer;
public:
	monotonic() = delete;
	monotonic(size_type n, A& a){
		auto pt = a.allocate(n);
		bP_ = std::shared_ptr<block<void_pointer, size_type>>{
			new block<void_pointer, size_type>{void_pointer(pt), n*sizeof(value_type)}, 
			[a, pt, n](auto p) mutable {a.deallocate(pt, n); delete p;}
		};
	}
	pointer allocate(size_type s, const void* hint = 0){
		using std::align;
		if(align(alignof(value_type), sizeof(value_type)*s, bP_->ptr, bP_->sz)){
            pointer result = pointer(bP_->ptr);
			bP_->ptr = char_ptr(bP_->ptr) + sizeof(value_type)*s; // *markP_ = char_ptr(*markP_) + sizeof(value_type)*s;
			bP_->sz -= sizeof(value_type)*s; // *szP_ -= sizeof(value_type)*s;
			return result;
		}else throw std::bad_alloc{};
		return nullptr;
	}
	void deallocate(pointer, size_type){/*perhaps TODO implement simple rewind*/}
	monotonic& operator=(monotonic const&) = delete;
	monotonic(monotonic const&) = default;
	template<class AA>
	monotonic(monotonic<AA> const& o) : bP_(o.bP_){}//markP_(o.markP_), szP_(o.szP_){}
	template<class AA> friend struct monotonic;
	monotonic(monotonic&&) = default;
	~monotonic() = default;
};

}}}

#ifdef _TEST_MPI3_DETAIL_MONOTONIC
#include<vector>
#include<cmath>
#include<iostream>
#include<numeric>
#include<cassert>

int main(){

	using boost::mpi3::detail::monotonic;

	std::allocator<double> stda;
	monotonic<std::allocator<double>> a(1024, stda);

	std::vector<double, monotonic<std::allocator<double>>> v1(64, a);
	std::vector<double, monotonic<std::allocator<double>>> vv = v1;
	
	std::iota(v1.begin(), v1.end(), 0);
	std::vector<double, monotonic<std::allocator<double>>> v2(256, a);
	std::vector<char, monotonic<std::allocator<char>>> v3(128, a);
	try{
		std::vector<double, monotonic<std::allocator<double>>> v4(1056, a);
	}catch(std::bad_alloc&){std::cerr << "no memory\n";}
	std::cout << v1[4] << '\n';
	assert(v1[13] == 13.);

}

#endif
#endif

