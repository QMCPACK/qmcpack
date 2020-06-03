#if COMPILATION_INSTRUCTIONS
(echo "#include\""$0"\"" > $0x.cpp) && mpic++ -O3 -std=c++14 -Wall -Wextra -Wpedantic -Wfatal-errors -D_TEST_MPI3_SHM_ALLOCATOR $0x.cpp -o $0x.x && time mpirun -n 3 $0x.x $@ && rm -f $0x.x $0x.cpp; exit
#endif
#ifndef MPI3_SHM_ALLOCATOR_HPP
#define MPI3_SHM_ALLOCATOR_HPP

//#include "../../mpi3/shared_window.hpp"
#include "../shm/pointer.hpp"

namespace boost{
namespace mpi3{
namespace shm{

//template<class... As>
//using allocator = mpi3::intranode::allocator<As...>;

//template<class T>
//using pointer = mpi3::intranode::array_ptr<T>;

template<class T = void>
struct allocator{
	template<class U> struct rebind{typedef allocator<U> other;};
	using value_type = T;
	using pointer = shm::pointer<T>;
	using const_pointer = shm::pointer<T const>;
	using size_type = mpi3::size_t;
	using difference_type = mpi3::ptrdiff_t;
private:
	mpi3::shared_communicator* comm_;
	allocator() = delete;
public:
	allocator(mpi3::shared_communicator* comm) : comm_{comm}{}
	template<class U> allocator(allocator<U> const& o) : comm_{o.comm_}{}
	~allocator() = default;

	pointer allocate(size_type n, const void* /*hint*/ = 0) const{
		pointer ret = 0;
	// this code using reset is to workaround a bug in outdated libc++ https://bugs.llvm.org/show_bug.cgi?id=42021
		ret.w_.reset(new shared_window<T>{
			comm_->make_shared_window<T>(comm_->root()?n:0)
		});
	//	ret.w_ = std::make_shared<shared_window<T>>(
	//		comm_->make_shared_window<T>(comm_->root()?n:0)
	//	);
		return ret;
	}
	void deallocate(pointer p, size_type){p.w_.reset();}
	allocator& operator=(allocator const& other) = default;
	bool operator==(allocator const& other) const{return comm_==other.comm_;}
	bool operator!=(allocator const& other) const{return not(other == *this);}
	template<class P = pointer, class... As>
	void construct(P p, As&&... as){
		if(comm_->root()) ::new((void*)p) typename std::pointer_traits<P>::element_type(std::forward<As>(as)...);
	}
	template<class P = pointer> void destroy(P p){
		using V = typename std::pointer_traits<P>::element_type;
		p->~V();
	}
};

}}
}

#ifdef _TEST_MPI3_SHM_ALLOCATOR

#include "../../mpi3/main.hpp"
#include "../../mpi3/shm/allocator.hpp"

namespace mpi3 = boost::mpi3;

#if 0
template<class A>
struct monotonic{
	A& a_;
	template<class U> struct rebind{
		using other = monotonic<typename A::template rebind<U>::other>;
	//	using other = typename std::iterator_traits<A>::template alloc_rebind<U>;
	};
	using value_type = typename A::value_type;
	using pointer = typename A::pointer;
	using size_type = typename A::size_type;
	pointer mark_;
	pointer end_;
	monotonic(monotonic::size_type s, A& a) : a_{a}{
		mark_ = a_.allocate(s);
		end_ = mark_ + s;
	}
	monotonic::pointer allocate(monotonic::size_type s, const void*/*hint*/= 0){
		auto ret = mark_;
		if(size_type(std::distance(mark_, end_)) < s) throw std::bad_alloc{};
		mark_ += s;
		return ret;
	}
	void deallocate(monotonic::pointer, monotonic::size_type){}
};
#endif

int mpi3::main(int, char*[], mpi3::communicator world){

	mpi3::shared_communicator node = world.split_shared();

	mpi3::shm::allocator<double> A1(std::addressof(node));
	mpi3::shm::pointer<double> data1 = A1.allocate(80);

	using ptr = decltype(data1);
	std::pointer_traits<ptr>::pointer pp = data1;
	double* dp = std::addressof(*data1);
	double* dp2 = mpi3::shm::pointer_traits<ptr>::to_address(data1);

	if(node.root()) data1[3] = 3.4;
	node.barrier();
	assert( *dp == *data1 );
	assert( *dp2 == *data1 );
	assert(data1);
	assert(!!data1);
	assert(not (data1 < data1));
	assert(data1[3] == 3.4);

	A1.deallocate(data1, 80);

	return 0;
}

#endif
#endif


