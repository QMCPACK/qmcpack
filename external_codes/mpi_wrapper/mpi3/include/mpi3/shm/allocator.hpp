//-*-indent-tabs-mode: t; c-basic-offset: 4; tab-width: 4;-*-
// © Alfredo A. Correa 2019-2020

#ifndef MPI3_SHM_ALLOCATOR_HPP
#define MPI3_SHM_ALLOCATOR_HPP

//#include "../../mpi3/shared_window.hpp"
#include "../shm/ptr.hpp"

namespace boost {
namespace mpi3 {
namespace shm {

template<class T = void>
struct allocator {
	template<class U> struct rebind{using other = allocator<U>;};
	using value_type = T;
	using pointer = shm::ptr<T>;
	using const_pointer = shm::ptr<T const>;
	using size_type = typename pointer::difference_type;
	using difference_type = typename pointer::difference_type;
private:
	mpi3::shared_communicator* comm_;
//	allocator() = delete;

 public:
	explicit allocator(mpi3::shared_communicator* comm = std::addressof(static_cast<mpi3::shared_communicator&>(mpi3::environment::get_self_instance()))) : comm_{comm} {}  // NOLINT(cppcoreguidelines-pro-type-static-cast-downcast) TODO(correaa) check
	// cppcheck-suppress noExplicitConstructor
	template<class U> allocator(allocator<U> const& o) : comm_{o.comm_} {}  // NOLINT(google-explicit-constructor,hicpp-explicit-conversions) common practice

	allocator(allocator const&) = default;  // cppcheck-suppress noExplicitConstructor ; bug in cppcheck 2.3
	allocator(allocator     &&) noexcept = default;  // cppcheck-suppress noExplicitConstructor ; bug in cppcheck 2.3

	~allocator() = default;
	pointer allocate(size_type n, const void* /*hint*/ = nullptr) const{
		pointer ret;
		ret.offset_ = 0;
		ret.wP_ = new shared_window<T>{comm_->make_shared_window<T>(comm_->root()?n:0)};  // NOLINT(cppcoreguidelines-owning-memory) TODO(correaa) use unique_ptr
		ret.wP_->lock_all();
		return ret;
	}
	void deallocate(pointer p, size_type /*size*/) {
		p.wP_->unlock_all();
		delete p.wP_;  // NOLINT(cppcoreguidelines-owning-memory) TODO(correaa) use unique_ptr
	}

	allocator& operator=(allocator const&) = default;
	allocator& operator=(allocator     &&) noexcept = default;

	bool operator==(allocator const& other) const{return comm_==other.comm_;}
	bool operator!=(allocator const& other) const{return not(other == *this);}
#if 0
	template<class... As>
	void construct(pointer p, As&&... as){
		if(comm_->root()) ::new((void*)raw_pointer_cast(p)) value_type(std::forward<As>(as)...);
	}
	template<class P = pointer> void destroy(P p){
		using V = typename std::pointer_traits<P>::element_type;
		p->~V();
	}
#endif
//	void destroy(pointer p){if(comm_->root()) raw_pointer_cast(p)->~value_type();}
	template<class TT, class Size, class TT1>
	auto alloc_uninitialized_fill_n(ptr<TT> it, Size n, TT1 const& t){
		if(comm_->root()) {std::uninitialized_fill_n(raw_pointer_cast(it), n, t);} // TODO(correaa) implement with with for_each in parallel
		it.wP_->fence();
		it.wP_->fence();
		comm_->barrier();
	}
	template<class InputIt, class Size, class TT>
	auto alloc_uninitialized_copy_n(InputIt first, Size count, ptr<TT> d_first){
	//	return alloc_uninitialized_copy_n(first, count, raw_pointer_cast(d_first));
		if(comm_->root()) {std::uninitialized_copy_n(first, count, raw_pointer_cast(d_first));} // TODO(correaa) implement with with for_each in parallel
		d_first.wP_->fence();
		comm_->barrier();
	}
//	template<class TT, class Size, class ForwardIt, std::enable_if_t<not is_a_ptr<ForwardIt>{}, int> = 0>
//	auto alloc_uninitialized_copy_n(ptr<TT> first, Size count, ForwardIt d_first){
//		first.w_->fence();
//		if(comm_->root()) std::uninitialized_copy_n(raw_pointer_cast(first), count, d_first); // TODO(correaa) implement with with for_each in parallel
//		comm_->barrier();
//	}
//	template<class TT1, class Size, class TT2>
//	auto alloc_uninitialized_copy_n(ptr<TT1> first, Size count, ptr<TT2> d_first){
//		first.w_->fence();
//		if(comm_->root()) std::uninitialized_copy_n(raw_pointer_cast(first), count, raw_pointer_cast(d_first)); // TODO(correaa) implement with with for_each in parallel
//		d_first.w_->fence();
//		comm_->barrier();
//	}
	template<class Ptr, class Size, class TT = typename std::pointer_traits<Ptr>::element_type>
	auto alloc_destroy_n(Ptr first, Size count) {
		first.wP_->fence();
		if(comm_->root()) {
			std::for_each_n(first, count, [](auto const& e) {raw_pointer_cast(&e)->~TT();});
			// for( ; count > 0; --count, ++first) {raw_pointer_cast(first)->~TT();}
		}
		first.wP_->fence();
		comm_->barrier();
		return first + count;
	}
};

}  // end namespace shm
}  // end namespace mpi3
}  // end namespace boost

//#if not __INCLUDE_LEVEL__ // TEST BELOW
//#include "../../mpi3/main.hpp"
//#include "../../mpi3/shm/allocator.hpp"

//namespace mpi3 = boost::mpi3;

//int mpi3::main(int, char*[], mpi3::communicator world){

//	mpi3::shared_communicator node = world.split_shared();

//	mpi3::shm::allocator<double> A1(&node);
//	mpi3::shm::ptr<double> data1 = A1.allocate(80);

//	using ptr = decltype(data1);
//	std::pointer_traits<ptr>::pointer pp = data1;
////	double* dp = std::addressof(*data1);
////	double* dp2 = mpi3::shm::pointer_traits<ptr>::to_address(data1);

//	if(node.root()) data1[3] = 3.4;
//	node.barrier();
//	assert( data1[3] == 3.4);

//	A1.deallocate(data1, 80);

//	return 0;
//}

//#endif
#endif
