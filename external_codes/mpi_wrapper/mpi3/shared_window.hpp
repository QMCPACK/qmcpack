#if COMPILATION_INSTRUCTIONS
(echo "#include\""$0"\"" > $0x.cpp) && mpic++ -O3 -std=c++14 -Wall -Wfatal-errors -D_TEST_MPI3_SHARED_WINDOW $0x.cpp -o $0x.x && time mpirun -n 3 $0x.x $@ && rm -f $0x.x $0x.cpp; exit
#endif
//  (C) Copyright Alfredo A. Correa 2018.
#ifndef MPI3_SHARED_WINDOW_HPP
#define MPI3_SHARED_WINDOW_HPP

#include "../mpi3/shared_communicator.hpp"
#include "../mpi3/dynamic_window.hpp"
#include "../mpi3/group.hpp"

#include <boost/interprocess/containers/vector.hpp>
#include <boost/interprocess/allocators/allocator.hpp>
#include <boost/interprocess/managed_shared_memory.hpp>

#include<mpi.h>

namespace boost{
namespace mpi3{

template<class T /*= void*/>
struct shared_window : window<T>{
//	shared_communicator& comm_;
	shared_window(shared_communicator& comm, mpi3::size_t n, int disp_unit = alignof(T)) : //sizeof(T)) : // here we assume that disp_unit is used for align
		window<T>()//, comm_{comm}
	{
		void* base_ptr = nullptr;
		int s = MPI_Win_allocate_shared(n*sizeof(T), disp_unit, MPI_INFO_NULL, &comm, &base_ptr, &this->impl_);
		if(s != MPI_SUCCESS) throw std::runtime_error{"cannot create shared window"};
	}
	shared_window(shared_communicator& comm, int disp_unit = alignof(T)) : 
		shared_window(comm, 0, disp_unit)
	{}
	shared_window(shared_window const&) = default;
	shared_window(shared_window&& other) : window<T>{std::move(other)}{}//, comm_{other.comm_}{}
	using query_t = std::tuple<mpi3::size_t, int, void*>;
	query_t query(int rank = MPI_PROC_NULL) const{
		query_t ret;
		MPI_Win_shared_query(this->impl_, rank, &std::get<0>(ret), &std::get<1>(ret), &std::get<2>(ret));
		return ret;
	}
	template<class TT = T>
	mpi3::size_t size(int rank = 0) const{
		return std::get<0>(query(rank))/sizeof(TT);
	}
	int disp_unit(int rank = 0) const{return std::get<1>(query(rank));}
	template<class TT = T>
	TT* base(int rank = 0) const{return static_cast<TT*>(std::get<2>(query(rank)));}
};

template<class T /*= char*/> 
shared_window<T> shared_communicator::make_shared_window(
	mpi3::size_t size
){
	return shared_window<T>(*this, size);
}

template<class T /*= char*/>
shared_window<T> shared_communicator::make_shared_window(){
	return shared_window<T>(*this);//, sizeof(T));
}

namespace intranode{

template<class T> struct array_ptr;

template<>
struct array_ptr<const void>{
	using T = const void;
	std::shared_ptr<shared_window<>> wSP_;
	std::ptrdiff_t offset = 0;
	array_ptr(std::nullptr_t = nullptr){}
	array_ptr(array_ptr const& other) = default;
	array_ptr& operator=(array_ptr const& other) = default;
	array_ptr& operator=(std::nullptr_t){wSP_.reset(); return *this;}
	bool operator==(std::nullptr_t) const{return (bool)wSP_;}
	bool operator!=(std::nullptr_t) const{return not operator==(nullptr);}
};

template<>
struct array_ptr<void>{
	using T = void;
	using element_type = T;
	std::shared_ptr<shared_window<>> wSP_;
	std::ptrdiff_t offset = 0;
	array_ptr(std::nullptr_t = nullptr){}
	array_ptr(array_ptr const& other) = default;
	array_ptr& operator=(array_ptr const& other) = default;
	array_ptr& operator=(std::nullptr_t){wSP_.reset(); return *this;}
	bool operator==(std::nullptr_t) const{return (bool)wSP_;}
	bool operator!=(std::nullptr_t) const{return not operator==(nullptr);}
};

template<class T>
struct array_ptr{
	using element_type = T;
	using difference_type = std::ptrdiff_t;
	using value_type = std::decay_t<T>; // std::remove_cv_t<T>; // T until C++20?
	using pointer = T*; // TODO self?
	using reference = T&; //TODO fancy_reference?
	using iterator_category = std::random_access_iterator_tag;
	std::shared_ptr<shared_window<value_type>> wSP_;
	std::ptrdiff_t offset = 0;
	array_ptr(){}
	array_ptr(std::nullptr_t){}
//	array_ptr(std::nullptr_t = nullptr) : offset(0){}
	array_ptr(array_ptr const& other) = default;
//	array_ptr(T* const& other = nullptr) : offset(0){}
//	array_ptr(T* const& other = nullptr) : offset(0){}
	array_ptr& operator=(array_ptr const& other) = default;
	array_ptr& operator=(std::nullptr_t){return *this;}
	~array_ptr() = default;
	T& operator*() const{return *((T*)(wSP_->base(0)) + offset);}
	T& operator[](int idx) const{return ((T*)(wSP_->base(0)) + offset)[idx];}
	T* operator->() const{return (T*)(wSP_->base(0)) + offset;}
	T* get() const{return wSP_->base(0) + offset;}
	explicit operator bool() const{return (bool)wSP_;}//.get();}
	bool operator==(std::nullptr_t) const{return (bool)wSP_;}
	bool operator!=(std::nullptr_t) const{return not operator==(nullptr);}
	operator array_ptr<T const>() const{
		array_ptr<T const> ret;
		ret.wSP_ = wSP_;
		return ret;
	}
	operator array_ptr<void const>() const{
		array_ptr<void const> ret;
		ret.wSP_ = wSP_;
		return ret;
	}
	array_ptr operator+(std::ptrdiff_t d) const{
		array_ptr ret(*this);
		ret += d;
		return ret;
	}
	std::ptrdiff_t operator-(array_ptr other) const{return offset-other.offset;}
	array_ptr& operator--(){--offset; return *this;}
	array_ptr& operator++(){++offset; return *this;}
	array_ptr& operator-=(std::ptrdiff_t d){offset -= d; return *this;}
	array_ptr& operator+=(std::ptrdiff_t d){offset += d; return *this;}
	bool operator==(array_ptr<T> const& other) const{
		return wSP_->base(0) == other.wSP_->base(0) and offset == other.offset;
	}
	bool operator!=(array_ptr<T> const& other) const{return not((*this)==other);}
	bool operator<(array_ptr<T> const& other) const{
		return wSP_->base(0) + offset < other.wSP_->base(0) + other.offset;
	}
	friend pointer to_address(array_ptr const& ap){return ap.wSP_->base(0) + ap.offset;}
};

template<class T, class F>
F for_each(array_ptr<T> first, array_ptr<T> last, F f){
	if(first == last) return f;
//	assert(first.wSP_->comm_ == last.wSP_->comm_);
	auto& comm = first.wSP_->comm_;
	// TODO do a partitioning std::for_each
	if(mpi3::group(*first.wSP_).root()) std::for_each(to_address(first), to_address(last), f);
//	if(first.wSP_->comm_.root()) std::for_each(to_address(first), to_address(last), f);
	comm.barrier();
	first.wSP_->fence();
	first.wSP_->fence();
}

template<typename T, typename Size, typename... Args>
array_ptr<T> uninitialized_fill_n(array_ptr<T> first, Size n, Args&&...args){
	if(n == 0) return first;
	if(mpi3::group(*first.wSP_).root()) std::uninitialized_fill_n(to_address(first), n, std::forward<Args>(args)...); // change to to_pointer
//	if(first.wSP_->comm_.root()) std::uninitialized_fill_n(to_address(first), n, std::forward<Args>(args)...); // change to to_pointer
	first.wSP_->fence();
	first.wSP_->fence();
//	first.wSP_->comm_.barrier();
	return first + n;
}
template<typename T, typename Size>
array_ptr<T> destroy_n(array_ptr<T> first, Size n){
//	if(first.wSP_->comm_.root()){
	if(mpi3::group(*first.wSP_).root()){
		auto first_ptr = to_address(first);
		for(; n > 0; (void) ++first_ptr, --n) first->~T();
	}
	first.wSP_->fence();
	first.wSP_->fence();
//	first.wSP_->comm_.barrier();
	return first + n;
}

//uninitialized_fill_n(
//			this->data_, this->num_elements(), 
//			typename array::element(std::forward<Args>(args)...)
//		)

template<class T> struct allocator{
	template<class U> struct rebind{typedef allocator<U> other;};
	using value_type = T;
	using pointer = array_ptr<T>;
	using const_pointer = array_ptr<T const>;
	using reference = T&;
	using const_reference = T const&;
	using size_type = std::size_t;
	using difference_type = std::ptrdiff_t;

	mpi3::shared_communicator& comm_;
	allocator(mpi3::shared_communicator& comm) : comm_(comm){}
	allocator() = delete;
	~allocator() = default;
	allocator(allocator const& other) : comm_(other.comm_){
	//	std::cout << "popd size " << other.comm_.size() << '\n';
	}
	template<class U>
	allocator(allocator<U> const& other) : comm_(other.comm_){}

//	template<class ConstVoidPtr = const void*>
	array_ptr<T> allocate(size_type n, const void* /*hint*/ = 0){
	/*	std::cerr << "allocating " << n << std::endl; 
		std::cerr << " from rank " << comm_.rank() << std::endl;
		std::cerr << "active1 " << bool(comm_) << std::endl;
		std::cerr << "active2 " << bool(&comm_ == MPI_COMM_NULL) << std::endl;
		std::cerr << "size " << comm_.size() << std::endl;
		std::cout << std::flush;*/
	//	comm_.barrier();
		array_ptr<T> ret = 0;
		if(n == 0){
			ret.wSP_ = std::make_shared<shared_window<T>>(
				comm_.make_shared_window<T>(0)
			);
			return ret;
		}
		ret.wSP_ = std::make_shared<shared_window<T>>(
			comm_.make_shared_window<T>(comm_.root()?n:0)
		//	comm_.allocate_shared(comm_.rank()==0?n*sizeof(T):1)
		);
		return ret;
	}
	void deallocate(array_ptr<T> ptr, size_type){ptr.wSP_.reset();}
	allocator& operator=(allocator const& other){
		assert( (*this)==other ); // TODO make comm a shared_ptr
		return *this;
	}
	bool operator==(allocator const& other) const{return comm_ == other.comm_;}
	bool operator!=(allocator const& other) const{return not(other == *this);}
	template<class U, class... Args>
	void construct(U* p, Args&&... args){
	//	std::cout << "construct: I am " << comm_.rank() << std::endl;
		::new((void*)p) U(std::forward<Args>(args)...);
	}
	template< class U >	void destroy(U* p){
	//	std::cout << "destroy: I am " << comm_.rank() << std::endl;
		p->~U();
	}
};

struct is_root{
	shared_communicator& comm_;
	template<class Alloc>
	is_root(Alloc& a) : comm_(a.comm_){}
	bool root(){return comm_.root();}
};

}

}}

#ifdef _TEST_MPI3_SHARED_WINDOW

#include "../mpi3/main.hpp"
#include "../mpi3/ostream.hpp"

namespace mpi3 = boost::mpi3; 

int mpi3::main(int, char*[], mpi3::communicator world){

	mpi3::ostream wout(world);

	double* p;
	double* b;
	wout << (p < b) << std::endl;
	wout << (p >= b) << std::endl;

	mpi3::shared_communicator node = world.split_shared();

	mpi3::shared_window<int> win = node.make_shared_window<int>(node.root()?node.size():0);

	assert(win.base() != nullptr);
	assert(win.size() == node.size());

	win.base()[node.rank()] = node.rank() + 1;
	node.barrier();
	for(int i = 0; i != node.size(); ++i) assert(win.base()[i] == i + 1);

	{
		mpi3::shared_window<int> win = node.make_shared_window<int>(0);
	}

	return 0;
}

#endif
#endif

