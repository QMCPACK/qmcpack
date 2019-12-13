# if COMPILATION_INSTRUCTIONS
(echo "#include\""$0"\"" > $0x.cpp) && mpic++ -O3 -std=c++14 -Wall -Wextra -Wpedantic -Wfatal-errors -D_TEST_MPI3_SHARED_WINDOW $0x.cpp -o $0x.x && time mpirun -n 3 $0x.x $@ && rm -f $0x.x $0x.cpp; exit
#endif
//  (C) Copyright Alfredo A. Correa 2018.
#ifndef MPI3_SHARED_WINDOW_HPP
#define MPI3_SHARED_WINDOW_HPP

#include "../mpi3/shared_communicator.hpp"
#include "../mpi3/dynamic_window.hpp"
#include "../mpi3/group.hpp"

//#include <boost/interprocess/containers/vector.hpp>
//#include <boost/interprocess/allocators/allocator.hpp>
//#include <boost/interprocess/managed_shared_memory.hpp>

#include<mpi.h>

namespace boost{
namespace mpi3{

template<class T>
struct shared_window : window<T>{
	shared_communicator& comm_;
	shared_window(shared_communicator& comm, mpi3::size_t n, int disp_unit = alignof(T)) : //sizeof(T)) : // here we assume that disp_unit is used for align
		window<T>(), comm_{comm}
	{
		void* base_ptr = nullptr;
		auto e = static_cast<enum error>(
			MPI_Win_allocate_shared(
				n*sizeof(T), disp_unit, 
				MPI_INFO_NULL, comm.get(), &base_ptr, &this->impl_
			)
		);
		if(e != mpi3::error::success) throw std::system_error{e, "cannot win_alloc"};
	}
	shared_window(shared_communicator& comm, int disp_unit = alignof(T)) : 
		shared_window(comm, 0, disp_unit)
	{}
	shared_window(shared_window const&) = default;
	shared_window(shared_window&& other) noexcept : window<T>{std::move(other)}, comm_{other.comm_}
	{}
	group get_group() const{
		group r; MPI3_CALL(MPI_Win_get_group)(this->impl_, &(&r)); return r;
	}
	shared_communicator& get_communicator() const{return comm_;}
	struct query_t{
		mpi3::size_t size;
		int disp_unit;
		void* base;
	};
	struct query_t query(int rank = MPI_PROC_NULL) const{
		query_t r;
		MPI3_CALL(MPI_Win_shared_query)(this->impl_, rank, &r.size, &r.disp_unit, &r.base);
		return r;
	}
	template<class TT = T>
	mpi3::size_t size(int rank = 0) const{return query(rank).size/sizeof(TT);}
	int disp_unit(int rank = 0) const{return query(rank).disp_unit;}
	template<class TT = T>
	TT* base(int rank = 0) const{return static_cast<TT*>(query(rank).base);}
};

template<class T /*= char*/> 
shared_window<T> shared_communicator::make_shared_window(mpi3::size_t size){
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
	explicit operator T*() const{return get();}
	explicit operator bool() const{return (bool)wSP_;}//.get();}
	bool operator==(std::nullptr_t) const{return not (bool)wSP_;}
	bool operator!=(std::nullptr_t) const{return not operator==(nullptr);}
	operator array_ptr<T const>() const{
		array_ptr<T const> ret;
		ret.wSP_ = wSP_;
		ret.offset = offset;
		return ret;
	}
	operator array_ptr<void const>() const{
		array_ptr<void const> ret;
		ret.wSP_ = wSP_;
		ret.offset = offset;
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
	static element_type* to_address(array_ptr p) noexcept{
		return p.wSP_->base(0) + p.offset;
	}
	friend pointer to_address(array_ptr const& p){return array_ptr::to_address(p);}
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

template<typename T, typename Size, typename TT>
array_ptr<T> uninitialized_fill_n(array_ptr<T> first, Size n, TT const& val){
	if(n == 0) return first;
	if(mpi3::group(*first.wSP_).root()) std::uninitialized_fill_n(to_address(first), n, val); // change to to_pointer
//	if(first.wSP_->comm_.root()) std::uninitialized_fill_n(to_address(first), n, std::forward<Args>(args)...); // change to to_pointer
	first.wSP_->fence();
	first.wSP_->fence();
//	first.wSP_->comm_.barrier();
	return first + n;
}

template<typename T, typename Size>
array_ptr<T> destroy_n(array_ptr<T> first, Size n){
	if(n == 0) return first;
	if(mpi3::group(*first.wSP_).root()) { 
		auto first_ptr = to_address(first);
                for(; n > 0; (void) ++first_ptr, --n) first->~T();
	}
	first.wSP_->fence();
	first.wSP_->fence();
	return first + n;
}

template<class It1, typename T, typename Size>
array_ptr<T> copy_n(It1 first, Size n, array_ptr<T> d_first){
//	if(n == 0) return d_first;
	d_first.wSP_->fence();
	using std::copy_n;
	if(mpi3::group(*d_first.wSP_).root()) copy_n(first, n, to_address(d_first));
	d_first.wSP_->fence();
	return d_first + n;
}

template<class It1, typename T>
array_ptr<T> copy(It1 first, It1 last, array_ptr<T> d_first){
	if(first == last) return d_first;
	first.wSP_->fence();
	using std::copy;
	if(mpi3::group(*d_first.wSP_).root()) copy(first, last, to_address(d_first));
	first.wSP_->fence();
	using std::distance;
	return d_first + distance(first, last);
}

template<class It1, class Size, typename T>
array_ptr<T> uninitialized_copy_n(It1 f, Size n, array_ptr<T> d){
	if(n == 0) return d;
	f.wSP_->fence();
	using std::uninitialized_copy_n;
	if(mpi3::group(*d.wSP_).root()) uninitialized_copy_n(f, n, to_address(d));
	f.wSP_->fence();
	return d + n;
}

template<class It1, typename T>
array_ptr<T> uninitialized_copy(It1 f, It1 l, array_ptr<T> d){
	if(f == l) return d;
	f.wSP_->fence();
	using std::uninitialized_copy;
	if(mpi3::group(*d.wSP_).root()) uninitialized_copy(f, l, to_address(d));
	f.wSP_->fence();
	using std::distance;
	return d + distance(f, l);
}

template<class T, class Size>
array_ptr<T> uninitialized_default_construct_n(array_ptr<T> f, Size n){
	if(n == 0) return f;
#if __cplusplus >= 201703L
	using std::uninitialized_default_construct_n;
#endif
	f.wSP_->fence();
	if(group(*f.wSP_).root())
		uninitialized_default_construct_n(to_address(f), n);
	f.wSP_->fence();
	return f + n;
}

template<class T, class Size>
array_ptr<T> uninitialized_value_construct_n(array_ptr<T> f, Size n){
	if(n == 0) return f;
#if __cplusplus >= 201703L
	using std::uninitialized_value_construct_n;
#endif
	f.wSP_->fence();
	if(group(*f.wSP_).root()) uninitialized_value_construct_n(to_address(f), n);
	f.wSP_->fence();
	return f + n;
}

template<class T = void> struct allocator{
	template<class U> struct rebind{typedef allocator<U> other;};
	using value_type = T;
	using pointer = array_ptr<T>;
	using const_pointer = array_ptr<T const>;
//	using reference = T&;
//	using const_reference = T const&;
	using size_type = mpi3::size_t; // std::size_t; 
	using difference_type = std::make_signed_t<size_type>;//std::ptrdiff_t;

	mpi3::shared_communicator& comm_;
	allocator() = delete;
	allocator(mpi3::shared_communicator& comm) : comm_(comm){}
	allocator(allocator const& other) : comm_(other.comm_){}
	~allocator() = default;
	template<class U>  allocator(allocator<U> const& o) : comm_(o.comm_){}

	array_ptr<T> allocate(size_type n, const void* /*hint*/ = 0){
		array_ptr<T> ret = 0;
		if(n == 0){
			//ret.wSP_ = std::make_shared<shared_window<T>>(
			//	comm_.make_shared_window<T>(0)
			//);
                        ret.wSP_.reset(new shared_window<T>{
                                comm_.make_shared_window<T>(comm_.root()?n:0)
                        });
			return ret;
		}
		//ret.wSP_ = std::make_shared<shared_window<T>>(
		//	comm_.make_shared_window<T>(comm_.root()?n:0)
		//);
                ret.wSP_.reset(new shared_window<T>{
                        comm_.make_shared_window<T>(comm_.root()?n:0)
                        });
		return ret;
	}
	void deallocate(array_ptr<T> ptr, size_type){ptr.wSP_.reset();}
	allocator& operator=(allocator const& other){
		assert( (*this)==other ); // TODO make comm a shared_ptr
		return *this;
	}
	bool operator==(allocator const& other) const{return comm_ == other.comm_;}
	bool operator!=(allocator const& other) const{return not(other == *this);}
	template<class U, class... As>
	void construct(U* p, As&&... as){::new((void*)p) U(std::forward<As>(as)...);}
	template< class U >	void destroy(U* p){p->~U();}
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

