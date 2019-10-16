#if COMPILATION_INSTRUCTIONS
(echo "#include\""$0"\"" > $0x.cpp) && mpic++ -O3 -std=c++17 `#-Wfatal-errors` -Wall -Wextra -Wpedantic -D_TEST_BOOST_MPI3_WINDOW $0x.cpp -o $0x.x && time mpirun -np 4 $0x.x $@ && rm -f $0x.x $0x.cpp; exit
#endif
//  (C) Copyright Alfredo A. Correa 2019
#ifndef BOOST_MPI3_WINDOW_HPP
#define BOOST_MPI3_WINDOW_HPP

#include "../mpi3/error.hpp"
#include "../mpi3/type.hpp"
#include "../mpi3/detail/datatype.hpp"
#include "../mpi3/detail/call.hpp"
#include "../mpi3/communicator.hpp"

#define OMPI_SKIP_MPICXX 1  // https://github.com/open-mpi/ompi/issues/5157
#include<mpi.h>

namespace boost{
namespace mpi3{

struct target{
	int rank;
	mpi3::size_t disp;
};

template<class T = void> class panel;

template<class T = void> struct window;

template<>
struct window<void>{
protected:
	MPI_Win impl_ = MPI_WIN_NULL;
public:
	MPI_Win& operator&(){return impl_;}
	MPI_Win const& operator&() const{return impl_;}
	void clear(){
		if(impl_ != MPI_WIN_NULL) MPI3_CALL(MPI_Win_free)(&impl_);
		assert(impl_ == MPI_WIN_NULL);
	}
protected:
	window() = default;
public:
	template<class T, class Size = mpi3::size_t>
	window(communicator const& c, T* b, Size n = 0){
		MPI3_CALL(MPI_Win_create)(b, n*sizeof(T), alignof(T), MPI_INFO_NULL, c.get(), &impl_);
		assert( alignof(T) == sizeof(T) ); // to see in what type it is different
	}
	window(window const&) = delete;// see text before §4.5 in Using Adv. MPI
	window(window&& o) noexcept : impl_{std::exchange(o.impl_, MPI_WIN_NULL)}{}
	window& operator=(window const&) = delete; // see cctor
	window& operator=(window&& other){// self assignment is undefined
		clear(); swap(*this, other); return *this;
	}
	friend void swap(window& a, window& b){std::swap(a.impl_, b.impl_);}
	~window(){clear();}
	template<typename It1, typename Size, typename V = typename std::iterator_traits<It1>::value_type>
	void accumulate_n(It1 first, Size count, int target_rank, int target_disp = 0){
		using detail::data;
		int target_count = count;
		int s = MPI_Accumulate(data(first), count, detail::basic_datatype<V>{}, target_rank, target_disp, target_count, detail::basic_datatype<V>{}, MPI_SUM, impl_); 
		if(s != MPI_SUCCESS) throw std::runtime_error("cannot accumulate_n");
	}
//	void attach(void* base, MPI_Aint size){MPI_Win_attach(impl_, base, size);}
//	void call_errhandler(int errorcode);
	void complete() const{MPI_Win_complete(impl_);}
//	void create_errhandler(...);
//	void create_keyval(...);
//	void delete_attr(...);
	void fence(int assert_mode = 0 /*MPI_MODE_NOCHECK*/){
		MPI_Win_fence(assert_mode, impl_);
	}
//	void free_keyval(...);
	void flush(int rank){MPI_Win_flush(rank, impl_);}
	void flush_all(){MPI_Win_flush_all(impl_);}
	void flush(){return flush_all();}
	void flush_local(int rank){MPI_Win_flush_local(rank, impl_);}
	void flush_local_all(){MPI_Win_flush_local_all(impl_);}
	void flush_local(){return flush_local_all();}
	void* base() const{
		void* base; int flag;
		MPI3_CALL(MPI_Win_get_attr)(impl_, MPI_WIN_BASE, &base, &flag);
		assert(flag);
		return base;
	}
	mpi3::size_t const& size() const{
		MPI_Aint* size_p; int flag;
		MPI3_CALL(MPI_Win_get_attr)(impl_, MPI_WIN_SIZE, &size_p, &flag);
		assert(flag);
		return *size_p;
	}
	int const& disp_unit() const{
		int* disp_unit_p; int flag;
		MPI3_CALL(MPI_Win_get_attr)(impl_, MPI_WIN_DISP_UNIT, &disp_unit_p, &flag);
		assert(flag);
		return *disp_unit_p;
	}
//	get_errhandler(...);
	group get_group() const{
		group ret;
		MPI3_CALL(MPI_Win_get_group)(impl_, &(&ret));
		return ret;
	}
//	group get_group(){use reinterpret_cast?}
//	... get_info
//	... get_name
	// lock arguments are reversed
	void lock(int rank, int lock_type = MPI_LOCK_EXCLUSIVE, int assert = MPI_MODE_NOCHECK){
		MPI3_CALL(MPI_Win_lock)(lock_type, rank, assert, impl_);
	}
	void lock_exclusive(int rank, int assert = MPI_MODE_NOCHECK){
		MPI3_CALL(MPI_Win_lock)(MPI_LOCK_EXCLUSIVE, rank, assert, impl_);
	}
	void lock_shared(int rank, int assert = MPI_MODE_NOCHECK){
		MPI3_CALL(MPI_Win_lock)(MPI_LOCK_SHARED, rank, assert, impl_);
	}
	void lock_all(int assert = MPI_MODE_NOCHECK){MPI_Win_lock_all(assert, impl_);}
//	void post(group const& g, int assert = MPI_MODE_NOCHECK) const{
//		MPI_Win_post(g.impl_, assert, impl_);
//	}
//	void set_attr(...)
//	void set_errhandler(...)
//	void set_info(...)
//	void set_name(...)
//	void shared_query(...) delegated to child class
//	void start(group const& g, int assert = MPI_MODE_NOCHECK){
//		MPI_Win_start(g.impl_, assert, impl_);
//	}
	void sync(){MPI_Win_sync(impl_);}
//	void test(...)
	void unlock(int rank) const{MPI_Win_unlock(rank, impl_);}
	void unlock_all(){MPI_Win_unlock_all(impl_);}
	void wait() const{MPI_Win_wait(impl_);}
//	void fetch_and_op(T const*  origin, T* target, int target_rank, int target_disp = 0) const{
//		MPI_Fetch_and_op(origin, target, detail::datatype<T>{}, target_rank, target_disp, , impl_);
//	}
//	template<class T, class Op, class datatype = detail::datatype<T>, >
//	void fetch_and_op(T const*  origin, T* target, int target_rank, int target_disp = 0) const{
//		MPI_Fetch_and_op(origin, target, datatype{}, target_rank, target_disp, , impl_);
//	}
//	void fetch_exchange(T const*  origin, T* target, int target_rank, int target_disp = 0) const{
//		MPI_Fetch_and_op(origin, target,detail::datatype<T>{}, target_rank, target_disp, MPI_REPLACE, impl_);
//	}
//	maybe this goes to a pointer impl
	template<class T>
	void fetch_sum_value(T const& origin, T& target, int target_rank, int target_disp=0) const{
		MPI3_CALL(MPI_Fetch_and_op)(&origin, &target, detail::basic_datatype<T>{}, target_rank, target_disp, MPI_SUM, impl_);
	}
	template<class T>
	void fetch_prod_value(T const& origin, T& target, int target_rank, int target_disp = 0) const{
		MPI3_CALL(MPI_Fetch_and_op)(&origin, &target, detail::basic_datatype<T>{}, target_rank, target_disp, MPI_PROD, impl_);
	}
	template<class T>
	void fetch_replace_value(T const&  origin, T& target, int target_rank, int target_disp = 0) const{
		MPI3_CALL(MPI_Fetch_and_op)(&origin, &target, detail::basic_datatype<T>{}, target_rank, target_disp, MPI_REPLACE, impl_);
	}
	template<class CI1, class CI2, class datatype = detail::basic_datatype<typename std::iterator_traits<CI1>::value_type> >
	void fetch_replace(CI1 it1, CI2 it2, int target_rank, int target_disp = 0) const{
		MPI3_CALL(MPI_Fetch_and_op)(std::addressof(*it1), std::addressof(*it2), datatype{}, target_rank, target_disp, MPI_REPLACE, impl_); 
	}
	template<class ContiguousIterator>
	void blocking_put_n(ContiguousIterator it, int count, int target_rank, int target_offset = 0){
		lock_shared(target_rank, 0);
		put_n(it, count, target_rank, target_offset);
		unlock(target_rank);
	}
	template<class ContiguousIterator>
	void put_n(ContiguousIterator it, std::size_t n, int target_rank, int target_disp = 0) const{
		using detail::data;
		MPI3_CALL(MPI_Put)(
			data(it), /* void* origin_address = a + i*/ 
			n, /*int origin_count = 1 */
			detail::basic_datatype<typename std::iterator_traits<ContiguousIterator>::value_type>{}, 
			target_rank, /*int target_rank = 1*/
			target_disp, /*int target_disp = i*/
			n, /*int target_count = 1*/
			detail::basic_datatype<typename std::iterator_traits<ContiguousIterator>::value_type>{}, 
			impl_
		);
	}
	template<class ContiguousIterator>
	void put(ContiguousIterator begin, ContiguousIterator end, int target_rank, int target_disp = 0) const{
		return put_n(begin, std::distance(begin, end), target_rank, target_disp);
	}
	
	template<class Value>
	void put_value(Value const& t, int target_rank, int target_disp = 0) const{
		put_n(&t, 1, target_rank, target_disp);
	}
	template<typename ContiguousIterator, typename Size>
	ContiguousIterator get_n(ContiguousIterator it, Size n, int target_rank, int target_disp = 0) const{
		using detail::data;
		MPI3_CALL(MPI_Get)(
			data(it), /* void* origin_address = b + i*/
			n, /*int origin_count = 1 */
			detail::basic_datatype<typename std::iterator_traits<ContiguousIterator>::value_type>{}, 
			target_rank, /*int target_rank = 1 */
			target_disp, /*int target_disp = size1 + i*/
			n, /*int target_count = 1 */
			detail::basic_datatype<typename std::iterator_traits<ContiguousIterator>::value_type>{}, 
			impl_
		);
		return it + n;
	}
	template<typename ContiguousIterator>
	ContiguousIterator get(ContiguousIterator it1, ContiguousIterator it2, int target_rank, int target_disp = 0) const{
		return get_n(it1, std::distance(it1, it2), target_rank, target_disp);
	}
	template<class Value>
	void get_value(Value& t, int target_rank, int target_disp = 0) const{
		get_n(&t, 1, target_rank, target_disp);
	}
	panel<> operator[](int rank) const;
};

template<class T>
struct window : window<void>{
protected:
	window() = default;
public:
	template<class Size = mpi3::size_t>
	window(communicator const& c, T* b, Size n = 0) : window<void>{c, b, n}{}
	T* base() const{return static_cast<T*>(window<void>::base());}
	mpi3::size_t size() const{return window<void>::size()/sizeof(T);}
};

template<class T>
class panel{
	window<T>& w_;
	int rank_;
	panel(window<T>& w, int rank) : w_(w), rank_(rank){}
//	friend window;
};

template<class T> struct reference;

template<class T>
struct shm_pointer : window<>{
//	T* ptr_ = nullptr;
	T* local_ptr(int rank) const{
		mpi3::size_t size;
		int disp_unit;
		void* baseptr;
		int i = MPI_Win_shared_query(window::impl_, rank, &size, &disp_unit, &baseptr);
		if(i != MPI_SUCCESS) throw std::runtime_error{"cannot query"};
		return static_cast<T*>(baseptr);
	}
	mpi3::size_t local_size(int rank) const{
		mpi3::size_t ret = -1;
		int disp_unit = -1;
		void* baseptr = nullptr;
		int s = MPI_Win_shared_query(window::impl_, rank, &ret, &disp_unit, &baseptr);
		if(s != MPI_SUCCESS) throw std::runtime_error("cannot get local size");
		assert(ret%disp_unit == 0);
		return ret/disp_unit;
	}
	reference<T> operator*() const;
};

//template<class T> reference<T> pointer<T>::operator*() const{
//	return {*this};
//}

#if 0
template<class T>
shm_pointer<T> communicator::allocate_shared(MPI_Aint size) const
{
	shm_pointer<T> ret;
	// this assumes that the communicator is contained in a node
	int i = MPI_Win_allocate_shared(
		size*sizeof(T), sizeof(T), MPI_INFO_NULL, impl_, 
		&ret.ptr_, //&static_cast<window&>(ret).impl_
		&ret.window::impl_
	);
	if(i!=0) assert(0);
	return ret;
}
#endif 

}}

#ifdef _TEST_BOOST_MPI3_WINDOW

#include "../mpi3/main.hpp"
#include<iostream>

namespace mpi3 = boost::mpi3; 
using std::cout;

int mpi3::main(int, char*[], mpi3::communicator world){

	std::vector<double> darr(world.rank()?0:20);
	std::iota(darr.begin(), darr.end(), 0);

	mpi3::window<double> win{world, darr.data(), darr.size()};
	if(world.rank() == 0){
		std::cout << win.size() << std::endl;
		assert( win.size() == 20 );
		assert( win.base()[13] == 13 );
	}else{
		assert(win.size() == 0);
		assert(win.base() == nullptr );
	}
	win.fence();
	if(world.rank() == 0){
		std::vector<double> a = {5., 6.};
		win.put(a.begin(), a.end(), 0);
	}
	mpi3::communicator{win.get_group()}.barrier();
	win.fence();
	std::vector<double> b(2);
	win.get(b.begin(), b.end(), 0);
	win.fence();
	assert( b[0] == 5. and b[1] == 6. );
	return 0;
}

#endif
#endif

