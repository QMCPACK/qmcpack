#if COMPILATION_INSTRUCTIONS
(echo "#include\""$0"\"" > $0x.cpp) && mpic++ -O3 -std=c++17 -Wall -Wextra -Wpedantic -D_TEST_MPI3_RMA_MEMORY $0x.cpp -o $0x.x && time mpirun -n 4 $0x.x $@ && rm -f $0x.x $0x.cpp; exit
#endif
// (C) Copyright 2019 Alfredo A. Correa
#ifndef MPI3_RMA_MEMORY_HPP
#define MPI3_RMA_MEMORY_HPP

#include "../../mpi3/window.hpp"

namespace boost{
namespace mpi3{
namespace rma{

template<class T = void> class ptr;

//template<> class ptr<void>{};

template<class T>
class ptr<T const>{
	window<T>* w_ = nullptr;
	mpi3::size_t global_offset_ = 0;
	ptr(std::nullptr_t) : w_{nullptr}, global_offset_{0}{}
	using decay_t = T;
	T operator*() const{
		T ret;
		w_->get_n(&ret, 1, global_offset_/w_->size(), global_offset_%w_->size());
		return ret;
	}
};

template<class T> class ref;

template<class T>
class ptr{
protected:
	template<class TT, class AA> friend class allocator;
	ptr(window<T>* w, mpi3::size_t global_offset = 0) : w_{w}, global_offset_{global_offset}{}
	window<T>* w_;
	mpi3::size_t global_offset_;
public:
//	int target() const{return global_offset_/w_->size();}
//	int offset() const{return global_offset_%w_->size();}
	auto address() const{
		auto d = std::div(global_offset_, w_->size());
		return mpi3::target{static_cast<int>(d.quot), d.rem};
	};
	using difference_type = mpi3::ptrdiff_t;
	ptr() = default;
	ptr(std::nullptr_t) : w_{nullptr}, global_offset_{0}{}
	ptr(ptr const&) = default;
	ref<T> operator*() const;
	ref<T> operator[](difference_type d) const{ptr tmp{*this}; tmp+=d; return *tmp;}
	ptr& operator+=(difference_type d){global_offset_ += d; return *this;}
	ptr operator+(difference_type d){ptr tmp{*this}; tmp+=d; return tmp;}
	operator void*() const{
		auto target = address();
		return w_->get_group().rank()==target.rank?(w_->base() + target.disp):nullptr;
	}
};

class mutex{
	window<void>* w_;
	mutex(mutex const&) = delete;
//	mutex& operator=(mutex const&) = deleted;
	
};

class shared_mutex{
};

template<class T>
class ref : public ptr<T>{
	ref<T>(ptr<T> const& p) : ptr<T>{p}{}
	friend class ptr<T>;
public:
	using decay_t = T;
	operator T() const&{
		T t;
		this->w_->fence(MPI_MODE_NOPUT | MPI_MODE_NOPRECEDE);
		auto target = this->address();
		this->w_->get_n(&t, 1, target.rank, target.disp);
		this->w_->fence(MPI_MODE_NOSUCCEED);
		return t;
	}
	ref<T> const& operator=(T const& t) const&{
		auto target = this->address();
	//	this->w_->fence(MPI_MODE_NOPRECEDE);
		this->w_->lock_exclusive(target.rank);
		if(this->w_->get_group().rank() == target.rank) 
			this->w_->put_n(&t, 1, target.rank, target.disp);
	//	this->w_->fence(MPI_MODE_NOSTORE | MPI_MODE_NOSUCCEED);
		this->w_->unlock(target.rank);
		return *this;
	}
//	ref<T> const& operator=(ref<T> const& other) const&{
//		if(this->w_ == other.w_) 
//	}
	template<class U> decltype(auto) 
	operator-=(U&& u) const&{T t{*this}; return operator=(t-=std::forward<U>(u));}
	ptr<T> const& operator&() const&{return *this;}
};

template<class T> ref<T> ptr<T>::operator*() const{return ref<T>{*this};}

template<class T = void, class Alloc = std::allocator<T> >
struct allocator{
	template<class U> struct rebind{typedef allocator<U> other;};
	using value_type = T;
	using pointer = rma::ptr<value_type>;
	using const_pointer = rma::ptr<value_type const>;
	using size_type = mpi3::size_t;
	using difference_type = mpi3::ptrdiff_t;
private:
	Alloc a_;
	mpi3::communicator const* c_;
	allocator() = delete;
public:
	allocator(mpi3::communicator const& c, Alloc const& a = {}) : a_{a}, c_{std::addressof(c)}{}
	template<class U> allocator(allocator<U> const& o) : a_{o.a_}, c_{o.c_}{}
	~allocator() = default;
	pointer allocate(size_type n, const void* /*hint*/ = 0){
		if(n==0) return pointer{nullptr};
		auto local_n = (n-1)/c_->size()+1;
		return {new mpi3::window<T>{*c_, a_.allocate(local_n), local_n}}; 
	}
	void deallocate(pointer p, size_type n){
		assert(n == p.w_->size());
		a_.deallocate(p.w_->base(), p.w_->size());
		delete p.w_;
	}
	allocator& operator=(allocator const& other) = default;
	bool operator==(allocator const& o) const{return a_==o.a_ and c_==o.c_;}
	bool operator!=(allocator const& o) const{return not(o == *this);}
	template<class... As>
	void construct(pointer p, As&&... as){
		void* v = static_cast<void*>(p);
		p.w_->fence();
		if(v){
			T t; a_.construct(static_cast<T*>(&t), std::forward<As>(as)...);
			p.w_->put_n(&t, 1, p.w_->get_group().rank(), p.global_offset_%p.w_->size());
		}
		p.w_->fence();
	}
	template<class P = pointer> void destroy(P p){
		void* v = static_cast<void*>(p);
		if(v) a_.destroy(v);
	}
};

}
#if 0
namespace shm{

template<class Ptr>
struct pointer_traits : std::pointer_traits<Ptr>{
	static auto to_address(Ptr const& p){
		return std::addressof(*p);
	}
};

template<class T>
struct pointer :
	std::pointer_traits<T*>,
	boost::dereferenceable<pointer<T>, T*>,
	boost::random_access_iteratable<pointer<T>, T*, std::ptrdiff_t, T&>
{
	template<class U> using rebind = pointer<U>;
//	template<class U> struct rebind{typedef pointer<U> other;};
	std::shared_ptr<mpi3::shared_window<std::decay_t<typename pointer::element_type>>> w_;
	typename pointer::difference_type offset_;
	pointer() = default;
	pointer(std::nullptr_t) : offset_(0){}
//	pointer(pointer const&) = default;
	pointer& operator=(pointer const& other) = default;
	template<class Other, typename = decltype(std::shared_ptr<mpi3::shared_window<typename pointer::element_type>>{std::declval<Other>().w_})> 
	pointer(Other&& o) : w_{o.w_}, offset_{o.offset_}{}
	pointer(pointer<std::remove_const_t<T>> const& o) : w_{o.w_}, offset_{o.offset_}{}
	pointer& operator=(std::nullptr_t){w_ = nullptr; offset_ = 0; return *this;}
	~pointer() = default;
	T& operator*() const{return *(static_cast<T*>(w_->base(0)) + offset_);}
	pointer& operator+=(typename pointer::difference_type d){offset_+=d; return *this;}
	pointer& operator++(){++offset_; return *this;}
//	pointer operator->() const{return wSP_->base(0) + offset_;}
//	reference operator[](difference_type d) const{return *((*this)+d);}
	explicit operator pointer() const{return w_->base(0) + offset_;}
	explicit operator bool() const{return bool{w_};}
	bool operator==(pointer const& o) const{assert(w_==o.w_); return offset_==o.offset_;}
	bool operator<(pointer const& o) const{assert(w_==o.w_); return offset_<o.offset_;}
	bool operator>(pointer const& o) const{assert(w_==o.w_); return offset_>o.offset_;}
	friend typename std::pointer_traits<T*>::pointer to_address(pointer const& p){
		p.w_->base(0) + p.offset_;
	}
};

template<class T> using ptr = pointer<T>;

template<class T, class F>
F for_each(pointer<T> f, pointer<T> l, F fun){ //TODO do a partitioning std::for_each
	auto& comm = f.wSP_->comm_; assert(comm == l.wSP_->comm_);
	using std::for_each;
	if(mpi3::group(*f.wSP_).root()) for_each(to_address(f), to_address(l), fun);
	f.wSP_->fence();
	f.wSP_->fence();
	return f;
}

template<class It1, typename T, typename Size>
pointer<T> copy_n(It1 f, Size n, pointer<T> d){
	d.wSP_->fence();
	using std::copy_n;
	if(mpi3::group(*d.wSP_).root()) copy_n(f, n, to_address(d));
	d.wSP_->fence();
	return d + n;
}

template<class It1, typename T>
pointer<T> copy(It1 f, It1 l, pointer<T> d){
	d.wSP_->fence();
	using std::copy;
	if(mpi3::group(*d.wSP_).root()) copy(f, l, to_address(d));
	d.wSP_->fence();
	using std::distance; using std::advance;
	advance(d, distance(f, l));
	return d;
}

template<typename T, typename Size, typename TT>
pointer<T> uninitialized_fill_n(pointer<T> f, Size n, TT const& val){
	using std::uninitialized_fill_n;
	if(mpi3::group(*f.wSP_).root()) uninitialized_fill_n(to_address(f), n, val);
	f.wSP_->fence();
	f.wSP_->fence();
	return f + n;
}

}
#endif

}}

#ifdef _TEST_MPI3_RMA_MEMORY

#include "../../mpi3/main.hpp"
#include "../../mpi3/ostream.hpp"

namespace mpi3 = boost::mpi3; 

int mpi3::main(int, char*[], mpi3::communicator world){

	mpi3::rma::allocator<int> a{world};
	auto p = a.allocate(8);
	for(int i = 0; i != 8; ++i) a.construct(&p[i], double(i));
	auto&& p2 = p[2];
	assert( p2 == 2 );
//	if(world.rank() == 1) 
	p[7] -= 7;
	world.barrier();
	assert(p[7] == 0);
	
	return 0;
}

#endif


#endif
