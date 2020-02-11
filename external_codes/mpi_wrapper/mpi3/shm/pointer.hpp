#if COMPILATION_INSTRUCTIONS
(echo "#include\""$0"\"" > $0.cpp) && mpic++ -O3 -std=c++14 -Wall `#-Wfatal-errors` -D_TEST_MPI3_SHM_POINTER $0.cpp -o $0x && time mpirun -n 3 $0x $@ && rm $0x $0.cpp; exit
#endif
// (C) Copyright 2019 Alfredo A. Correa
#ifndef MPI3_SHM_POINTER_HPP
#define MPI3_SHM_POINTER_HPP

#include "../../mpi3/shared_communicator.hpp"
#include "../../mpi3/shared_window.hpp"

#include<boost/operators.hpp> // dereferenceable, random_access_iteratable

namespace boost{
namespace mpi3{
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
//	explicit operator pointer() const{return w_->base(0) + offset_;}
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

}}}

#ifdef _TEST_MPI3_SHM_POINTER

#include "../../mpi3/main.hpp"
#include "../../mpi3/ostream.hpp"

namespace mpi3 = boost::mpi3; 

int mpi3::main(int, char*[], mpi3::communicator world){
	using p = mpi3::shm::pointer<double>;
	using cp = std::pointer_traits<p>::template rebind<double const>;//::other;
//	whatis<cp>();
	p P;
	cp CP = P;
	return 0;
}

#endif


#endif
