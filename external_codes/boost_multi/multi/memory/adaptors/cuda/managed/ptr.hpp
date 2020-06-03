#ifdef COMPILATION_INSTRUCTIONS
nvcc -D_TEST_MULTI_MEMORY_ADAPTORS_CUDA_MANAGED_PTR -x c++ $0 -o $0x&&$0x&&rm $0x; exit
#endif

#ifndef BOOST_MULTI_MEMORY_ADAPTORS_CUDA_MANAGED_PTR_HPP
#define BOOST_MULTI_MEMORY_ADAPTORS_CUDA_MANAGED_PTR_HPP

#include<cstddef> // nullptr_t
#include<iterator> // random_access_iterator_tag

#include<type_traits> // is_const

#include "../../cuda/ptr.hpp"

#include "../../../../detail/memory.hpp"

#ifndef _DISABLE_CUDA_SLOW
#ifdef NDEBUG
#define SLOW deprecated("because it implies a slow element access to GPU memory")
#else
#define SLOW
#endif
#else
#define SLOW
#endif

#ifndef HD
#ifdef __CUDA_ARCH__
#define HD __host__ __device__
#else
#define HD
#endif
#endif

namespace boost{
namespace serialization{
	template<class T> class array_wrapper;
	template<class T, class S> const array_wrapper<T> make_array(T* t, S s);
}}

namespace boost{namespace multi{
namespace memory{namespace cuda{

namespace managed{

template<typename T, typename Ptr = T*> struct ptr;

template<typename RawPtr>
struct ptr<void const, RawPtr>{
	using T = void const;
	using raw_pointer = RawPtr;
	raw_pointer rp_;
	template<typename, typename> friend struct ptr;
	template<class TT> friend ptr<TT> const_pointer_cast(ptr<TT const> const&);
	ptr(raw_pointer rp) : rp_{rp}{}
public:
	ptr() = default;
	ptr(ptr const&) = default;
	ptr(std::nullptr_t n) : rp_{n}{}
	template<class Other, typename = decltype(raw_pointer{std::declval<Other const&>().rp_})>
	ptr(Other const& o) : rp_{o.rp_}{}
	ptr& operator=(ptr const&) = default;

	using pointer = ptr<T>;
	using element_type = typename std::pointer_traits<raw_pointer>::element_type;
	using difference_type = void;//typename std::pointer_traits<impl_t>::difference_type;
	explicit operator bool() const{return rp_;}
//	explicit operator raw_pointer&()&{return rp_;}
	bool operator==(ptr const& other) const{return rp_==other.rp_;}
	bool operator!=(ptr const& other) const{return rp_!=other.rp_;}
	friend ptr to_address(ptr const& p){return p;}
	void operator*() const = delete;
	template<class U> using rebind = ptr<U, typename std::pointer_traits<raw_pointer>::template rebind<U>>;
};

template<typename RawPtr>
struct ptr<void, RawPtr>{
protected:
	using T = void;
	using raw_pointer = RawPtr;
	raw_pointer rp_;
private:
	ptr(ptr<void const> const& p) : rp_{const_cast<void*>(p.rp_)}{}
	template<class TT> friend ptr<TT> const_pointer_cast(ptr<TT const> const&);
	template<class, class> friend struct ptr;
	template<class> friend class allocator;
public:
	template<class Other> ptr(ptr<Other> const& p) : rp_{p.rp_}{}
	explicit ptr(raw_pointer rp) : rp_{rp}{}
	ptr() = default;
	ptr(ptr const& p) = default;
	ptr(std::nullptr_t n) : rp_{n}{}
	template<class Other, typename = decltype(raw_pointer{std::declval<Other const&>().impl_})>
	ptr(Other const& o) : rp_{o.rp_}{}
	ptr& operator=(ptr const&) = default;
	bool operator==(ptr const& other) const{return rp_==other.rp_;}
	bool operator!=(ptr const& other) const{return rp_!=other.rp_;}
	operator cuda::ptr<void>(){return {rp_};}
	using pointer = ptr<T>;
	using element_type    = typename std::pointer_traits<raw_pointer>::element_type;
	using difference_type = typename std::pointer_traits<raw_pointer>::difference_type;
	template<class U> using rebind = ptr<U, typename std::pointer_traits<raw_pointer>::template rebind<U>>;

	explicit operator bool() const{return rp_;}
	explicit operator raw_pointer&()&{return rp_;}
	friend ptr to_address(ptr const& p){return p;}
	void operator*() = delete;
};

template<class T> class allocator;

template<typename T, typename RawPtr>
struct ptr{
	using raw_pointer = RawPtr;
protected:
	friend struct cuda::ptr<T, RawPtr>; // to allow automatic conversions
	raw_pointer rp_;
	template<class TT> friend class allocator;
	template<typename, typename> friend struct ptr;
//	template<class TT, typename = typename std::enable_if<not std::is_const<TT>{}>::type> 
//	ptr(ptr<TT const> const& p) : rp_{const_cast<T*>(p.impl_)}{}
	template<class TT> friend ptr<TT> const_pointer_cast(ptr<TT const> const&);
public:
	template<class U> using rebind = ptr<U, typename std::pointer_traits<RawPtr>::template rebind<U>>;
//	explicit ptr(cuda::ptr<T, RawPtr> const& other) : rp_{other.rp_}{}
	template<class Other, typename = std::enable_if_t<std::is_convertible<std::decay_t<decltype(std::declval<ptr<Other>>().rp_)>, raw_pointer>{}>>
	/*explicit(false)*/ ptr(ptr<Other> const& o) HD : rp_{static_cast<raw_pointer>(o.rp_)}{}
	template<class Other, typename = std::enable_if_t<not std::is_convertible<std::decay_t<decltype(std::declval<ptr<Other>>().rp_)>, raw_pointer>{}>>
	explicit/*(true)*/ ptr(ptr<Other> const& o, void** = 0) HD : rp_{static_cast<raw_pointer>(o.rp_)}{}
//	template<class Other, typename = std::enable_if_t<std::is_convertible<std::decay_t<decltype(std::declval<ptr<Other>>().rp_)>, raw_pointer>{}>>
//	ptr(ptr<Other> const& o) HD : rp_{static_cast<raw_pointer>(o.rp_)}{}
//	template<class Other, typename = std::enable_if_t<not std::is_convertible<std::decay_t<decltype(std::declval<ptr<Other>>().rp_)>, raw_pointer>{}>>
//	explicit ptr(ptr<Other> const& o, void** = 0) HD : rp_{static_cast<raw_pointer>(o.rp_)}{}
	explicit ptr(cuda::ptr<T, raw_pointer> const& other) : ptr{other.rp_}{
		assert(other.rp_!=nullptr or Cuda::pointer::attributes(other.rp_).type == cudaMemoryTypeManaged);
	}
	explicit ptr(raw_pointer p) HD : rp_{p}{}//Cuda::pointer::is_device(p);}
	ptr() = default;
	ptr(ptr const&) = default;
	ptr(std::nullptr_t n) : rp_{n}{}
	ptr& operator=(ptr const&) = default;
	bool operator==(ptr const& other) const{return rp_==other.rp_;}
	bool operator!=(ptr const& other) const{return rp_!=other.rp_;}

	using element_type = typename std::pointer_traits<raw_pointer>::element_type;
	using difference_type = typename std::pointer_traits<raw_pointer>::difference_type;
	using value_type = T;
	using pointer = ptr<T>;
	using iterator_category = typename std::iterator_traits<raw_pointer>::iterator_category; //	using iterator_concept  = typename std::iterator_traits<impl_t>::iterator_concept;
	explicit operator bool() const{return rp_;}
	bool operator not() const{return !rp_;}
	operator raw_pointer()const&{return rp_;}
	operator ptr<void>() const{return ptr<void>{rp_};}
//	template<class PM>
//	decltype(auto) operator->*(PM pm) const{return *ptr<std::decay_t<decltype(rp_->*pm)>, decltype(&(rp_->*pm))>{&(rp_->*pm)};}
	explicit operator typename std::pointer_traits<raw_pointer>::template rebind<void>() const{return typename std::pointer_traits<raw_pointer>::template rebind<void>{rp_};}
	explicit operator typename std::pointer_traits<raw_pointer>::template rebind<void const>() const{return typename std::pointer_traits<raw_pointer>::template rebind<void const>{rp_};}
	ptr& operator++(){++rp_; return *this;}
	ptr& operator--(){--rp_; return *this;}
	ptr  operator++(int){auto tmp = *this; ++(*this); return tmp;}
	ptr  operator--(int){auto tmp = *this; --(*this); return tmp;}
	ptr& operator+=(typename ptr::difference_type n) HD{rp_+=n; return *this;}
	ptr& operator-=(typename ptr::difference_type n) HD{rp_-=n; return *this;}
	ptr operator+(typename ptr::difference_type n) const HD{return ptr{rp_ + n};}
	ptr operator-(typename ptr::difference_type n) const HD{return (*this) + (-n);}
	using reference = typename std::pointer_traits<raw_pointer>::element_type&;//ref<element_type>;
//	[[SLOW]] 
//	[[deprecated]] 
	reference operator*() const HD{return *rp_;}
	HD reference operator[](difference_type n){return *((*this)+n);}
	friend inline ptr to_address(ptr const& p){return p;}
	typename ptr::difference_type operator-(ptr const& other) const{return rp_-other.rp_;}
	friend raw_pointer raw_pointer_cast(ptr const& self){return self.rp_;}
	friend cuda::ptr<T, RawPtr> cuda_pointer_cast(ptr const& self){return cuda::ptr<T, RawPtr>{self.rp_};}
	operator cuda::ptr<T, RawPtr>() const{return cuda::ptr<T, RawPtr>{rp_};}
	friend allocator<std::decay_t<T>> get_allocator(ptr const&){return {};}
	using default_allocator_type = allocator<std::decay_t<T>>;
	default_allocator_type default_allocator() const{return {};}

};

template<class T, class S> const boost::serialization::array_wrapper<T> make_array(ptr<T> t, S s){
	using boost::serialization::make_array;
	return make_array(raw_pointer_cast(t), s);
}


}

}}
}}

#undef SLOW

#ifdef _TEST_MULTI_MEMORY_ADAPTORS_CUDA_MANAGED_PTR

#include "../../cuda/managed/clib.hpp" // cuda::malloc
#include "../../cuda/managed/malloc.hpp"

#include<memory>
#include<cstring>
#include<iostream>

namespace multi = boost::multi;
namespace cuda = multi::memory::cuda;

void add_one(double& d){d += 1.;}
template<class T>
void add_one(T&& t){std::forward<T>(t) += 1.;}

// * Functions with a __global__ qualifier, which run on the device but are called by the host, cannot use pass by reference. 
//__global__ void set_5(cuda::ptr<double> const& p){
//__global__ void set_5(cuda::ptr<double> p){*p = 5.;}
//__global__ void check_5(cuda::ptr<double> p){assert(*p == 5.);}

double const* g(){double* p{nullptr}; return p;}

cuda::managed::ptr<double const> f(){
	return cuda::managed::ptr<double>{nullptr};
}

cuda::managed::ptr<double> ff(){
	return cuda::managed::ptr<double>{cuda::ptr<double>{nullptr}};
}

int main(){
	f();
	using T = double; static_assert( sizeof(cuda::managed::ptr<T>) == sizeof(T*) );
	std::size_t const n = 100;
	{
		auto p = static_cast<cuda::managed::ptr<T>>(cuda::managed::malloc(n*sizeof(T)));
		cuda::managed::ptr<void> vp = p;
		T* rp = p;
		void* vrp = p;
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
		*p = 99.; 
		if(*p != 99.) assert(0);
		if(*p == 11.) assert(0);
#pragma GCC diagnostic pop
		cuda::managed::free(p);
	}
	{
		auto p = static_cast<cuda::managed::ptr<T>>(cuda::managed::malloc(n*sizeof(T)));
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
		double* ppp = p; *ppp = 3.14;
		assert( *p == 3.14 );
#pragma GCC diagnostic pop
		cuda::managed::ptr<T> P = nullptr;
	}
	{
		cuda::managed::ptr<double> p = nullptr;
		cuda::managed::ptr<double const> pc = nullptr; 
		pc = static_cast<cuda::managed::ptr<double const>>(p);
		double* dp = cuda::managed::ptr<double>{nullptr};
		auto f = [](double const*){};
		f(p);
		cuda::ptr<double> pp = p;
//		std::reinterpret_pointer_cast<double*>(pp);
	//	cuda::managed::ptr<double> ppp{pp};
	}
	{
		auto p = static_cast<cuda::managed::ptr<T>>(cuda::managed::malloc(n*sizeof(T)));
		cuda::ptr<T> cp = p;
		cuda::managed::ptr<T> mcp{cp};
	}
	{
		static_assert(std::is_same<std::pointer_traits<cuda::managed::ptr<double>>::rebind<double const>, cuda::managed::ptr<double const>>{}, "!");
	}
	std::cout << "Finish" << std::endl;
}
#endif
#endif


