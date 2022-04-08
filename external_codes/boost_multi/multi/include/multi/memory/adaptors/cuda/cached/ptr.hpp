#ifdef COMPILATION// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;-*-
$CXXX $CXXFLAGS $0 -o $0.$X `pkg-config --cflags --libs cudart-11.0`&&$0.$X&&rm $0.$X;exit
#endif
// © Alfredo A. Correa 2019-2020

#ifndef BOOST_MULTI_MEMORY_ADAPTORS_CUDA_CACHED_PTR_HPP
#define BOOST_MULTI_MEMORY_ADAPTORS_CUDA_CACHED_PTR_HPP

#include<cstddef> // nullptr_t
#include<iterator> // random_access_iterator_tag

#include<type_traits> // is_const

#include "../../cuda/ptr.hpp"

#include "../../../../detail/memory.hpp"

#include<cuda_runtime.h> // cudaDeviceSynchronize

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

namespace cached{

template<typename T, typename Ptr = T*> struct ptr;

template<typename RawPtr>
struct ptr<void const, RawPtr> : cuda::ptr<void const, RawPtr> {
	using T = void const;
	using raw_pointer = RawPtr;
//	raw_pointer rp_;
	template<typename, typename> friend struct ptr;
	template<class TT> friend ptr<TT> const_pointer_cast(ptr<TT const> const&);
	explicit ptr(raw_pointer rp) : cuda::ptr<void const, RawPtr>{rp} {}

 public:
	ptr() = default;
	ptr(ptr const&) = default;

	// cppcheck-suppress noExplicitConstructor ; initialized from nullptr
	ptr(std::nullptr_t n) : cuda::ptr<void const, RawPtr>{n} {}

	template<class Other, typename = decltype(raw_pointer{std::declval<Other const&>().rp_})>
	// cppcheck-suppress noExplicitConstructor ; any pointer is convertible to void pointer
	ptr(Other const& o) : cuda::ptr<void const, RawPtr>{o.rp_} {}

	ptr& operator=(ptr const&) = default;

	using pointer = ptr<T>;
	using element_type = typename std::pointer_traits<raw_pointer>::element_type;
	using difference_type = void;//typename std::pointer_traits<impl_t>::difference_type;
//	explicit operator bool() const{return rp_;}
//	explicit operator raw_pointer&()&{return rp_;}
	friend constexpr bool operator==(ptr const& self, ptr const& other) {return self.rp_ == other.rp_;}
	friend constexpr bool operator!=(ptr const& self, ptr const& other) {return self.rp_ != other.rp_;}

	void operator*() const = delete;
	template<class U> using rebind = ptr<U, typename std::pointer_traits<raw_pointer>::template rebind<U>>;
//	friend raw_pointer raw_pointer_cast(ptr const& self) {return self.rp_;}
};

template<typename RawPtr>
struct ptr<void, RawPtr> : cuda::ptr<void, RawPtr> {
	using pointer = ptr;
	using element_type    = void;
	using difference_type = typename std::pointer_traits<RawPtr>::difference_type;

 protected:
	using raw_pointer = RawPtr;
//	raw_pointer rp_;

 private:
	ptr(ptr<void const> const& p) : cuda::ptr<void, RawPtr>{const_cast<void*>(p.rp_)} {}
	template<class TT> friend ptr<TT> const_pointer_cast(ptr<TT const> const&);
	template<class, class> friend struct ptr;
	template<class TT, class DP> friend class allocator;

 public:
	template<class Other> ptr(ptr<Other> const& p) : cuda::ptr<void, RawPtr>{p.rp_} {}
	explicit ptr(raw_pointer rp) : cuda::ptr<void, RawPtr>{rp} {}
	ptr() = default;
	ptr(ptr const& p) = default;

	// cppcheck-suppress noExplicitConstructor ; initialized from nullptr
	ptr(std::nullptr_t n) : cuda::ptr<void, RawPtr>{n} {}

	template<class Other, typename = decltype(raw_pointer{std::declval<Other const&>().impl_})>
	// cppcheck-suppress noExplicitConstructor ; any pointer is convertible to void pointer
	ptr(Other const& o) : cuda::ptr<void, RawPtr>{o.rp_}{}

	ptr& operator=(ptr const&) = default;

	friend constexpr bool operator==(ptr const& self, ptr const& other){return self.rp_==other.rp_;}
	friend constexpr bool operator!=(ptr const& self, ptr const& other){return self.rp_!=other.rp_;}

	template<class U> using rebind = ptr<U, typename std::pointer_traits<raw_pointer>::template rebind<U>>;

//	explicit operator bool() const {return this->rp_;}
	explicit operator raw_pointer&()& {return this->rp_;}

	void operator*() = delete;
	friend raw_pointer raw_pointer_cast(ptr const& self){return self.rp_;}
};

template<class T, class PrefetchDevice = std::integral_constant<int, -99> > class allocator;

template<typename T, typename RawPtr>
struct ptr : cuda::ptr<T, RawPtr> {
	using raw_pointer = RawPtr;
//	raw_pointer rp_;

 protected:
	friend struct cuda::ptr<T, RawPtr>; // to allow automatic conversions
	template<class TT, class DP> friend class allocator;
	template<typename, typename> friend struct ptr;
//	template<class TT, typename = typename std::enable_if<not std::is_const<TT>{}>::type> 
//	ptr(ptr<TT const> const& p) : rp_{const_cast<T*>(p.impl_)}{}
	template<class TT> friend ptr<TT> const_pointer_cast(ptr<TT const> const&);

 public:
	template<class U> using rebind = ptr<U, typename std::pointer_traits<RawPtr>::template rebind<U>>;
//	explicit ptr(cuda::ptr<T, RawPtr> const& other) : rp_{other.rp_}{}

	template<class Other, typename = std::enable_if_t<std::is_convertible<std::decay_t<decltype(std::declval<ptr<Other>>().rp_)>, raw_pointer>{}>>
	// cppcheck-suppress noExplicitConstructor ; propagate implicit of underlying pointer
	constexpr /*explicit(false)*/ ptr(ptr<Other> const& o) : cuda::ptr<T, RawPtr>{static_cast<raw_pointer>(o.rp_)} {}

	template<class Other, typename = std::enable_if_t<not std::is_convertible<std::decay_t<decltype(std::declval<ptr<Other>>().rp_)>, raw_pointer>{}>, typename = decltype(static_cast<raw_pointer>(std::declval<ptr<Other>>().rp_))>
	constexpr explicit/*(true)*/ ptr(ptr<Other> const& o, void** = 0) : cuda::ptr<T, RawPtr>{static_cast<raw_pointer>(o.rp_)} {}

	constexpr explicit           ptr(void* vp) : cuda::ptr<T, RawPtr>{static_cast<raw_pointer>(vp)} {}
//	template<class Other, typename = std::enable_if_t<std::is_convertible<std::decay_t<decltype(std::declval<ptr<Other>>().rp_)>, raw_pointer>{}>>
//	ptr(ptr<Other> const& o) HD : rp_{static_cast<raw_pointer>(o.rp_)}{}
//	template<class Other, typename = std::enable_if_t<not std::is_convertible<std::decay_t<decltype(std::declval<ptr<Other>>().rp_)>, raw_pointer>{}>>
//	explicit ptr(ptr<Other> const& o, void** = 0) HD : rp_{static_cast<raw_pointer>(o.rp_)}{}
	explicit ptr(cuda::ptr<T, raw_pointer> const& other) : ptr{other.rp_} {
	//	assert(other.rp_!=nullptr or Cuda::pointer::type(other.rp_) == cudaMemoryTypeCached);
	}
	constexpr explicit ptr(raw_pointer p) : cuda::ptr<T, RawPtr>{p} {}
	ptr() = default;

	// cppcheck-suppress noExplicitConstructor ; bug in cppcheck 2.3
	ptr(ptr const&) = default;

	// cppcheck-suppress noExplicitConstructor ; initialize from nullptr
	constexpr ptr(std::nullptr_t n) : cuda::ptr<T, RawPtr>{n} {}

	ptr& operator=(ptr const&) = default;
	friend constexpr bool operator==(ptr const& s, ptr const& o) {return s.rp_==o.rp_;}
	friend constexpr bool operator!=(ptr const& s, ptr const& o) {return s.rp_!=o.rp_;}

	using element_type = typename std::pointer_traits<raw_pointer>::element_type;
	using difference_type = typename std::pointer_traits<raw_pointer>::difference_type;
	using value_type = T;
	using pointer = ptr<T>;
	using iterator_category = typename std::iterator_traits<raw_pointer>::iterator_category; //	using iterator_concept  = typename std::iterator_traits<impl_t>::iterator_concept;
	explicit constexpr operator bool() const {return this->rp_;}
//	bool operator not() const{return !rp_;}
	constexpr 
#ifndef MULTI_ALLOW_IMPLICIT_CPU_CONVERSION
	explicit
#endif
	operator raw_pointer() const& {return this->rp_;} // do not =delete
	constexpr operator ptr<void>() const {return ptr<void>{this->rp_};}
//	template<class PM>
//	decltype(auto) operator->*(PM pm) const{return *ptr<std::decay_t<decltype(rp_->*pm)>, decltype(&(rp_->*pm))>{&(rp_->*pm)};}
	explicit constexpr operator typename std::pointer_traits<raw_pointer>::template rebind<void>() const{return typename std::pointer_traits<raw_pointer>::template rebind<void>{this->rp_};}
	explicit operator typename std::pointer_traits<raw_pointer>::template rebind<void const>() const{return typename std::pointer_traits<raw_pointer>::template rebind<void const>{this->rp_};}

	constexpr ptr& operator++() {++(this->rp_); return *this;} // remove
	constexpr ptr& operator--() {--(this->rp_); return *this;} // remove

	ptr  operator++(int) {auto tmp = *this; ++(*this); return tmp;} // remove
	ptr  operator--(int) {auto tmp = *this; --(*this); return tmp;} // remove

	constexpr ptr& operator+=(typename ptr::difference_type n) {(this->rp_)+=n; return *this;} // remove
	constexpr ptr& operator-=(typename ptr::difference_type n) HD {(this->rp_)-=n; return *this;} // remove

	constexpr ptr operator+(typename ptr::difference_type n) const {return ptr{(this->rp_) + n};} // remove
	constexpr ptr operator-(typename ptr::difference_type n) const {return (*this) + (-n);} // remove

	using reference = typename std::pointer_traits<raw_pointer>::element_type&;//ref<element_type>;
	constexpr reference operator*() const {return *(this->rp_);}
	constexpr reference operator[](difference_type n) const {return *(this->rp_ +n);}

	constexpr typename ptr::difference_type operator-(ptr const& other) const {return (this->rp_)-other.rp_;}
	constexpr raw_pointer raw_pointer_cast() const& {return this->rp_;} // remove
	friend raw_pointer raw_pointer_cast(ptr const& self) {return self.rp_;}
	friend cuda::ptr<T, RawPtr> cuda_pointer_cast(ptr const& self) {return cuda::ptr<T, RawPtr>{self.rp_};}
//	constexpr operator cuda::ptr<T, RawPtr>() const{return cuda::ptr<T, RawPtr>{this->rp_};}
	friend constexpr allocator<std::decay_t<T>> get_allocator(ptr const&) {return {};} // do not =delete
	using default_allocator_type = allocator<std::decay_t<T>>;
	default_allocator_type default_allocator() const {return {};}

	template<class T1, class... A1, class Size, class T2, class... A2>//, std::enable_if_t<std::is_trivially_assignable<T2&, T1>{}, int> =0>
	static auto copy_n(
		cached::ptr<T1, A1...> first, Size count,
		cached::ptr<T2, A2...> result
	) {
		return adl_copy_n(cuda::ptr<T1>(first), count, cuda::ptr<T2>(result)), result + count;
	}
public:
	friend allocator<std::decay_t<T>> default_allocator_of(ptr const&){return {};}

	template <typename ToPointer>//, typename FromElement>
	friend constexpr ToPointer
	reinterpret_pointer_cast(ptr self) {
		using to_element = typename std::pointer_traits<ToPointer>::element_type;
		return ToPointer(reinterpret_cast<to_element*>(self.raw_pointer_cast()));
	}
};

template<class T, class S> const boost::serialization::array_wrapper<T> make_array(ptr<T> t, S s) {
	using boost::serialization::make_array;
	return make_array(raw_pointer_cast(t), s);
}

}

}}
}}

#undef SLOW

#if defined(__INCLUDE_LEVEL__) and not __INCLUDE_LEVEL__

#include "../../cuda/cached/clib.hpp" // cuda::malloc
#include "../../cuda/cached/malloc.hpp"

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

cuda::cached::ptr<double const> f(){
	return cuda::cached::ptr<double>{nullptr};
}

cuda::cached::ptr<double> ff(){
	return cuda::cached::ptr<double>{cuda::ptr<double>{nullptr}};
}

std::string full_overload(double*){return "cpu";}
std::string full_overload(cuda::ptr<double>){return "gpu";}
std::string full_overload(cuda::cached::ptr<double>){return "mng";}

std::string cpugpu_overload(double*){return "cpu";}
std::string cpugpu_overload(cuda::ptr<double>){return "gpu";}

std::string cpuonly_overload(double*){return "cpu";}

std::string gpuonly_overload(cuda::ptr<double>){return "gpu";}

template<class T> void what(T&&) = delete;

int main(){


	f();
	using T = double; static_assert( sizeof(cuda::cached::ptr<T>) == sizeof(T*) , "!");
	std::size_t const n = 100;
	{
		auto p = static_cast<cuda::cached::ptr<T>>(cuda::cached::malloc(n*sizeof(T)));
	//	cuda::cached::ptr<void> vp = p;
	//	T* rp = p;
	//	void* vrp = p;
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
		*p = 99.; 
		if(*p != 99.) assert(0);
		if(*p == 11.) assert(0);
#pragma GCC diagnostic pop
		cuda::cached::free(p);		
	}
	{
		double d = 1.;
		assert( full_overload(&d) == "cpu" );
		assert( cpugpu_overload(&d) == "cpu" );
		assert( cpugpu_overload(&d) == "cpu" );

		cuda::ptr<double> p = nullptr;
		assert( full_overload(p) == "gpu" );
		assert( cpugpu_overload(p) == "gpu" );
		assert( gpuonly_overload(p) == "gpu" );

		cuda::cached::ptr<double> pm = nullptr;
		assert( full_overload(pm) == "mng" );
		assert( cpugpu_overload(pm) == "gpu" );
		assert( cpuonly_overload(pm) == "cpu" );
		assert( gpuonly_overload(pm) == "gpu" );
	}
	{
		auto p = static_cast<cuda::cached::ptr<T>>(cuda::cached::malloc(n*sizeof(T)));
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
		double* ppp = p; *ppp = 3.14;
		assert( *p == 3.14 );
#pragma GCC diagnostic pop
	//	cuda::cached::ptr<T> P = nullptr;
	}
	{
		cuda::cached::ptr<double> p = nullptr;
		cuda::cached::ptr<double const> pc = nullptr;
		assert( p == pc );
		pc = static_cast<cuda::cached::ptr<double const>>(p);
	//	double* dp = cuda::cached::ptr<double>{nullptr};
		auto f = [](double const*){};
		f(p);
	//	cuda::ptr<double> pp = p;
//		std::reinterpret_pointer_cast<double*>(pp);
	//	cuda::cached::ptr<double> ppp{pp};
	}
	{
		static_assert(std::is_convertible<cuda::cached::ptr<double>, double*>{});
	}
	{
		auto p = static_cast<cuda::cached::ptr<T>>(cuda::cached::malloc(n*sizeof(T)));
		cuda::ptr<T> cp = p;
		cuda::cached::ptr<T> mcp{cp};
	}
	{
		static_assert(std::is_same<std::pointer_traits<cuda::cached::ptr<double>>::rebind<double const>, cuda::cached::ptr<double const>>{}, "!");
	}
	std::cout << "Finish" << std::endl;
}
#endif
#endif


