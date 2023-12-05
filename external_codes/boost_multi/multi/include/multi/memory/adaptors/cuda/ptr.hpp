//#ifdef COMPILATION// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4-*-
//$CXXX $CXXFLAGS $0 -o $0.$X `pkg-config --cflags --libs cudart-11.0` -lboost_unit_test_framework&&$0.$X&&rm $0.$X;exit
//#endif

#ifndef BOOST_MULTI_MEMORY_ADAPTORS_CUDA_PTR_HPP
#define BOOST_MULTI_MEMORY_ADAPTORS_CUDA_PTR_HPP

#include "../../adaptors/cuda/clib.hpp"
#include "../../adaptors/cuda/error.hpp"
#include "../../../array_ref.hpp"
#include "../../../complex.hpp" // adl_conj

#include "../../../config/DEPRECATED.hpp"

#include<cassert> // debug
#include<utility> // exchange

#include<thrust/system/cuda/pointer.h>

#ifndef _DISABLE_CUDA_SLOW
	#ifdef NDEBUG
		#define SLOW DEPRECATED("WARNING: slow memory operation")
	#else
		#define SLOW
	#endif
#else
	#define SLOW
#endif

#define CUDA_SLOW(ExpR) NO_DEPRECATED(ExpR)

namespace boost {namespace multi {
namespace memory {namespace cuda {

template<class T> struct ref;

template<typename T, typename Ptr = T*> struct ptr;

namespace managed{template<typename T, typename RawPtr> struct ptr;}

template<typename RawPtr>
struct ptr<void const, RawPtr> {
	using pointer = ptr;
	using element_type = void const;
//	using difference_type = void;//typename std::pointer_traits<impl_t>::difference_type;

	operator ::thrust::cuda::pointer<void const>() const {return ::thrust::cuda::pointer<void const>{rp_};}

 protected:
	using raw_pointer = RawPtr;
	template<class, class> friend struct managed::ptr;
	raw_pointer rp_;
	template<typename, typename> friend struct ptr;
	template<class TT> friend ptr<TT> const_pointer_cast(ptr<TT const> const&);
	explicit ptr(raw_pointer rp) : rp_{rp} {}

 public:
	ptr() = default;
	ptr(ptr const&) = default;

	// cppcheck-suppress noExplicitConstructor ; initialize from nullptr
	ptr(std::nullptr_t n) : rp_{n} {}

	template<class Other, typename = decltype(raw_pointer{std::declval<Other const&>().rp_})>
	// cppcheck-suppress noExplicitConstructor ; any other pointer can be converted to void const pointer
	ptr(Other const& o) : rp_{o.rp_} {}
	ptr& operator=(ptr const&) = default;

	explicit operator bool() const {return rp_;}

	friend constexpr bool operator==(ptr const& s, ptr const& o) {return s.rp_==o.rp_;}
	friend constexpr bool operator!=(ptr const& s, ptr const& o) {return s.rp_!=o.rp_;}

	friend constexpr raw_pointer to_address(ptr const& self) {return self.rp_;}
	friend raw_pointer raw_pointer_cast(ptr const& self) {return self.rp_;}
};

template<class T> class allocator;

template<typename RawPtr>
struct ptr<void, RawPtr> {
	operator ::thrust::cuda_cub::pointer<void>() const {return ::thrust::cuda_cub::pointer<void>{rp_};}

 protected:
	using T = void;
	using raw_pointer = RawPtr;
	using raw_pointer_traits = std::pointer_traits<raw_pointer>;
	static_assert(std::is_same<void, typename raw_pointer_traits::element_type>{}, "!");
	raw_pointer rp_;
	friend ptr<void> malloc(size_t);
	friend void free();
	friend ptr<void> memset(ptr<void> dest, int ch, std::size_t byte_count);
	template<class, class> friend struct managed::ptr;

 protected:
	template<class TT> friend ptr<TT> const_pointer_cast(ptr<TT const> const&);
	template<class, class> friend struct ptr;
	explicit ptr(raw_pointer rp) : rp_{rp} {}
	operator raw_pointer() const{return rp_;}
	friend ptr<void> malloc(std::size_t);
	friend void free(ptr<void>);

 public:
	ptr() = default;
	ptr(ptr const& other) : rp_{other.rp_}{}//= default;

	// cppcheck-suppress noExplicitConstructor ; initialize from nullptr
	ptr(std::nullptr_t n) : rp_{n}{}

	template<class Other, typename = decltype(raw_pointer{std::declval<Other const&>().rp_})>
	// cppcheck-suppress noExplicitConstructor ; any pointer can be converted to void pointer
	ptr(Other const& o) : rp_{o.rp_}{}

	ptr& operator=(ptr const&) = default;
	friend constexpr bool operator==(ptr const& s, ptr const& o){return s.rp_==o.rp_;}
	friend constexpr bool operator!=(ptr const& s, ptr const& o){return s.rp_!=o.rp_;}
	using pointer = ptr<T>;
	using element_type    = typename std::pointer_traits<raw_pointer>::element_type;
	using difference_type = typename std::pointer_traits<raw_pointer>::difference_type;
	template<class U> using rebind = ptr<U, typename std::pointer_traits<raw_pointer>::template rebind<U>>;
//	using default_allocator_type = typename cuda::allocator<typename std::iterator_traits<raw_pointer>::value_type>;
	explicit operator bool() const{return rp_;}
//	explicit operator raw_pointer&()&{return impl_;}

	friend constexpr raw_pointer to_address(ptr const& p) {return p.rp_;}
};

template<typename T, typename RawPtr>
struct ptr {
	operator ::thrust::cuda_cub::pointer<T>() const {return ::thrust::cuda_cub::pointer<T>{rp_};}

	using raw_pointer = RawPtr;
	using default_allocator_type = typename cuda::allocator<std::decay_t<T>>;
	raw_pointer rp_ = {};

	static_assert( not std::is_same<raw_pointer, void*>{} , "!");

 protected:
	using raw_pointer_traits = typename std::pointer_traits<raw_pointer>;
	template<class TT> friend class allocator;

	template<typename, typename> friend struct ptr;
	template<typename> friend struct ref;

	template<class TT> friend ptr<TT> const_pointer_cast(ptr<TT const> const&);
	friend struct managed::ptr<T, RawPtr>;

 public:
	template<class U> using rebind = ptr<U, typename std::pointer_traits<raw_pointer>::template rebind<U>>;

	template<class Other, typename = std::enable_if_t<std::is_convertible<std::decay_t<decltype(std::declval<ptr<Other>>().rp_)>, raw_pointer>{} and not std::is_same<Other, T>{} >>
	// cppcheck-suppress noExplicitConstructor ;
	HD constexpr /*explicit(false)*/ ptr(ptr<Other> const& o) : rp_{static_cast<raw_pointer>(o.rp_)} {}
	template<class Other, typename = std::enable_if_t<not std::is_convertible<std::decay_t<decltype(std::declval<ptr<Other>>().rp_)>, raw_pointer>{} and not std::is_same<Other, T>{}>, typename = decltype(static_cast<raw_pointer>(std::declval<ptr<Other>>().rp_))>
	HD constexpr explicit/*(true)*/ ptr(ptr<Other> const& o, void** = 0) : rp_{static_cast<raw_pointer>(o.rp_)} {}
	HD constexpr explicit           ptr(raw_pointer rp) : rp_{rp} {}

	template<class TT> friend auto reinterpret_pointer_cast(ptr p)
	->decltype(ptr<TT>{reinterpret_cast<TT*>(std::declval<raw_pointer>())}){
		return ptr<TT>{reinterpret_cast<TT*>(p.rp_)};}

	template<class Other, typename = decltype(static_cast<raw_pointer>(std::declval<Other const&>().rp_))>
	HD constexpr explicit ptr(Other const& o) : rp_{static_cast<raw_pointer>(o.rp_)}{}
	ptr() = default;

	// cppcheck-suppress noExplicitConstructor ; bug in cppcheck 2.3
	ptr(ptr const&) = default;

	// cppcheck-suppress noExplicitConstructor ; initialize from nullptr
	constexpr ptr(std::nullptr_t nu) : rp_{nu} {}

	ptr& operator=(ptr const&) = default;

	friend constexpr bool operator==(ptr const& s, ptr const& o) {return s.rp_==o.rp_;}
	friend constexpr bool operator!=(ptr const& s, ptr const& o) {return s.rp_!=o.rp_;}

	using element_type    = typename raw_pointer_traits::element_type;
	using difference_type = typename raw_pointer_traits::difference_type;
	using size_type       = difference_type;
	using value_type      = T;

	using pointer = ptr<T, RawPtr>;
	using iterator_category = typename std::iterator_traits<raw_pointer>::iterator_category;
	explicit constexpr operator bool() const {return rp_;}
	explicit constexpr operator void const*() const {return rp_;}
	template<class TT=T, typename = decltype(static_cast<TT*>(raw_pointer{}))>
	explicit constexpr operator TT*() const {return static_cast<TT*>(rp_);}
	ptr& operator++() {
		static_assert(not std::is_same<raw_pointer, void*>{}, "!");
		++rp_;
		return *this;
	}
	ptr& operator--() {--rp_; return *this;}

	ptr  operator++(int) {auto tmp = *this; ++(*this); return tmp;}
	ptr  operator--(int) {auto tmp = *this; --(*this); return tmp;}

	constexpr ptr& operator+=(difference_type n) {rp_+=n; return *this;}
	constexpr ptr& operator-=(difference_type n) {rp_-=n; return *this;}

	constexpr ptr operator+(difference_type n) const {
	//	static_cast(not std::is_same<raw_pointer, void*>{} , "!");
		return ptr{rp_ + n};
	}
	constexpr ptr operator-(difference_type n) const {return ptr{rp_ - n};}

	using reference = ref<element_type>;

	[[deprecated("slow")]] constexpr auto operator*() const {return reference{*this};}
	constexpr auto operator[](difference_type n) const {return reference{*((*this)+n)};}

	constexpr difference_type operator-(ptr const& o) const {return rp_-o.rp_;}
	operator ptr<void>() {return ptr<void>{rp_};}
	HD auto get() const {return rp_;}

	friend constexpr raw_pointer to_address(ptr const& p) {return p.rp_;}  // TODO(correaa) consider returning T* , from https://en.cppreference.com/w/cpp/memory/to_address
	explicit constexpr operator raw_pointer() const {return rp_;}
	constexpr raw_pointer raw_pointer_cast() const {return this->rp_;}
	friend constexpr raw_pointer raw_pointer_cast(ptr const& self) {return self.rp_;}

	template<class PM>
	constexpr auto operator->*(PM&& pm) const
	->decltype(ref<std::decay_t<decltype(rp_->*std::forward<PM>(pm))>>{ptr<std::decay_t<decltype(rp_->*std::forward<PM>(pm))>>{&(rp_->*std::forward<PM>(pm))}}) {
		return ref<std::decay_t<decltype(rp_->*std::forward<PM>(pm))>>{ptr<std::decay_t<decltype(rp_->*std::forward<PM>(pm))>>{&(rp_->*std::forward<PM>(pm))}}; }

 public:
	friend allocator<std::decay_t<T>> get_allocator(ptr const&) {return {};}
	friend allocator<std::decay_t<T>> default_allocator_of(ptr const&) {return {};}
};

template<class T>
DEPRECATED("experimental function, it might be removed soon https://gitlab.com/correaa/boost-multi/-/issues/91")
T* raw_pointer_cast(T* p) {return p;}

template<class T> allocator<T> get_allocator(ptr<T> const&){return {};}

template<
	class InputIt, class Size, class... T, class ForwardIt = ptr<T...>,
	typename InputV = typename std::pointer_traits<InputIt>::element_type,
	typename ForwardV = typename std::pointer_traits<ForwardIt>::element_type,
	std::enable_if_t<std::is_trivially_constructible<ForwardV, InputV>{}, int> =0
>
ForwardIt uninitialized_copy_n(InputIt f, Size n, ptr<T...> d) {
	memcpy(d, f, n*sizeof(ForwardV));
	return d + n;
}

template<class It, typename Size, class T2, class Q2, typename T = typename std::iterator_traits<It>::value_type, typename = std::enable_if_t<std::is_trivially_constructible<T2, T>{}>>
auto uninitialized_copy_n(It first, Size count, boost::multi::iterator<T2, 1, ptr<Q2>> result)
//->decltype(memcpy2D(base(result), sizeof(T2)*stride(result), first, sizeof(T)*stride(first), sizeof(T), count), result + count){
{	return memcpy2D(base(result), sizeof(T2)*stride(result), base(first), sizeof(T)*stride(first), sizeof(T), count), result + count;}

template<class It, typename Size, class T2, class Q2, typename T = typename std::iterator_traits<It>::value_type, typename = std::enable_if_t<std::is_trivially_constructible<T2, T>{}>>
auto uninitialized_move_n(It first, Size count, boost::multi::iterator<T2, 1, ptr<Q2>> result)
//->decltype(memcpy2D(base(result), sizeof(T2)*stride(result), first, sizeof(T)*stride(first), sizeof(T), count), result + count){
{	return memcpy2D(base(result), sizeof(T2)*stride(result), base(first), sizeof(T)*stride(first), sizeof(T), count), result + count;}

template<class... T1, class Size, class... T2, class Element = typename std::pointer_traits<ptr<T1...>>::element_type>
auto uninitialized_move_n(ptr<T1...> first, Size n, ptr<T2...> dest) {
	assert(( std::is_trivially_constructible<Element, Element>{} ));
	memcpy(dest, first, n*sizeof(Element));
	return dest + n;
}


// copy_n

template<class Q1, class L1, class Size, class Q2, class L2>
auto copy_n(
	boost::multi::elements_iterator_t<    Q1*, L1>   first, Size count,
	boost::multi::elements_iterator_t<ptr<Q2>, L2> d_first
)-> boost::multi::elements_iterator_t<ptr<Q2>, L2> {
	copy_n(
	                                                                               first,      count,
		static_cast<boost::multi::elements_iterator_t<::thrust::cuda::pointer<Q2>, L2>>(d_first)
	);
	return d_first + count;
}

template<class Q1, class L1, class Size, class Q2, class L2>
auto copy_n(
	boost::multi::elements_iterator_t<ptr<Q1>, L1>   first, Size count,
	boost::multi::elements_iterator_t<    Q2*, L2> d_first
)-> boost::multi::elements_iterator_t<    Q2*, L2> {
	copy_n(
		static_cast<boost::multi::elements_iterator_t<::thrust::cuda::pointer<Q1>, L1>>(  first), count,
		                                                                                d_first
	);
	return d_first + count;
}

template<class Q1, class L1, class Q2, class L2>
auto copy(
	boost::multi::elements_iterator_t<    Q1*, L1>   first,
	boost::multi::elements_iterator_t<    Q1*, L1>   last ,
	boost::multi::elements_iterator_t<ptr<Q2>, L2> d_first
)-> boost::multi::elements_iterator_t<ptr<Q2>, L2> {
	return copy_n(first, last - first, d_first);
}

template<class Q1, class L1, class Q2, class L2>
auto copy(
	boost::multi::elements_iterator_t<ptr<Q1>, L1>   first,
	boost::multi::elements_iterator_t<ptr<Q1>, L1>   last ,
	boost::multi::elements_iterator_t<    Q2*, L2> d_first
)-> boost::multi::elements_iterator_t<    Q2*, L2> {
	return copy_n(first, last - first, d_first);
}


//template<
//	multi::dimensionality_type D,
//	class T1, class Q1,
//	class Size,
//	class T2, class Q2
//>
//auto copy_n(
//	boost::multi::array_iterator<T1, D, ptr<Q1>>   first , Size count,
//	boost::multi::array_iterator<T2, D,     Q2*> d_first
//)-> boost::multi::array_iterator<T2, D,     Q2*> {
//	copy_n(
//		static_cast<boost::multi::array_iterator<T1, D, ::thrust::cuda::pointer<Q1>>>(  first), count,
//		                                                                              d_first
//	);
//	return d_first + count;
//}

//template<
//	multi::dimensionality_type D,
//	class T1, class Q1,
//	class Size,
//	class T2, class Q2
//>
//auto copy_n(
//	boost::multi::array_iterator<T1, D, ptr<Q1>>   first , Size count,
//	boost::multi::array_iterator<T2, D, ptr<Q2>> d_first
//)-> boost::multi::array_iterator<T2, D, ptr<Q2>> {
//	copy_n(
//		static_cast<boost::multi::array_iterator<T1, D, ::thrust::cuda::pointer<Q1>>>(  first), count,
//		static_cast<boost::multi::array_iterator<T2, D, ::thrust::cuda::pointer<Q2>>>(d_first)
//	);
//	return d_first + count;
//}

// copy

template<
	multi::dimensionality_type D,
	class T1, class Q1,
	class T2, class Q2
>
auto copy(
	boost::multi::array_iterator<T1, D,     Q1*>   first,
	boost::multi::array_iterator<T1, D,     Q1*>   last ,
	boost::multi::array_iterator<T2, D, ptr<Q2>> d_first
)-> boost::multi::array_iterator<T2, D, ptr<Q2>> {
	copy_n(first, last - first, static_cast<boost::multi::array_iterator<T2, D, ::thrust::cuda::pointer<Q2>>>(d_first));
	return d_first + (last - first);
}

template<
	multi::dimensionality_type D,
	class T1, class Q1,
	class T2, class Q2
>
auto copy(
	boost::multi::array_iterator<T1, D, ptr<Q1>>   first,
	boost::multi::array_iterator<T1, D, ptr<Q1>>   last ,
	boost::multi::array_iterator<T2, D,     Q2*> d_first
)-> boost::multi::array_iterator<T2, D,     Q2*> {
	return copy_n(first, last - first, d_first);
}

template<
	multi::dimensionality_type D,
	class T1, class Q1,
	class T2, class Q2
>
auto copy(
	boost::multi::array_iterator<T1, D, ptr<Q1>>   first,
	boost::multi::array_iterator<T1, D, ptr<Q1>>   last ,
	boost::multi::array_iterator<T2, D, ptr<Q2>> d_first
)-> boost::multi::array_iterator<T2, D, ptr<Q2>> {
	return copy_n(first, last - first, d_first);
}

// uninitialized_copy_n

template<
	multi::dimensionality_type D,
	class T1, class Q1,
	class Size,
	class T2, class Q2
>
auto uninitialized_copy_n(
	boost::multi::array_iterator<T1, D,     Q1*>   first , Size count,
	boost::multi::array_iterator<T2, D, ptr<Q2>> d_first
)-> boost::multi::array_iterator<T2, D, ptr<Q2>> {
	uninitialized_copy_n(
		                                                                                first , count,
		static_cast<boost::multi::array_iterator<T2, D, ::thrust::cuda::pointer<Q2>>>(d_first)
	);
	return d_first + count;
}

template<
	multi::dimensionality_type D,
	class T1, class Q1,
	class Size,
	class T2, class Q2
>
auto uninitialized_copy_n(
	boost::multi::array_iterator<T1, D, ptr<Q1>>   first , Size count,
	boost::multi::array_iterator<T2, D,     Q2*> d_first
)-> boost::multi::array_iterator<T2, D,     Q2*> {
	uninitialized_copy_n(
		static_cast<boost::multi::array_iterator<T1, D, ::thrust::cuda::pointer<Q1>>>(  first), count,
		                                                                              d_first
	);
	return d_first + count;
}

template<
	multi::dimensionality_type D,
	class T1, class Q1,
	class Size,
	class T2, class Q2
>
auto uninitialized_copy_n(
	boost::multi::array_iterator<T1, D, ptr<Q1>>   first , Size count,
	boost::multi::array_iterator<T2, D, ptr<Q2>> d_first
)-> boost::multi::array_iterator<T2, D, ptr<Q2>> {
	uninitialized_copy_n(
		static_cast<boost::multi::array_iterator<T1, D, ::thrust::cuda::pointer<Q1>>>(  first), count,
		static_cast<boost::multi::array_iterator<T2, D, ::thrust::cuda::pointer<Q2>>>(d_first)
	);
	return d_first + count;
}

// uninitalized copy

template<
	multi::dimensionality_type D,
	class T1, class Q1,
	class T2, class Q2
>
auto uninitialized_copy(
	boost::multi::array_iterator<T1, D,     Q1*>   first,
	boost::multi::array_iterator<T1, D,     Q1*>   last ,
	boost::multi::array_iterator<T2, D, ptr<Q2>> d_first
)-> boost::multi::array_iterator<T2, D, ptr<Q2>> {
	return uninitialized_copy_n(first, last - first, static_cast<boost::multi::array_iterator<T2, D, ::thrust::cuda::pointer<Q2>>>(d_first));
}

template<
	multi::dimensionality_type D,
	class T1, class Q1,
	class T2, class Q2
>
auto uninitialized_copy(
	boost::multi::array_iterator<T1, D, ptr<Q1>>   first,
	boost::multi::array_iterator<T1, D, ptr<Q1>>   last ,
	boost::multi::array_iterator<T2, D,     Q2*> d_first
)-> boost::multi::array_iterator<T2, D,     Q2*> {
	return uninitialized_copy_n(first, last - first, d_first);
}

template<
	multi::dimensionality_type D,
	class T1, class Q1,
	class T2, class Q2
>
auto uninitialized_copy(
	boost::multi::array_iterator<T1, D, ptr<Q1>>   first,
	boost::multi::array_iterator<T1, D, ptr<Q1>>   last ,
	boost::multi::array_iterator<T2, D, ptr<Q2>> d_first
)-> boost::multi::array_iterator<T2, D, ptr<Q2>> {
	return uninitialized_copy_n(first, last - first, d_first);
}

template<
	class Alloc, class InputIt, class Size, class... T, class ForwardIt = ptr<T...>,
	typename InputV = typename std::pointer_traits<InputIt>::element_type, 
	typename ForwardV = typename std::pointer_traits<ForwardIt>::element_type
//	, typename = std::enable_if_t<std::is_constructible<ForwardV, InputV>{}>
>
ForwardIt alloc_uninitialized_copy_n(Alloc&, InputIt f, Size n, ptr<T...> d){ 
	if(std::is_trivially_constructible<ForwardV, InputV>{}) {
		memcpy(d, f, n*sizeof(ForwardV)); // TODO, this is not correct whe InputIt is not a pointer
		return d + n;
	} else {assert(0);}
	return d;
}

template<class Alloc, class InputIt, typename Size, class... T, class ForwardIt = ptr<T...>>
ForwardIt alloc_uninitialized_move_n(Alloc& a, InputIt f, Size n, ptr<T...> d) {
	return alloc_uninitialized_copy_n(a, f, n, d);
}

template<class T> 
ptr<T> const_pointer_cast(ptr<T const> const& p){return ptr<T>{p.impl_};}

template<class TTT>
static std::true_type is_ref_aux(ref<TTT> const&);
std::false_type is_ref_aux(...);

template<class TTT> struct is_ref : decltype(is_ref_aux(std::declval<TTT>())){};

template<class T>
struct ref {
	using value_type = T;
	using reference = value_type&;
	using pointer = ptr<T>;
	using raw_reference = value_type&;

 private:
	pointer pimpl_;
	constexpr explicit ref(pointer const& p) : pimpl_{p}{}
	template<class TT> friend struct ref;

 public:
	constexpr explicit ref(T& t) : pimpl_{&t} {}
	template<class Other, typename = decltype(multi::implicit_cast<pointer>(std::declval<ref<Other>>().pimpl_))>
	/*explicit(false)*/ constexpr ref(ref<Other>&& o) /*HD*/ : pimpl_{multi::implicit_cast<pointer>(std::move(o).pimpl_)} {}
	template<class Other, typename = std::enable_if_t<not std::is_convertible<std::decay_t<decltype(std::declval<ptr<Other>>())>, pointer>{}>>
	explicit/*(true) */ constexpr ref(ref<Other> const& o, void** = 0) /*HD*/ : pimpl_{static_cast<pointer>(o)} {}
	template<class TT, class PP> friend struct ptr;

	pointer operator&() &      __host__ __device__ {return pimpl_;}
	pointer operator&() const& __host__ __device__ {return pimpl_;}
	pointer operator&() &&     __host__ __device__ {return pimpl_;}

	struct skeleton_t {
		std::array<char, sizeof(T)> buff; T* p_;
		SLOW
		explicit skeleton_t(T* p) /*HD*/ : p_{p} {
			#if __CUDA_ARCH__
			#else
			{cudaError_t s = cudaMemcpy(buff.data(), p_, buff.size(), cudaMemcpyDeviceToHost); (void)s; assert(s == cudaSuccess);}
			#endif
		}
		operator T&()&& /*HD*/{return reinterpret_cast<T&>(buff);}
		void conditional_copyback_if_not(std::false_type) const /*HD*/{
			#if __CUDA_ARCH__
		//	*p_ = reinterpret_cast<T const&>(
			#else
			{cudaError_t s = cudaMemcpy(p_, buff.data(), buff.size(), cudaMemcpyHostToDevice); (void)s; assert(s == cudaSuccess);}
			#endif
		}
		void conditional_copyback_if_not(std::true_type) const /*HD*/{
			#if __CUDA_ARCH__
		//	*p_ = reinterpret_cast<T const&>(
			#else
		//	[[maybe_unused]] 
			cudaError_t s = cudaMemcpy(p_, buff.data(), buff.size(), cudaMemcpyHostToDevice);
			(void)s; assert(s == cudaSuccess);
			#endif
		}
		~skeleton_t() /*HD*/{conditional_copyback_if_not(std::is_const<T>{});}
	};
	skeleton_t skeleton()&& /*HD*/{return skeleton_t{raw_pointer_cast(pimpl_.rp_)};}

 public:
	constexpr ref(ref&& r) : pimpl_{std::move(r.pimpl_).rp_} {}

 private:
	ref& move_assign(ref&& other, std::true_type) & {
		cudaError_t s = cudaMemcpy(pimpl_.rp_, other.rp_, sizeof(T), cudaMemcpyDeviceToDevice); (void)s; assert(s == cudaSuccess);
		return *this;
	}
	ref& move_assign(ref&& other, std::false_type) & {
		cudaError_t s = cudaMemcpy(pimpl_.rp_, other.rp_, sizeof(T), cudaMemcpyDeviceToDevice); (void)s; assert(s == cudaSuccess);
		return *this;
	}

 public:
	template<class TT, std::enable_if_t<std::is_trivially_assignable<T&, TT>{}, int> =0>
	[[deprecated]]
	ref&& operator=(ref<TT> const& other) && {
		cudaError_t s = cudaMemcpy(pimpl_.rp_, other.pimpl_.rp_, sizeof(T), cudaMemcpyDeviceToDevice); assert(s==cudaSuccess); (void)s;
		return std::move(*this);
	}
#if __CUDA__
#ifdef __NVCC__
  #ifndef __CUDA_ARCH__
	template<class TT>
	[[deprecated]]
    __host__ auto operator=(TT&& t) &&
	->decltype(*pimpl_.rp_ = std::forward<TT>(t), std::move(*this)){
		assert(0); return std::move(*this);}
  #else
	template<class TT>
    __device__ auto operator=(TT&& t) &&
	->decltype(*pimpl_.rp_ = std::forward<TT>(t), std::move(*this)){
		return *pimpl_.rp_ = std::forward<TT>(t), std::move(*this);}
  #endif
#else
	template<class TT>
	[[deprecated("because it implies slow memory access, suround code with CUDA_SLOW")]]
	__host__ auto operator=(TT&& t) &&
	->decltype(*pimpl_.rp_ = std::forward<TT>(t), std::move(*this)){	//	assert(0);
		static_assert(std::is_trivially_assignable<T&, TT&&>{}, "!");
		cudaError_t s=cudaMemcpy(pimpl_.rp_, std::addressof(t), sizeof(T), cudaMemcpyHostToDevice);assert(s==cudaSuccess);(void)s;
		return std::move(*this);
	}
	template<class TT>
	__device__ ref&& operator=(TT&& t) &&{*pimpl_.rp_ = std::forward<TT>(t); return std::move(*this);}
#endif
#else
	template<class TT, class=std::enable_if_t<std::is_trivially_assignable<T&, TT>{}> >
	SLOW
	ref&& operator=(TT const& t) &&{
		static_assert(std::is_trivially_assignable<T&, TT&>{});
		cudaError_t s=cudaMemcpy(pimpl_.rp_, std::addressof(t), sizeof(T), cudaMemcpyHostToDevice);assert(s==cudaSuccess);(void)s;
		return std::move(*this);
	}
#endif

#if defined(__clang__)
#if defined(__CUDA__) //&& !defined(__CUDA_ARCH__)
	operator T()&& __device__{return *(pimpl_.rp_);}
	operator T()&& __host__  {static_assert( std::is_trivially_copyable<std::decay_t<T>>{}, "!" );
		typename std::aligned_storage<sizeof(T), alignof(T)>::type ret;
		{cudaError_t s=cudaMemcpy((void*)&ret, pimpl_.rp_, sizeof(T), cudaMemcpyDeviceToHost); assert(s == cudaSuccess); (void)s;}
		return *reinterpret_cast<T*>(&ret);
	}
	operator T() const& __host__  {static_assert( std::is_trivially_copyable<std::decay_t<T>>{}, "!" );
		typename std::aligned_storage<sizeof(T), alignof(T)>::type ret;
		{cudaError_t s=cudaMemcpy((void*)&ret, pimpl_.rp_, sizeof(T), cudaMemcpyDeviceToHost); assert(s == cudaSuccess); (void)s;}
		return *reinterpret_cast<T*>(&ret);
	}
#else
	SLOW operator T() && {
		std::array<char, sizeof(T)> buff;
		cudaError_t s = cudaMemcpy(buff.data(), pimpl_.rp_, buff.size(), cudaMemcpyDeviceToHost);
		switch(s) {
			case cudaSuccess                    : break;
			case cudaErrorInvalidValue          : throw std::runtime_error{"cudaErrorInvalidValue"};
			case cudaErrorInvalidMemcpyDirection: throw std::runtime_error{"cudaErrorInvalidMemcpyDirection"};
			default                             : throw std::runtime_error{"unknown error"};
		}
		return std::move(reinterpret_cast<T&>(buff));
	}
	SLOW operator T() const&{
		std::array<char, sizeof(T)> buff;  // char buff[sizeof(T)];
		cudaError_t s = cudaMemcpy(buff.data(), pimpl_.rp_, buff.size(), cudaMemcpyDeviceToHost);
		switch(s) {
			case cudaSuccess                    : break;
			case cudaErrorInvalidValue          : throw std::runtime_error{"cudaErrorInvalidValue"};
			case cudaErrorInvalidMemcpyDirection: throw std::runtime_error{"cudaErrorInvalidMemcpyDirection"};
			default                             : throw std::runtime_error{"unknown error"};
		}
		return reinterpret_cast<T&>(buff);
	}
#endif
#else // no clang
#if __CUDA_ARCH__
	operator T()&& __device__{return *(pimpl_.rp_);}
#else
	SLOW 
	operator T() && __host__ {
		std::array<char, sizeof(T)> buff;  // char buff[sizeof(T)];
		{cudaError_t s = cudaMemcpy(buff.data(), pimpl_.rp_, buff.size(), cudaMemcpyDeviceToHost); assert(s == cudaSuccess); (void)s;}
		return std::move(reinterpret_cast<T&>(buff));
	}
#endif
#if defined(__clang__)
	[[SLOW]] operator T() const& __host__{
		std::array<char, sizeof(T)> buff;  // char buff[sizeof(T)];
		{
		//	cudaError_t s = cudaMemcpy(buff, this->rp_, sizeof(T), cudaMemcpyDeviceToHost);
			auto e = static_cast<Cuda::error>(cudaMemcpy(buff.data(), pimpl_.rp_, buff.size(), cudaMemcpyDeviceToHost));
			if(e != Cuda::error::success) {throw std::system_error(e, " when trying to memcpy for element access");}
		}
		return std::move(reinterpret_cast<T&>(buff));
	}
	operator T() const& __device__{return *(pimpl_.rp_);}
#else //no clang
#if __CUDA_ARCH__
	operator T() const& __device__{return *(pimpl_.rp_);}
#else
	SLOW
	operator T() const& __host__{
		std::array<char, sizeof(T)> buff;  // char buff[sizeof(T)];
		{
			auto e = static_cast<Cuda::error>(cudaMemcpy(buff.data(), pimpl_.rp_, buff.size(), cudaMemcpyDeviceToHost));
			if(e != Cuda::error::success) {throw std::system_error(e, " when trying to memcpy for element access");}
		}
		return std::move(reinterpret_cast<T&>(buff));
	}
#endif
#endif
	#endif

#ifndef _MULTI_MEMORY_CUDA_DISABLE_ELEMENT_ACCESS
	bool operator!=(ref const& other) const&{return not(*this == other);}
	template<class Other>
	bool operator!=(ref<Other>&& other)&&{
		std::array<char, sizeof(T)> buff1;  // char buff1[sizeof(T)];
		{cudaError_t s1 = cudaMemcpy(buff1.data(), this->impl_, buff1.size(), cudaMemcpyDeviceToHost); assert(s1 == cudaSuccess); (void)s1;}
		std::array<char, sizeof(Other)> buff2;
		{cudaError_t s2 = cudaMemcpy(buff2.data(), other.impl_, buff2.size(), cudaMemcpyDeviceToHost); assert(s2 == cudaSuccess); (void)s2;}
		return reinterpret_cast<T const&>(buff1) != reinterpret_cast<Other const&>(buff2);
	}
#else
//	bool operator==(ref const& other) const = delete;
#endif
#if 1

#if defined(__clang__)
#if defined(__CUDA__) && defined(__CUDA_ARCH__)
	template<class Other, typename = std::enable_if_t<not is_ref<std::decay_t<Other>>{}> > 
	friend auto operator==(ref&& self, Other&& other) __host__{
//#if __CUDA_ARCH__
//		return std::forward<Other>(other)==*(this->rp_);
//		return *(self->rp_) == std::forward<Other>(other);
//#else
		return std::move(self).operator T() == std::forward<Other>(other);
	//	return static_cast<T>(std::move(self)) == std::forward<Other>(other);
//#endif
	}
	template<class Other, typename = std::enable_if_t<not is_ref<Other>{}> > 
	friend auto operator==(ref&& self, Other&& other) __device__ {
		return *(self->rp_) == std::forward<Other>(other);
	}
#else
	template<class Other, typename = std::enable_if_t<not is_ref<Other>{}> > 
	friend auto operator==(ref&& self, Other&& other) __host__ {
		return std::move(self).operator T() == std::forward<Other>(other);
	}
#endif
#else // no clang
	template<class Other, typename = std::enable_if_t<not is_ref<Other>{}> > 
	friend auto operator==(ref&& self, Other&& other) /*HD*/{
//	#if __CUDA_ARCH__
//		return *(self.pimpl_.rp_) == std::forward<Other>(other);
//	#else
		return static_cast<T>(std::move(self)) == std::forward<Other>(other);
//	#endif
	}
#endif

	friend __host__ __device__ decltype(auto) raw_reference_cast(ref&& r) {return *raw_pointer_cast(&r);}
	friend __host__ __device__          auto  raw_value_cast(ref&& r) {return std::move(r).operator T();}
	auto raw_value_cast() && {return std::move(*this).operator T();}

	template<class Other, typename = std::enable_if_t<not is_ref<Other>{}> >
	friend constexpr bool operator==(Other&& other, ref const& self) {
#if __CUDA_ARCH__
//		return std::forward<Other>(other)==*(this->rp_);
		return std::forward<Other>(other)==*(self.pimpl_);
#else
		return std::forward<Other>(other)== self.operator T();//static_cast<T>(std::move(self));
#endif
	}
	template<class Other, typename = std::enable_if_t<not std::is_same<T, Other>{}> >
	SLOW
	bool operator==(ref<Other>&& other) && {
		std::array<char, sizeof(T)> buff1;  // char buff1[sizeof(T)];
	//	cuda::memcpy(buff1, ref::rp_, sizeof(T));
		{cudaError_t s1 = cudaMemcpy(buff1.data(), pimpl_.rp_, buff1.size(), cudaMemcpyDeviceToHost); assert(s1 == cudaSuccess); (void)s1;}
		std::array<char, sizeof(T)> buff2;  // char buff2[sizeof(Other)];
		{cudaError_t s2 = cudaMemcpy(buff2.data(), raw_pointer_cast(&other), buff2.size(), cudaMemcpyDeviceToHost); assert(s2 == cudaSuccess); (void)s2;}
		return reinterpret_cast<T const&>(buff1) == reinterpret_cast<Other const&>(buff2);
	}
#if 1
	SLOW
	bool operator==(ref const& other) &&{
		std::array<char, sizeof(T)> buff1;  // char buff1[sizeof(T)];
		{cudaError_t s1 = cudaMemcpy(buff1.data(), pimpl_.rp_, buff1.size(), cudaMemcpyDeviceToHost); assert(s1 == cudaSuccess); (void)s1;}
		std::array<char, sizeof(T)> buff2;  // char buff2[sizeof(T)];
		{cudaError_t s2 = cudaMemcpy(buff2.data(), other.pimpl_.rp_, buff2.size(), cudaMemcpyDeviceToHost); assert(s2 == cudaSuccess); (void)s2;}
		return reinterpret_cast<T const&>(buff1) == reinterpret_cast<T const&>(buff2);
	}
#endif
#endif

#if __CUDA_ARCH__
	template<class O, typename = decltype(std::declval<T&>() += std::declval<O&&>())> __device__ ref& operator+=(O&& o) && {*(pimpl_.rp_) += std::forward<O>(o); return std::move(*this);}
	template<class O, typename = decltype(std::declval<T&>() -= std::declval<O&&>())> __device__ ref& operator-=(O&& o) && {*(pimpl_.rp_) -= std::forward<O>(o); return std::move(*this);}
#else
	template<class O, typename = decltype(std::declval<T&>() += std::declval<O&&>())> __host__ SLOW ref&& operator+=(O&& o) && {
		std::array<char, sizeof(T)> buff;
		{cudaError_t s = cudaMemcpy(buff.data(), pimpl_.rp_, buff.size(), cudaMemcpyDeviceToHost); assert(s == cudaSuccess); (void)s;}
		reinterpret_cast<T&>(buff) += std::forward<O>(o);
		{cudaError_t s = cudaMemcpy(pimpl_.rp_, buff.data(), buff.size(), cudaMemcpyHostToDevice); assert(s == cudaSuccess); (void)s;}
		return std::move(*this);
	}
	template<class O, typename = decltype(std::declval<T&>() -= std::declval<O&&>())> __host__ SLOW ref&& operator-=(O&& o) && {
		std::array<char, sizeof(T)> buff;
		{cudaError_t s = cudaMemcpy(buff.data(), pimpl_.rp_, buff.size(), cudaMemcpyDeviceToHost); assert(s == cudaSuccess); (void)s;}
		reinterpret_cast<T&>(buff) -= std::forward<O>(o);
		{cudaError_t s = cudaMemcpy(pimpl_.rp_, buff.data(), buff.size(), cudaMemcpyHostToDevice); assert(s == cudaSuccess); (void)s;}
		return std::move(*this);
	}
#endif

 private:
	template<class Ref>
	void swap(Ref&& b) &&{
		T tmp = std::move(*this);
	BEGIN_CUDA_SLOW
		*this = std::forward<Ref>(b);
		b = std::move(tmp);
	END_CUDA_SLOW
	}

 public:
	template<class Ref>
#if __NVCC__
	__attribute__((deprecated))
#else
	[[deprecated("WARNING: slow cuda memory operation")]]
#endif
	friend void swap(ref&& a, Ref&& b){a.swap(std::forward<Ref>(b));}
	template<class Ref> DEPRECATED("WARNING: slow cuda memory operation")
	friend void swap(Ref&& a, ref&& b){std::move(b).swap(a);}
	DEPRECATED("WARNING: slow cuda memory operation")
	friend void swap(ref&& a, ref&& b){std::move(a).swap(std::move(b));}
	ref<T>&& operator++()&&{++(std::move(*this).skeleton()); return std::move(*this);}
	ref<T>&& operator--()&&{--(std::move(*this).skeleton()); return std::move(*this);}
//	template<class Self, typename = std::enable_if_t<std::is_base_of<ref, Self>{}>>
//	friend auto conj(Self&& self, ref* = 0){
//		using std::conj;
//		return conj(std::forward<Self>(self).skeleton());
//	}
	template<class RRef, std::enable_if_t<std::is_same<RRef, ref>{}, int> =0>
	friend /*std::decay_t<T>*/ auto conj(RRef const& self){
		return adl_conj(self.operator T());
	}
	template<class RRef, std::enable_if_t<std::is_same<RRef, ref>{}, int> =0>
	friend auto /*std::decay_t<T>*/ imag(RRef const& self){
		return adl_imag(self.operator T());
	}
	template<class RRef, std::enable_if_t<std::is_same<RRef, ref>{}, int> =0>
	friend auto /*std::decay_t<T>*/ real(RRef const& self){
		return adl_real(self.operator T());
	}
};

}}}}

namespace thrust {
template<class T, class P> P raw_pointer_cast(boost::multi::memory::cuda::ptr<T, P> const& p) __host__ __device__ {
	return p.raw_pointer_cast();
}
}
#undef SLOW

//namespace boost {
//namespace multi {

//template<
//	class T1, class Q1,
//	class Size,
//	class T2, class P2//, class E2 = typename std::pointer_traits<P2>::element_type //, typename TP2 = decltype(ptr<E2>{std::declval<P2>()})
//> struct adl_custom_copy<
//	boost::multi::array_iterator<T1, 1L, Q1*>, boost::multi::array_iterator<T1, 1L, Q1*>,
//	boost::multi::array_iterator<T2, 1L, P2 >
//> {
//	auto copy(
//		boost::multi::array_iterator<T1, 1L, Q1*>   first , boost::multi::array_iterator<T1, 1L, Q1*> last,
//		boost::multi::array_iterator<T2, 1L, P2 > d_first
//	)-> boost::multi::array_iterator<T2, 1L, P2 > {
//		return copy_n(first, last - first, d_first);
//	}
//}

//}}

#if not __INCLUDE_LEVEL__ // def _TEST_MULTI_MEMORY_ADAPTORS_CUDA_PTR

#define BOOST_TEST_MODULE "C++ Unit Tests for Multi CUDA pointers"
#define BOOST_TEST_DYN_LINK
#include<boost/test/unit_test.hpp>

#include "../cuda/malloc.hpp"

#include "../../../adaptors/blas/numeric.hpp"

namespace multi = boost::multi;
namespace cuda = multi::memory::cuda;

template<class T> __device__ void WHAT(T&&) = delete;

#if __CUDA_ARCH__
__device__ void f(cuda::ptr<double>){
//	printf("%f", *p);
//	printf("%f", static_cast<double const&>(*p));
}
#endif

BOOST_AUTO_TEST_CASE(multi_memory_cuda_ptr){

//	static_assert( not std::is_convertible<std::complex<double>*, multi::memory::cuda::ptr<std::complex<double>>>{}, "!" );
//	static_assert( not std::is_convertible<multi::memory::cuda::ptr<std::complex<double>>, std::complex<double>*>{}, "!" );

	multi::memory::cuda::ptr<std::complex<double>> xxx = nullptr;
	std::complex<double>* ppp = raw_pointer_cast(xxx); (void)ppp;
	{
		auto ppp2 = static_cast<multi::memory::cuda::ptr<std::complex<double>>>(cuda::malloc(1*sizeof(std::complex<double>)));
		std::complex<double> const dd{*ppp2};
		assert( dd == std::complex<double>{0} );
	}
	using T = double; 
	static_assert( sizeof(cuda::ptr<T>) == sizeof(T*), "!");
	std::size_t const n = 100;
	{
		using cuda::ptr;
		auto p = static_cast<ptr<T>>(cuda::malloc(n*sizeof(T)));
CUDA_SLOW( 
		*p = 99.; 
)
		{
			ptr<T const> pc = p;
			BOOST_REQUIRE( *p == *pc );
		}
		BOOST_REQUIRE( CUDA_SLOW( *p == 99. ) );
		BOOST_REQUIRE( *p != 11. );
		cuda::free(p);

		cuda::ptr<T> P = nullptr; 
		BOOST_REQUIRE( P == nullptr );
		ptr<void> pv = p; (void)pv;
	}
//	what<typename cuda::ptr<T>::rebind<T const>>();
//	what<typename std::pointer_traits<cuda::ptr<T>>::rebind<T const>>();
	static_assert( std::is_same<typename std::pointer_traits<cuda::ptr<T>>::rebind<T const>, cuda::ptr<T const>>{} , "!");	
}

BOOST_AUTO_TEST_CASE(ptr_conversion){
	cuda::ptr<double> p = nullptr;
	cuda::ptr<double const> pc = p; (void)pc;
	static_assert(not std::is_convertible<cuda::ptr<double>, double*>{});
}

template<class T> struct Complex_{T real; T imag;};

BOOST_AUTO_TEST_CASE(multi_memory_cuda_ptr_member_pointer){
	
	Complex_<double> c{10.,20.};
//	double Complex_<double>::* 
	Complex_<double>* p = &c;
	auto pm = &Complex_<double>::imag;
	BOOST_REQUIRE( p->*pm == 20. );
	BOOST_REQUIRE( *p.*pm == 20. );

//	cuda::ptr<Complex_<double>> pcu;
//	pcu->*pm;
}



#endif
#endif

