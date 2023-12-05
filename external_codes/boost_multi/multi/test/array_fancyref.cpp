// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;autowrap:nil;-*-
// Copyright 2018-2022 Alfredo A. Correa

#define BOOST_TEST_MODULE "C++ Unit Tests for Multi constructors"
#include<boost/test/unit_test.hpp>

#include "multi/array.hpp"

namespace fancy {

template<class T> class ref;

template<class T = void> class ptr {  // NOLINT(cppcoreguidelines-special-member-functions,hicpp-special-member-functions)
	static double value;

 public:
	using difference_type = std::ptrdiff_t;
	using value_type = std::decay_t<T>;
	using pointer = T*;
	using reference = ref<T>;
	using iterator_category = std::random_access_iterator_tag;

	ptr() = default;//ptr(ptr const=default; ptr& operator=(ptr const&)=default;
	explicit ptr(std::nullptr_t) {}
	template<class Other> constexpr explicit ptr(ptr<Other> const& /*other*/) {}
	constexpr ptr(ptr const& /*other*/) {}  // NOLINT(hicpp-use-equals-default,modernize-use-equals-default)

	// vvv it is important that these two functions are device or device host functions
	// NOLINTNEXTLINE(fuchsia-overloaded-operator, fuchsia-trailing-return): this class simulates pointer
	constexpr auto operator*() const -> reference {return reference{};}
	// NOLINTNEXTLINE(fuchsia-overloaded-operator, fuchsia-trailing-return): this class simulates pointer
	constexpr auto operator+(difference_type /*unused*/) const -> ptr {return *this;}
	// NOLINTNEXTLINE(fuchsia-overloaded-operator, fuchsia-trailing-return): this class simulates pointer

	auto operator+=(difference_type /*difference*/) -> ptr& {return *this;}
	// NOLINTNEXTLINE(fuchsia-overloaded-operator, fuchsia-trailing-return): this class simulates pointer
	auto operator++() -> ptr& {return operator+=(1);}
	// NOLINTNEXTLINE(fuchsia-overloaded-operator, fuchsia-trailing-return): this class simulates pointer
	friend auto operator-(ptr const& /*a*/, ptr const& /*b*/) -> difference_type {return 0;}
	// NOLINTNEXTLINE(fuchsia-overloaded-operator, fuchsia-trailing-return): this class simulates pointer
	auto operator==(ptr const& /*other*/) const -> bool {return true;}
	// NOLINTNEXTLINE(fuchsia-overloaded-operator, fuchsia-trailing-return): this class simulates pointer
	auto operator!=(ptr const& /*other*/) const -> bool {return false;}
//	explicit operator T*() const{return &value;}
	// NOLINTNEXTLINE(fuchsia-overloaded-operator, fuchsia-trailing-return): this class simulates pointer
	auto operator->() const -> ptr const& {return *this;}
	// NOLINTNEXTLINE(fuchsia-trailing-return): this class simulates pointer
	friend auto to_address(ptr const& pointer) -> ptr {return pointer;}
	explicit operator bool() {return false;}
//	operator double*() const{return &value;}
	friend auto get_allocator(ptr const& /*self*/){return std::allocator<value_type>{};}
};

template<> double ptr<double>::value = 42.;
template<> double ptr<double const>::value = 42.;

template<class T> class ref {
	friend class ptr<T>;
	friend class ref<T const>;
	ref() = default;

 public:
//  explicit ref(ref<std::remove_const_t<T>> const& other) : p_{other.p_} {}
	~ref() = default;
	auto operator=(ref const& other) -> ref& = delete;//{
//		if(this == &other) {return *this;}
//		*p_ = *other.p_; return *this;
//	}
//  ref(ref const&) = delete;
	constexpr ref(ref const& /*other*/) = delete;
	constexpr ref(ref&& /*other*/) noexcept {}  // this is needed by nvcc, needs to be a device function for nvcc 11.2 and lower

	auto operator=(ref     && other) noexcept -> ref& = delete;// {*p_ = std::move(*other.p_); return *this;}
	constexpr operator T const&() const& {return ptr<T>::value;}  // NOLINT(google-explicit-constructor,hicpp-explicit-conversions)
	// NOLINTNEXTLINE(fuchsia-overloaded-operator): this class simulates a reference
	auto operator==(ref const& /*other*/) const {return true;}
	// NOLINTNEXTLINE(fuchsia-overloaded-operator): this class simulates a reference
	auto operator!=(ref const& /*other*/) const {return false;}
	using decay_t = std::decay_t<T>;
};

template<class T> struct allocator {
	using pointer = ptr<T>;
	using value_type = T;
	auto allocate(std::size_t /*size*/) {return pointer{};}
	void deallocate(pointer /*base*/, std::size_t /*size*/) {}
//	std::true_type operator==(allocator const&){return {};}
	allocator() = default;
	template<class T2> explicit allocator(allocator<T2> const& /*other*/) {}
	template<class... Args>
	void construct(pointer /*location*/, Args&&... /*args*/) {}
	void destroy(pointer /*location*/) {}
};

// all these are optional, depending on the level of specialization needed
template<class Ptr, class T, class Size>
auto copy_n(Ptr /*first*/, Size /*count*/, ptr<T> result) {
//	std::cerr<< "called Pointer-based copy_n(Ptr, n, fancy::ptr)" <<std::endl;
	return result;
}
template<class Ptr, class T, class Size>
auto copy_n(ptr<T> /*first*/, Size /*count*/, Ptr result) {
//	std::cerr<< "called Pointer-based copy_n(fancy::ptr, n, Ptr)" <<std::endl;
	return result;
}
template<class T1, class T2, class Size>
auto copy_n(ptr<T1> /*first*/, Size /*count*/, ptr<T2> result) {
//	std::cerr<< "called Pointer-based copy_n(fancy::ptr, n, fancy::ptr)" <<std::endl;
	return result;
}

} // namespace fancy

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//  multi-fancy glue, where should this be?
//  In boost/multi/adaptors/MyFancyApaptor.hpp if anything, or in user code if it is very specialized
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

namespace boost::multi {

template<class It, class T>
auto copy(It first, It last, fancy::ptr<T> dest) {
	return copy(first, last, multi::array_iterator<T, 1, fancy::ptr<T>>{dest});
//	std::cerr << "1D copy(it1D, it1D, it1D) with strides " << stride(first) << " " << stride(dest) << std::endl;
//	return dest;
}

template<class It, class T>  // custom copy 1D (aka strided copy)
auto copy(It/*first*/, It/*last*/, multi::array_iterator<T, 1, fancy::ptr<T>> dest) {
//	std::cerr << "1D copy(it1D, it1D, it1D) with strides " << stride(first) << " " << stride(dest) << std::endl;
	return dest;
}

template<class It, class T> // custom copy 2D (aka double strided copy)
auto copy(It/*first*/, It/*last*/, multi::array_iterator<T, 2, fancy::ptr<T>> dest) {
//	std::cerr<<"2D copy(It, It, it2D) with strides 1"<< first.stride() <<" "<< dest.stride() <<std::endl;
	return dest;
}

//template<class Alloc, class It, class T> // custom copy 2D (aka double strided copy)
//auto uninitialized_copy(Alloc&, It first, It last, multi::array_iterator<T, 2, fancy::ptr<T>> const& dest){
//	std::cerr << "2D uninitialized_copy(...) calls raw copy 2D" << std::endl;
//	return copy(first, last, dest);
//}

} // namespace boost::multi

////////////////////////////////////////////////////////////////////////////////
// user code
////////////////////////////////////////////////////////////////////////////////

BOOST_AUTO_TEST_CASE(multi_fancy) {
	namespace multi = boost::multi;

	multi::array<double, 2, fancy::allocator<double>> arr({5, 5});
	BOOST_REQUIRE( arr.size() == 5 );
	BOOST_REQUIRE( arr[1][1] == arr[2][2] );
}
