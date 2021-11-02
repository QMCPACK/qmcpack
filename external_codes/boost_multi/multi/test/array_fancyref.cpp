// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;autowrap:nil;-*-
// © Alfredo A. Correa 2018-2021

#define BOOST_TEST_MODULE "C++ Unit Tests for Multi constructors"
#define BOOST_TEST_DYN_LINK
#include<boost/test/unit_test.hpp>

#include "../array.hpp"

namespace fancy {

template<class T> class ref;

template<class T = void> class ptr{
	static double value;

 public:
	using difference_type = std::ptrdiff_t;
	using value_type = std::decay_t<T>;
	using pointer = T*;
	using reference = ref<T>;
	using iterator_category = std::random_access_iterator_tag;
	ptr() = default;//ptr(ptr const=default; ptr& operator=(ptr const&)=default;
	explicit ptr(std::nullptr_t) {}
	template<class Other> explicit ptr(ptr<Other> const& /*unused*/) {}
	// NOLINTNEXTLINE(fuchsia-overloaded-operator, fuchsia-trailing-return): this class simulates pointer
	auto operator*() const -> reference {return reference{&value};}
	// NOLINTNEXTLINE(fuchsia-overloaded-operator, fuchsia-trailing-return): this class simulates pointer
	auto operator+(difference_type /*unused*/) const -> ptr {return *this;}
	// NOLINTNEXTLINE(fuchsia-overloaded-operator, fuchsia-trailing-return): this class simulates pointer
	auto operator+=(difference_type /*unused*/) -> ptr& {return *this;}
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
	friend auto to_address(ptr const& p) -> ptr {return p;}
	explicit operator bool() {return false;}
//	operator double*() const{return &value;}
	friend auto get_allocator(ptr const& /*self*/){return std::allocator<value_type>{};}
};
template<> double ptr<double>::value = 42.;
template<> double ptr<double const>::value = 42.;

template<class T> class ref{
	T* p_;
	explicit ref(T* p) : p_{p} {}
	friend class ptr<T>;
	friend struct ref<T const>;

 public:
	explicit ref(ref<std::remove_const_t<T>> const& other) : p_{other.p_} {}
	// NOLINTNEXTLINE(fuchsia-overloaded-operator): this class simulates a reference
	auto operator==(ref const& /*other*/) const {return true;}
	// NOLINTNEXTLINE(fuchsia-overloaded-operator): this class simulates a reference
	auto operator!=(ref const& /*other*/) const {return false;}
	using decay_t = std::decay_t<T>;
};

template<class T> struct allocator{
	using pointer = ptr<T>;
	using value_type = T;
	auto allocate(std::size_t /*size*/) {return pointer{};}
	void deallocate(pointer /*base*/, std::size_t /*size*/) {}
//	std::true_type operator==(allocator const&){return {};}
	allocator() = default;
	template<class T2> explicit allocator(allocator<T2> const& /*other*/) {}
	template<class... Args>
	void construct(pointer /*location*/, Args&&... /*args*/) {}
	void destroy(pointer /*location*/){}
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

namespace boost {
namespace multi {

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

} // namespace multi
} // namespace boost

////////////////////////////////////////////////////////////////////////////////
// user code
////////////////////////////////////////////////////////////////////////////////

//namespace multi = boost::multi;

BOOST_AUTO_TEST_CASE(multi_fancy) {
//	multi::array<double, 2, fancy::allocator<double>> A({5, 5});
//	BOOST_REQUIRE( A.size() == 5 );
//	BOOST_REQUIRE( A[1][1] == A[2][2] );
}

