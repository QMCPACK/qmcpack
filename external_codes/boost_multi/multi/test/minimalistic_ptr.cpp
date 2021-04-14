#ifdef COMPILATION// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4-*-
$CXX $0 -o $0x -lboost_unit_test_framework&&$0x&&rm $0x;exit
#endif
// © Alfredo A. Correa 2018-2020

#define BOOST_TEST_MODULE "C++ Unit Tests for Multi minimalistic pointer"
#define BOOST_TEST_DYN_LINK
#include<boost/test/unit_test.hpp>

#include<array>

#include "../array.hpp"

namespace multi = boost::multi;

namespace minimalistic{

template<class T> class ptr : public std::iterator_traits<T*>{ // minimalistic pointer
	T* impl_;
public:
	constexpr explicit ptr(T* impl) : impl_{impl}{}
	using typename std::iterator_traits<T*>::reference;
	using typename std::iterator_traits<T*>::difference_type;
	// NOLINTNEXTLINE(fuchsia-overloaded-operator, fuchsia-trailing-return): operator* used because this class simulates a pointer, trailing return helps
	constexpr auto operator*() const -> reference {return *impl_;} 
	// NOLINTNEXTLINE(fuchsia-overloaded-operator, cppcoreguidelines-pro-bounds-pointer-arithmetic): operator+ is overloaded to simulate a pointer
	constexpr auto operator+(difference_type n) const {return ptr{impl_ + n};}
//	T& operator[](difference_type n) const{return impl_[n];} // optional
	using default_allocator_type = std::allocator<T>;
	template<class> friend class ptr2;
};

template<class T> class ptr2 : public std::iterator_traits<T*>{ // minimalistic pointer
	T* impl_;
public:
	constexpr explicit ptr2(T* impl) : impl_{impl}{}
	constexpr explicit ptr2(ptr<T> p) : impl_{p.impl_}{} 
	using typename std::iterator_traits<T*>::reference;
	using typename std::iterator_traits<T*>::difference_type;
	// NOLINTNEXTLINE(fuchsia-overloaded-operator, fuchsia-trailing-return): operator* used because this class simulates a pointer, trailing return helps
	constexpr auto operator*() const -> reference {return *impl_;}
	// NOLINTNEXTLINE(fuchsia-overloaded-operator, cppcoreguidelines-pro-bounds-pointer-arithmetic): operator+ is overloaded to simulate a pointer
	constexpr auto operator+(typename ptr2::difference_type n) const{return ptr2{impl_ + n};}
//	T& operator[](std::ptrdiff_t n) const{return impl_[n];} // optional
	using default_allocator_type = std::allocator<T>;
};

} // namespace minimalistic

BOOST_AUTO_TEST_CASE(test_minimalistic_ptr){

	std::array<double, 100> buffer{};

	multi::array_ptr<double, 2, minimalistic::ptr<double> > CCP(minimalistic::ptr<double>{buffer.data()}, {10, 10});
	(*CCP)[2]; // requires operator+ 
	(*CCP)[1][1]; // requires operator*
	(*CCP)[1][1] = 9;
	BOOST_REQUIRE( &(*CCP)[1][1] == &buffer[11] );

	auto&& CC2 = CCP->template static_array_cast<double, minimalistic::ptr2<double>>();
	BOOST_REQUIRE( &CC2[1][1] == &(*CCP)[1][1] );

}

