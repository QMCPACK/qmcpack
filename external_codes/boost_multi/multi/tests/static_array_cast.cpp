#ifdef COMPILATION// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;-*-
$CXX $0 -o $0x -lboost_unit_test_framework&&$0x&&rm $0x;exit
#endif
// © Alfredo A. Correa 2019-2020

#include "../array.hpp"

#include<complex>
#include<iostream>

#define BOOST_TEST_MODULE "C++ Unit Tests for Multi comparisons"
#define BOOST_TEST_DYN_LINK
#include<boost/test/unit_test.hpp>

using std::cout; using std::cerr;
namespace multi = boost::multi;

template<class It, class F> class involuter;

template<class Ref, class Involution>
class involuted{
protected:
	Ref r_; // [[no_unique_address]] 
	Involution f_;
public:
	using decay_type =std::decay_t<decltype(std::declval<Involution>()(std::declval<Ref>()))>;
	explicit involuted(Ref r, Involution f = {}) : r_{std::forward<Ref>(r)}, f_{f}{}
	involuted& operator=(involuted const& other)=delete;//{r_ = other.r_; return *this;}
public:
	involuted(involuted const&) = default;
public:
	involuted(involuted&&) = default; // for C++14
	operator decay_type() const&{return f_(r_);}
	decltype(auto) operator&()&&{return involuter<decltype(&std::declval<Ref>()), Involution>{&r_, f_};}
//	template<class DecayType>
//	auto operator=(DecayType&& other)&&
//	->decltype(r_=f_(std::forward<DecayType>(other)), *this){
//		return r_=f_(std::forward<DecayType>(other)), *this;}
	template<class DecayType>
	auto operator=(DecayType&& other)&
	->decltype(r_=f_(std::forward<DecayType>(other)), *this){
		return r_=f_(std::forward<DecayType>(other)), *this;}
//	template<class OtherRef>
//	auto operator=(involuted<OtherRef, Involution> const& o)&
//	->decltype(r_=f_==o.f_?std::forward<decltype(o.r_)>(o.r_):f_(o), *this){
//		return r_=f_==o.f_?std::forward<decltype(o.r_)>(o.r_):f_(o), *this;}
//	template<class DecayType>
//	auto operator==(DecayType&& other) const&
//	->decltype(this->operator decay_type()==other){
//		return this->operator decay_type()==other;}
//	template<class Any> friend auto operator<<(Any&& a, involuted const& self)->decltype(a << std::declval<decay_type>()){return a << self.operator decay_type();}
};

#if __cpp_deduction_guides
template<class T, class F> involuted(T&&, F)->involuted<T const, F>;
//template<class T, class F> involuted(T&, F)->involuted<T&, F>;
//template<class T, class F> involuted(T const&, F)->involuted<T const&, F>;
#endif

template<class It, class F>
class involuter{//: public std::iterator_traits<It>{
	It it_; // [[no_unique_address]] 
	F f_;
	template<class, class> friend class involuter;
public:
	using pointer = void;
	using element_type = typename std::pointer_traits<It>::element_type;
	using difference_type = typename std::pointer_traits<It>::difference_type;
	template <class U> using rebind = involuter<typename std::pointer_traits<It>::template rebind<U>, F>;
	using reference = involuted<typename std::iterator_traits<It>::reference, F>;
	using value_type = typename std::iterator_traits<It>::value_type;
	using iterator_category = typename std::iterator_traits<It>::iterator_category;

//	using rebind_const = involuter<typename multi::iterator_traits<It>::rebind_const, F>;
	explicit involuter(It it, F f = {}) : it_{std::move(it)}, f_{std::move(f)}{}
	involuter(involuter const& other) = default;
	template<class Other> involuter(involuter<Other, F> const& o) : it_{o.it_}, f_{o.f_}{}
//	using reference = involuted<typename std::iterator_traits<It>::reference, F>;
	auto operator*() const{return reference{*it_, f_};}
	bool operator==(involuter const& o) const{return it_==o.it_;}
	bool operator!=(involuter const& o) const{return it_!=o.it_;}
	involuter& operator+=(typename involuter::difference_type n){it_+=n; return *this;}
	auto operator+(typename involuter::difference_type n) const{return involuter{it_+n, f_};}
	decltype(auto) operator->() const{
		return involuter<typename std::iterator_traits<It>::pointer, F>{&*it_, f_};
	}
};

template<class Ref> using negated = involuted<Ref, std::negate<>>;
template<class It>  using negater = involuter<It, std::negate<>>;

BOOST_AUTO_TEST_CASE(static_array_cast){
{
	double a = 5;
	double& b = a; assert( b == 5 );
	auto&& c = involuted<double&, std::negate<>>(a);
	BOOST_REQUIRE( c == -5. );
	c = 10.;
	BOOST_REQUIRE( a = -10. );

	auto m5 = involuted<double, std::negate<>>(5.);
	BOOST_REQUIRE( m5 == -5. );
}
{
	multi::array<double, 1> A = { 0,  1,  2,  3,  4};
	auto&& A_ref = multi::static_array_cast<double, double const*>(A);
	BOOST_REQUIRE( A_ref[2] == A[2] );
	BOOST_REQUIRE( A[2] == A_ref[2] );

	BOOST_REQUIRE( std::equal(begin(A_ref), end(A_ref), begin(A)) );
	BOOST_REQUIRE( A_ref == A );
	BOOST_REQUIRE( A == A_ref );
}
{
	multi::array<double, 2> A = {
		{ 0, 1, 2, 3, 4},
		{ 5, 6, 7, 8, 9}
	};
	auto&& A_ref = multi::static_array_cast<double, double const*>(A);
	BOOST_REQUIRE( A_ref[1][1] == A[1][1] );
	BOOST_REQUIRE( std::equal(begin(A_ref[1]), end(A_ref[1]), begin(A[1])) );
	BOOST_REQUIRE( A_ref[1] == A[1] );
	BOOST_REQUIRE( std::equal(begin(A_ref), end(A_ref), begin(A)) );
	BOOST_REQUIRE( A_ref == A );
	BOOST_REQUIRE( A == A_ref );
}
{
	multi::array<double, 1> A = { 0,  1,  2,  3,  4};
	multi::array<double, 1> mA = { -0,  -1,  -2,  -3, -4};
	auto&& mA_ref = multi::static_array_cast<double, involuter<double*, std::negate<>>>(A);
	BOOST_REQUIRE( mA_ref[2] == mA[2] );
	BOOST_REQUIRE( mA[2] == mA_ref[2] );
	BOOST_REQUIRE( std::equal(begin(mA_ref), end(mA_ref), begin(mA)) );
	BOOST_REQUIRE( mA_ref == mA );
	BOOST_REQUIRE( mA == mA_ref );
}
{
	multi::array<double, 2> A = {
		{ 0,  1,  2,  3,  4}, 
		{ 5,  6,  7,  8,  9}, 
		{10, 11, 12, 13, 14}, 
		{15, 16, 17, 18, 19}
	};
	multi::array<double, 2> mA = {
		{ -0,  -1,  -2,  -3, -4}, 
		{ -5,  -6,  -7,  -8, -9}, 
		{-10, -11, -12, -13, -14}, 
		{-15, -16, -17, -18, -19}
	};
	auto&& mA_ref = multi::static_array_cast<double, negater<double*>>(A);
	BOOST_REQUIRE( mA_ref[1][1] == mA[1][1] );
	BOOST_REQUIRE( mA[1][1] == mA_ref[1][1] );
	BOOST_REQUIRE( std::equal(begin(mA[1]), end(mA[1]), begin(mA_ref[1])) );
	BOOST_REQUIRE( mA[1] == mA_ref[1] );
	BOOST_REQUIRE( mA_ref[1] == mA[1] );
	BOOST_REQUIRE( std::equal(begin(mA), end(mA), begin(mA_ref)) );	
	assert( mA_ref == mA );
	assert( mA == mA_ref );
}

}

