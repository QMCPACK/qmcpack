#ifdef COMPILATION// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;-*-
$CXXX $CXXFLAGS $0 -o $0.$X -lboost_unit_test_framework&&$0.$X&&rm $0.$X;exit
#endif
// © Alfredo A. Correa 2019-2020

#ifndef MULTI_ADAPTORS_BLAS_NUMERIC_HPP
#define MULTI_ADAPTORS_BLAS_NUMERIC_HPP

#include "../../memory/pointer_traits.hpp"
#include "../../array_ref.hpp"
#include "../../complex.hpp"

#include "numeric/is_complex.hpp"

namespace boost{
namespace multi{
namespace blas{

template<class T> struct Complex_{T real; T imag;};

template<class A, typename T=typename std::decay_t<A>::element_type::value_type, typename C_=Complex_<T>>
auto real_aux(A&& a, std::complex<T> const&)
->decltype(std::forward<A>(a).template reinterpret_array_cast<C_>().template member_cast<T>(&C_::real)){
	return std::forward<A>(a).template reinterpret_array_cast<C_>().template member_cast<T>(&C_::real);}
	
template<class A, typename Complex, typename T=typename std::decay_t<A>::element_type::value_type, typename C_=Complex_<T>, class=std::enable_if_t<blas::numeric::is_complex_of<Complex, T>::value>>
auto real_aux(A&& a, Complex const&)
->decltype(std::forward<A>(a).template reinterpret_array_cast<C_>().template member_cast<T>(&C_::real)){
	return std::forward<A>(a).template reinterpret_array_cast<C_>().template member_cast<T>(&C_::real);}

template<class A, class T>
auto real_aux(A&& a, T const&)
->decltype(std::forward<A>(a).template member_cast<T>(&T::real)){
	return std::forward<A>(a).template member_cast<T>(&T::real);}

template<class A>
auto real(A&& a)
->decltype(real_aux(std::forward<A>(a), typename std::decay_t<A>::element_type{})){
	return real_aux(std::forward<A>(a), typename std::decay_t<A>::element_type{});}

template<class A, typename T=typename std::decay_t<A>::element_type::value_type, typename C_=Complex_<T>>
auto imag_aux(A&& a, std::complex<T> const&)
->decltype(std::forward<A>(a).template reinterpret_array_cast<C_>().template member_cast<T>(&C_::imag)){
	return std::forward<A>(a).template reinterpret_array_cast<C_>().template member_cast<T>(&C_::imag);}
	
template<class A, class Complex, typename T=typename std::decay_t<A>::element_type::value_type, typename C_=Complex_<T>, class=std::enable_if_t<blas::numeric::is_complex_of<Complex, T>::value>
>
auto imag_aux(A&& a, Complex const&)
->decltype(std::forward<A>(a).template reinterpret_array_cast<C_>().template member_cast<T>(&C_::imag)){
	return std::forward<A>(a).template reinterpret_array_cast<C_>().template member_cast<T>(&C_::imag);}

template<class A, class T>
auto imag_aux(A&& a, T const&)
->decltype(std::forward<A>(a).template member_cast<T>(&T::imag)){
	return std::forward<A>(a).template member_cast<T>(&T::imag);}

template<class A>
auto imag(A&& a)
->decltype(imag_aux(std::forward<A>(a), typename std::decay_t<A>::element_type{})){
	return imag_aux(std::forward<A>(a), typename std::decay_t<A>::element_type{});}

template<class Ref, class Involution> class involuted;

template<class It, class F, class Reference = involuted<typename std::iterator_traits<It>::reference, F> > class involuter;

template<class Ref, class Involution>
class involuted{
protected:
	Ref r_; // [[no_unique_address]] 
	Involution f_;
public:
	using decay_type =std::decay_t<decltype(std::declval<Involution>()(std::declval<Ref>()))>;
	explicit constexpr involuted(Ref r, Involution f = {}) : r_{std::forward<Ref>(r)}, f_{f}{}
	involuted& operator=(involuted const& other)=delete;//{r_ = other.r_; return *this;}
public:
	involuted(involuted const&) = delete;
	involuted(involuted&&) = default; // for C++14
	decay_type decay() const&{return f_(r_);}
	constexpr operator decay_type() const&{return f_(r_);}
	constexpr operator decay_type() &&{return f_(r_);}
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
	template<class DecayType>
	auto operator==(DecayType&& other) const
	->decltype(this->operator decay_type()==other){
		return this->operator decay_type()==other;}
	template<class DecayType>
	auto operator!=(DecayType&& other) const
	->decltype(this->operator decay_type()!=other){
		return this->operator decay_type()!=other;}
	template<class DecayType, std::enable_if_t<not std::is_base_of<involuted, DecayType>{}, int> =0>
	friend auto operator==(DecayType&& other, involuted const& self){
		return other == self.operator decay_type();}
	template<class DecayType, std::enable_if_t<not std::is_base_of<involuted, DecayType>{}, int> =0>
	friend auto operator!=(DecayType&& other, involuted const& self){
		return other != self.operator decay_type();}
//	auto imag() const{return static_cast<decay_type>(*this).imag();}
	template<class Any> friend Any& operator<<(Any&& a, involuted const& self)
//	->decltype(a << self.operator decay_type())
	{
		return a << self.operator decay_type();}
	auto conj() const&{return adl_conj(operator decay_type());}
	template<class T = void*>
	friend auto imag(involuted const& self, T = nullptr)
	->decltype(adl_imag(std::declval<decay_type>())){
		return adl_imag(self.operator decay_type());}
};

#if defined(__cpp_deduction_guides)
template<class T, class F> involuted(T&&, F)->involuted<T const, F>;
//template<class T, class F> involuted(T&, F)->involuted<T&, F>;
//template<class T, class F> involuted(T const&, F)->involuted<T const&, F>;
#endif

//template<class It, class F>
//class involuter;

template<class It, class F>
auto get_allocator(involuter<It, F> const& s);

template<class It, class F>
auto default_allocator_of(involuter<It, F> const& iv){
	return default_allocator_of(iv.it_);
}

template<class It, class F, class Reference>
class involuter{// : public std::iterator_traits<It>{
	It it_; // [[no_unique_address]] 
	F f_;
	template<class, class, class> friend class involuter;
public:
	using difference_type = typename std::iterator_traits<It>::difference_type;
	using value_type 	  = typename std::iterator_traits<It>::value_type;
	using pointer         = involuter<It, F>;//svoid; // typename std::iterator_traits<It>::pointer
	using reference 	  = Reference;
	using iterator_category = typename std::iterator_traits<It>::iterator_category;
	using element_type 	  = typename std::pointer_traits<It>::element_type;
	template<class U> using rebind = involuter<typename std::pointer_traits<It>::template rebind<U>, F>;

	involuter() = default;
	explicit involuter(It it, F f = {}) : it_{std::move(it)}, f_{std::move(f)}{}
	involuter(involuter const& other) = default;
//	template<class Other, > constexpr involuter(Other const& other) : it_{other.it_}, f_{other.f_}{}

	template<class Other, typename = decltype(_implicit_cast<It>(typename Other::underlying_type{}))> 
	constexpr involuter(Other const& o) : it_{o.it_}, f_{o.f_}{}
	template<class Other, typename = decltype(_explicit_cast<It>(typename Other::underlying_type{}))> 
	explicit constexpr involuter(Other const& o, int = 0) : it_{o.it_}, f_{o.f_}{}

	constexpr auto operator*() const {return reference{*it_, f_};}
	bool operator==(involuter const& o) const{return it_==o.it_;}
	bool operator!=(involuter const& o) const{return it_!=o.it_;}
	constexpr involuter& operator+=(typename involuter::difference_type n){it_+=n; return *this;}
	constexpr auto operator+(typename involuter::difference_type n) const{return involuter{it_+n, f_};}
//	decltype(auto) operator->() const{
//		return &const_cast<reference&>(reinterpret_cast<reference const&>(*this));
//		return reference{*it_, f_};
//		return involuter<typename std::iterator_traits<It>::pointer, F>{&*it_, f_};
//	}
	auto operator-(involuter const& other) const{return it_-other.it_;}
	explicit operator bool() const{return it_;}
	using underlying_type = It;
	friend constexpr underlying_type underlying(involuter const& self){return self.it_;}
	constexpr explicit operator It() const {return underlying(*this);}
	template<class Itt, class FF> friend auto get_allocator(involuter<Itt, FF> const&);
	friend auto default_allocator_of(involuter const& inv){
		using multi::default_allocator_of;
		return default_allocator_of(inv.it_);
	}
	using default_allocator_type = typename multi::pointer_traits<It>::default_allocator_type;
	friend auto get_allocator(involuter const& inv){
		using boost::multi::get_allocator;
		return get_allocator(inv.it_);
	}
};

template<class It, class F>
auto get_allocator(involuter<It, F> const& inv){
	using multi::get_allocator;
	return get_allocator(inv.it_);
}

template<class Ref> using negated = involuted<Ref, std::negate<>>;
template<class It>  using negater = involuter<It, std::negate<>>;

#if 1
struct conjugate{
	template<class T>
	decltype(auto) operator()(T&& a) const{
	//	using std::conj; /*for doubles?*/ 
	//	using std::conj;
	//	std::complex<double> A = static_cast<std::complex<double>>(a);
		return multi::adl_conj(std::forward<T>(a)); // this is needed by icc
	}
};
#endif

#if 0
namespace detail{
template<class Ref> struct conjugated : involuted<Ref, conjugate>{
	using involuted<Ref, conjugate>::involuted;
	template<class Other>
	conjugated(conjugated<Other> const& other) : involuted<Ref, conjugate>{static_cast<involuted<Ref, conjugate> const&>(other)}{}
	auto real() const{return static_cast<typename conjugated::decay_type>(*this).real();}
	auto imag() const{return static_cast<typename conjugated::decay_type>(*this).imag();}
	friend auto imag(conjugated const& self){return self.imag();}
	friend auto real(conjugated const& self){return self.real();}
public:
	decltype(auto) operator->() const{return this;}
//	friend auto conj(conjugated const& self){
//		return conjugate{}(static_cast<typename conjugated::decay_type>(self));
//	}
};
}
#endif

template<class Ref> using conjugated = involuted<Ref, conjugate>;

template<class It> using conjugater = involuter<It, conjugate>;//, conjugated<typename std::iterator_traits<It>::reference> >;

template<class It> auto make_conjugater(It it){return conjugater<It>{it};}
template<class It> It make_conjugater(conjugater<It> it){return underlying(it);}

template<class T> auto imag(involuted<T, conjugate> const& inv){return inv.decay().imag();}
template<class T> auto real(involuted<T, conjugate> const& inv){return inv.decay().real();}

template<class T> auto has_imag_fun_aux(T const& t)->decltype(imag(t), std::true_type {});
                  auto has_imag_fun_aux(...       )->decltype(         std::false_type{});
template<class T> struct has_imag_fun : decltype(has_imag_fun_aux(std::declval<T>())){};


template<class T> auto has_imag_mem_aux(T const& t)->decltype(t.imag(), std::true_type {});
                  auto has_imag_mem_aux(...       )->decltype(         std::false_type{});
template<class T> struct has_imag_mem : decltype(has_imag_mem_aux(std::declval<T>())){};

template<class T> struct has_imag : std::integral_constant<bool, (has_imag_fun<T>{} or has_imag_mem<T>{})>{};

template<class A = void> struct is_complex_array{
	template<class T> static auto _(T const& t) -> has_imag<T>;
	constexpr operator bool() const{return decltype(_(*base(std::declval<A>()))){};}
	template<class AA> constexpr auto operator()(AA&&){return _(*base(std::declval<A>()));}
};

template<class A = void> struct is_conjugated{
	template<class It> static std::true_type  _(conjugater<It> a);
	                   static std::false_type _(...             );
	constexpr operator bool() const{return decltype(_(base(std::declval<A>()))){};}
	template<class AA> constexpr auto operator()(AA&&){return _(base(std::declval<A>()));}
};

template<class A, class D = std::decay_t<A>, typename Elem=typename D::element_type, typename Ptr=typename D::element_ptr,
	std::enable_if_t<not is_complex_array<A>{}, int> =0>
A&& conj(A&& a){
//	return multi::static_array_cast<Elem, conjugater<Ptr>>(a);
	return std::forward<A>(a);
}

template<class A, class D = std::decay_t<A>, typename Elem=typename D::element_type, typename Ptr=typename D::element_ptr,
	std::enable_if_t<not is_conjugated<A>{} and is_complex_array<A>{}, int> =0>
decltype(auto) conj(A&& a){
//	return multi::static_array_cast<Elem, conjugater<Ptr>>(a);
	return std::forward<A>(a).template static_array_cast<Elem, conjugater<Ptr>>();
}

template<class A, class D = std::decay_t<A>, typename Elem=typename D::element_type, typename Ptr=typename D::element_ptr::underlying_type,
	std::enable_if_t<    is_conjugated<A>{}, int> =0>
auto conj(A&& a)
->decltype(std::forward<A>(a).template static_array_cast<Elem, Ptr>()){
	return std::forward<A>(a).template static_array_cast<Elem, Ptr>();}
//	return multi::static_array_cast<Elem, Ptr>(a);}
//	return multi::static_array_cast<Elem, Ptr>(a);}

}

template<class It, class F, class Reference>
auto default_allocator_of(blas::involuter<It, F, Reference> it){
	return multi::default_allocator_of(underlying(it));
}

}
}

namespace std{
	template<> struct is_convertible<boost::multi::blas::Complex_<double>*, std::complex<double>*> : std::true_type{};
	template<class T> struct is_convertible<boost::multi::blas::Complex_<double>*, T*> : boost::multi::blas::numeric::is_complex_of<T, double>{};
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

#if not __INCLUDE_LEVEL__ // _TEST_MULTI_ADAPTORS_BLAS_NUMERIC

#define BOOST_TEST_MODULE "C++ Unit Tests for Multi BLAS numeric"
#define BOOST_TEST_DYN_LINK
#include<boost/test/unit_test.hpp>

//#include "../blas/gemm.hpp"

#include "../../array.hpp"
#include "../../utility.hpp"

//#include<thrust/complex.h>

//#include "../../adaptors/cuda.hpp"

#include<cassert>
#include<iostream>

namespace multi = boost::multi;

template<class M> decltype(auto) print(M const& C){
	using std::cout;
	using boost::multi::size;
	for(int i = 0; i != size(C); ++i){
		for(int j = 0; j != size(C[i]); ++j) cout<< C[i][j] <<' ';
		cout<<std::endl;
	}
	return cout<<std::endl;
}

BOOST_AUTO_TEST_CASE(multi_blas_numeric_real_imag_part){

	using complex = std::complex<double>; complex const I{0, 1};

	multi::array<double, 2> A = {
		{1., 3., 4.}, 
		{9., 7., 1.}
	};
	multi::array<complex, 2> Acplx = A;

	multi::array<complex, 2> B = {
		{1. - 3.*I, 6. + 2.*I},
		{8. + 2.*I, 2. + 4.*I},
		{2. - 1.*I, 1. + 1.*I}
	};

	multi::array<double, 2> Breal = {
		{1., 6.},
		{8., 2.},
		{2., 1.}
	};
	multi::array<double, 2> Bimag = {
		{-3., +2.},
		{+2., +4.},
		{-1., +1.}
	};

	using multi::blas::real;
	using multi::blas::imag;

	BOOST_REQUIRE( Breal == real(B) );
	BOOST_REQUIRE( real(B) == Breal );
	BOOST_REQUIRE( imag(B) == Bimag );

	BOOST_REQUIRE( B[1][0] == 8. + 2.*I );
	BOOST_REQUIRE( B[1][0].imag() == 2. );
// 	using multi::blas::hermitized;
//	BOOST_REQUIRE( hermitized(B)[0][1] == 8. - 2.*I );
//	BOOST_REQUIRE( imag(hermitized(B)[0][1]) == -2. );
	
}

BOOST_AUTO_TEST_CASE(multi_blas_numeric_real_conjugated){

	using complex = std::complex<double>; complex const I{0, 1};

	multi::array<complex, 2> B = {
		{1. - 3.*I, 6. + 2.*I},
		{8. + 2.*I, 2. + 4.*I},
		{2. - 1.*I, 1. + 1.*I}
	};
	BOOST_REQUIRE( B[0][0] == 1. - 3.*I );

	multi::array<complex, 2> const Bconst = {
		{1. - 3.*I, 6. + 2.*I},
		{8. + 2.*I, 2. + 4.*I},
		{2. - 1.*I, 1. + 1.*I}
	};
	BOOST_REQUIRE( Bconst[0][0] == 1. - 3.*I );

	auto BdataC = multi::blas::make_conjugater(B.data_elements());
	auto BconstdataC = multi::blas::make_conjugater(Bconst.data_elements());
	decltype(BconstdataC) ppp = BdataC;
	ppp = BdataC;

	BOOST_REQUIRE( *BdataC == 1. + 3.*I );

//	static_assert(    multi::blas::is_complex_array<multi::array<thrust::complex<double>, 2>>{}, "!");
	static_assert(    multi::blas::is_complex_array<decltype(B)>{}, "!");
	static_assert(not multi::blas::is_conjugated<decltype(B)>{}, "!");

	auto&& Bconj = multi::blas::conj(B);
	static_assert(multi::blas::is_conjugated<decltype(Bconj)>{}, "!");

	BOOST_REQUIRE( Bconj[0][0] == 1. + 3.*I );
	BOOST_TEST_REQUIRE( imag(*base(Bconj)) == +3 );

//	BOOST_TEST_REQUIRE( base(Bconj)->imag() == +3 );
	BOOST_REQUIRE( rotated(Bconj)[1][0] == Bconj[0][1] );

//	BOOST_REQUIRE( base(Bconj) == -3.*I );
	static_assert(multi::blas::is_complex_array<decltype(Bconj)>{}, "!");

	BOOST_REQUIRE( conj(Bconj) == B );
#if 0
	BOOST_REQUIRE( base(conj(Bconj)) == base(B) );
	BOOST_REQUIRE( base(conj(Bconj))->imag() == -3. );
//	BOOST_REQUIRE( base(conjugated(Bconj))->imag() == -3. );
#endif

	BOOST_REQUIRE( multi::blas::conj(B)[1][0] == std::conj(B[1][0]) );
}

#if 0

	namespace cuda = multi::cuda;
	{
		cuda::array<complex, 2> Bgpu = B;
		using multi::blas::imag;
		BOOST_REQUIRE( imag(Bgpu)[1][1] == imag(B)[1][1] );
		BOOST_REQUIRE( real(Bgpu)[1][1] == real(B)[1][1] );
	}
	{
		cuda::managed::array<complex, 2> Bgpu = B;
		using multi::blas::imag;
		BOOST_REQUIRE( imag(Bgpu)[1][1] == imag(B)[1][1] );
		BOOST_REQUIRE( real(Bgpu)[1][1] == real(B)[1][1] );
	}

	multi::array_ref<double, 2> rB(reinterpret_cast<double*>(data_elements(B)), {size(B), 2*size(*begin(B))});

	auto&& Bconj = multi::static_array_cast<complex, multi::blas::detail::conjugater<complex*>>(B);
	assert( size(Bconj) == size(B) );
	assert( conj(B[1][2]) == Bconj[1][2] );

//	auto&& BH = multi::blas::hermitized(B);
//	assert( BH[1][2] == conj(B[2][1]) );
//	std::cout << BH[1][2] << " " << B[2][1] << std::endl;

//	auto&& BH1 = multi::static_array_cast<complex, multi::blas::detail::conjugater<complex*>>(rotated(B));
//	auto&& BH2 = rotated(multi::static_array_cast<complex, multi::blas::detail::conjugater<complex*>>(B));

//	what( BH1, BH2 );
//	using multi::blas::imag;

//	assert( real(A)[1][2] == 1. );
//	assert( imag(A)[1][2] == -3. );

//	print(A) <<"--\n";
//	print(real(A)) <<"--\n";
//	print(imag(A)) <<"--\n";

	multi::array<complex, 2> C({2, 2});
	multi::array_ref<double, 2> rC(reinterpret_cast<double*>(data_elements(C)), {size(C), 2*size(*begin(C))});

//	gemm('T', 'T', 1., A, B, 0., C);
//	gemm('T', 'T', 1., A, B, 0., C);
//	gemm('T', 'T', 1., real(A), B, 0., C);
}
#endif
#endif
#endif

