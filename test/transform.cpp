// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;autowrap:nil;-*-
// Copyright 2019-2022 Alfredo A. Correa

#define BOOST_TEST_MODULE "C++ Unit Tests for Multi transformed array"
#include<boost/test/unit_test.hpp>

#include "multi/array.hpp"

#include<complex>

namespace test {
	constexpr struct neg_t {
		template<class T>
		constexpr auto operator()(T const& value) const -> decltype(-value) {return -value;}
	} neg;
} // end namespace test

namespace test {

template<class Involution, class It> class involuter;

template<class Involution, class Ref>
class involuted {
	Ref r_;
	template<class Involuted, std::enable_if_t<std::is_base_of<involuted, std::decay_t<Involuted>>{}, int> =0>
	friend auto underlying(Involuted&& self) ->decltype(self.r_) {return self.r_;}

 public:
	using decay_type = std::decay_t<decltype(std::declval<Involution>()(std::declval<Ref>()))>;
	constexpr involuted(Involution /*stateless*/, Ref ref) : r_{std::forward<Ref>(ref)} {}
	auto operator=(decay_type const& other) -> involuted& {  // NOLINT(fuchsia-trailing-return) simulate reference
		r_ = Involution{}(other);
		return *this;
	}
	constexpr explicit operator decay_type() const {return Involution{}(r_);}
	// NOLINTNEXTLINE(google-runtime-operator): simulated reference
	constexpr auto operator&()&& {return involuter<Involution, decltype(&std::declval<Ref>())>{Involution{}, &r_};}  // NOLINT(runtime/operator)
	// NOLINTNEXTLINE(google-runtime-operator): simulated reference
	constexpr auto operator&() & {return involuter<Involution, decltype(&std::declval<Ref>())>{Involution{}, &r_};}  // NOLINT(runtime/operator)
	// NOLINTNEXTLINE(google-runtime-operator): simulated reference
	constexpr auto operator&() const& {return involuter<Involution, decltype(&std::declval<decay_type const&>())>{Involution{}, &r_};}  // NOLINT(runtime/operator)

	auto operator==(involuted  const& other) const {return r_ == other.r_;}
	auto operator!=(involuted  const& other) const {return r_ == other.r_;}

	auto operator==(decay_type const& other) const {return Involution{}(r_) == other;}
	auto operator!=(decay_type const& other) const {return Involution{}(r_) != other;}
};

template<class Involution, class It>
class involuter {
	It it_;
	template<class, class> friend class involuter;

 public:
	using pointer           = involuter<Involution, typename std::iterator_traits<It>::pointer>;
	using element_type      = typename std::pointer_traits<It>::element_type;
	using difference_type   = typename std::pointer_traits<It>::difference_type;
	template<class U>
	using rebind            = involuter<Involution, typename std::pointer_traits<It>::template rebind<U>>;
	using reference         = involuted<Involution, typename std::iterator_traits<It>::reference>;
	using value_type        = typename std::iterator_traits<It>::value_type;
	using iterator_category = typename std::iterator_traits<It>::iterator_category;

	constexpr explicit involuter(It it) : it_{std::move(it)} {}
	constexpr involuter(Involution /*stateless*/, It it) : it_{std::move(it)} {}// f_{std::move(f)}{}
	template<class Other>
	explicit involuter(involuter<Involution, Other> const& other) : it_{other.it_} {}

	constexpr auto operator*() const {return reference{Involution{}, *it_};}
	constexpr auto operator->() const {return pointer{&*it_};}

	constexpr auto operator==(involuter const& other) const {return it_ == other.it_;}
	constexpr auto operator!=(involuter const& other) const {return it_ != other.it_;}

	constexpr auto operator+=(difference_type n) -> involuter& {it_ += n; return *this;}
	constexpr auto operator-=(difference_type n) -> involuter& {it_ -= n; return *this;}

	constexpr auto operator+(difference_type n) const {return involuter{it_ + n};}
	constexpr auto operator-(difference_type n) const {return involuter{it_ - n};}
};

template<class Ref> using negated = involuted<std::negate<>, Ref>;
template<class It>  using negater = involuter<std::negate<>, It >;

class basic_conjugate_t {
	template<int N> struct prio : std::conditional_t<N!=0, prio<N-1>, std::true_type>{};
	template<class T> static auto _(prio<0>/**/, T const& value) DECLRETURN(std::conj(value))
	template<class T> static auto _(prio<1>/**/, T const& value) DECLRETURN(     conj(value))
	template<class T> static auto _(prio<2>/**/, T const& value) DECLRETURN(  T::conj(value))
	template<class T> static auto _(prio<3>/**/, T const& value) DECLRETURN(   value.conj( ))

 public:
	template<class T> static auto _(T const& value) DECLRETURN(_(prio<3>{}, value))
};

template<class T = void>
struct conjugate : private basic_conjugate_t {
	constexpr auto operator()(T const& arg) const DECLRETURN(_(arg))
};

template<>
struct conjugate<> : private basic_conjugate_t {
	template<class T>
	constexpr auto operator()(T const& arg) const DECLRETURN(_(arg))
};

#if defined(__NVCC__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wsubobject-linkage"
#endif
template<class ComplexRef> struct conjd : test::involuted<conjugate<>, ComplexRef>{
	explicit conjd(ComplexRef ref) : test::involuted<conjugate<>, ComplexRef>(conjugate<>{}, ref) {}
	auto real() const {return underlying(*this).real();}
	auto imag() const {return negated<decltype(underlying(std::declval<test::involuted<conjugate<>, ComplexRef> const&>()).imag())>{std::negate<>{}, underlying(*this).imag()};}
	friend auto real(conjd const& self) -> decltype(auto) {using std::real; return real(static_cast<typename conjd::decay_type>(self));}
	friend auto imag(conjd const& self) -> decltype(auto) {using std::imag; return imag(static_cast<typename conjd::decay_type>(self));}
};
#if defined(__NVCC__)
#pragma GCC diagnostic pop
#endif

#if defined(__cpp_deduction_guides)
template<class T> conjd(T&&)->conjd<T>;
#endif

template<class Complex> using conjr = test::involuter<conjugate<>, Complex>;

template<class P = std::complex<double>*>
class indirect_real {
	P impl_;

 public:
	explicit indirect_real(P const& ptr) : impl_{ptr} {}
	auto operator+(std::ptrdiff_t n) const {return indirect_real{impl_ + n};}
	// NOLINTNEXTLINE(cppcoreguidelines-pro-type-reinterpret-cast): extra real part as reference
	auto operator*() const -> decltype(auto) {return reinterpret_cast<std::array<double, 2>&>(*impl_)[0];}

	using difference_type = std::ptrdiff_t;
	using value_type = typename std::iterator_traits<P>::value_type;
	using pointer = void;
	using reference = void;
	using iterator_category = typename std::iterator_traits<P>::iterator_category;
};

}  // namespace test

BOOST_AUTO_TEST_CASE(transformed_array) {
	namespace multi = boost::multi;
	{
		using complex = std::complex<double>;
		complex cee{1., 2.};

		auto&& zee = test::conjd<complex&>{cee};
		BOOST_REQUIRE(( zee == complex{1., -2.} ));

		BOOST_REQUIRE( real(zee)  ==  1. );
		BOOST_REQUIRE( imag(zee)  == -2. );
		BOOST_REQUIRE( zee.real() ==  1. );
		BOOST_REQUIRE( zee.imag() == -2. );
	}
	{
		double doub = 5;

		auto&& negd_a = test::involuted<test::neg_t, double&>(test::neg, doub);
		BOOST_REQUIRE( negd_a == -5. );

		negd_a = 10.;
		BOOST_REQUIRE( negd_a == 10. );
		BOOST_REQUIRE( doub = -10. );
	}
	{
		multi::array<double, 1> arr = { 0,  1,  2,  3,  4};
		auto&& ref = arr.static_array_cast<double, double const*>();
		BOOST_REQUIRE( ref[2] == arr[2] );
	}

	{
		multi::array<double, 1> arr = { +0.0, +1.0, +2.0, +3.0, +4.0};
		multi::array<double, 1> neg = { -0.0, -1.0, -2.0, -3.0, -4.0};
		auto&& negd_arr = arr.static_array_cast<double, test::negater<double*>>();
		BOOST_REQUIRE( negd_arr[2] == neg[2] );
	}
	{
		multi::array<double, 2> arr = {
			{ +0.0,  +1.0,  +2.0,  +3.0,  +4.0},
			{ +5.0,  +6.0,  +7.0,  +8.0,  +9.0},
			{+10.0, +11.0, +12.0, +13.0, +14.0},
			{+15.0, +16.0, +17.0, +18.0, +19.0}
		};
		multi::array<double, 2> neg = {
			{ -0.0,  -1.0,  -2.0,  -3.0,  -4.0},
			{ -5.0,  -6.0,  -7.0,  -8.0,  -9.0},
			{-10.0, -11.0, -12.0, -13.0, -14.0},
			{-15.0, -16.0, -17.0, -18.0, -19.0}
		};
		auto&& negd_arr = arr.static_array_cast<double, test::negater<double*>>();
		BOOST_REQUIRE( negd_arr[1][1] == neg[1][1] );
	}
	{
	#if defined(__cpp_deduction_guides)
		double zee[4][5] {  // NOLINT(cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays) : testing legacy types
			{ 0,  1,  2,  3,  4},
			{ 5,  6,  7,  8,  9},
			{10, 11, 12, 13, 14},
			{15, 16, 17, 18, 19}
		};
		auto&& d2DC = multi::make_array_ref(test::involuter<decltype(test::neg), double*>{test::neg, &zee[0][0]}, {4, 5});

		d2DC[1][1] = -66.;
		BOOST_REQUIRE( zee[1][1] == 66 );
	#endif
		{
			using complex = std::complex<double>;
			multi::array<complex, 2> d2D = {
				{ { 0., 3.}, { 1., 9.}, { 2., 4.}, { 3., 0.}, { 4., 0.} },
				{ { 5., 0.}, { 6., 3.}, { 7., 5.}, { 8., 0.}, { 9., 0.} },
				{ { 1., 4.}, { 9., 1.}, {12., 0.}, {13., 0.}, {14., 0.} },
				{ {15., 0.}, {16., 0.}, {17., 0.}, {18., 0.}, {19., 0.} }
			};

			auto&& d2Dreal = d2D.reinterpret_array_cast<double>();
			BOOST_REQUIRE( d2Dreal[2][1] == 9. );

			d2Dreal[2][1] = 12.;
			BOOST_REQUIRE( d2D[2][1] == complex(12., 1.) );

			auto&& d2DrealT = rotated(d2D).reinterpret_array_cast<double>();
			BOOST_REQUIRE( d2DrealT[2][1] == 7. );

			multi::array<double, 2> d2Dreal_copy = d2D.template reinterpret_array_cast<double>();
			BOOST_REQUIRE( d2Dreal_copy == d2Dreal );
		}
		{
			using complex = std::complex<double>;
			constexpr auto const I = complex{0., 1.};  // NOLINT(readability-identifier-length) imaginary unit
			multi::array<complex, 2> arr = {
				{ 1. + 3.*I, 3.- 2.*I, 4.+ 1.*I},
				{ 9. + 1.*I, 7.- 8.*I, 1.- 3.*I}
			};
			auto conjd_arr = arr.static_array_cast<complex, test::conjr<complex*>>();
			BOOST_REQUIRE( conjd_arr[1][2] == conj(arr[1][2]) );
		}
	}
}
