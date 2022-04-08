// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;autowrap:nil;-*-
// © Alfredo Correa 2019-2021

#define BOOST_TEST_MODULE "C++ Unit Tests for Multi transformed array"
#define BOOST_TEST_DYN_LINK
#include<boost/test/unit_test.hpp>

#include "../array.hpp"

//#include<boost/iterator/transform_iterator.hpp> //might need -DBOOST_RESULT_OF_USE_DECLTYPE

#include<complex>

namespace test {
	constexpr struct neg_t {
		template<class T>
		constexpr auto operator()(T const& x) const -> decltype(-x) {return -x;}
	} neg;
} // end namespace test

namespace test {

template<class Involution, class It> class involuter;

template<class Involution, class Ref>
class involuted{
	Ref r_;
	template<class Involuted, std::enable_if_t<std::is_base_of<involuted, std::decay_t<Involuted>>{}, int> =0>
	friend auto underlying(Involuted&& self) ->decltype(self.r_) {return self.r_;}

 public:
	using decay_type = std::decay_t<decltype(std::declval<Involution>()(std::declval<Ref>()))>;
	constexpr involuted(Involution /*stateless*/, Ref r) : r_{std::forward<Ref>(r)} {}
#if __cplusplus <= 201402L
	involuted(involuted&&) noexcept = default;
	~involuted() = default;
	involuted(involuted const&) = default;
	auto operator=(involuted const&) -> involuted& = default; // NOLINT(fuchsia-trailing-return): simulate reference
	auto operator=(involuted&&     ) noexcept -> involuted& = default; // NOLINT(fuchsia-trailing-return): simulate reference
#endif
	auto operator=(decay_type const& other) -> involuted& { //  NOLINT(fuchsia-trailing-return): simulate reference
		r_ = Involution{}(other);
		return *this;
	}
	constexpr explicit operator decay_type() const{return Involution{}(r_);}
	// NOLINTNEXTLINE(google-runtime-operator): simulated reference
	constexpr auto operator&()&& {return involuter<Involution, decltype(&std::declval<Ref>())>{Involution{}, &r_};} // NOLINT(runtime/operator)
	// NOLINTNEXTLINE(google-runtime-operator): simulated reference
	constexpr auto operator&() & {return involuter<Involution, decltype(&std::declval<Ref>())>{Involution{}, &r_};}  // NOLINT(runtime/operator)
	// NOLINTNEXTLINE(google-runtime-operator): simulated reference
	constexpr auto operator&() const&{return involuter<Involution, decltype(&std::declval<decay_type const&>())>{Involution{}, &r_};} // NOLINT(runtime/operator)

	auto operator==(involuted  const& other) const {return r_ == other.r_;}
	auto operator!=(involuted  const& other) const {return r_ == other.r_;}

	auto operator==(decay_type const& other) const {return Involution{}(r_) == other;}
	auto operator!=(decay_type const& other) const {return Involution{}(r_) != other;}
};

template<class Involution, class It>
class involuter{
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
	explicit involuter(involuter<Involution, Other> const& o) : it_{o.it_} {}
	constexpr auto operator*() const {return reference{Involution{}, *it_};}
	constexpr auto operator==(involuter const& o) const {return it_==o.it_;}
	constexpr auto operator!=(involuter const& o) const {return it_!=o.it_;}
	constexpr auto operator+=(difference_type n) -> involuter& {it_+=n; return *this;}
	constexpr auto operator-=(difference_type n) -> involuter& {it_-=n; return *this;}
	constexpr auto operator+(difference_type n) const {return involuter{it_+n};}
	constexpr auto operator-(difference_type n) const {return involuter{it_-n};}
	constexpr auto operator->() const {return pointer{&*it_};}
};

//#if defined(__cpp_deduction_guides)
//template<class T, class F> involuted(F, T&&)->involuted<F, T const>;
//#endif

template<class Ref> using negated = involuted<std::negate<>, Ref>;
template<class It>  using negater = involuter<std::negate<>, It >;

class basic_conjugate_t{
	template<int N> struct prio : std::conditional_t<N!=0, prio<N-1>, std::true_type>{};
	template<class T> static auto _(prio<0>/**/, T const& t) DECLRETURN(std::conj(t))
	template<class T> static auto _(prio<1>/**/, T const& t) DECLRETURN(     conj(t))
	template<class T> static auto _(prio<2>/**/, T const& t) DECLRETURN(  T::conj(t))
	template<class T> static auto _(prio<3>/**/, T const& t) DECLRETURN(   t.conj( ))
 public:
	template<class T> static auto _(T const& t) DECLRETURN(_(prio<3>{}, t))
} basic_conjugate;

template<class T = void>
struct conjugate : private basic_conjugate_t{
	constexpr auto operator()(T const& arg) const DECLRETURN(_(arg))
};

template<>
struct conjugate<> : private basic_conjugate_t{
	template<class T>
	constexpr auto operator()(T const& arg) const DECLRETURN(_(arg))
};

//MAYBE_UNUSED static auto conj = [](std::complex<double> const& a){return std::conj(std::forward<decltype(a)>(a));};

#if defined(__NVCC__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wsubobject-linkage"
#endif
template<class ComplexRef> struct conjd : test::involuted<conjugate<>, ComplexRef>{
	explicit conjd(ComplexRef r) : test::involuted<conjugate<>, ComplexRef>(conjugate<>{}, r) {}
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
class indirect_real{// : std::iterator_traits<typename std::pointer_traits<P>::element_type::value_type*>{
	P impl_;
 public:
	explicit indirect_real(P const& p) : impl_{p} {}
	auto operator+(std::ptrdiff_t n) const {return indirect_real{impl_ + n};}
	// NOLINTNEXTLINE(cppcoreguidelines-pro-type-reinterpret-cast): extra real part as reference
	auto operator*() const -> decltype(auto) {return reinterpret_cast<std::array<double, 2>&>(*impl_)[0];}

	using difference_type = std::ptrdiff_t;
	using value_type = typename std::iterator_traits<P>::value_type;
	using pointer = void;
	using reference = void;
	using iterator_category = typename std::iterator_traits<P>::iterator_category;
};

} // namespace test

BOOST_AUTO_TEST_CASE(transformed_array) {
	namespace multi = boost::multi;
	{
		using complex = std::complex<double>;
		complex c{1., 2.};

		auto&& z = test::conjd<complex&>{c};
		BOOST_REQUIRE(( z == complex{1., -2.} ));

		BOOST_REQUIRE( real(z)  ==  1. );
		BOOST_REQUIRE( imag(z)  == -2. );
		BOOST_REQUIRE( z.real() ==  1. );
		BOOST_REQUIRE( z.imag() == -2. );
	}
	{
		double a = 5;

		auto&& c = test::involuted<test::neg_t, double&>(test::neg, a);
		BOOST_REQUIRE( c == -5. );

		c = 10.;
		BOOST_REQUIRE( c == 10. );
		BOOST_REQUIRE( a = -10. );
	}
	{
		multi::array<double, 1> A = { 0,  1,  2,  3,  4};
		auto&& A_ref = A.static_array_cast<double, double const*>();
		BOOST_REQUIRE( A_ref[2] == A[2] );
	}

	{
		multi::array<double, 1> A = { 0,  1,  2,  3,  4};
		multi::array<double, 1> mA = { -0,  -1,  -2,  -3, -4};
		auto&& mA_ref = A.static_array_cast<double, test::negater<double*>>();
		BOOST_REQUIRE( mA_ref[2] == mA[2] );
	}
	{
		multi::array<double, 2> A = {
			{ 0,  1,  2,  3,  4},
			{ 5,  6,  7,  8,  9},
			{10, 11, 12, 13, 14},
			{15, 16, 17, 18, 19}
		};
		multi::array<double, 2> mA = {
			{ -0,  -1,  -2,  -3,  -4},
			{ -5,  -6,  -7,  -8,  -9},
			{-10, -11, -12, -13, -14},
			{-15, -16, -17, -18, -19}
		};
		auto&& mA_ref = A.static_array_cast<double, test::negater<double*>>();
		BOOST_REQUIRE( mA_ref[1][1] == mA[1][1] );
	}
	 {
	#if defined(__cpp_deduction_guides)
		double Z[4][5] {  // NOLINT(cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays) : testing legacy types
			{ 0,  1,  2,  3,  4},
			{ 5,  6,  7,  8,  9},
			{10, 11, 12, 13, 14},
			{15, 16, 17, 18, 19}
		};
		auto&& d2DC = multi::make_array_ref(test::involuter<decltype(test::neg), double*>{test::neg, &Z[0][0]}, {4, 5});

		d2DC[1][1] = -66.;
		BOOST_REQUIRE( Z[1][1] == 66 );
	#endif
		{
			using complex = std::complex<double>;
			multi::array<complex, 2> d2D = {
				{ {0., 3.}, { 1., 9.}, { 2., 4.},  3.,  4.},
				{  5.     , { 6., 3.}, { 7., 5.},  8.,  9.},
				{ {1., 4.}, { 9., 1.},  12.     , 13., 14.},
				{  15.    ,  16.    ,   17.     , 18., 19.}
			};

			auto&& d2Dreal = d2D.reinterpret_array_cast<double>();
			BOOST_REQUIRE( d2Dreal[2][1] == 9. );

			d2Dreal[2][1] = 12.;
			BOOST_REQUIRE( d2D[2][1] == complex(12., 1.) );

			auto&& d2DrealT = rotated(d2D).reinterpret_array_cast<double>();
			BOOST_REQUIRE( d2DrealT[2][1] == 7. );

			multi::array<double, 2> d2real_copy = d2D.template reinterpret_array_cast<double>();//d2Dreal;
		}
		 {
			using complex = std::complex<double>;
			constexpr auto const I = complex{0., 1.};
			multi::array<complex, 2> A = {
				{ 1. + 3.*I, 3.- 2.*I, 4.+ 1.*I},
				{ 9. + 1.*I, 7.- 8.*I, 1.- 3.*I}
			};
			auto Aconj = A.static_array_cast<complex, test::conjr<complex*>>();
			BOOST_REQUIRE( Aconj[1][2] == conj(A[1][2]) );
		}
	}
}

