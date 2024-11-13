// Copyright 2019-2024 Alfredo A. Correa
// Copyright 2024 Matt Borland
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/multi/array.hpp>  // for array, subarray, static_array

#include <boost/core/lightweight_test.hpp>

#include <array>        // for array  // IWYU pragma: keep
#include <complex>      // for complex, operator*, operator+
#include <cstddef>      // for ptrdiff_t
#include <cstdlib>      // for abs
#include <functional>   // for negate  // IWYU pragma: keep
#include <iterator>     // for iterator_traits
#include <memory>       // for pointer_traits
#include <type_traits>  // for decay_t, conditional_t, true_type
#include <utility>      // for move, declval

#define BOOST_MULTI_DECLRETURN(ExpR) \
	->decltype(ExpR) { return ExpR; }  // NOLINT(cppcoreguidelines-macro-usage) saves a lot of typing

namespace test {

struct neg_t {
	template<class T>
	constexpr auto operator()(T const& value) const -> decltype(-value) { return -value; }
};
constexpr inline neg_t neg;

}  // end namespace test

namespace test {

template<class Involution, class Ref>
class involuted {
	Ref         r_;  // NOLINT(cppcoreguidelines-avoid-const-or-ref-data-members)
	friend auto underlying(involuted& self) -> decltype(auto) { return self.r_; }
	friend auto underlying(involuted&& self) -> decltype(auto) { return std::move(self).r_; }
	friend auto underlying(involuted const& self) -> decltype(auto) { return self.r_; }

 public:
	using decay_type = std::decay_t<decltype(std::declval<Involution>()(std::declval<Ref>()))>;
	constexpr involuted(Involution /*stateless*/, Ref ref) : r_{ref} {}
	involuted(involuted const&)     = default;
	involuted(involuted&&) noexcept = default;

	auto operator=(involuted&&) -> involuted&      = delete;
	auto operator=(involuted const&) -> involuted& = default;
	auto operator=(decay_type const& other) -> involuted& {  // NOLINT(fuchsia-trailing-return) simulate reference
		r_ = Involution{}(other);
		return *this;
	}
	constexpr explicit operator decay_type() const { return Involution{}(r_); }
	// NOLINTNEXTLINE(google-runtime-operator): simulated reference
	// constexpr auto operator&() && { return involuter<Involution, decltype(&std::declval<Ref>())>{Involution{}, &r_}; }  // NOLINT(runtime/operator)
	// NOLINTNEXTLINE(google-runtime-operator): simulated reference
	// constexpr auto operator&() & { return involuter<Involution, decltype(&std::declval<Ref>())>{Involution{}, &r_}; }  // NOLINT(runtime/operator)
	// NOLINTNEXTLINE(google-runtime-operator): simulated reference
	// constexpr auto operator&() const& {  // NOLINT(runtime/operator)
	//  return involuter<Involution, decltype(&std::declval<decay_type const&>())>{Involution{}, &r_};
	// }

	auto operator==(involuted const& other) const { return r_ == other.r_; }
	auto operator!=(involuted const& other) const { return r_ == other.r_; }

#if defined(__clang__)
	#pragma clang diagnostic push
	#pragma clang diagnostic ignored "-Wfloat-equal"
#elif defined(__GNUC__)
	#pragma GCC diagnostic push
	#pragma GCC diagnostic ignored "-Wfloat-equal"
#endif
	auto operator==(decay_type const& other) const { return Involution{}(r_) == other; }
	auto operator!=(decay_type const& other) const { return Involution{}(r_) != other; }
#if defined(__clang__)
	#pragma clang diagnostic pop
#elif defined(__GNUC__)
	#pragma GCC diagnostic pop
#endif
	~involuted() = default;
};

template<class Involution, class It>
class involuter {
	It it_;
	template<class, class> friend class involuter;

 public:
	using pointer         = involuter<Involution, typename std::iterator_traits<It>::pointer>;
	using element_type    = typename std::pointer_traits<It>::element_type;
	using difference_type = typename std::pointer_traits<It>::difference_type;
	template<class U>
	using rebind            = involuter<Involution, typename std::pointer_traits<It>::template rebind<U>>;
	using reference         = involuted<Involution, typename std::iterator_traits<It>::reference>;
	using value_type        = typename std::iterator_traits<It>::value_type;
	using iterator_category = typename std::iterator_traits<It>::iterator_category;

	constexpr explicit involuter(It it) : it_{std::move(it)} {}
	constexpr involuter(Involution /*stateless*/, It it) : it_{std::move(it)} {}  // f_{std::move(f)}{}
	template<class Other>
	explicit involuter(involuter<Involution, Other> const& other) : it_{other.it_} {}

	constexpr auto operator*() const { return reference{Involution{}, *it_}; }
	constexpr auto operator->() const { return pointer{&*it_}; }

	constexpr auto operator==(involuter const& other) const { return it_ == other.it_; }
	constexpr auto operator!=(involuter const& other) const { return it_ != other.it_; }

	constexpr auto operator+=(difference_type n) -> involuter& {
		it_ += n;
		return *this;
	}
	constexpr auto operator-=(difference_type n) -> involuter& {
		it_ -= n;
		return *this;
	}

#if defined(__clang__)
	#pragma clang diagnostic push
	#pragma clang diagnostic ignored "-Wunknown-warning-option"
	#pragma clang diagnostic ignored "-Wunsafe-buffer-usage"
#endif

	constexpr auto operator+(difference_type n) const { return involuter{it_ + n}; }
	constexpr auto operator-(difference_type n) const { return involuter{it_ - n}; }

#if defined(__clang__)
	#pragma clang diagnostic pop
#endif
};

template<class Ref> using negated = involuted<std::negate<>, Ref>;
template<class It> using negater  = involuter<std::negate<>, It>;

class basic_conjugate_t {
	// clang-format off
	template<int N> struct prio : std::conditional_t<N != 0, prio<N - 1>, std::true_type> {};

	template<class T> static auto _(prio<0> /**/, T const& value) BOOST_MULTI_DECLRETURN( std::conj(value))
	template<class T> static auto _(prio<1> /**/, T const& value) BOOST_MULTI_DECLRETURN(      conj(value))
	template<class T> static auto _(prio<2> /**/, T const& value) BOOST_MULTI_DECLRETURN(   T::conj(value))
	template<class T> static auto _(prio<3> /**/, T const& value) BOOST_MULTI_DECLRETURN(value.conj()     )

 public:
	template<class T>
	static auto _(T const& value) BOOST_MULTI_DECLRETURN(_(prio<3>{}, value))
	// clang-format on
};

template<class T = void>
struct conjugate : private basic_conjugate_t {
	constexpr auto operator()(T const& arg) const BOOST_MULTI_DECLRETURN(_(arg))
};

template<>
struct conjugate<> : private basic_conjugate_t {
	template<class T>
	constexpr auto operator()(T const& arg) const BOOST_MULTI_DECLRETURN(_(arg))
};

#if defined(__NVCC__)
	#pragma GCC diagnostic push
	#pragma GCC diagnostic ignored "-Wsubobject-linkage"
#endif
template<class ComplexRef> struct conjd : test::involuted<conjugate<>, ComplexRef> {
	explicit conjd(ComplexRef ref) : test::involuted<conjugate<>, ComplexRef>(conjugate<>{}, ref) {}
	auto real() const { return underlying(*this).real(); }
	auto imag() const {
		return negated<decltype(underlying(std::declval<test::involuted<conjugate<>, ComplexRef> const&>()).imag())>{std::negate<>{}, underlying(*this).imag()};
	}
	friend auto real(conjd const& self) -> decltype(auto) {
		using std::real;
		return real(static_cast<typename conjd::decay_type>(self));
	}
	friend auto imag(conjd const& self) -> decltype(auto) {
		using std::imag;
		return imag(static_cast<typename conjd::decay_type>(self));
	}
};
#if defined(__NVCC__)
	#pragma GCC diagnostic pop
#endif

#if defined(__cpp_deduction_guides)
template<class T> conjd(T&&) -> conjd<T>;
#endif

template<class Complex> using conjr = test::involuter<conjugate<>, Complex>;

template<class P = std::complex<double>*>
class indirect_real {
	P impl_;

 public:
	explicit indirect_real(P const& ptr) : impl_{ptr} {}
	auto operator+(std::ptrdiff_t n) const { return indirect_real{impl_ + n}; }
	// NOLINTNEXTLINE(cppcoreguidelines-pro-type-reinterpret-cast): extra real part as reference
	auto operator*() const -> decltype(auto) { return reinterpret_cast<std::array<double, 2>&>(*impl_)[0]; }

	using difference_type   = std::ptrdiff_t;
	using value_type        = typename std::iterator_traits<P>::value_type;
	using pointer           = void;
	using reference         = void;
	using iterator_category = typename std::iterator_traits<P>::iterator_category;
};

}  // namespace test

#define BOOST_AUTO_TEST_CASE(CasenamE) /**/

auto main() -> int {  // NOLINT(readability-function-cognitive-complexity,bugprone-exception-escape)
	BOOST_AUTO_TEST_CASE(transformed_array) {
		namespace multi = boost::multi;
		{
			using complex = std::complex<double>;
			complex cee{1.0, 2.0};

			auto&& zee = test::conjd<complex&>{cee};
			BOOST_TEST(( zee == complex{1.0, -2.0} ));

			BOOST_TEST( std::abs( real(zee) - +1.0 ) < 1E-6);
			BOOST_TEST( std::abs( imag(zee) - -2.0 ) < 1E-6);
		}
		{
			int val = 50;

			auto&& negd_a = test::involuted<test::neg_t, int&>(test::neg, val);
			BOOST_TEST( negd_a == -50 );

			negd_a = 100;
			BOOST_TEST( negd_a == 100 );
			BOOST_TEST( val == -100 );
		}
		{
			multi::array<double, 1> arr = {0.0, 1.0, 2.0, 3.0, 4.0};
			auto&&                  ref = arr.static_array_cast<double, double const*>();
			BOOST_TEST( std::abs(ref[2] - arr[2]) < 1E-6 );
		}
		{
			multi::array<double, 1> const arr      = {+0.0, +1.0, +2.0, +3.0, +4.0};
			multi::array<double, 1>       neg      = {-0.0, -1.0, -2.0, -3.0, -4.0};
			auto&&                        negd_arr = arr.static_array_cast<double, test::negater<double*>>();
			BOOST_TEST( negd_arr[2] == neg[2] );
		}
		{
			multi::array<double, 2> const arr = {
				{ +0.0,  +1.0,  +2.0,  +3.0,  +4.0},
				{ +5.0,  +6.0,  +7.0,  +8.0,  +9.0},
				{+10.0, +11.0, +12.0, +13.0, +14.0},
				{+15.0, +16.0, +17.0, +18.0, +19.0},
			};
			multi::array<double, 2> neg = {
				{ -0.0,  -1.0,  -2.0,  -3.0,  -4.0},
				{ -5.0,  -6.0,  -7.0,  -8.0,  -9.0},
				{-10.0, -11.0, -12.0, -13.0, -14.0},
				{-15.0, -16.0, -17.0, -18.0, -19.0},
			};
			// auto&& negd_arr = arr.static_array_cast<double, test::negater<double*>>();  // not compile, ok, read only
			auto&& negd_arr = arr.static_array_cast<double, test::negater<double const*>>();
			BOOST_TEST( negd_arr[1][1] == neg[1][1] );
			BOOST_TEST( negd_arr[1][1] == -6.0 );
			// negd_arr2[1][1] = 3.0;  // can't compile, ok, read-only
		}
		{
#if defined(__cpp_deduction_guides)
	#if defined(__clang__)
		#pragma clang diagnostic push
		#pragma clang diagnostic ignored "-Wunknown-warning-option"
		#pragma clang diagnostic ignored "-Wunsafe-buffer-usage"
	#endif

			// NOLINTNEXTLINE(cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays) : testing legacy types
			double zee[4][5]{
				{ 0.0,  1.0,  2.0,  3.0,  4.0},
				{ 5.0,  6.0,  7.0,  8.0,  9.0},
				{10.0, 11.0, 12.0, 13.0, 14.0},
				{15.0, 16.0, 17.0, 18.0, 19.0},
			};
			auto&& d2DC = multi::make_array_ref(test::involuter<decltype(test::neg), double*>{test::neg, &zee[0][0]}, {4, 5});

			d2DC[1][1] = -66.0;

			BOOST_TEST( std::abs( zee[1][1] - 66.0) < 1E-6 );

	#if defined(__clang__)
		#pragma clang diagnostic pop
	#endif

#endif
			{
				using complex = std::complex<double>;

				multi::array<complex, 2> d2D = {
					{ {0.0, 3.0},  {1.0, 9.0},  {2.0, 4.0},  {3.0, 0.0},  {4.0, 0.0}},
					{ {5.0, 0.0},  {6.0, 3.0},  {7.0, 5.0},  {8.0, 0.0},  {9.0, 0.0}},
					{ {1.0, 4.0},  {9.0, 1.0}, {12.0, 0.0}, {13.0, 0.0}, {14.0, 0.0}},
					{{15.0, 0.0}, {16.0, 0.0}, {17.0, 0.0}, {18.0, 0.0}, {19.0, 0.0}},
				};

				auto&& d2Dreal = d2D.reinterpret_array_cast<double>();
				BOOST_TEST( std::abs( d2Dreal[2][1] - 9.0) < 1E-6 );

				d2Dreal[2][1] = 12.0;
				BOOST_TEST( d2D[2][1] == complex(12.0, 1.0) );

				auto&& d2DrealT = d2D.rotated().reinterpret_array_cast<double>();
				BOOST_TEST( std::abs( d2DrealT[2][1] -  7.0) < 1E-6);

				multi::array<double, 2> const d2Dreal_copy{d2D.template reinterpret_array_cast<double>()};
				BOOST_TEST( d2Dreal_copy == d2Dreal );
			}
			{
				using complex = std::complex<double>;

				auto const I = complex{0.0, 1.0};  // NOLINT(readability-identifier-length) imaginary unit

				multi::array<complex, 2> arr = {
					{1.0 + 3.0 * I, 3.0 - 2.0 * I, 4.0 + 1.0 * I},
					{9.0 + 1.0 * I, 7.0 - 8.0 * I, 1.0 - 3.0 * I},
				};
				auto conjd_arr = arr.static_array_cast<complex, test::conjr<complex*>>();
				BOOST_TEST( conjd_arr[1][2] == conj(arr[1][2]) );
			}
		}
	}

#if !defined(__NVCC__) && defined(__GNU_MINOR__) && (__GNUC_MINOR__ > 7)
	BOOST_AUTO_TEST_CASE(transformed_to_string) {
		namespace multi = boost::multi;

		multi::array<int, 2> const AA = {
			{1, 2},
			{3, 4},
		};
		multi::array<std::string, 2> BB = AA.element_transformed([](int ee) noexcept { return std::to_string(ee); });

		BOOST_TEST( BB[1][1] == "4" );
	}
#endif

	/* accumulate rows */
	{
		namespace multi = boost::multi;

		multi::array<int, 2> Arr({100, 200}, 1);

		multi::array<int, 1> vv({200}, 0);
		for(auto i : Arr.extension()) {  // NOLINT(altera-unroll-loops)
		    for(auto j : vv.extension()) {  // NOLINT(altera-unroll-loops)
		        vv[j] += Arr[i][j];
		    }
		}

		// auto const v = std::reduce(
		//     A.begin(), A.end(),
		//     multi::array<int, 1>({200}, 0),
		//     [](auto&& acc, auto&& row) {
		//         std::transform(
		//             acc.begin(), vv.end(), row.begin(), acc.begin(),
		//             [](auto const& e1, auto const& e2) {
		//                 return e1 + e2;
		//             }
		//         );
		//         return acc;
		//     }
		// );

		// auto const v = std::transform_reduce(
		//     A.begin(), A.end(),
		//     multi::array<int, 1>({200}, 0),
		//     [](auto&& vv, auto&& row) {
		//         std::transform(
		//             vv.begin(), vv.end(), row.begin(), vv.begin(),
		//             [](auto const& e1, auto const& e2) {
		//                 return e1 + e2;
		//             }
		//         );
		//         return vv;
		//     },
		//     [](auto&& vv) {return vv;}
		// );

		// multi::array<int, 1> v({200}, 0);

		// std::transform(
		//     (~A).begin(), (~A).end(), v.begin(),
		//     [](auto const& col) { return std::reduce(col.begin(), col.end()); }
		// );

		BOOST_TEST(vv[0] == 100);
	}

	return boost::report_errors();
}

#undef BOOST_MULTI_DECLRETURN
