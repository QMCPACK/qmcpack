// Copyright 2019-2024 Alfredo A. Correa
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#ifndef BOOST_MULTI_ADAPTORS_BLAS_NUMERIC_HPP
#define BOOST_MULTI_ADAPTORS_BLAS_NUMERIC_HPP
#pragma once

// #include <boost/multi/adaptors/complex.hpp>

#include <boost/multi/array_ref.hpp>

#include <boost/multi/adaptors/complex/adl.hpp>

#include <boost/multi/adaptors/blas/complex_traits.hpp>

// #include <boost/multi/detail/pointer_traits.hpp>

#include <boost/multi/adaptors/blas/numeric/is_complex.hpp>

// #include <boost/multi/adaptors/complex.hpp>

#include <functional>                                        // for negate
#include <iterator>                                          // for iterator...
#include <memory>                                            // for pointer_...
#include <type_traits>                                       // for decay_t
#include <utility>                                           // for declval

#if defined(__NVCC__)
#define BOOST_MULTI_HD __host__ __device__
#else
#define BOOST_MULTI_HD
#endif

namespace boost {
namespace multi::blas {

template<class T> struct complex_dummy {
	T real;
	T imag;
};

template<
	class A, typename Complex = typename std::decay_t<A>::element, typename T = typename multi::blas::complex_traits<Complex>::real_type,
	class = std::enable_if_t<blas::numeric::is_complex_of<Complex, T>::value>>  // NOLINT(modernize-use-constraints) TODO(correaa) for C++20
auto real(A&& array)
	-> decltype(std::forward<A>(array).template reinterpret_array_cast<complex_dummy<T>>().template member_cast<T>(&complex_dummy<T>::real)) {
	return std::forward<A>(array).template reinterpret_array_cast<complex_dummy<T>>().template member_cast<T>(&complex_dummy<T>::real);
}

template<
	class A, class Complex = typename std::decay_t<A>::element_type, typename T = typename complex_traits<Complex>::real_type,
	class = std::enable_if_t<blas::numeric::is_complex_of<Complex, T>::value>>  // NOLINT(modernize-use-constraints) TODO(correaa) for C++20
auto imag(A&& array)
	-> decltype(std::forward<A>(array).template reinterpret_array_cast<complex_dummy<T>>().template member_cast<T>(&complex_dummy<T>::imag)) {
	return std::forward<A>(array).template reinterpret_array_cast<complex_dummy<T>>().template member_cast<T>(&complex_dummy<T>::imag);
}

template<class ComplexArr, class ComplexElem = typename std::decay_t<ComplexArr>::element, typename RealElem = typename ComplexElem::value_type,
	class = std::enable_if_t<blas::numeric::is_complex_of<ComplexElem, RealElem>::value>>  // NOLINT(modernize-use-constraints) TODO(correaa) for C++20
auto real_doubled(ComplexArr&& array) {  // produces a real view of complex array with the last dimension duplicated and with interleaved real imaginary parts
	return std::forward<ComplexArr>(array).template reinterpret_array_cast<RealElem>(2).rotated().flatted().unrotated();
}

template<class Ref, class Involution> class involuted;

template<class It, class F, class Reference = involuted<typename std::iterator_traits<It>::reference, F>> class involuter;  // IWYU pragma: keep  // bug in iwyu 0.22/18.1.8?

template<class Ref, class Involution>
class involuted {
	Ref        r_;  // [[no_unique_address]]  // NOLINT(cppcoreguidelines-avoid-const-or-ref-data-members)
	Involution f_;

 public:
	using decay_type = std::decay_t<decltype(std::declval<Involution>()(std::declval<Ref>()))>;

	constexpr explicit involuted(Ref& ref, Involution fun) : r_{ref}, f_{fun} {}  // r_{std::forward<Ref>(ref)}, f_{fun} {}
	constexpr explicit involuted(Ref& ref) : r_{ref}, f_{} {}

	~involuted() = default;

	involuted(involuted const&)     = delete;
	involuted(involuted&&) noexcept = default;

	auto operator=(involuted const& other) -> involuted&     = delete;
	auto operator=(involuted&& other) noexcept -> involuted& = default;

	constexpr auto decay() const& -> decay_type { return f_(r_); }

	constexpr explicit operator decay_type() & { return f_(r_); }
	constexpr explicit operator decay_type() const& { return f_(r_); }
	constexpr /*plct*/ operator decay_type() && { return f_(r_); }  // NOLINT(google-explicit-constructor,hicpp-explicit-conversions) //NOSONAR to allow terse syntax

	// constexpr auto operator*(decay_type const& other) const { return f_(r_) * other; }
	constexpr friend auto operator*(involuted const& self, decay_type const& other) { return self.f_(self.r_) * other; }

	template<class DecayType, class = decltype(std::declval<Ref&>() = (std::declval<Involution&>())(std::declval<DecayType&&>()))>
	constexpr auto operator=(DecayType&& other) & -> involuted& {
		r_ = f_(std::forward<DecayType>(other));
		return *this;
	}

	template<class DecayType, class = decltype(std::declval<Ref&>() = (std::declval<Involution&>())(std::declval<DecayType&&>()))>
	constexpr auto operator=(DecayType&& other) && -> involuted& {
		r_ = f_(std::forward<DecayType>(other));
		return *this;
	}

	template<class DecayType>
	friend constexpr auto operator==(involuted const& self, DecayType const& other)
		-> decltype(std::declval<decay_type>() == other) {
		return self.operator decay_type() == other;
	}
	template<class DecayType>
	friend constexpr auto operator!=(involuted const& self, DecayType const& other)
		-> decltype(std::declval<decay_type>() != other) {
		return self.operator decay_type() != other;
	}

	friend constexpr auto operator==(decay_type const& other, involuted const& self) -> bool {
		return other == self.operator decay_type();
	}

	friend constexpr auto operator!=(decay_type const& other, involuted const& self) -> bool {
		return other != self.operator decay_type();
	}

	template<class DecayType,
		std::enable_if_t<!std::is_base_of_v<involuted, DecayType>, int> =0>  // NOLINT(modernize-use-constraints) TODO(correaa) for C++20
	friend constexpr auto operator==(DecayType const& other, involuted const& self) {
		return other == self.operator decay_type();
	}
	template<class DecayType,
		std::enable_if_t<!std::is_base_of_v<involuted, DecayType>, int> =0>  // NOLINT(modernize-use-constraints) TODO(correaa) for C++20
	friend constexpr auto operator!=(DecayType const& other, involuted const& self) {
		return other != self.operator decay_type();
	}

	template<class Sink> friend constexpr auto operator<<(Sink&& sink, involuted const& self) -> Sink& {
		return std::forward<Sink>(sink) << self.operator decay_type();
	}

	constexpr auto conj() const& { return ::boost::multi::adl_conj(this->operator decay_type()); }

	template<class T = void*>
	friend constexpr auto imag(involuted const& self) {
		//->decltype(imag(std::declval<decay_type>())) {
		return self.operator decay_type().imag();
	}
};

#if defined(__cpp_deduction_guides)
template<class T, class F> involuted(T&&, F) -> involuted<T const, F>;
#endif

template<class It, class F>
auto default_allocator_of(involuter<It, F> const& iv) {
	return default_allocator_of(iv.it_);
}

template<class It, class F, class Reference>
class involuter {
	It it_;
	F  f_;  // [[no_unique_address]]
	template<class, class, class> friend class involuter;

 public:
	using difference_type          = typename std::iterator_traits<It>::difference_type;
	using value_type               = typename std::iterator_traits<It>::value_type;
	using pointer                  = involuter<It, F>;  // svoid; // typename std::iterator_traits<It>::pointer
	using reference                = Reference;
	using iterator_category        = typename std::iterator_traits<It>::iterator_category;
	using element_type             = typename std::pointer_traits<It>::element_type;
	template<class U> using rebind = involuter<typename std::pointer_traits<It>::template rebind<U>, F>;

	involuter() = default;

	BOOST_MULTI_HD constexpr explicit involuter(It it) : it_{std::move(it)}, f_{} {}
	BOOST_MULTI_HD constexpr explicit involuter(It it, F fun) : it_{std::move(it)}, f_{std::move(fun)} {}

	template<class Other, decltype(detail::implicit_cast<It>(typename Other::underlying_type{}))* = nullptr>
	// cppcheck-suppress noExplicitConstructor
	BOOST_MULTI_HD constexpr /*implct*/ involuter(Other const& other) : it_{other.it_}, f_{other.f_} {}  // NOLINT(google-explicit-constructor,hicpp-explicit-conversions) // NOSONAR inherit implicit conversion of underlying type
	template<class Other, decltype(detail::explicit_cast<It>(typename Other::underlying_type{}))* = nullptr>
	BOOST_MULTI_HD constexpr explicit involuter(Other const& other) : it_{other.it_}, f_{other.f_} {}

	constexpr auto operator*() const { return reference{*it_, f_}; }
	constexpr auto operator[](difference_type n) const { return reference{*(it_ + n), f_}; }

	// auto operator==(involuter const& other) const -> bool { return it_ == other.it_; }
	// auto operator!=(involuter const& other) const -> bool { return it_ != other.it_; }
	friend auto operator==(involuter const& slf, involuter const& thr) { return slf.it_ == thr.it_; }
	friend auto operator!=(involuter const& slf, involuter const& thr) { return slf.it_ != thr.it_; }

	constexpr auto operator+=(difference_type n) -> involuter& {
		it_ += n;
		return *this;
	}
	constexpr auto operator-=(difference_type n) -> involuter& {
		it_ -= n;
		return *this;
	}

	template<class = void>  // workaround for nvcc
	constexpr friend auto operator+(involuter lhs, difference_type n) { return lhs += n; }
	template<class = void>  // workaround for nvcc
	constexpr friend auto operator-(involuter lhs, difference_type n) { return lhs -= n; }

	auto operator-(involuter const& other) const { return it_ - other.it_; }

	explicit operator bool() const { return it_; }
	using underlying_type = It;
	friend /*constexpr*/ auto underlying(involuter const& self) -> underlying_type { return self.it_; }
	constexpr explicit operator It() const { return underlying(*this); }

	friend auto default_allocator_of(involuter const& inv) {
		using multi::default_allocator_of;
		return default_allocator_of(inv.it_);
	}

	using default_allocator_type = typename multi::pointer_traits<It>::default_allocator_type;

	friend auto get_allocator(involuter const& inv) {
		using boost::multi::get_allocator;
		return get_allocator(inv.it_);
	}
};

template<class Ref> using negated = involuted<Ref, std::negate<>>;
template<class It> using negater  = involuter<It, std::negate<>>;

struct conjugate {
	template<class Complex>
	constexpr auto operator()(Complex const& zee) const {
		//  using std::conj;  // for doubles?
		return conj(zee);
	}

#if defined(__CUDACC__)
	template<class Complex>
	constexpr auto operator()(::thrust::tagged_reference<Complex, ::thrust::cuda_cub::tag> zee) const {
		return conj(static_cast<Complex>(zee));
	}
#endif
#if defined(__HIPCC__)
	template<class Complex>
	constexpr auto operator()(::thrust::tagged_reference<Complex, ::thrust::hip::tag> zee) const {
		return conj(static_cast<Complex>(zee));
	}
#endif
};

template<class Ref> using conjugated = involuted<Ref, conjugate>;

template<class It> using conjugater = involuter<It, conjugate>;

template<class It> auto make_conjugater(It it) { return conjugater<It>{it}; }
template<class It> auto make_conjugater(conjugater<It> it) -> It { return underlying(it); }

template<class T> auto imag(involuted<T, conjugate> const& inv) { return inv.decay().imag(); }
template<class T> auto real(involuted<T, conjugate> const& inv) { return inv.decay().real(); }

template<class T> auto has_imag_fun_aux(T const& value) -> decltype((void)imag(value), std::true_type{});
inline auto            has_imag_fun_aux(...) -> decltype(std::false_type{});
template<class T> struct has_imag_fun : decltype(has_imag_fun_aux(std::declval<T>())) {};  // NOLINT(cppcoreguidelines-pro-type-vararg,hicpp-vararg)

template<class T> auto has_imag_mem_aux(T const& value) -> decltype((void)value.imag(), std::true_type{});
inline auto            has_imag_mem_aux(...) -> decltype(std::false_type{});
template<class T> struct has_imag_mem : decltype(has_imag_mem_aux(std::declval<T>())) {};  // NOLINT(cppcoreguidelines-pro-type-vararg,hicpp-vararg)

template<class T> struct has_imag : std::integral_constant<bool, (has_imag_fun<T>{} || has_imag_mem<T>{})> {};

template<class A = void>
struct is_complex_array : has_imag<std::decay_t<typename std::pointer_traits<std::decay_t<decltype(std::declval<A>().base())>>::element_type>> {};

template<class V> struct is_complex : has_imag<V> {};

template<class It>
auto        is_conjugated_aux(conjugater<It> const& /*self*/) -> std::true_type;
inline auto is_conjugated_aux(...) -> std::false_type;

template<class A = void> struct is_conjugated : decltype(is_conjugated_aux((std::declval<A>()).base())) {  // NOLINT(cppcoreguidelines-pro-type-vararg,hicpp-vararg)
	template<class AA> constexpr auto operator()(AA&& /*unused*/) { return is_conjugated_aux((std::declval<A>()).base()); }  // NOLINT(cppcoreguidelines-missing-std-forward)
};

template<class A, class D = std::decay_t<A>, typename Elem = typename D::element_type, typename Ptr = typename D::element_ptr,
	std::enable_if_t<!is_complex_array<A>{}, int> =0>  // NOLINT(modernize-use-constraints) TODO(correaa) for C++20
auto conj(A&& array) -> A&& {
	return std::forward<A>(array);
}

template<
	class A, class D = std::decay_t<A>, typename Elem = typename D::element_type,
	typename Ptr = std::decay_t<decltype(std::declval<A&&>().base())>,
	std::enable_if_t<!is_conjugated<A>{} && is_complex_array<A>{}, int> =0>  // NOLINT(modernize-use-constraints) TODO(correaa) for C++20
auto conj(A&& array) -> decltype(auto) {
	return std::forward<A>(array).template static_array_cast<Elem, conjugater<Ptr>>();
}

template<class A, class D = std::decay_t<A>, typename Elem = typename D::element_type,
         typename Ptr = typename decltype(std::declval<A&&>().base())::underlying_type,
		 std::enable_if_t<is_conjugated<A>{}, int> =0>  // NOLINT(modernize-use-constraints) TODO(correaa) for C++20
auto conj(A&& array)
	-> decltype(std::forward<A>(array).template static_array_cast<Elem, Ptr>()) {
	return std::forward<A>(array).template static_array_cast<Elem, Ptr>();
}

}  // end namespace multi::blas

template<class It, class F, class Reference>
auto default_allocator_of(multi::blas::involuter<It, F, Reference> it) {
	return multi::default_allocator_of(it.underlying());
}

}  // end namespace boost

#undef BOOST_MULTI_HD

#endif
