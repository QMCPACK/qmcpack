// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;autowrap:nil;-*-
// Copyright 2019-2022 Alfredo A. Correa

#ifndef MULTI_ADAPTORS_BLAS_NUMERIC_HPP
#define MULTI_ADAPTORS_BLAS_NUMERIC_HPP

#include "../../array_ref.hpp"
#include "../../complex.hpp"

#include "../../memory/pointer_traits.hpp"

#include "numeric/is_complex.hpp"

namespace boost {
namespace multi::blas {

template<class T> struct complex_dummy {T real; T imag;};

template<
	class A, typename Complex = typename std::decay_t<A>::element, typename T=typename Complex::value_type,
	class=std::enable_if_t<blas::numeric::is_complex_of<Complex, T>::value>
>
auto real(A&& array)
->decltype(std::forward<A>(array).template reinterpret_array_cast<complex_dummy<T>>().template member_cast<T>(&complex_dummy<T>::real)){
	return std::forward<A>(array).template reinterpret_array_cast<complex_dummy<T>>().template member_cast<T>(&complex_dummy<T>::real);}

template<
	class A, class Complex = typename std::decay_t<A>::element_type, typename T=typename Complex::value_type,
	class=std::enable_if_t<blas::numeric::is_complex_of<Complex, T>::value>
>
auto imag(A&& array)
->decltype(std::forward<A>(array).template reinterpret_array_cast<complex_dummy<T>>().template member_cast<T>(&complex_dummy<T>::imag)){
	return std::forward<A>(array).template reinterpret_array_cast<complex_dummy<T>>().template member_cast<T>(&complex_dummy<T>::imag);}

template<class ComplexArr, class ComplexElem = typename std::decay_t<ComplexArr>::element, typename RealElem = typename ComplexElem::value_type,
	class=std::enable_if_t<blas::numeric::is_complex_of<ComplexElem, RealElem>::value>
>
auto real_doubled(ComplexArr&& array) {  // produces a real view of complex array with the last dimension duplicated and with interleaved real imaginary parts
	return std::forward<ComplexArr>(array).template reinterpret_array_cast<RealElem>(2).rotated().flatted().unrotated();
}

template<class Ref, class Involution> class involuted;

template<class It, class F, class Reference = involuted<typename std::iterator_traits<It>::reference, F> > class involuter;

template<class Ref, class Involution>
class involuted {
	Ref r_;  // [[no_unique_address]] 
	Involution f_;

public:
	using decay_type =std::decay_t<decltype(std::declval<Involution>()(std::declval<Ref>()))>;

	constexpr explicit involuted(Ref ref, Involution fun) : r_{std::forward<Ref>(ref)}, f_{fun}{}
	constexpr explicit involuted(Ref ref) : r_{std::forward<Ref>(ref)}, f_{}{}

	auto operator=(involuted const& other) -> involuted& = delete;

	~involuted() = default;
	involuted(involuted const&) = delete;
	involuted(involuted&&) noexcept = default; // for C++14
	auto operator=(involuted&& other) noexcept -> involuted&{
		r_ = std::move(other.r_);
		return *this;
	}

	constexpr auto decay() const& -> decay_type{return f_(r_);}

	constexpr explicit operator decay_type()      &{return f_(r_);}
	constexpr explicit operator decay_type() const&{return f_(r_);}
	constexpr explicit operator decay_type()     &&{return f_(r_);}

	constexpr auto operator*(decay_type const& other) const{return f_(r_)*other;}
	constexpr auto operator&()&& -> decltype(auto){ // NOLINT(google-runtime-operator) : reference-like object
		return involuter<decltype(&std::declval<Ref>()), Involution>{&r_, f_};
	}

	template<class DecayType, class = decltype(std::declval<Ref&>() = (std::declval<Involution&>())(std::declval<DecayType&&>()))>
	constexpr auto operator=(DecayType&& other)& -> involuted&{
		r_=f_(std::forward<DecayType>(other));
		return *this;
	}

	template<class DecayType, class = decltype(std::declval<Ref&>() = (std::declval<Involution&>())(std::declval<DecayType&&>()))>
	constexpr auto operator=(DecayType&& other)&& -> involuted&{
		r_=f_(std::forward<DecayType>(other));
		return *this;
	}

	template<class DecayType>
	constexpr auto operator==(DecayType&& other) const
	->decltype(this->operator decay_type()==other){
		return this->operator decay_type()==other;}
	template<class DecayType>
	constexpr auto operator!=(DecayType&& other) const
	->decltype(this->operator decay_type()!=other){
		return this->operator decay_type()!=other;}

	friend constexpr auto operator==(decay_type const& other, involuted const& self){
		return other == self.operator decay_type();}

	template<class DecayType, std::enable_if_t<not std::is_base_of<involuted, DecayType>{}, int> =0>
	friend constexpr auto operator==(DecayType&& other, involuted const& self) {
		return other == self.operator decay_type();
	}
	template<class DecayType, std::enable_if_t<not std::is_base_of<involuted, DecayType>{}, int> =0>
	friend constexpr auto operator!=(DecayType&& other, involuted const& self) {
		return other != self.operator decay_type();\
	}
//	auto imag() const{return static_cast<decay_type>(*this).imag();}
	template<class Sink> friend constexpr auto operator<<(Sink&& sink, involuted const& self) -> Sink& {
		return sink<< self.operator decay_type();
	}
	constexpr auto conj() const& {return adl_conj(operator decay_type());}

	template<class T = void*>
	friend constexpr auto imag(involuted const& self)
	->decltype(adl_imag(std::declval<decay_type>())) {
		return adl_imag(self.operator decay_type()); }
};

#if defined(__cpp_deduction_guides)
template<class T, class F> involuted(T&&, F) -> involuted<T const, F>;
//template<class T, class F> involuted(T&, F)->involuted<T&, F>;
//template<class T, class F> involuted(T const&, F)->involuted<T const&, F>;
#endif

template<class It, class F, class Reference>
class involuter;

template<class It, class F>
auto default_allocator_of(involuter<It, F> const& iv) {
	return default_allocator_of(iv.it_);
}

template<class It, class F, class Reference>
class involuter {
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
//	~involuter() = default;

	HD constexpr explicit involuter(It it)      : it_{std::move(it)}, f_{} {}
	HD constexpr explicit involuter(It it, F fun) : it_{std::move(it)}, f_{std::move(fun)} {}

//	involuter(involuter const& other) = default;

	template<class Other, decltype(multi::implicit_cast<It>(typename Other::underlying_type{}))* = nullptr>
	// cppcheck-suppress noExplicitConstructor
	HD constexpr/*implct*/involuter(Other const& other) : it_{other.it_}, f_{other.f_}{} // NOLINT(google-explicit-constructor,hicpp-explicit-conversions) : inherit implicit conversion of underlying type
	template<class Other, decltype(explicit_cast<It>(typename Other::underlying_type{}))* = nullptr>
	HD constexpr explicit involuter(Other const& other) : it_{other.it_}, f_{other.f_}{}

	constexpr auto operator*() const {return reference{*it_, f_};}
	constexpr auto operator[](difference_type n) const {return reference{*(it_ + n), f_};}

	auto operator==(involuter const& other) const -> bool {return it_ == other.it_;}
	auto operator!=(involuter const& other) const -> bool {return it_ != other.it_;}

	constexpr auto operator+=(difference_type n) -> involuter& {it_ += n; return *this;}
	constexpr auto operator-=(difference_type n) -> involuter& {it_ -= n; return *this;}

	constexpr auto operator+(difference_type n) const {return involuter{it_ + n, f_};}
	constexpr auto operator-(difference_type n) const {return involuter{it_ - n, f_};}

	auto operator-(involuter const& other) const{return it_ - other.it_;}

	explicit operator bool() const{return it_;}
	using underlying_type = It;
	friend /*constexpr*/ auto underlying(involuter const& self) -> underlying_type{return self.it_;}
	constexpr explicit operator It() const {return underlying(*this);}
//	friend auto get_allocator(involuter const& self){return get_allocator(self.it_);}
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

template<class Ref> using negated = involuted<Ref, std::negate<>>;
template<class It>  using negater = involuter<It, std::negate<>>;

struct conjugate {
	template<class Complex>
	auto operator()(Complex&& zee) const -> decltype(auto) {
	//  using std::conj;  /*for doubles?*/
		return multi::adl_conj(std::forward<Complex>(zee));  // this is needed by icc
	}
};

template<class Ref> using conjugated = involuted<Ref, conjugate>;

template<class It> using conjugater = involuter<It, conjugate>;//, conjugated<typename std::iterator_traits<It>::reference> >;

template<class It> auto make_conjugater(It it){return conjugater<It>{it};}
template<class It> auto make_conjugater(conjugater<It> it) -> It {return underlying(it);}

template<class T> auto imag(involuted<T, conjugate> const& inv) {return inv.decay().imag();}
template<class T> auto real(involuted<T, conjugate> const& inv) {return inv.decay().real();}

template<class T> auto has_imag_fun_aux(T const& value) -> decltype(imag(value), std::true_type {});
                  auto has_imag_fun_aux(...           ) -> decltype(             std::false_type{});
template<class T> struct has_imag_fun : decltype(has_imag_fun_aux(std::declval<T>())){};


template<class T> auto has_imag_mem_aux(T const& value) -> decltype(value.imag(), std::true_type {});
                  auto has_imag_mem_aux(...           ) -> decltype(         std::false_type{});
template<class T> struct has_imag_mem : decltype(has_imag_mem_aux(std::declval<T>())){};

template<class T> struct has_imag : std::integral_constant<bool, (has_imag_fun<T>{} or has_imag_mem<T>{})>{};

template<class A = void> 
struct is_complex_array : has_imag<std::decay_t<decltype(*base(std::declval<A>()))>> {};
//	template<class T> static auto _(T const& t) -> has_imag<T>;
//	constexpr explicit operator bool()      &{return decltype(_(*base(std::declval<A>()))){};}
//	constexpr explicit operator bool()     &&{return decltype(_(*base(std::declval<A>()))){};}
//	constexpr operator bool() const&{return decltype(_(*base(std::declval<A>()))){};}
//	static constexpr bool value = decltype(_(*base(std::declval<A>()))){};
//	template<class AA> constexpr auto operator()(AA&& /*unused*/){return _(*base(std::declval<A>()));}
//};

template<class V> struct is_complex : has_imag<V> {};

template<class It> 
auto is_conjugated_aux(conjugater<It> /*self*/) -> std::true_type ;
auto is_conjugated_aux(...                    ) -> std::false_type;

template<class A = void> struct is_conjugated : decltype(is_conjugated_aux(base(std::declval<A>()))) {
	template<class AA> constexpr auto operator()(AA&& /*unused*/) {return is_conjugated_aux(base(std::declval<A>()));}
};

template<class A, class D = std::decay_t<A>, typename Elem=typename D::element_type, typename Ptr=typename D::element_ptr,
	std::enable_if_t<not is_complex_array<A>{}, int> =0>
auto conj(A&& array) -> A&& {
	return std::forward<A>(array);
}

template<class A, class D = std::decay_t<A>, typename Elem=typename D::element_type, typename Ptr=typename D::element_ptr,
	std::enable_if_t<not is_conjugated<A>{} and is_complex_array<A>{}, int> =0>
auto conj(A&& array) -> decltype(auto) {
	return std::forward<A>(array).template static_array_cast<Elem, conjugater<Ptr>>();
}

template<class A, class D = std::decay_t<A>, typename Elem=typename D::element_type, typename Ptr=typename D::element_ptr::underlying_type,
	std::enable_if_t<    is_conjugated<A>{}, int> =0>
auto conj(A&& array)
->decltype(std::forward<A>(array).template static_array_cast<Elem, Ptr>()) {
	return std::forward<A>(array).template static_array_cast<Elem, Ptr>(); }

} // end namespace multi::blas

template<class It, class F, class Reference>
auto default_allocator_of(multi::blas::involuter<It, F, Reference> it) {
	return multi::default_allocator_of(it.underlying());
}

} // end namespace boost

#endif
