// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;-*-
// Copyright 2019-2022 Alfredo A. Correa

#ifndef MULTI_ADAPTORS_BLAS_HERK_HPP
#define MULTI_ADAPTORS_BLAS_HERK_HPP

#include "../blas/copy.hpp" 
#include "../blas/core.hpp"
#include "../blas/filling.hpp"
#include "../blas/operations.hpp"
#include "../blas/side.hpp"
#include "../blas/syrk.hpp" // fallback to real case

namespace boost::multi::blas {

template<class A, std::enable_if_t<not is_conjugated<A>{}, int> =0> 
auto base_aux(A&& array)
->decltype(base(array)) {
	return base(array); }

template<class A, std::enable_if_t<    is_conjugated<A>{}, int> =0>
auto base_aux(A&& array)
->decltype(underlying(base(array))) {
	return underlying(base(array)); }

using core::herk;

template<class AA, class BB, class A2D, class C2D, class = typename A2D::element_ptr, std::enable_if_t<is_complex_array<C2D>{}, int> =0>
auto herk(filling c_side, AA alpha, A2D const& a, BB beta, C2D&& c) -> C2D&& {  // NOLINT(readability-function-cognitive-complexity,readability-identifier-length) 74, BLAS naming
	assert( a.size() == c.size() ); // NOLINT(cppcoreguidelines-pro-bounds-array-to-pointer-decay,hicpp-no-array-decay)
	assert( c.size() == rotated(c).size() ); // NOLINT(cppcoreguidelines-pro-bounds-array-to-pointer-decay,hicpp-no-array-decay)
	if(c.is_empty()) {return std::forward<C2D>(c);}
	if constexpr(is_conjugated<C2D>{}) {
		herk(flip(c_side), alpha, a, beta, hermitized(c));
	} else {
		auto base_a = base_aux(a);
		auto base_c = base_aux(c); //  static_assert( not is_conjugated<C2D>{}, "!" );
		if constexpr(is_conjugated<A2D>{}) {
		//	auto& ctxt = *blas::default_context_of(underlying(a.base()));
			// if you get an error here might be due to lack of inclusion of a header file with the backend appropriate for your type of iterator
			if     (stride(a)==1 and stride(c)!=1) {herk(c_side==filling::upper?'L':'U', 'N', size(c), size(rotated(a)), &alpha, base_a, stride(rotated(a)), &beta, base_c, stride(c));}
			else if(stride(a)==1 and stride(c)==1) {
				if(size(a)==1)                     {herk(c_side==filling::upper?'L':'U', 'N', size(c), size(rotated(a)), &alpha, base_a, stride(rotated(a)), &beta, base_c, stride(c));}
				else                               {assert(0);} // NOLINT(cppcoreguidelines-pro-bounds-array-to-pointer-decay,hicpp-no-array-decay)
			}
			else if(stride(a)!=1 and stride(c)==1) {herk(c_side==filling::upper?'U':'L', 'C', size(c), size(rotated(a)), &alpha, base_a, stride(        a ), &beta, base_c, stride(rotated(c)));}
			else if(stride(a)!=1 and stride(c)!=1) {herk(c_side==filling::upper?'L':'U', 'C', size(c), size(rotated(a)), &alpha, base_a, stride(        a ), &beta, base_c, stride(        c ));}
			else                                   {assert(0);} // NOLINT(cppcoreguidelines-pro-bounds-array-to-pointer-decay,hicpp-no-array-decay)
		} else {
		//	auto& ctxt = *blas::default_context_of(           a.base() );
			if     (stride(a)!=1 and stride(c)!=1) {herk(c_side==filling::upper?'L':'U', 'C', size(c), size(rotated(a)), &alpha, base_a, stride(        a ), &beta, base_c, stride(c));}
			else if(stride(a)!=1 and stride(c)==1) {
				if(size(a)==1)                     {herk(c_side==filling::upper?'L':'U', 'N', size(c), size(rotated(a)), &alpha, base_a, stride(rotated(a)), &beta, base_c, stride(rotated(c)));}
				else                               {assert(0);} // NOLINT(cppcoreguidelines-pro-bounds-array-to-pointer-decay,hicpp-no-array-decay)
			}
			else if(stride(a)==1 and stride(c)!=1) {assert(0);} // NOLINT(cppcoreguidelines-pro-bounds-array-to-pointer-decay,hicpp-no-array-decay)
			else if(stride(a)==1 and stride(c)==1) {herk(c_side==filling::upper?'U':'L', 'N', size(c), size(rotated(a)), &alpha, base_a, stride(rotated(a)), &beta, base_c, stride(rotated(c)));}
		//  else                                   {assert(0);} // NOLINT(cppcoreguidelines-pro-bounds-array-to-pointer-decay,hicpp-no-array-decay)
		}
	}
	return std::forward<C2D>(c);
}

template<class AA, class BB, class A2D, class C2D, class = typename A2D::element_ptr, std::enable_if_t<not is_complex_array<C2D>{}, int> =0>
auto herk(filling c_side, AA alpha, A2D const& a, BB beta, C2D&& c)  // NOLINT(readability-identifier-length) BLAS naming
->decltype(syrk(c_side, alpha, a, beta, std::forward<C2D>(c))) {
	return syrk(c_side, alpha, a, beta, std::forward<C2D>(c)); }

template<class AA, class A2D, class C2D, class = typename A2D::element_ptr>
auto herk(filling c_side, AA alpha, A2D const& a, C2D&& c)  // NOLINT(readability-identifier-length) BLAS naming
->decltype(herk(c_side, alpha, a, 0., std::forward<C2D>(c))) {
	return herk(c_side, alpha, a, 0., std::forward<C2D>(c)); }

template<typename AA, class A2D, class C2D>
auto herk(AA alpha, A2D const& a, C2D&& c)  // NOLINT(readability-identifier-length) BLAS naming
->decltype(herk(filling::lower, alpha, a, herk(filling::upper, alpha, a, std::forward<C2D>(c)))) {
	return herk(filling::lower, alpha, a, herk(filling::upper, alpha, a, std::forward<C2D>(c))); }

template<class A2D, class C2D>
auto herk(A2D const& a, C2D&& c)  // NOLINT(readability-identifier-length) BLAS naming
->decltype(herk(1., a, std::forward<C2D>(c))) {
	return herk(1., a, std::forward<C2D>(c)); }

template<class AA, class A2D, class Ret = typename A2D::decay_type>
[[nodiscard]]  // ("when argument is read-only")
auto herk(AA alpha, A2D const& a) {  // NOLINT(readability-identifier-length) BLAS naming
	return herk(alpha, a, Ret({size(a), size(a)}));//Ret({size(a), size(a)}));//, get_allocator(a)));
}

template<class T> struct numeric_limits : std::numeric_limits<T> {};
template<class T> struct numeric_limits<std::complex<T>> : std::numeric_limits<std::complex<T>> {
	static auto quiet_NaN() -> std::complex<T> {auto nana = numeric_limits<T>::quiet_NaN(); return {nana, nana};}  // NOLINT(readability-identifier-naming) conventional std name
};

template<class AA, class A2D, class Ret = typename A2D::decay_type>
[[nodiscard]]  // ("because argument is read-only")]]
auto herk(filling cs, AA alpha, A2D const& a)  // NOLINT(readability-identifier-length) BLAS naming
->std::decay_t<
decltype(  herk(cs, alpha, a, Ret({size(a), size(a)}, 0., get_allocator(a))))> {
	return herk(cs, alpha, a, Ret({size(a), size(a)},
#ifdef NDEBUG
		numeric_limits<typename Ret::element_type>::quiet_NaN(),
#endif
		get_allocator(a)
	));
}

template<class A2D> auto herk(filling s, A2D const& a)  // NOLINT(readability-identifier-length) BLAS naming
->decltype(herk(s, 1., a)) {
	return herk(s, 1., a); }

template<class A2D> auto herk(A2D const& array) {
		return herk(1., array);
}

}  // end namespace boost::multi::blas
#endif
