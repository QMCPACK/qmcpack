// Copyright 2019-2024 Alfredo A. Correa

#ifndef BOOST_MULTI_ADAPTORS_BLAS_HERK_HPP
#define BOOST_MULTI_ADAPTORS_BLAS_HERK_HPP
#pragma once

#include <boost/multi/adaptors/blas/copy.hpp>
#include <boost/multi/adaptors/blas/core.hpp>
#include <boost/multi/adaptors/blas/filling.hpp>
#include <boost/multi/adaptors/blas/operations.hpp>
#include <boost/multi/adaptors/blas/side.hpp>
#include <boost/multi/adaptors/blas/syrk.hpp>  // fallback to real case

// IWYU pragma: no_include "boost/multi/adaptors/blas/traits.hpp"      // for blas

namespace boost::multi::blas {

template<class A,
	std::enable_if_t<! is_conjugated<A>{}, int> =0>   // NOLINT(modernize-use-constraints) TODO(correaa) for C++20
auto base_aux(A&& array)
->decltype((std::forward<A>(array)).base()) {
	return (std::forward<A>(array)).base(); }

template<class A,
	std::enable_if_t<    is_conjugated<A>{}, int> =0>  // NOLINT(modernize-use-constraints) TODO(correaa) for C++20
auto base_aux(A&& array)
->decltype(underlying((std::forward<A>(array)).base())) {
	return underlying((std::forward<A>(array)).base()); }

using core::herk;

template<class AA, class BB, class A2DCursor, class Size, class C2DCursor>
auto herk_nm(filling c_side, AA alpha, A2DCursor a_home, Size n, Size m, BB beta, C2DCursor c_home) {  // NOLINT(readability-function-cognitive-complexity,readability-identifier-length) BLAS naming
	auto a_base = base_aux(a_home);  // NOLINT(llvm-qualified-auto,readability-qualified-auto) TODO(correaa)
	auto c_base = base_aux(c_home);  // NOLINT(llvm-qualified-auto,readability-qualified-auto) TODO(correaa)

	if constexpr(is_conjugated<A2DCursor>::value) {
		assert(0);  // TODO(correaa) implement
	} else {
		if     (a_home.stride()!=1 && c_home.stride()!=1) { 
														    herk(c_side==filling::upper?'L':'U', 'C', n, m, &alpha, a_base, a_home.template stride<0>(), &beta, c_base, c_home.template stride<0>() );
		} else if(a_home.stride()!=1 && c_home.stride()==1) {
			if(n==1)                                      { herk(c_side==filling::upper?'L':'U', 'N', n, m, &alpha, a_base, a_home.template stride<1>(), &beta, c_base, c_home.template stride<1>() ); }
			else                                          { assert(0);} // NOLINT(cppcoreguidelines-pro-bounds-array-to-pointer-decay,hicpp-no-array-decay)
		}
		else if(a_home.stride()==1 && c_home.stride()!=1) { assert(0);} // NOLINT(cppcoreguidelines-pro-bounds-array-to-pointer-decay,hicpp-no-array-decay)
		else if(a_home.stride()==1 && c_home.stride()==1) { herk(c_side==filling::upper?'U':'L', 'N', n, m, &alpha, a_base, a_home.template stride<1>(), &beta, c_base, c_home.template stride<1>() );}
	}

	return c_home;
}

template<class Scalar, class A2DCursor, class Size, class DecayType>
class herk_range {
	// ContextPtr ctxtp_;
	multi::blas::filling cs_;
	Scalar scale_;
	A2DCursor a_home_;
	Size a_rows_;
	Size a_cols_;

 public:
	herk_range(herk_range const&) = delete;
	herk_range(herk_range&&) = delete;
	auto operator=(herk_range const&) -> herk_range& = delete;
	auto operator=(herk_range&&) -> herk_range& = delete;
	~herk_range() = default;

	herk_range(multi::blas::filling cs, Scalar scale, A2DCursor a_home, Size a_rows, Size a_cols) : cs_{cs}, scale_{scale}, a_home_{a_home}, a_rows_{a_rows}, a_cols_{a_cols} {}  // NOLINT(bugprone-easily-swappable-parameters) TODO(correaa)

//  herk_range(ContextPtr ctxtp, Scalar s, A2D&& a)  // NOLINT(bugprone-easily-swappable-parameters,readability-identifier-length) BLAS naming
//  : ctxtp_{ctxtp}
//  , s_{s}, a_{std::forward<2D>(a_first)}, a_end_{std::move(a_last)}
//  {}

	struct iterator {
		herk_range const* self_;
		Size index_;
	};

//  // using iterator = herk_iterator<ContextPtr, Scalar, ItA>;
//  using decay_type = DecayType;
//  using size_type = typename decay_type::size_type;

	auto begin() const { return iterator{this, 0}; }
	auto end() const { return iterator{this, a_rows_}; }

	auto size() const { return a_rows_;}
	auto extensions() const { return multi::extensions_t<2>{a_rows_, a_rows_}; }

	// template<class ItOut>
	// friend auto copy(iterator const& first, iterator const& last, ItOut d_first) {
	//  auto const& a = multi::ref(first.self_, last);
	//  blas::herk_nm(first.self_->cs_, a.home(), a.size(), (~a).size(), 0.0, d_first);  // NOLINT(fuchsia-default-arguments-calls)
	// }

//  // auto operator+() const -> decay_type {return *this;} // TODO(correaa) : investigate why return decay_type{*this} doesn't work
//  // template<class Arr>
//  // friend auto operator+=(Arr&& a, gemm_range const& self) -> Arr&& {  // NOLINT(readability-identifier-length) BLAS naming
//  //  blas::gemm_n(self.ctxtp_, self.s_, self.a_begin_, self.a_end_ - self.a_begin_, self.b_begin_, 1., a.begin());
//  //  return std::forward<Arr>(a);
//  // }
//  // friend auto operator*(Scalar factor, gemm_range const& self) {
//  //  return gemm_range{self.ctxtp_, factor*self.s_, self.a_begin_, self.a_end_, self.b_begin_};
//  // }
};

template<class AA, class BB, class A2D, class C2D, class = typename A2D::element_ptr,
	std::enable_if_t<is_complex_array<C2D>{}, int> =0>  // NOLINT(modernize-use-constraints) TODO(correaa) for C++20
auto herk(filling c_side, AA alpha, A2D const& a, BB beta, C2D&& c) -> C2D&& {  // NOLINT(readability-function-cognitive-complexity,readability-identifier-length) 74, BLAS naming
	assert( a.size() == c.size() ); // NOLINT(cppcoreguidelines-pro-bounds-array-to-pointer-decay,hicpp-no-array-decay)
	assert( c.size() == c.rotated().size() ); // NOLINT(cppcoreguidelines-pro-bounds-array-to-pointer-decay,hicpp-no-array-decay)
	if(c.is_empty()) {return std::forward<C2D>(c);}
	if constexpr(is_conjugated<C2D>{}) {
		herk(flip(c_side), alpha, a, beta, hermitized(c));
		return std::forward<C2D>(c);
	}

	auto base_a = base_aux(a);  // NOLINT(llvm-qualified-auto,readability-qualified-auto) TODO(correaa)
	auto base_c = base_aux(c);  // NOLINT(llvm-qualified-auto,readability-qualified-auto) TODO(correaa)
	if constexpr(is_conjugated<A2D>{}) {
	//  auto& ctxt = *blas::default_context_of(underlying(a.base()));
		// if you get an error here might be due to lack of inclusion of a header file with the backend appropriate for your type of iterator
		if     (stride(a)==1 && c.stride()!=1) {herk(c_side==filling::upper?'L':'U', 'N', c.size(), a.rotated().size(), &alpha, base_a, a.rotated().stride(), &beta, base_c, c.stride());}
		else if(stride(a)==1 && c.stride()==1) {
			if(size(a)==1)                     {herk(c_side==filling::upper?'L':'U', 'N', c.size(), a.rotated().size(), &alpha, base_a, a.rotated().stride(), &beta, base_c, c.stride());}
			else                               {assert(0);} // NOLINT(cppcoreguidelines-pro-bounds-array-to-pointer-decay,hicpp-no-array-decay)
		}
		else if(stride(a)!=1 && c.stride()==1) { herk(c_side==filling::upper?'U':'L', 'C', c.size(), a.rotated().size(), &alpha, base_a, stride(        a ), &beta, base_c, c.rotated().stride());}
		else if(stride(a)!=1 && c.stride()!=1) { herk(c_side==filling::upper?'L':'U', 'C', c.size(), a.rotated().size(), &alpha, base_a, stride(        a ), &beta, base_c, c          .stride());}
		else                                  { assert(0);} // NOLINT(cppcoreguidelines-pro-bounds-array-to-pointer-decay,hicpp-no-array-decay)
	} else {
	//  auto& ctxt = *blas::default_context_of(           a.base() );
		if     (stride(a)!=1 && c.stride()!=1) { herk(c_side==filling::upper?'L':'U', 'C', c.size(), a.rotated().size(), &alpha, base_a, stride(        a ), &beta, base_c, c.stride());}
		else if(stride(a)!=1 && c.stride()==1) {
			if(size(a)==1)                    { herk(c_side==filling::upper?'L':'U', 'N', c.size(), a.rotated().size(), &alpha, base_a, a.rotated().stride(), &beta, base_c, c.rotated().stride());}
			else                              { assert(0);} // NOLINT(cppcoreguidelines-pro-bounds-array-to-pointer-decay,hicpp-no-array-decay)
		}
		else if(stride(a)==1 && c.stride()!=1) {assert(0);} // NOLINT(cppcoreguidelines-pro-bounds-array-to-pointer-decay,hicpp-no-array-decay)
		else if(stride(a)==1 && c.stride()==1) {herk(c_side==filling::upper?'U':'L', 'N', c.size(), a.rotated().size(), &alpha, base_a, a.rotated().stride(), &beta, base_c, c.rotated().stride());}
	//  else                                   {assert(0);} // NOLINT(cppcoreguidelines-pro-bounds-array-to-pointer-decay,hicpp-no-array-decay)
	}

	return std::forward<C2D>(c);
}

template<class AA, class BB, class A2D, class C2D, class = typename A2D::element_ptr,
	std::enable_if_t<! is_complex_array<C2D>{}, int> =0>  // NOLINT(modernize-use-constraints) TODO(correaa) for C++20
auto herk(filling c_side, AA alpha, A2D const& a, BB beta, C2D&& C)  // NOLINT(readability-identifier-length) BLAS naming
->decltype(syrk(c_side, alpha, a, beta, std::forward<C2D>(C))) {
	return syrk(c_side, alpha, a, beta, std::forward<C2D>(C)); }

template<class AA, class A2D, class C2D, class = typename A2D::element_ptr>
auto herk(filling c_side, AA alpha, A2D const& a, C2D&& C)  // NOLINT(readability-identifier-length) BLAS naming
->decltype(herk(c_side, alpha, a, 0., std::forward<C2D>(C))) {
	return herk(c_side, alpha, a, 0., std::forward<C2D>(C)); }

template<typename AA, class A2D, class C2D>
auto herk(AA alpha, A2D const& a, C2D&& C)  // NOLINT(readability-identifier-length) BLAS naming
->decltype(herk(filling::lower, alpha, a, herk(filling::upper, alpha, a, std::forward<C2D>(C)))) {
	return herk(filling::lower, alpha, a, herk(filling::upper, alpha, a, std::forward<C2D>(C))); }

template<class A2D, class C2D>
auto herk(A2D const& A, C2D&& C)  // NOLINT(readability-identifier-length) BLAS naming
->decltype(herk(1.0, A, std::forward<C2D>(C))) {
	return herk(1.0, A, std::forward<C2D>(C)); }

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
{
	return herk_range<AA, typename A2D::const_cursor, typename A2D::size_type, Ret>(cs, alpha, a.home(), std::get<0>(a.sizes()), std::get<1>(a.sizes()) );
}
// ->std::decay_t<
// decltype(  herk(cs, alpha, a, Ret({size(a), size(a)}, 0.0, get_allocator(a))))> {
//  return herk(cs, alpha, a, Ret({size(a), size(a)},
// #ifdef NDEBUG
//      numeric_limits<typename Ret::element_type>::quiet_NaN(),
// #endif
//      get_allocator(a)
//  ));
// }

template<class A2D> auto herk(filling f, A2D const& a)  // NOLINT(readability-identifier-length) BLAS naming
->decltype(herk(f, 1.0, a)) {
	return herk(f, 1.0, a); }

template<class A2D> auto herk(A2D const& array) {
		return herk(1.0, array);
}

}  // end namespace boost::multi::blas
#endif
