// Copyright 2019-2024 Alfredo A. Correa
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#ifndef BOOST_MULTI_ADAPTORS_BLAS_DOT_HPP
#define BOOST_MULTI_ADAPTORS_BLAS_DOT_HPP

#include <boost/multi/adaptors/blas/core.hpp>
#include <boost/multi/adaptors/blas/numeric.hpp>     // for is_complex
#include <boost/multi/adaptors/blas/operations.hpp>  // for blas::C

namespace boost::multi::blas {

using core::dot;
using core::dotu;
using core::dotc;

template<class Context, class XIt, class Size, class YIt, class RPtr>
auto dot_n(Context&& ctxt, XIt x_first, Size count, YIt y_first, RPtr rp) {
	if constexpr(! is_complex<typename XIt::value_type>{}) {
                                                                           std::forward<Context>(ctxt)->dot (count,            x_first.base() , stride(x_first), y_first.base(), stride(y_first), rp);
	} else {
		if      constexpr(!is_conjugated<XIt>{} && !is_conjugated<YIt>{}) {std::forward<Context>(ctxt)->dotu(count,            x_first.base() , stride(x_first), y_first.base(), stride(y_first), rp);}
		else if constexpr(!is_conjugated<XIt>{} &&  is_conjugated<YIt>{}) {std::forward<Context>(ctxt)->dotc(count, underlying(y_first.base()), stride(y_first), x_first.base(), stride(x_first), rp);}
		else if constexpr( is_conjugated<XIt>{} && !is_conjugated<YIt>{}) {std::forward<Context>(ctxt)->dotc(count, underlying(x_first.base()), stride(x_first), y_first.base(), stride(y_first), rp);}
		else if constexpr( is_conjugated<XIt>{} &&  is_conjugated<YIt>{}) {static_assert(!sizeof(XIt*), "not implemented in blas");}
	}

	struct{XIt x_last; YIt y_last;} ret{x_first + count, y_first + count};
	return ret;
}

template<class XIt, class Size, class YIt, class RPtr>
auto dot_n(XIt x_first, Size count, YIt y_first, RPtr rp) {  //->decltype(dot_n(blas::context{}, x_first, count, y_first, rp)){
	if constexpr(is_conjugated<XIt>{}) {
		auto ctxtp = blas::default_context_of(underlying(x_first.base()));
		return dot_n(ctxtp, x_first, count, y_first, rp);
	} else {
		auto ctxtp = blas::default_context_of(x_first.base());
		return dot_n(ctxtp, x_first, count, y_first, rp);
	}
}

template<class Context, class X1D, class Y1D, class R,
	std::enable_if_t<! multi::has_base<std::decay_t<R>>::value, int> =0>  // NOLINT(modernize-use-constraints) TODO(correaa) for C++20
auto dot(Context ctxt, X1D const& x, Y1D const& y, R&& res) -> R&& {  // NOLINT(readability-identifier-length) res = \sum_i x_i y_i
	assert( size(x) == size(y) ); // NOLINT(cppcoreguidelines-pro-bounds-array-to-pointer-decay,hicpp-no-array-decay)
	return blas::dot_n(ctxt, begin(x), size(x), begin(y), &res), std::forward<R>(res);
}

template<class Context, class X1D, class Y1D, class R,
	std::enable_if_t<multi::has_base<std::decay_t<R>>::value, int> =0>  // NOLINT(modernize-use-constraints) TODO(correaa) for C++20
auto dot(Context ctxt, X1D const& x, Y1D const& y, R&& res) -> R&& {  // NOLINT(readability-identifier-length) res = \sum_i x_i y_i
	assert( size(x) == size(y) ); // NOLINT(cppcoreguidelines-pro-bounds-array-to-pointer-decay,hicpp-no-array-decay)
	return blas::dot_n(ctxt, begin(x), size(x), begin(y), res.base()), std::forward<R>(res);
}

template<class X1D, class Y1D, class R>
auto dot(X1D const& x, Y1D const& y, R&& res) -> R&& {  // NOLINT(readability-identifier-length) res = \sum_i x_i y_i
	assert( size(x) == size(y) ); // NOLINT(cppcoreguidelines-pro-bounds-array-to-pointer-decay,hicpp-no-array-decay)
	if constexpr(is_conjugated<X1D>{}) {
		auto ctxtp = blas::default_context_of(underlying(x.base()));
		return blas::dot(ctxtp, x, y, std::forward<R>(res));
	} else {
		auto ctxtp = blas::default_context_of(           x.base() );
		return blas::dot(ctxtp, x, y, std::forward<R>(res));
	}
}

template<class ContextPtr, class ItX, class Size, class ItY>
class dot_ptr {
	ContextPtr ctxt_;
	ItX  x_first_;
	Size count_;
	ItY  y_first_;

 protected:
	dot_ptr(ContextPtr ctxt, ItX x_first, Size count, ItY y_first) : ctxt_{ctxt}, x_first_{x_first}, count_{count}, y_first_{y_first} {}

 public:
	constexpr explicit operator bool() const {return true;}

	template<class ItOut, class Size2>
	friend constexpr auto copy_n(dot_ptr first, Size2 count, ItOut d_first)
	->decltype(blas::dot_n(std::declval<ContextPtr>(), std::declval<ItX>(), Size{}      , std::declval<ItY>(), d_first), d_first + count) {
		assert(count == 1); // NOLINT(cppcoreguidelines-pro-bounds-array-to-pointer-decay,hicpp-no-array-decay)
		return blas::dot_n(first.ctxt_               , first.x_first_     , first.count_, first.y_first_     , d_first), d_first + count;
	}

	template<class ItOut, class Size2>
	friend constexpr auto uninitialized_copy_n(dot_ptr first, Size2 count, ItOut d_first)
	->decltype(blas::dot_n(std::declval<ContextPtr>(), std::declval<ItX>(), Size{}      , std::declval<ItY>(), d_first), d_first + count) {assert(count == 1);
		return blas::dot_n(first.ctxt_               , first.x_first_     , first.count_, first.y_first_     , d_first), d_first + count; }
};

template<class ContextPtr, class X, class Y, class Ptr = dot_ptr<ContextPtr, typename X::const_iterator, typename X::size_type, typename Y::const_iterator>>
struct dot_ref : private Ptr {
	using decay_type = decltype(std::declval<typename X::value_type>()*std::declval<typename Y::value_type>());
	dot_ref(ContextPtr ctxt, X const& x, Y const& y) : Ptr{ctxt, begin(x), size(x), begin(y)} {  // NOLINT(readability-identifier-length) BLAS naming
		assert(( size(x) == size(y) ));  // NOLINT(cppcoreguidelines-pro-bounds-array-to-pointer-decay,hicpp-no-array-decay)
	}

	constexpr auto operator&() const& -> Ptr const& {return *this;}  // NOLINT(google-runtime-operator) reference type  // NOSONAR

	auto decay() const -> decay_type {decay_type ret; copy_n(operator&(), 1, &ret); return ret;}  // NOLINT(fuchsia-default-arguments-calls)
	operator decay_type()       const {return decay();}  // NOLINT(google-explicit-constructor,hicpp-explicit-conversions,hicpp-explicit-conversion) //NOSONAR to allow terse syntax
#if ! defined(__CUDACC__) || ! defined(__INTEL_COMPILER)
	friend auto operator*(decay_type const& lhs, dot_ref const& self) {return lhs*self.decay();}
#endif
	auto operator+() const -> decay_type {return decay();}

	// friend auto operator==(dot_ref const& self, dot_ref const& other) -> bool {return self.decay() == other.decay();}
	// friend auto operator!=(dot_ref const& self, dot_ref const& other) -> bool {return self.decay() != other.decay();}

	friend auto operator==(dot_ref const& self, dot_ref const& other) -> bool {
		return self.decay() == other.decay();
	}
	friend auto operator!=(dot_ref const& self, dot_ref const& other) -> bool {
		return self.decay() != other.decay();
	}

	template<class Other>  // , class = std::enable_if_t<!std::is_same_v<dot_ref, Other> > >
	friend auto operator==(dot_ref const& self, Other const& other) -> bool {
		return self.decay() == other;
	}
	template<class Other>  // , class = std::enable_if_t<!std::is_same_v<dot_ref, Other> > >
	friend auto operator!=(dot_ref const& self, Other const& other) -> bool {
		return self.decay() != other;
	}
	template<class Other,
		class = std::enable_if_t<!std::is_same_v<dot_ref, Other> > >  // NOLINT(modernize-use-constraints) TODO(correaa) for C++20
	friend auto operator==(Other const& other, dot_ref const& self) -> bool {
		return other == self.decay();
	}
	template<class Other,
		class = std::enable_if_t<!std::is_same_v<dot_ref, Other> > >  // NOLINT(modernize-use-constraints) TODO(correaa) for C++20
	friend auto operator!=(Other const& other, dot_ref const& self) -> bool {
		return other != self.decay();
	}

};

template<class Context, class X, class Y> [[nodiscard]]
auto dot(Context ctxt, X const& x, Y const& y) {  // NOLINT(readability-identifier-length) BLAS naming
	return dot_ref<Context, X, Y>{ctxt, x, y};
}

#if defined __NVCC__   // in place of global -Xcudafe \"--diag_suppress=implicit_return_from_non_void_function\"
	#ifdef __NVCC_DIAG_PRAGMA_SUPPORT__
		#pragma nv_diagnostic push
		#pragma nv_diag_suppress = implicit_return_from_non_void_function
	#else
		#pragma    diagnostic push
		#pragma    diag_suppress = implicit_return_from_non_void_function
	#endif
#elif defined __NVCOMPILER
	#pragma    diagnostic push
	#pragma    diag_suppress = implicit_return_from_non_void_function
#endif

template<class X, class Y> [[nodiscard]]
auto dot(X const& x, Y const& y) {  // NOLINT(readability-identifier-length) BLAS naming
	if constexpr(is_conjugated<X>{}) {
		auto ctxtp = blas::default_context_of(underlying(x.base()));
		return blas::dot(ctxtp, x, y);
	} else {
		auto ctxtp = blas::default_context_of(x.base());
		return blas::dot(ctxtp, x, y);
	}
}
#if defined __NVCC__
	#ifdef __NVCC_DIAG_PRAGMA_SUPPORT__
		#pragma nv_diagnostic pop
	#else
		#pragma    diagnostic pop
	#endif
#elif defined __NVCOMPILER
	#pragma    diagnostic pop
#endif

namespace operators {
	template<class X1D, class Y1D> [[nodiscard]]
	auto operator,(X1D const& x, Y1D const& y)  // NOLINT(readability-identifier-length) BLAS naming
	->decltype(dot(x, y)) {
		return dot(x, y); }
}  // end namespace operators

}  // end namespace boost::multi::blas

#endif
