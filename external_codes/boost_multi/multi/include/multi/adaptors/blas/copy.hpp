 // -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;autowrap:nil;-*-
// Copyright 2020-2022 Alfredo A. Correa

#ifndef MULTI_ADAPTORS_BLAS_COPY_HPP
#define MULTI_ADAPTORS_BLAS_COPY_HPP

#include "../blas/core.hpp"
#include "../blas/operations.hpp"

#include "../../config/NODISCARD.hpp"

#include<type_traits>

namespace boost::multi::blas {

using core::copy;

template<class It, typename Size, class OutIt>
auto copy_n(It first, Size n, OutIt d_first)
->decltype(copy(n, first.base(), first.stride(), d_first.base(), d_first.stride()), d_first + n) {
	return copy(n, first.base(), first.stride(), d_first.base(), d_first.stride()), d_first + n; }

template<class Context, class It, typename Size, class OutIt, class=std::enable_if_t<blas::is_context<Context>{}>>
auto copy_n(Context&& ctxt, It first, Size n, OutIt d_first)
->decltype(copy(std::forward<Context>(ctxt), n, first.base(), first.stride(), d_first.base(), d_first.stride()), d_first + n) {
	return copy(std::forward<Context>(ctxt), n, first.base(), first.stride(), d_first.base(), d_first.stride()), d_first + n; }

template<class It, class OutIt>
auto copy(It first, It last, OutIt d_first)
->decltype(copy_n(first, last - first, d_first)) {
	return copy_n(first, last - first, d_first); }

template<class Context, class It, class OutIt, class=std::enable_if_t<blas::is_context<Context>{}>>
auto copy(Context&& ctxt, It first, It last, OutIt d_first)
->decltype(copy_n(std::forward<Context>(ctxt), first, last - first, d_first)) {
	return copy_n(std::forward<Context>(ctxt), first, last - first, d_first); }

template<class X1D, class Y1D>
auto copy(X1D const& x, Y1D&& y)  // NOLINT(readability-identifier-length) BLAS naming
->decltype(blas::copy_n(x.begin(), size(x), y.begin()), std::forward<Y1D>(y)) {
	assert( (x.size() == y.size()) );  // NOLINT(cppcoreguidelines-pro-bounds-array-to-pointer-decay,hicpp-no-array-decay) : assert
	return blas::copy_n(x.begin(), x.size(), y.begin()), std::forward<Y1D>(y); }

template<class Context, class X1D, class Y1D>
auto copy(Context&& ctxt, X1D const& x, Y1D&& y)  // NOLINT(readability-identifier-length) BLAS naming
->decltype(blas::copy_n(std::forward<Context>(ctxt), x.begin(), size(x), y.begin()), std::forward<Y1D>(y)) {
	assert(x.size()==y.size());  // NOLINT(cppcoreguidelines-pro-bounds-array-to-pointer-decay,hicpp-no-array-decay) : normal in a constexpr assert
	return blas::copy_n(std::forward<Context>(ctxt), x.begin(), x.size(), y.begin()), std::forward<Y1D>(y);
}

template<class ContextPtr, class It1D>
class copy_iterator{
	ContextPtr ctxt = {};
	It1D it_;
public:
	using difference_type = typename std::iterator_traits<It1D>::difference_type;
	using value_type 	  = typename std::iterator_traits<It1D>::value_type;
	using pointer 	      = void;
	using reference 	  = void;
	using iterator_category = std::output_iterator_tag;
	using iterator_type   = It1D;
	using context_type    = ContextPtr;
	constexpr explicit copy_iterator(It1D it) : it_{it}{}
	constexpr          copy_iterator(ContextPtr cp, It1D it) : ctxt{cp}, it_{it}{}
	constexpr auto base() const -> iterator_type{return it_;}
	template<class It1DOut>
	friend constexpr auto copy_n(copy_iterator first, difference_type count, It1DOut result) -> It1DOut{
		return blas::copy_n(first.ctxt, first.base(), count, result);
	}
	template<class It1DOut> 
	friend constexpr auto copy(copy_iterator first, copy_iterator last, It1DOut d_first) -> It1DOut{
		return copy_n(first, distance(first, last), d_first);
	}
	template<class It1DOut>
	friend constexpr auto uninitialized_copy(copy_iterator first, copy_iterator last, It1DOut d_first) -> It1DOut{
		return copy_n(first, distance(first, last), d_first);
	}
	friend constexpr auto distance(copy_iterator const& self, copy_iterator const& other) -> difference_type{
		assert(stride(other.it_) == stride(self.it_)); // NOLINT(cppcoreguidelines-pro-bounds-array-to-pointer-decay,hicpp-no-array-decay) : normal in a constexpr assert
		return other.it_ - self.it_;
	}
	constexpr auto operator*() const -> value_type{return *it_;}
};

template<class ContextPtr, class It1D, class DecayType = void, class DiffType = typename std::iterator_traits<It1D>::difference_type>
class copy_range {
	ContextPtr ctxp_ = {};
	It1D begin_, end_;

 public:
	using difference_type = DiffType;
	using iterator        = copy_iterator<ContextPtr, It1D>;
	using decay_type      = DecayType;

	copy_range(copy_range const&) = delete;
	copy_range(copy_range&&) noexcept = default;
	auto operator=(copy_range const&) -> copy_range& = delete;
	auto operator=(copy_range&&) -> copy_range& = delete;
	~copy_range() = default;

	constexpr copy_range(It1D first, It1D last) : begin_{first}, end_{last} {}
	constexpr copy_range(ContextPtr ctxp, It1D first, It1D last) : ctxp_{ctxp}, begin_{first}, end_{last} {}
	constexpr auto size()  const -> difference_type {return end_ - begin_;}
	constexpr auto begin() const {return iterator{ctxp_, begin_};}
	constexpr auto end()   const {return iterator{ctxp_, end_  };}
	constexpr auto extensions() const -> typename decay_type::extensions_type {return {multi::iextension{size()}};}
	template<class Other, class=decltype(Other(std::declval<iterator>(), std::declval<iterator>()))>
	operator Other() const{return Other(begin(), end());} // NOLINT(google-explicit-constructor,hicpp-explicit-conversions) : to allow terse syntax
	friend auto operator+(copy_range const& self) {return self.operator decay_type();}
};

template<class DecayType, class It> NODISCARD()
auto copy(It const& first, It const& last)
->decltype(copy_range<void*, It, DecayType>{first, last}){
	return copy_range<void*, It, DecayType>{first, last};}

template<class DecayType, class Context, class It> NODISCARD()
auto copy(Context&& ctxt, It const& first, It const& last)
->decltype(copy_range<Context, It, DecayType>{ctxt, first, last}){
	return copy_range<Context, It, DecayType>{ctxt, first, last};}

template<class A> NODISCARD()
auto copy(A const& array)  // need to specify templates (instead of deduced for intel)
->decltype(copy<typename A::decay_type, typename A::const_iterator>(array.begin(), array.end())){
	return copy<typename A::decay_type, typename A::const_iterator>(array.begin(), array.end());}

template<class Context, class A, class=std::enable_if_t<blas::is_context<Context>{}>> NODISCARD()
auto copy(Context&& ctxt, A const& a)  // NOLINT(readability-identifier-length) conventional name in BLAS
->decltype(copy<typename A::decay_type, Context, typename A::const_iterator>(std::forward<Context>(ctxt), a.begin(), a.end())){
	return copy<typename A::decay_type, Context, typename A::const_iterator>(std::forward<Context>(ctxt), a.begin(), a.end());}

} // end namespace boost::multi::blas

#endif
