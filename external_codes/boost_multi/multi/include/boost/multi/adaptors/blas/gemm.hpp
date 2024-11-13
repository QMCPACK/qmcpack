// Copyright 2019-2024 Alfredo A. Correa
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#ifndef BOOST_MULTI_ADAPTORS_BLAS_GEMM_HPP
#define BOOST_MULTI_ADAPTORS_BLAS_GEMM_HPP

#include <boost/multi/adaptors/blas/core.hpp>
// #include <boost/multi/adaptors/blas/gemv.hpp>
#include <boost/multi/adaptors/blas/numeric.hpp>
// #include <boost/multi/adaptors/blas/operations.hpp>

#include <boost/multi/array_ref.hpp>              // for base, size, begin

#include <cassert>                               // for assert
#include <cstddef>                                // for nullptr_t
#include <exception>                              // for exception
#include <iterator>                               // for iterator_traits
#include <stdexcept>                              // for logic_error
#include <string>                                 // for to_string, operator""s
#include <type_traits>                            // for enable_if_t, integr...
#include <utility>                                // for forward, declval

namespace boost::multi::blas {

using core::gemm;

template<class It>
auto xbase_aux(It const& it, std::true_type const& /*true */)
->decltype(underlying(base(it))) {
	return underlying(base(it)); }

template<class It>
auto xbase_aux(It const& it, std::false_type const& /*false*/)
->decltype(base(it)) {
	return base(it); }

template<class It>
auto xbase(It const& it)
->decltype(xbase_aux(it, std::integral_constant<bool, is_conjugated<It>{}>{})) {
	return xbase_aux(it, std::integral_constant<bool, is_conjugated<It>{}>{}); }

#define CTXT std::forward<Context>(ctxt)

template<class Context, class It2DA, class Size, class It2DB, class It2DC,
	std::enable_if_t<(!is_conjugated<It2DA>{} && !is_conjugated<It2DB>{}), int> = 0>  // NOLINT(modernize-use-constraints) TODO(correaa) for C++20
auto gemm_n(Context&& ctxt, typename It2DA::element alpha, It2DA a_first, Size a_count, It2DB b_first, typename It2DA::element beta, It2DC c_first) // NOLINT(readability-function-cognitive-complexity) : 125
{
	assert( (*b_first).size() == (*c_first).size() );          // NOLINT(cppcoreguidelines-pro-bounds-array-to-pointer-decay,hicpp-no-array-decay)
	assert( a_first.stride()==1 || (*a_first).stride()==1 ); // NOLINT(cppcoreguidelines-pro-bounds-array-to-pointer-decay,hicpp-no-array-decay)
	assert( b_first.stride()==1 || (*b_first).stride()==1 ); // NOLINT(cppcoreguidelines-pro-bounds-array-to-pointer-decay,hicpp-no-array-decay)
	assert( c_first.stride()==1 || (*c_first).stride()==1 ); // NOLINT(cppcoreguidelines-pro-bounds-array-to-pointer-decay,hicpp-no-array-decay)

	if(a_count == 0) { return c_first; }

	if      ((*a_first).stride()==1 && (*b_first).stride()==1 && (*c_first).stride()==1) {
		if     ( a_count==1 && (*b_first).size()==1 ) {CTXT->gemm('N', 'N', (*b_first).size(), a_count, (*a_first).size(), &alpha, b_first.base(), (*b_first).size(), a_first.base(), (*a_first).size()  , &beta, c_first.base(), (*c_first).size()  );}
		else if( a_count==1                        ) {CTXT->gemm('N', 'N', (*b_first).size(), a_count, (*a_first).size(), &alpha, b_first.base(), b_first. stride(), a_first.base(), (*a_first).size()  , &beta, c_first.base(), (*c_first).size()  );}
		else                                         {CTXT->gemm('N', 'N', (*b_first).size(), a_count, (*a_first).size(), &alpha, b_first.base(), b_first. stride(), a_first.base(), a_first. stride(), &beta, c_first.base(), c_first. stride());}
	}else if((*a_first).stride()==1 && (*b_first).stride()==1 && c_first. stride()==1) {
		if  (a_count==1)                            {CTXT->gemm('T', 'T', a_count, (*b_first).size(), (*a_first).size(), &alpha, a_first.base(), a_first. stride(), b_first.base(), (*b_first).size()  , &beta, c_first.base(), (*a_first).size()  );}
		else                                        {CTXT->gemm('T', 'T', a_count, (*b_first).size(), (*a_first).size(), &alpha, a_first.base(), a_first. stride(), b_first.base(), b_first. stride(), &beta, c_first.base(), (*c_first).stride());}
	}else if(a_first. stride()==1 && (*b_first).stride()==1 && (*c_first).stride()==1) { 
		if  (a_count==1)                            {CTXT->gemm('N', 'T', (*c_first).size(), a_count, (*a_first).size(), &alpha, b_first.base(), b_first. stride(), a_first.base(), (*a_first).stride(), &beta, c_first.base(), a_count         );}
		else                                        {CTXT->gemm('N', 'T', (*c_first).size(), a_count, (*a_first).size(), &alpha, b_first.base(), b_first. stride(), a_first.base(), (*a_first).stride(), &beta, c_first.base(), c_first.stride());}
	}else if(a_first. stride()==1 && (*b_first).stride()==1 && c_first. stride()==1) {
		if  (a_count==1)                            {CTXT->gemm('N', 'T', a_count, (*b_first).size(), (*a_first).size(), &alpha, a_first.base(), (*a_first).stride(), b_first.base(), (*a_first).size()  , &beta, c_first.base(), (*b_first).size()  );}
		else                                        {CTXT->gemm('N', 'T', a_count, (*b_first).size(), (*a_first).size(), &alpha, a_first.base(), (*a_first).stride(), b_first.base(), b_first. stride(), &beta, c_first.base(), (*c_first).stride());}
	}else if((*a_first).stride()==1 && b_first.stride()==1 && c_first. stride()==1) {
		if     (a_count==1 && (*b_first).size()==1)  {CTXT->gemm('N', 'N', (*c_first).size(), a_count, (*a_first).size(), &alpha, b_first.base(), (*b_first).size()  , a_first.base(), (*a_first).size()  , &beta, c_first.base(), (*c_first).stride());}
		else if(a_count==1)                         {CTXT->gemm('N', 'T', (*c_first).size(), a_count, (*a_first).size(), &alpha, b_first.base(), (*b_first).stride(), a_first.base(), (*a_first).size()  , &beta, c_first.base(), (*c_first).stride());}
		else if((*a_first).size() == 1 && (*b_first).size() == 1)
		                                            {CTXT->gemm('N', 'N', (*c_first).size(), a_count, (*a_first).size(), &alpha, b_first.base(), (*b_first).stride(), a_first.base(), a_first. stride(), &beta, c_first.base(), (*c_first).stride());}
		else                                        {CTXT->gemm('N', 'T', (*c_first).size(), a_count, (*a_first).size(), &alpha, b_first.base(), (*b_first).stride(), a_first.base(), a_first. stride(), &beta, c_first.base(), (*c_first).stride());}
	}else if((*a_first).stride()==1 && b_first. stride()==1 && (*c_first).stride()==1) {
		if  (a_count==1)                            {CTXT->gemm('T', 'N', a_count, (*c_first).size(), (*a_first).size(), &alpha, b_first.base(), (*b_first).stride(), a_first.base(), (*a_first).size(), &beta, c_first.base(), c_first. stride());}
		else                                        {CTXT->gemm('T', 'N', (*c_first).size(), a_count, (*a_first).size(), &alpha, b_first.base(), (*b_first).stride(), a_first.base(), a_first. stride(), &beta, c_first.base(), c_first. stride());}
	}else if(a_first. stride()==1 && b_first.stride( )==1 && c_first. stride()==1) {
		if  ((*b_first).size()==1)                   {CTXT->gemm('N', 'N', a_count, (*b_first).size(), (*a_first).size(), &alpha, a_first.base(), (*a_first).stride(), b_first.base(), (*b_first).stride(), &beta, c_first.base(), a_count          );}
		else                                        {CTXT->gemm('N', 'N', a_count, (*b_first).size(), (*a_first).size(), &alpha, a_first.base(), (*a_first).stride(), b_first.base(), (*b_first).stride(), &beta, c_first.base(), (*c_first).stride());}
	}else if(a_first. stride()==1 && b_first.stride( )==1 && (*c_first).stride()==1) {          
	                                                {CTXT->gemm('T', 'T', (*b_first).size(), a_count, (*a_first).size(), &alpha, b_first.base(), (*b_first).stride(), a_first.base(), (*a_first).stride(), &beta, c_first.base(), c_first. stride());}
	} else {assert(0);}  // NOLINT(cppcoreguidelines-pro-bounds-array-to-pointer-decay,hicpp-no-array-decay)

	return c_first + a_count;
}

template<class Context, class It2DA, class Size, class It2DB, class It2DC,
	std::enable_if_t<(!is_conjugated<It2DA>{} && is_conjugated<It2DB>{}), int> =0>  // NOLINT(modernize-use-constraints) TODO(correaa) for C++20
auto gemm_n(Context&& ctxt, typename It2DA::element alpha, It2DA a_first, Size a_count, It2DB b_first, typename It2DA::element beta, It2DC c_first) // NOLINT(readability-function-cognitive-complexity) : 125
{
	assert( (*b_first).size() == (*c_first).size() );          // NOLINT(cppcoreguidelines-pro-bounds-array-to-pointer-decay,hicpp-no-array-decay)
	assert( a_first.stride()==1 || (*a_first).stride()==1 ); // NOLINT(cppcoreguidelines-pro-bounds-array-to-pointer-decay,hicpp-no-array-decay)
	assert( b_first.stride()==1 || (*b_first).stride()==1 ); // NOLINT(cppcoreguidelines-pro-bounds-array-to-pointer-decay,hicpp-no-array-decay)
	assert( c_first.stride()==1 || (*c_first).stride()==1 ); // NOLINT(cppcoreguidelines-pro-bounds-array-to-pointer-decay,hicpp-no-array-decay)

	if(a_count == 0) { return c_first; }

	if      ((*a_first).stride()==1 && (*b_first).stride()==1 && (*c_first).stride()==1) {
	                            {CTXT->gemm('C', 'N', (*c_first).size(), a_count, (*a_first).size(), &alpha, underlying(b_first.base()), (*b_first).stride(), a_first.base(), (*a_first).size()  , &beta, c_first.base(), c_first.stride());}
	}else if((*a_first).stride()==1 && b_first. stride()==1 && (*c_first).stride()==1){
		if  (a_count==1)        {CTXT->gemm('C', 'N', a_count, (*c_first).size(), (*a_first).size(), &alpha, underlying(b_first.base()), (*b_first).stride(), a_first.base(), (*a_first).size()  , &beta, c_first.base(), c_first.stride());}
		else                    {CTXT->gemm('C', 'N', (*c_first).size(), a_count, (*a_first).size(), &alpha, underlying(b_first.base()), (*b_first).stride(), a_first.base(), a_first.stride(), &beta, c_first.base(), c_first.stride());}
	}else if((*a_first).stride()==1 && b_first. stride()==1 && c_first. stride()==1){
								{CTXT->gemm('C', 'N', (*c_first).size(), a_count, (*a_first).size(), &alpha, underlying(b_first.base()), (*b_first).stride(), a_first.base(), a_first. stride(), &beta, c_first.base(), (*c_first).stride());}
	}else if(a_first. stride()==1 && b_first. stride()==1 && c_first. stride()==1){
								{CTXT->gemm('C', 'T', (*c_first).size(), a_count, (*a_first).size(), &alpha, underlying(b_first.base()), (*b_first).stride(), a_first.base(), (*a_first).stride(), &beta, c_first.base(), (*c_first).stride());}
	}else if(a_first. stride()==1 && b_first. stride()==1 && (*c_first).stride()==1){
								{CTXT->gemm('C', 'T', a_count, (*c_first).size(), (*a_first).size(), &alpha, underlying(b_first.base()), (*b_first).stride(), a_first.base(), (*a_first).stride(), &beta, c_first.base(), c_first. stride());}
	}else{assert(0);}  // NOLINT(cppcoreguidelines-pro-bounds-array-to-pointer-decay,hicpp-no-array-decay)

	return c_first + a_count;
}

template<class Context, class It2DA, class Size, class It2DB, class It2DC,
	std::enable_if_t<(is_conjugated<It2DA>{} && !is_conjugated<It2DB>{}), int> =0>  // NOLINT(modernize-use-constraints) TODO(correaa) for C++20
auto gemm_n(Context&& ctxt, typename It2DA::element alpha, It2DA a_first, Size a_count, It2DB b_first, typename It2DA::element beta, It2DC c_first) // NOLINT(readability-function-cognitive-complexity) : 125
{
	assert( (*b_first).size() == (*c_first).size() );          // NOLINT(cppcoreguidelines-pro-bounds-array-to-pointer-decay,hicpp-no-array-decay)
	assert( a_first.stride()==1 || (*a_first).stride()==1 ); // NOLINT(cppcoreguidelines-pro-bounds-array-to-pointer-decay,hicpp-no-array-decay)
	assert( b_first.stride()==1 || (*b_first).stride()==1 ); // NOLINT(cppcoreguidelines-pro-bounds-array-to-pointer-decay,hicpp-no-array-decay)
	assert( c_first.stride()==1 || (*c_first).stride()==1 ); // NOLINT(cppcoreguidelines-pro-bounds-array-to-pointer-decay,hicpp-no-array-decay)

	if(a_count == 0) { return c_first; }

	if      (a_first. stride()==1 && (*b_first).stride()==1 && (*c_first).stride()==1){
		if  (a_count==1)        {CTXT->gemm('N', 'C', (*c_first).size(), a_count, (*a_first).size(), &alpha, b_first.base(), b_first. stride(), underlying(a_first.base()), (*a_first).stride(), &beta, base(c_first), (*a_first).size()); }
		else                    {CTXT->gemm('N', 'C', (*c_first).size(), a_count, (*a_first).size(), &alpha, b_first.base(), b_first. stride(), underlying(a_first.base()), (*a_first).stride(), &beta, base(c_first), c_first.stride() ); }
	} else                      {throw std::logic_error{"not BLAS-implemented"};}

	return c_first + a_count;
}

template<class Context, class It2DA, class Size, class It2DB, class It2DC,
	std::enable_if_t<(is_conjugated<It2DA>{} && is_conjugated<It2DB>{}), int> =0>  // NOLINT(modernize-use-constraints) TODO(correaa) for C++20
auto gemm_n(Context&& ctxt, typename It2DA::element alpha, It2DA a_first, Size a_count, It2DB b_first, typename It2DA::element beta, It2DC c_first) // NOLINT(readability-function-cognitive-complexity) : 125
{
	assert( (*b_first).size() == (*c_first).size() );          // NOLINT(cppcoreguidelines-pro-bounds-array-to-pointer-decay,hicpp-no-array-decay)
	assert( a_first.stride()==1 || (*a_first).stride()==1 ); // NOLINT(cppcoreguidelines-pro-bounds-array-to-pointer-decay,hicpp-no-array-decay)
	assert( b_first.stride()==1 || (*b_first).stride()==1 ); // NOLINT(cppcoreguidelines-pro-bounds-array-to-pointer-decay,hicpp-no-array-decay)
	assert( c_first.stride()==1 || (*c_first).stride()==1 ); // NOLINT(cppcoreguidelines-pro-bounds-array-to-pointer-decay,hicpp-no-array-decay)

	if(a_count == 0) { return c_first; }
	if      (a_first. stride()==1 && b_first. stride()==1 && (*c_first).stride()==1){
	                            {CTXT->gemm('C', 'C', a_count, (*c_first).size(), (*a_first).size(), &alpha, underlying(base(b_first)), (*b_first).stride(), underlying(base(a_first)), (*a_first).stride(), &beta, base(c_first), c_first. stride());}
	} else                      {throw std::logic_error{"not BLAS-implemented"};}
	return c_first + a_count;
}

#undef CTXT

template<class It2DA, class Size, class It2DB, class It2DC, class Context = blas::context*> // TODO(correaa) automatic deduction of context
auto gemm_n(typename It2DA::element alpha, It2DA a_first, Size a_count, It2DB b_first, typename It2DA::element beta, It2DC c_first)
->decltype(gemm_n(Context{}, alpha, a_first, a_count, b_first, beta, c_first)) {
	return gemm_n(Context{}, alpha, a_first, a_count, b_first, beta, c_first); }

template<class Context, class A, class B, class C>
auto gemm(Context&& ctx, typename A::element alpha, A const& a, B const& b, typename A::element beta, C&& c) -> C&& {  // NOLINT(readability-identifier-length) BLAS naming
	assert( size( a) == size( c) );
	if(! a.is_empty()) {assert( size(~a) == size( b) );}
	if constexpr(is_conjugated<C>{}) {blas::gemm  (std::forward<Context>(ctx), conj(alpha), conj(a),           conj(b) , conj(beta), conj(c) );}
	else                             {blas::gemm_n(std::forward<Context>(ctx),      alpha , begin(a), size(a), begin(b),      beta , begin(c));}
	return std::forward<C>(c);
}

template<class A, class B, class C>
auto gemm(typename A::element alpha, A const& a, B const& b, typename A::element beta, C&& c) -> C&& {  // NOLINT(readability-identifier-length) BLAS naming
	if constexpr(is_conjugated<A>{}) {
		auto ctxt = blas::default_context_of(underlying(a.base()));
		return gemm(ctxt, alpha, a, b, beta, std::forward<C>(c));
	} else {
		auto ctxt = blas::default_context_of(a.base());
		return gemm(ctxt, alpha, a, b, beta, std::forward<C>(c));
	}
}

// template<class ContextPtr, class Scalar, class ItA, class ItB, class DecayType> class gemm_range;

template<class Ext>
class gemm_reference {  // TODO(correaa) implement this in terms of gemv_range?
	Ext exts_;

 public:
	explicit gemm_reference(Ext exts) : exts_{std::move(exts)} {}
	auto extensions() const {return exts_;}
	friend auto extensions(gemm_reference const& self) {return self.extensions();}
};

template<class ContextPtr, class Scalar, class ItA, class ItB>
class gemm_iterator {
	ContextPtr ctxtp_;
	Scalar s_;
	ItA a_it_;
	ItB b_begin_;
	gemm_iterator(ContextPtr ctxtp, Scalar s, ItA a_it, ItB b_begin) : ctxtp_{ctxtp}, s_{s}, a_it_{std::move(a_it)}, b_begin_{std::move(b_begin)} {}  // NOLINT(readability-identifier-length) BLAS naming
	template<class ContextPtr2, class Scalar2, class ItA2, class ItB2, class DecayType2>
	friend class gemm_range;

 public:
	gemm_iterator(gemm_iterator const&) = default;
	gemm_iterator(gemm_iterator&&) noexcept = default;
	~gemm_iterator() = default;

	auto operator=(gemm_iterator&&) -> gemm_iterator& = delete;
	auto operator=(gemm_iterator const&) -> gemm_iterator& = delete;

	using difference_type = typename std::iterator_traits<ItA>::difference_type;
	using value_type = typename std::iterator_traits<ItA>::value_type;
	using pointer = std::nullptr_t;
	using reference = gemm_reference<decltype((*b_begin_).extensions())>;
	using iterator_category = std::random_access_iterator_tag;

	auto operator+=(difference_type n) -> gemm_iterator& {a_it_ += n; return *this;}
	auto operator-=(difference_type n) -> gemm_iterator& {a_it_ -= n; return *this;}

	auto operator++() -> gemm_iterator& { return operator+=(1); }  // required by random access concept requires even if not used explicitly
	auto operator--() -> gemm_iterator& { return operator-=(1); }

	friend auto operator+(gemm_iterator ret, difference_type n) { return ret += n; }

	friend auto operator-(gemm_iterator const& a, gemm_iterator const& b) -> difference_type {  // NOLINT(readability-identifier-length) BLAS naming
		assert(a.b_begin_ == b.b_begin_);  // NOLINT(cppcoreguidelines-pro-bounds-array-to-pointer-decay,hicpp-no-array-decay)
		return a.a_it_ - b.a_it_;
	}
	friend auto operator==(gemm_iterator const& self, gemm_iterator const& other) -> bool {return self.a_it_ == other.a_it_;}
	friend auto operator!=(gemm_iterator const& self, gemm_iterator const& other) -> bool {return self.a_it_ != other.a_it_;}

	template<class ItOut>
	friend auto copy_n(gemm_iterator const& first, difference_type count, ItOut d_first)
	->decltype(blas::gemm_n(std::declval<ContextPtr>(), std::declval<typename ItA::element>()       , std::declval<ItA>(), count, std::declval<ItB>(), 0.0, d_first)) try {  // std::complex NOLINT(fuchsia-default-arguments-calls)
		return blas::gemm_n(first.ctxtp_              , static_cast<typename ItA::element>(first.s_), first.a_it_        , count, first.b_begin_     , 0.0, d_first);  // NOLINT(fuchsia-default-arguments-calls)
	} catch(std::exception const& e) {
		throw std::logic_error(
			"in `copy_n`\nCouldn't decay product of arrays of size "+ std::to_string(count) +"x"+ std::to_string((*first.a_it_).size()) + " and " + // NOLINT(cppcoreguidelines-pro-bounds-array-to-pointer-decay,hicpp-no-array-decay)
			std::to_string((*first.a_it_).size())+ "x" +std::to_string((*first.b_begin_).size()) + " into " + std::to_string(count) +"x" + std::to_string((*first.b_begin_).size()) +
			"\nbecause\n" + e.what()
		);
	}

	template<class ItOut>
	friend auto copy(gemm_iterator const& first, gemm_iterator const& last, ItOut d_first) {assert(first.s_ == last.s_);
		return copy_n(first, last - first, d_first);
	}

	template<class ItOut>
	friend auto uninitialized_copy_n(gemm_iterator const& first, difference_type count, ItOut d_first) {
		return copy_n(first, count, d_first);
	}

	template<class ItOut>
	friend auto uninitialized_copy(gemm_iterator const& first, gemm_iterator const& last, ItOut d_first) {
		assert( first.s_ == last.s_ ); // NOLINT(cppcoreguidelines-pro-bounds-array-to-pointer-decay,hicpp-no-array-decay)
		return uninitialized_copy_n(first, last - first, d_first);
	}

	auto operator*() const {return reference{(*b_begin_).extensions()};}
};

template<class ContextPtr, class Scalar, class ItA, class ItB, class DecayType>
class gemm_range {
	ContextPtr ctxtp_;
	Scalar s_;
	ItA a_begin_;
	ItA a_end_;
	ItB b_begin_;

 public:
	gemm_range(gemm_range const&) = delete;
	gemm_range(gemm_range&&) = delete;
	auto operator=(gemm_range const&) -> gemm_range& = delete;
	auto operator=(gemm_range&&) -> gemm_range& = delete;
	~gemm_range() = default;

	gemm_range(ContextPtr ctxtp, Scalar s, ItA a_first, ItA a_last, ItB b_first)  // NOLINT(bugprone-easily-swappable-parameters,readability-identifier-length) BLAS naming
	: ctxtp_{ctxtp}
	, s_{s}, a_begin_{std::move(a_first)}, a_end_{std::move(a_last)}
	, b_begin_{std::move(b_first)}
	{}

	using iterator = gemm_iterator<ContextPtr, Scalar, ItA, ItB>;
	using decay_type = DecayType;
	using size_type = typename decay_type::size_type;

	       auto begin()          const& -> iterator {return {ctxtp_, s_, a_begin_, b_begin_};}
	       auto end()            const& -> iterator {return {ctxtp_, s_, a_end_  , b_begin_};}
	friend auto begin(gemm_range const& self) {return self.begin();}
	friend auto end  (gemm_range const& self) {return self.end  ();}

	auto size() const -> size_type {return a_end_ - a_begin_;}

	auto extensions() const -> typename decay_type::extensions_type {return size()*(*b_begin_).extensions();}
	friend auto extensions(gemm_range const& self) {return self.extensions();}

	auto operator+() const -> decay_type {return *this;} // TODO(correaa) : investigate why return decay_type{*this} doesn't work
	template<class Arr>
	friend auto operator+=(Arr&& a, gemm_range const& self) -> Arr&& {  // NOLINT(readability-identifier-length) BLAS naming
		blas::gemm_n(self.ctxtp_, self.s_, self.a_begin_, self.a_end_ - self.a_begin_, self.b_begin_, 1., a.begin());
		return std::forward<Arr>(a);
	}
	friend auto operator*(Scalar factor, gemm_range const& self) {
		return gemm_range{self.ctxtp_, factor*self.s_, self.a_begin_, self.a_end_, self.b_begin_};
	}
};

template<class ContextPtr, class Scalar, class A2D, class B2D,
	class = std::enable_if_t<is_context<decltype(*std::declval<ContextPtr>())>{}> >  // NOLINT(modernize-use-constraints) TODO(correaa) for C++20
auto gemm(ContextPtr ctxtp, Scalar s, A2D const& a, B2D const& b)  // NOLINT(readability-identifier-length) BLAS naming
->gemm_range<ContextPtr, Scalar, typename A2D::const_iterator, typename B2D::const_iterator, typename A2D::decay_type/*B2D*/>
{
	return
		gemm_range<ContextPtr, Scalar, typename A2D::const_iterator, typename B2D::const_iterator, typename A2D::decay_type/*B2D*/>
			(ctxtp, s, a.begin(), a.end(), b.begin())
		;
}

#if defined __NVCC__
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
template<class A2D, class B2D, class Scalar = typename A2D::element_type, class = decltype(Scalar{0.0})>
auto gemm(Scalar s, A2D const& a, B2D const& b) {  // NOLINT(readability-identifier-length) conventional BLAS naming
	if constexpr(is_conjugated<A2D>{}) {
		auto ctxtp = blas::default_context_of(underlying(a.base()));
		return blas::gemm(ctxtp, s, a, b);
	} else {
		auto ctxtp = blas::default_context_of(a.base());
		return blas::gemm(ctxtp, s, a, b);
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
	template<class A2D, class B2D,
		std::enable_if_t<(A2D::dimensionality == 2) && (B2D::dimensionality == 2),int> =0>  // NOLINT(modernize-use-constraints) TODO(correaa) for C++20
	auto operator*(A2D const& A, B2D const& B)  // NOLINT(readability-identifier-length) conventional BLAS names
	->decltype(blas::gemm(1.0, A, B)) {
		return blas::gemm(1.0, A, B); }
}  // end namespace operators

}  // end namespace boost::multi::blas
#endif
