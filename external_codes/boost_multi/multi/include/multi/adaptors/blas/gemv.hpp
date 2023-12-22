// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;autowrap:nil;-*-
// Copyright 2019-2022 Alfredo A. Correa

#ifndef MULTI_ADAPTORS_BLAS_GEMV_HPP
#define MULTI_ADAPTORS_BLAS_GEMV_HPP

#include "../blas/core.hpp"
#include "../blas/dot.hpp"

#include "./../../detail/../utility.hpp"

namespace boost::multi::blas {

using core::gemv;

template<class Context, class MIt, class Size, class XIt, class YIt>
auto gemv_n(Context ctxt, typename MIt::element a, MIt m_first, Size count, XIt x_first, typename MIt::element b, YIt y_first) {  // NOLINT(readability-identifier-length) BLAS naming
	assert(m_first->stride()==1 or m_first.stride()==1); // blas doesn't implement this case
	assert( x_first.base() != y_first.base() );

	if constexpr(not is_conjugated<MIt>::value) {
		if     (m_first .stride()==1) {ctxt->gemv('N', count, m_first->size(), &a, m_first.base()            , m_first->stride(), x_first.base(), x_first.stride(), &b, y_first.base(), y_first.stride());}
		else if(m_first->stride()==1) {ctxt->gemv('T', m_first->size(), count, &a, m_first.base()            , m_first. stride(), x_first.base(), x_first.stride(), &b, y_first.base(), y_first.stride());}
		else                          {throw std::logic_error{"not BLAS-implemented"};}  // LCOV_EXCL_LINE
	} else {
		if     (m_first->stride()==1) {ctxt->gemv('C', m_first->size(), count, &a, underlying(m_first.base()), m_first. stride(), x_first.base(), x_first.stride(), &b, y_first.base(), y_first.stride());}
	//  else if(m_first. stride()==1) {assert(0);} // not implemented in blas (use cblas?)
		else                          {throw std::logic_error{"not BLAS-implemented"};}  // LCOV_EXCL_LINE
	}

	struct {
		MIt m_last;
		YIt y_last;
	} ret{m_first + count, y_first + count};

	return ret;
}

template<class A, class MIt, class Size, class XIt, class B, class YIt>
auto gemv_n(A a, MIt m_first, Size count, XIt x_first, B b, YIt y_first) {  // NOLINT(readability-identifier-length) BLAS naming
	blas::context ctxt;
	return gemv_n(&ctxt, a, m_first, count, x_first, b, y_first);
}

template<class Ctxt, class A, class M, class V, class B, class W>
auto gemv(Ctxt ctxt, A const& a, M const& m, V const& v, B const& b, W&& w) -> W&& {  // NOLINT(readability-identifier-length) BLAS naming
	assert(size( m) == size(w) );
	assert(size(~m) == size(v) );
	// gemv_n(a, begin(m), size(m), begin(v), b, begin(w));
	gemv_n(ctxt, a, begin(m), size(m), begin(v), b, begin(w));  // NOLINT(fuchsia-default-arguments-calls)
	return std::forward<W>(w);


	// assert(size( m) == size(w) );
	// assert(size(~m) == size(v) );
	// gemv_n(a, begin(m), size(m), begin(v), b, begin(w));
	// return std::forward<W>(w);

	// if constexpr(is_conjugated<M>{}) {
	//  auto ctxtp = blas::default_context_of(underlying(m.base()));
	//  return blas::gemv(ctxtp, a, m, v, b, w);
	// } else {
	//  auto ctxtp = blas::default_context_of(m.base());
	//  return blas::gemv(ctxtp, a, m, v, b, w);
	// }
}


template<class A, class M, class V, class B, class W>
auto gemv(A const& a, M const& m, V const& v, B const& b, W&& w) -> W&& {  // NOLINT(readability-identifier-length) BLAS naming
	assert(size( m) == size(w) );
	// assert(size(~m) == size(v) );
	// gemv_n(a, begin(m), size(m), begin(v), b, begin(w));
	// return std::forward<W>(w);

	if constexpr(is_conjugated<M>{}) {
		auto ctxtp = blas::default_context_of(underlying(m.base()));
		return blas::gemv(ctxtp, a, m, v, b, std::forward<W>(w));
	} else {
		auto ctxtp = blas::default_context_of(m.base());
		return blas::gemv(ctxtp, a, m, v, b, std::forward<W>(w));
	}
}

template<class Scalar, class It2D, class It1D, class Context>
class gemv_iterator {
	Scalar alpha_ = 1.0;
	It2D m_it_;
	It1D v_first_;
	Context ctxt_;

 public:
	using difference_type = typename std::iterator_traits<It2D>::difference_type;
	using value_type = typename std::iterator_traits<It1D>::value_type;
	using pointer = void;
	using reference = void;
	using iterator_category = std::random_access_iterator_tag;

	friend auto operator-(gemv_iterator const& self, gemv_iterator const& other) -> difference_type {
		assert(self.v_first_ == other.v_first_);
		return self.m_it_ - other.m_it_;
	}
	template<class It1DOut>
	friend auto copy_n(gemv_iterator first, difference_type count, It1DOut result){
		if constexpr(std::is_same_v<Context, void>) {blas::gemv_n(             first.alpha_, first.m_it_, count, first.v_first_, 0., result);}  // NOLINT(fuchsia-default-arguments-calls)
		else                                        {blas::gemv_n(first.ctxt_, first.alpha_, first.m_it_, count, first.v_first_, 0., result);}  // NOLINT(fuchsia-default-arguments-calls)
		return result + count;
	}
	template<class It1DOut>
	friend auto copy(gemv_iterator first, gemv_iterator last, It1DOut result){return copy_n(first, last - first, result);}
	template<class It1DOut>
	friend auto uninitialized_copy(gemv_iterator first, gemv_iterator last, It1DOut result) {
		// static_assert(std::is_trivially_default_constructible_v<typename It1DOut::value_type>);
		return copy(first, last, result);
	}
	gemv_iterator(Scalar alpha, It2D m_it, It1D v_first, Context ctxt)
	: alpha_{alpha}, m_it_{std::move(m_it)}, v_first_{std::move(v_first)}, ctxt_{ctxt} {}
	auto operator*() const { return value_type{0.0}; }  // could be std::complex NOLINT(fuchsia-default-arguments-calls)
};

template<class Scalar, class It2D, class It1D, class DecayType, class Context>
class gemv_range {
	Scalar alpha_ = 1.0;
	It2D m_begin_;
	It2D m_end_;
	It1D v_first_;
	Context ctxt_;

 public:
	gemv_range(gemv_range&&) noexcept = default;
	gemv_range(gemv_range const&) = delete;
	~gemv_range() = default;
	auto operator=(gemv_range const&) = delete;
	auto operator=(gemv_range&&) = delete;

	gemv_range(Scalar alpha, It2D m_first, It2D m_last, It1D v_first)  // NOLINT(bugprone-easily-swappable-parameters)
		: alpha_{alpha}, m_begin_{std::move(m_first)}, m_end_{std::move(m_last)}, v_first_{std::move(v_first)} {
		assert(m_begin_.stride() == m_end_.stride());
	}
	gemv_range(Context ctxt, Scalar alpha, It2D m_first, It2D m_last, It1D v_first)  // NOLINT(bugprone-easily-swappable-parameters)
	: alpha_{alpha}
	, m_begin_{std::move(m_first)}, m_end_{std::move(m_last)}
	, v_first_{std::move(v_first)}
	, ctxt_{std::forward<Context>(ctxt)} {
		assert(m_begin_.stride() == m_end_.stride());
	}
	using iterator = gemv_iterator<Scalar, It2D, It1D, Context>;
	using decay_type = DecayType;

	auto begin() const -> iterator{return {alpha_, m_begin_, v_first_, ctxt_};}
	auto end()   const -> iterator{return {alpha_, m_end_  , v_first_, ctxt_};}

	auto size() const -> size_type{return end() - begin();}
	auto extensions() const -> typename decay_type::extensions_type{return typename decay_type::extensions_type{{0, size()}};}
	auto decay() const{return decay_type{*this};}

	friend auto operator+(gemv_range const& self) {return self.decay();}
	template<class V>
	friend auto operator+=(V&& v, gemv_range const& s) -> V&& {  // NOLINT(readability-identifier-length) BLAS naming
		// if constexpr(std::is_same<Context, void*>{}) {blas::gemv_n(         s.alpha_, s.m_begin_, s.m_end_ - s.m_begin_, s.v_first_, 1.0, v.begin());}
		// else                                         
		{blas::gemv_n(s.ctxt_, s.alpha_, s.m_begin_, s.m_end_ - s.m_begin_, s.v_first_, 1., v.begin());}
		return std::forward<V>(v);
	}
};

template<class Context, class Scalar, class M, class V>
auto gemv(Context ctxt, Scalar s, M const& m, V const& v) {  // NOLINT(readability-identifier-length) BLAS naming
	assert(size(~m) == size(v));
	return gemv_range<Scalar, typename M::const_iterator, typename V::const_iterator, typename V::decay_type, Context>(ctxt, s, m.begin(), m.end(), v.begin());
}

template<class Scalar, class M, class V>
auto gemv(Scalar s, M const& m, V const& v)  // NOLINT(readability-identifier-length) BLAS naming
{
	if constexpr(is_conjugated<M>{}) {
		auto ctxtp = blas::default_context_of(underlying(m.base()));
		return blas::gemv(ctxtp, s, m, v);
	} else {
		auto ctxtp = blas::default_context_of(m.base());
		return blas::gemv(ctxtp, s, m, v);
	}
}

template<class T, class Matrix>
struct scaled_matrix {
	T aa_;
	Matrix const& A_;  // NOLINT(readability-identifier-length,cppcoreguidelines-avoid-const-or-ref-data-members) BLAS naming
	
	template<class Vector>
	friend auto operator%(scaled_matrix const& aaA, Vector const& x) {  // NOLINT(readability-identifier-length) BLAS naming
		return blas::gemv(aaA.aa_, aaA.A_, x);
	}
};

namespace operators {
	template<class M, class V>
	auto operator%(M const& m, V const& v)  // NOLINT(readability-identifier-length) BLAS naming
	->decltype(+blas::gemv(1.0, m, v)) {
		return +blas::gemv(1.0, m, v); }

	template<class Matrix, std::enable_if_t<Matrix::dimensionality == 2, int> =0>
	auto operator*(typename Matrix::element_type aa, Matrix const& A) {  // NOLINT(readability-identifier-length) BLAS naming
		return scaled_matrix<typename Matrix::element_type, Matrix const&>{aa, A};
	}

} // end namespace operators

} // end namespace boost::multi::blas

#endif
