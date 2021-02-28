// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;autowrap:nil;-*-
// © Alfredo A. Correa 2019-2020

#ifndef MULTI_ADAPTORS_BLAS_GEMV_HPP
#define MULTI_ADAPTORS_BLAS_GEMV_HPP

#include "../blas/core.hpp"

#include "../blas/dot.hpp"

#include "./../../detail/../utility.hpp"

namespace boost{
namespace multi::blas{

using core::gemv;

template<class Context, class A, class MIt, class Size, class XIt, class B, class YIt>
auto gemv_n(Context&& ctxt, A a, MIt m_first, Size count, XIt x_first, B b, YIt y_first){
	assert(m_first->stride()==1 or m_first.stride()==1); // blas doesn't implement this case
	assert( x_first.base() != y_first.base() ); 
	if constexpr(not is_conjugated<MIt>{}){
		assert( y_first.base() !=            m_first.base()  );
		;;;; if(m_first .stride()==1) std::forward<Context>(ctxt).gemv('N', count, m_first->size(), a, m_first.base()            , m_first->stride(), x_first.base(), x_first.stride(), b, y_first.base(), y_first.stride());
		else if(m_first->stride()==1) std::forward<Context>(ctxt).gemv('T', m_first->size(), count, a, m_first.base()            , m_first. stride(), x_first.base(), x_first.stride(), b, y_first.base(), y_first.stride());
		else                          assert(0);
	}else{
		assert( y_first.base() != underlying(m_first.base()) );
		;;;; if(m_first->stride()==1) std::forward<Context>(ctxt).gemv('C', m_first->size(), count, a, underlying(m_first.base()), m_first. stride(), x_first.base(), x_first.stride(), b, y_first.base(), y_first.stride());
		else if(m_first. stride()==1) assert(0); // not implemented in blas (use cblas?)
		else                          assert(0); // not implemented in blas
	}
	struct{
		MIt m_last;
		YIt y_last;
	} ret{m_first + count, y_first + count};
	return ret;
}

template<class A, class MIt, class Size, class XIt, class B, class YIt>
auto gemv_n(A a, MIt m_first, Size count, XIt x_first, B b, YIt y_first){
	return gemv_n(blas::context{}, a, m_first, count, x_first, b, y_first);
}

template<class A, class M, class V, class B, class W>
W&& gemv(A const& a, M const& m, V const& v, B const& b, W&& w){
	assert(size( m) == size(w) );
	assert(size(~m) == size(v) );
	gemv_n(a, begin(m), size(m), begin(v), b, begin(w));
	return std::forward<W>(w);
}

template<class Scalar, class It2D, class It1D, class Context>
class gemv_iterator{
	Scalar alpha_ = 1.;
	It2D m_it_; 
	It1D v_first_;
	Context ctxt_;
public:
	using difference_type = typename std::iterator_traits<It2D>::difference_type;
	using value_type = typename std::iterator_traits<It1D>::value_type;
	using pointer = void;
	using reference = void;
	using iterator_category = std::random_access_iterator_tag;
//	using iterator_category = std::output_iterator_tag;
//	friend difference_type distance(gemv_iterator const& a, gemv_iterator const& b){assert(a.v_first_ == b.v_first_);
//		return b.m_it_ - a.m_it_;
//	}
	friend difference_type operator-(gemv_iterator const& a, gemv_iterator const& b){assert(a.v_first_ == b.v_first_);
		return a.m_it_ - b.m_it_;
	}
	template<class It1DOut>
	friend auto copy_n(gemv_iterator first, difference_type count, It1DOut result){
		if constexpr
			(std::is_same<Context, void>{}) blas::gemv_n(             first.alpha_, first.m_it_, count, first.v_first_, 0., result);
		else                                blas::gemv_n(first.ctxt_, first.alpha_, first.m_it_, count, first.v_first_, 0., result);
		return result + count;
	}
	template<class It1DOut>
	friend auto copy(gemv_iterator first, gemv_iterator last, It1DOut result){return copy_n(first, last - first, result);}
	template<class It1DOut>
	friend auto uninitialized_copy(gemv_iterator first, gemv_iterator last, It1DOut result){
		static_assert(std::is_trivially_default_constructible<typename It1DOut::value_type>{});
		return copy(first, last, result);
	}
	gemv_iterator(Scalar alpha, It2D m_it, It1D v_first, Context ctxt) 
		: alpha_{alpha}, m_it_{m_it}, v_first_{v_first}, ctxt_{ctxt}{}
	value_type operator*() const{return 0.;}
};

template<class Scalar, class It2D, class It1D, class DecayType, class Context>
class gemv_range{
	Scalar alpha_ = 1.;
	It2D m_begin_;
	It2D m_end_;
	It1D v_first_;
	Context ctxt_ = {};
public:
	gemv_range(gemv_range const&) = delete;
	gemv_range(Scalar alpha, It2D m_first, It2D m_last, It1D v_first) 
		: alpha_{alpha}, m_begin_{m_first}, m_end_{m_last}, v_first_{v_first}{
		assert(m_begin_.stride() == m_end_.stride());
	}
	gemv_range(Context&& ctxt, Scalar alpha, It2D m_first, It2D m_last, It1D v_first) 
		: alpha_{alpha}, m_begin_{m_first}, m_end_{m_last}, v_first_{v_first}, ctxt_{std::forward<Context>(ctxt)}{
		assert(m_begin_.stride() == m_end_.stride());
	}
	using iterator = gemv_iterator<Scalar, It2D, It1D, Context>;
	using decay_type = DecayType;
	iterator begin() const{return {alpha_, m_begin_, v_first_, ctxt_};}
	iterator end()   const{return {alpha_, m_end_  , v_first_, ctxt_};}
	size_type size() const{return end() - begin();}
	typename decay_type::extensions_type extensions() const{return {{0, size()}};}
	decay_type decay() const{
		decay_type ret; 
		ret = *this;
		return ret;
	}
	friend auto operator+(gemv_range const& self){return self.decay();}
	template<class V>
	friend V&& operator+=(V&& v, gemv_range const& s){
		if constexpr
			(std::is_same<Context, void*>{}) blas::gemv_n(         s.alpha_, s.m_begin_, s.m_end_ - s.m_begin_, s.v_first_, 1., v.begin());
		else                                 blas::gemv_n(s.ctxt_, s.alpha_, s.m_begin_, s.m_end_ - s.m_begin_, s.v_first_, 1., v.begin());
		return std::forward<V>(v);
	}
};

template<class Scalar, class M, class V>
auto gemv(Scalar s, M const& m, V const& v)
{//->decltype(gemv_range{s, m, v}){
	assert(size(~m) == size(v));
	return gemv_range<Scalar, typename M::const_iterator, typename V::const_iterator, typename V::decay_type, blas::context>(s, m.begin(), m.end(), v.begin());}

template<class Context, class Scalar, class M, class V>
auto gemv(Context&& ctxt, Scalar s, M const& m, V const& v)
//->decltype(gemv_ranges, m, v})
{	assert(size(~m) == size(v));
	return gemv_range<Scalar, typename M::const_iterator, typename V::const_iterator, typename V::decay_type, Context&&>(std::forward<Context>(ctxt), s, m.begin(), m.end(), v.begin());}

namespace operators{
	template<class M, class V>
	auto operator%(M const& m, V const& v)
	->decltype(+blas::gemv(1., m, v)){
		return +blas::gemv(1., m, v);}
}

}
}

#endif

