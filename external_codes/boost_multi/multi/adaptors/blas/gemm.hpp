// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;autowrap:nil;-*-
// © Alfredo A. Correa 2019-2021

#ifndef MULTI_ADAPTORS_BLAS_GEMM_HPP
#define MULTI_ADAPTORS_BLAS_GEMM_HPP

#include "../blas/core.hpp"
#include "../blas/gemv.hpp"
#include "../blas/numeric.hpp"
#include "../blas/operations.hpp"

namespace boost{
namespace multi{
namespace blas{

using core::gemm;

template<class It>
auto xbase_aux(It const& it, std::true_type const&)
->decltype(underlying(base(it))){
	return underlying(base(it));}

template<class It>
auto xbase_aux(It const& it, std::false_type const&)
->decltype(base(it)){
	return base(it);}

template<class It>
auto xbase(It const& it)
->decltype(xbase_aux(it, std::integral_constant<bool, is_conjugated<It>{}>{})){
	return xbase_aux(it, std::integral_constant<bool, is_conjugated<It>{}>{});}

template<class Context, class It2DA, class Size, class It2DB, class It2DC>
auto gemm_n(Context&& ctxt, typename It2DA::element alpha, It2DA a_first, Size a_count, It2DB b_first, typename It2DA::element beta, It2DC c_first)
//->decltype(std::forward<Context>(ctxt).gemm('N', 'N', b_first->size(), a_count, a_first->size(), &alpha, xbase(b_first), b_first->size()  , xbase(a_first), a_first->size(), &beta, c_first.base(), c_first->size()  ), It2DC{})
try{
	assert( b_first->size() == c_first->size() );
	assert( a_first.stride()==1 or a_first->stride()==1 );
	assert( b_first.stride()==1 or b_first->stride()==1 );
	assert( c_first.stride()==1 or c_first->stride()==1 );

	if(a_count != 0){
		#define CTXT std::forward<Context>(ctxt)
		;;;;; if constexpr(!is_conjugated<It2DA>{} and !is_conjugated<It2DB>{}){
			;;;;; if(a_first->stride()==1 and b_first->stride()==1 and c_first->stride()==1){
				;;;; if( a_count==1 and b_first->size()==1 ){CTXT.gemm('N', 'N', b_first->size(), a_count, a_first->size(), &alpha, base(b_first), b_first->size()  , base(a_first), a_first->size() , &beta, base(c_first), c_first->size()  );}
				else if( a_count==1                        ){CTXT.gemm('N', 'N', b_first->size(), a_count, a_first->size(), &alpha, base(b_first), b_first. stride(), base(a_first), a_first->size()  , &beta, base(c_first), c_first->size()  );}
				else                                        {CTXT.gemm('N', 'N', b_first->size(), a_count, a_first->size(), &alpha, base(b_first), b_first. stride(), base(a_first), a_first. stride(), &beta, base(c_first), c_first. stride());}
			}else if(a_first->stride()==1 and b_first->stride()==1 and c_first. stride()==1){
				if  (a_count==1)        {CTXT.gemm('T', 'T', a_count, b_first->size(), a_first->size(), &alpha, base(a_first), a_first. stride(), base(b_first), b_first->size() , &beta, base(c_first), a_first->size()  );}
				else                    {CTXT.gemm('T', 'T', a_count, b_first->size(), a_first->size(), &alpha, base(a_first), a_first. stride(), base(b_first), b_first.stride(), &beta, base(c_first), c_first->stride());}
			}else if(a_first. stride()==1 and b_first->stride()==1 and c_first->stride()==1){
				if  (a_count==1)        {CTXT.gemm('N', 'T', c_first->size(), a_count, a_first->size(), &alpha, base(b_first), b_first. stride(), base(a_first), a_first->stride(), &beta, base(c_first), a_first->size()  );}
				else                    {CTXT.gemm('N', 'T', c_first->size(), a_count, a_first->size(), &alpha, base(b_first), b_first. stride(), base(a_first), a_first->stride(), &beta, base(c_first), c_first.stride());}
			}else if(a_first. stride()==1 and b_first->stride()==1 and c_first. stride()==1){
				if  (a_count==1)        {CTXT.gemm('N', 'T', a_count, b_first->size(), a_first->size(), &alpha, base(a_first), a_first->stride(), base(b_first), a_first->size()  , &beta, base(c_first), b_first->size()  );}
				else                    {CTXT.gemm('N', 'T', a_count, b_first->size(), a_first->size(), &alpha, base(a_first), a_first->stride(), base(b_first), b_first. stride(), &beta, base(c_first), c_first->stride());}
			}else if(a_first->stride()==1 and b_first.stride()==1 and c_first. stride()==1){
				;;;; if(a_count==1 and b_first->size()){CTXT.gemm('N', 'N', c_first->size(), a_count, a_first->size(), &alpha, base(b_first), b_first->size()  , base(a_first), a_first->size()  , &beta, base(c_first), c_first->stride());}
				else if(a_count==1)                    {CTXT.gemm('N', 'T', c_first->size(), a_count, a_first->size(), &alpha, base(b_first), b_first->stride(), base(a_first), a_first->size()  , &beta, base(c_first), c_first->stride());}
				else                                   {CTXT.gemm('N', 'T', c_first->size(), a_count, a_first->size(), &alpha, base(b_first), b_first->stride(), base(a_first), a_first.stride() , &beta, base(c_first), c_first->stride());}
			}else if(a_first->stride()==1 and b_first. stride()==1 and c_first->stride()==1){
				if  (a_count==1)        {CTXT.gemm('T', 'N', a_count, c_first->size(), a_first->size(), &alpha, base(b_first), b_first->stride(), base(a_first), a_first->size()  , &beta, base(c_first), c_first.stride());}
				else                    {CTXT.gemm('T', 'N', c_first->size(), a_count, a_first->size(), &alpha, base(b_first), b_first->stride(), base(a_first), a_first.stride(), &beta, base(c_first), c_first.stride());}
			}else if(a_first. stride()==1 and b_first.stride( )==1 and c_first. stride()==1){
				                        {CTXT.gemm('N', 'N', c_first->size(), a_count, a_first->size(), &alpha, base(a_first), a_first->stride(), base(b_first), b_first->stride(), &beta, base(c_first), c_first->stride());}
			}else if(a_first. stride()==1 and b_first.stride( )==1 and c_first->stride()==1){
				                        {CTXT.gemm('T', 'T', a_count, c_first->size(), a_first->size(), &alpha, base(b_first), b_first->stride(), base(a_first), a_first->stride(), &beta, base(c_first), c_first. stride());}
			}else assert(0);
		}else if constexpr(!is_conjugated<It2DA>{} and  is_conjugated<It2DB>{}){
			;;;;; if(a_first->stride()==1 and b_first->stride()==1 and c_first->stride()==1){
				if(b_first->size()==1)  {CTXT.gemm('C', 'N', c_first->size(), a_count, a_first->size(), &alpha, underlying(base(b_first)), b_first->stride(), base(a_first), a_first->size()  , &beta, base(c_first), c_first.stride());}
				else                    {CTXT.gemm('C', 'N', c_first->size(), a_count, a_first->size(), &alpha, underlying(base(b_first)), b_first->stride(), base(a_first), a_first->size()  , &beta, base(c_first), c_first.stride());}
			}else if(a_first->stride()==1 and b_first. stride()==1 and c_first->stride()==1){
				if  (a_count==1)        {CTXT.gemm('C', 'N', a_count, c_first->size(), a_first->size(), &alpha, underlying(base(b_first)), b_first->stride(), base(a_first), a_first->size()  , &beta, base(c_first), c_first.stride());}
				else                    {CTXT.gemm('C', 'N', c_first->size(), a_count, a_first->size(), &alpha, underlying(base(b_first)), b_first->stride(), base(a_first), a_first.stride(), &beta, base(c_first), c_first.stride());}
			}else if(a_first->stride()==1 and b_first. stride()==1 and c_first. stride()==1){
				                        {CTXT.gemm('C', 'N', c_first->size(), a_count, a_first->size(), &alpha, underlying(base(b_first)), b_first->stride(), base(a_first), a_first. stride(), &beta, base(c_first), c_first->stride());}
			}else if(a_first. stride()==1 and b_first. stride()==1 and c_first. stride()==1){
				                        {CTXT.gemm('C', 'T', c_first->size(), a_count, a_first->size(), &alpha, underlying(base(b_first)), b_first->stride(), base(a_first), a_first->stride(), &beta, base(c_first), c_first->stride());}
			}else if(a_first. stride()==1 and b_first. stride()==1 and c_first->stride()==1){
				                        {CTXT.gemm('C', 'T', a_count, c_first->size(), a_first->size(), &alpha, underlying(base(b_first)), b_first->stride(), base(a_first), a_first->stride(), &beta, base(c_first), c_first. stride());}
			}else assert(0);
		}else if constexpr( is_conjugated<It2DA>{} and !is_conjugated<It2DB>{}){
			;;;;; if(a_first. stride()==1 and b_first->stride()==1 and c_first->stride()==1){
				if  (a_count==1)        {CTXT.gemm('N', 'C', c_first->size(), a_count, a_first->size(), &alpha, base(b_first), b_first. stride(), underlying(base(a_first)), a_first->stride(), &beta, base(c_first), a_first->size()  );}
				else                    {CTXT.gemm('N', 'C', c_first->size(), a_count, a_first->size(), &alpha, base(b_first), b_first. stride(), underlying(base(a_first)), a_first->stride(), &beta, base(c_first), c_first.stride());}
			}else assert(0);
		}else if constexpr( is_conjugated<It2DA>{} and  is_conjugated<It2DB>{}){
			;;;;; if(a_first. stride()==1 and b_first. stride()==1 and c_first->stride()==1){
				                        {CTXT.gemm('C', 'C', a_count, c_first->size(), a_first->size(), &alpha, underlying(base(b_first)), b_first->stride(), underlying(base(a_first)), a_first->stride(), &beta, base(c_first), c_first. stride());}
			}else assert(0);
		}
		#undef CTXT
	}
	return c_first + a_count;
}catch(std::logic_error& e){
	using std::to_string;
	throw std::logic_error{
		"couldn't do "+std::string(__PRETTY_FUNCTION__)+" of layout a_count="+std::to_string(a_count)
		+" a_strides="+to_string(a_first.stride())+","+to_string(a_first->stride())+" a->size="+to_string(a_first->size())
		+" b_strides="+to_string(b_first.stride())+","+to_string(b_first->stride())+" b->size="+to_string(b_first->size())
		+" c_strides="+to_string(c_first.stride())+","+to_string(c_first->stride())+" c->size="+to_string(c_first->size())
		+" because " + e.what()
	};
}

template<class It2DA, class Size, class It2DB, class It2DC, class Context = blas::context> // TODO automatic deduction of context
auto gemm_n(typename It2DA::element alpha, It2DA a_first, Size a_count, It2DB b_first, typename It2DA::element beta, It2DC c_first)
->decltype(gemm_n(Context{}, alpha, a_first, a_count, b_first, beta, c_first)){
	return gemm_n(Context{}, alpha, a_first, a_count, b_first, beta, c_first);}

template<class Context, class A, class B, class C>
C&& gemm(Context&& ctx, typename A::element alpha, A const& a, B const& b, typename A::element beta, C&& c){
	assert( size( a) == size( c) );
	if(not a.is_empty()) assert( size(~a) == size( b) );
	if constexpr(is_conjugated<C>{}){blas::gemm  (std::forward<Context>(ctx), conj(alpha), conj(a),           conj(b) , conj(beta), conj(c) );}
	else                            {blas::gemm_n(std::forward<Context>(ctx),      alpha , begin(a), size(a), begin(b),      beta , begin(c));}
	return std::forward<C>(c);
}

template<class A, class B, class C>
C&& gemm(typename A::element alpha, A const& a, B const& b, typename A::element beta, C&& c){
	return gemm(blas::context{}, alpha, a, b, beta, std::forward<C>(c));
}

template<class ContextPtr, class Scalar, class ItA, class ItB, class DecayType>
class gemm_range;

template<class Ext>
struct gemm_reference{ // TODO implement this in terms of gemv_range
	Ext x;
	Ext const& extensions() const{return x;}
	friend Ext const& extensions(gemm_reference const& self){return self.extensions();}
};

template<class ContextPtr, class Scalar, class ItA, class ItB>
class gemm_iterator{
	ContextPtr ctxtp_;
	Scalar s_;
	ItA a_it_;
	ItB b_begin_;
	gemm_iterator(ContextPtr ctxtp, Scalar s, ItA a_it, ItB b_begin) : ctxtp_{ctxtp}, s_{s}, a_it_{a_it}, b_begin_{b_begin}{}
	template<class ContextPtr2, class Scalar2, class ItA2, class ItB2, class DecayType2>
	friend class gemm_range;
public:
	gemm_iterator(gemm_iterator const&) = default;
	using difference_type = typename std::iterator_traits<ItA>::difference_type;
	using value_type = typename std::iterator_traits<ItA>::value_type;
	using pointer = void*;
	using reference = gemm_reference<decltype(b_begin_->extensions())>;
	using iterator_category = std::random_access_iterator_tag; // using iterator_category = std::input_iterator_tag;
	
	static_assert( std::is_base_of<std::random_access_iterator_tag, typename std::iterator_traits<gemm_iterator>::iterator_category>{} );
	
	gemm_iterator& operator+=(difference_type n){a_it_ += n; return *this;}
	gemm_iterator& operator-=(difference_type n){a_it_ -= n; return *this;}

	gemm_iterator& operator++(){return operator+=(1);} // required by random access concept requires even if not used explicitly
	gemm_iterator& operator--(){return operator-=(1);}
	
	auto operator+(difference_type n) const{gemm_iterator ret{*this}; ret+=n; return ret;}

	friend difference_type operator-(gemm_iterator const& a, gemm_iterator const& b){assert(a.b_begin_ == b.b_begin_);
		return a.a_it_ - b.a_it_;
	}
	friend bool operator==(gemm_iterator const& a, gemm_iterator const& b){return a.a_it_ == b.a_it_;}
	friend bool operator!=(gemm_iterator const& a, gemm_iterator const& b){return a.a_it_ != b.a_it_;}

	template<class ItOut> 
	friend auto copy_n(gemm_iterator const& first, difference_type count, ItOut d_first)
	->decltype(blas::gemm_n(*std::declval<ContextPtr>(), std::declval<Scalar>(), std::declval<ItA>(), count, std::declval<ItB>(), 0., d_first)) try{
		return blas::gemm_n(*first.ctxtp_              , first.s_              , first.a_it_        , count, first.b_begin_     , 0., d_first);
	}catch(std::exception const& e){
		throw std::logic_error(
			"in " + std::string(__PRETTY_FUNCTION__) + "\nCouldn't decay product of arrays of size " + std::to_string(count) +"x"+ std::to_string(first.a_it_->size()) + " and " + 
			std::to_string(first.a_it_->size())+ "x" +std::to_string(first.b_begin_->size()) + " into " + std::to_string(count) +"x" + std::to_string(first.b_begin_->size()) +
			"\nbecause\n"+e.what()
		);
	}
	
	template<class ItOut>
	friend auto copy(gemm_iterator const& first, gemm_iterator const& last, ItOut d_first){assert(first.s_ == last.s_);
		return copy_n(first, last - first, d_first);
	}

	template<class ItOut>
	friend auto uninitialized_copy_n(gemm_iterator const& first, difference_type count, ItOut d_first){
		return copy_n(first, count, d_first);
	}

	template<class ItOut>
	friend auto uninitialized_copy(gemm_iterator const& first, gemm_iterator const& last, ItOut d_first){assert( first.s_ == last.s_ );
		return uninitialized_copy_n(first, last - first, d_first);}

	reference operator*() const{return {b_begin_->extensions()};}
};

template<class ContextPtr, class Scalar, class ItA, class ItB, class DecayType>
class gemm_range{
	ContextPtr ctxtp_;
	Scalar s_;
	ItA a_begin_;
	ItA a_end_;
	ItB b_begin_;
public:
	gemm_range(gemm_range const&) = delete;
	gemm_range(ContextPtr ctxtp, Scalar s, ItA a_first, ItA a_last, ItB b_first) : ctxtp_{ctxtp}, s_{s}, a_begin_{a_first}, a_end_{a_last}, b_begin_{b_first}{}
	using iterator = gemm_iterator<ContextPtr, Scalar, ItA, ItB>;
	using decay_type = DecayType;
	using size_type = typename decay_type::size_type;
	iterator begin() const{return {ctxtp_, s_, a_begin_, b_begin_};}
	iterator end()   const{return {ctxtp_, s_, a_end_  , b_begin_};}
	friend auto begin(gemm_range const& self){return self.begin();}
	friend auto end  (gemm_range const& self){return self.end  ();}
	size_type size() const{return a_end_ - a_begin_;}
	typename decay_type::extensions_type extensions() const{return size()*b_begin_->extensions();}
	friend auto extensions(gemm_range const& self){return self.extensions();}
//	operator decay_type() const{return decay_type(*this);} // do not use curly { }
	decay_type operator+() const{return *this;}
	template<class Arr>
	friend Arr&& operator+=(Arr&& a, gemm_range const& gr){
		blas::gemm_n(*gr.ctxtp_, gr.s_, gr.a_begin_, gr.a_end_ - gr.a_begin_, gr.b_begin_, 1., a.begin());
		return std::forward<Arr>(a);
	}
};

template<class ContextPtr, class Scalar, class A2D, class B2D, class=std::enable_if_t<is_context<decltype(*ContextPtr{})>{}> >
gemm_range<ContextPtr, Scalar, typename A2D::const_iterator, typename B2D::const_iterator, typename A2D::decay_type/*B2D*/> 
gemm(ContextPtr ctxtp, Scalar s, A2D const& a, B2D const& b){
	return {ctxtp, s, begin(a), end(a), begin(b)};
}

//#pragma warning (disable:1011)
//#pragma diag_suppress 0117 //"implicit_return_from_non_void_function"
//#pragma diag_suppress 940 //"implicit_return_from_non_void_function"
// -Xcudafe "--diag_suppress=implicit_return_from_non_void_function"
template<               class Scalar, class A2D, class B2D> 
auto gemm(                Scalar s, A2D const& a, B2D const& b){
	if constexpr(is_conjugated<A2D>{}){
		auto ctxtp = blas::default_context_of(underlying(a.base()));
		return blas::gemm(ctxtp, s, a, b);
	}else{
		auto ctxtp = blas::default_context_of(a.base());
		return blas::gemm(ctxtp, s, a, b);
	}
}

namespace operators{
	template<class A2D, class B2D> 
	auto operator*(A2D const& A, B2D const& B)
	->decltype(+blas::gemm(1., A, B)){
		return +blas::gemm(1., A, B);}
}

}}}

#endif

