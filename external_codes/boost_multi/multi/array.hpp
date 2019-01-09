#ifdef COMPILATION_INSTRUCTIONS
(echo "#include\""$0"\"" > $0x.cpp) && c++ `#-DNDEBUG` -std=c++17 -Wall -Wextra `#-Wfatal-errors` -D_TEST_BOOST_MULTI_ARRAY $0x.cpp -o $0x.x && time $0x.x $@ && rm -f $0x.x $0x.cpp; exit
#endif
#ifndef BOOST_MULTI_ARRAY_HPP
#define BOOST_MULTI_ARRAY_HPP

#include "./array_ref.hpp"

#include "detail/memory.hpp"

#include "utility.hpp"

namespace boost{
namespace multi{

template<int N> struct priority : priority<N-1>{};
template<> struct priority<0>{};

namespace detail{
    using std::begin;
    template<class C> auto maybestd_begin(C&& c) 
    ->decltype(begin(std::forward<C>(c))){
        return begin(std::forward<C>(c));}
	using std::end;
    template<class C> auto maybestd_end(C&& c) 
    ->decltype(end(std::forward<C>(c))){
        return end(std::forward<C>(c));}	
}

template<class C> auto maybestd_begin(C&& c)
->decltype(detail::maybestd_begin(std::forward<C>(c))){
    return detail::maybestd_begin(std::forward<C>(c));}
template<class C> auto maybestd_end(C&& c)
->decltype(detail::maybestd_end(std::forward<C>(c))){
    return detail::maybestd_end(std::forward<C>(c));}

template<class T, dimensionality_type D, class Alloc>
struct array : 
	private Alloc, 
	array_ref<T, D, typename std::allocator_traits<Alloc>::pointer>, 
	boost::multi::random_iterable<array<T, D, Alloc>>
{
	using allocator_type = Alloc;
	using ref = array_ref<T, D, typename std::allocator_traits<Alloc>::pointer>;
//	using add_lvalue_reference = array_ref<T, D, typename std::allocator_traits<Alloc>::pointer&>;
//	using value_type = typename std::conditional<
//		D == 1,
//		typename array::element,
//		array<typename array::element, D-1, allocator_type>
//	>::type;
private:
//	allocator_type allocator_;
	using alloc_traits = std::allocator_traits<allocator_type>;
	using typename ref::reference;
public:
	using typename ref::value_type;
	using ref::operator<;
	array() noexcept(noexcept(Alloc())) : Alloc{}, ref{{}, nullptr}{}           //1a
	explicit array(Alloc const& a) : Alloc{a}, ref{{}, nullptr}{}               //1b
	array(typename array::extensions_type const& e, typename array::element const& el, Alloc const& a = {})  //2
	:	Alloc{a} , ref{e, allocate(typename array::layout_t{e}.num_elements())}{
		uninitialized_fill(el);
	}
	array(typename array::index_extension n, value_type const& v, Alloc const& a = {})
	: 	Alloc{a}, 
		ref{
			std::tuple_cat(std::make_tuple(n), multi::extensions(v)),
			allocate(layout_t<D>{std::tuple_cat(std::make_tuple(n), multi::extensions(v))}.num_elements())
		}
	{
		for(auto it = this->begin(); it != this->end(); ++it) *it = v;
	}
	array(typename array::size_type n, value_type const& v, Alloc const& a)
	: 	array(typename array::index_extension(n), v, a){}
	array(typename array::size_type n, value_type const& v)
	: 	array(typename array::index_extension(n), v){}

	explicit array(typename array::extensions_type const& e, Alloc const& a={}) //3
	:	Alloc{a}, ref{e, allocate(typename array::layout_t{e}.num_elements())}{
		uninitialized_value_construct();
	}
	template<class It> static auto distance(It a, It b){using std::distance; return distance(a, b);}
	template<class It, typename=decltype(std::distance(std::declval<It>(), std::declval<It>()), *std::declval<It>())>      
	array(It first, It last, allocator_type const& a = {}) :                    //4
		Alloc{a}, 
		ref{
			std::tuple_cat(std::make_tuple(index_extension{index{distance(first, last)}}), multi::extensions(*first)),
			allocate(layout_t<D>{std::tuple_cat(std::make_tuple(index_extension{index{distance(first, last)}}), multi::extensions(*first))}.num_elements())
		}
	{
		using std::next;
		using std::all_of;
		if(first!=last) assert( all_of(next(first), last, [x=multi::extensions(*first)](auto& e){return extensions(e)==x;}) );
		multi::uninitialized_copy<D>(first, last, ref::begin());
	}
	template<class Array, typename=std::enable_if_t<!std::is_base_of<array, Array>{}>, typename=std::enable_if_t<std::rank<std::decay_t<Array>>{}==D> >
	array(Array&& other, allocator_type const& a = {})
	:	Alloc{a}, ref{extensions(other), allocate(num_elements(other))}{
		using std::begin; using std::end;
		multi::uninitialized_copy<D>(begin(other), end(other), ref::begin());
	}
	array(array const& other)                                                   // 5a
	:	Alloc{other}, ref{extensions(other), allocate(other.num_elements())}{
		uninitialized_copy(other.data());
	}
	array(array const& other, allocator_type const& a)                          // 5b
	:	Alloc{a}, ref{extensions(other), allocate(other.num_elements())}{
		uninitialized_copy(other.data());
	}
	array(array&& other) noexcept                                              //6a
	:	Alloc{static_cast<Alloc const&>(other)},
		ref{extensions(other), std::exchange(other.base_, nullptr)}
	{
		other.layout_t<D>::operator=({});
	}
	array(array&& other, allocator_type const& a)                             //6b
	:	Alloc{a},
		ref(other.extensions(), std::exchange(other.base_, nullptr))
	{
		//TODO
		other.layout_t<D>::operator=({});
	}
	array(std::initializer_list<value_type> il, allocator_type const& a={}) 
	:	array(il.begin(), il.end(), a){}

	template<class A>
	array& operator=(A&& a){
		auto ext = extensions(a);
		if(ext==array::extensions()){
			ref::operator=(std::forward<A>(a));
		}else{
			clear();
			this->layout_t<D>::operator=(layout_t<D>{extensions(a)});
			this->base_ = allocate(this->num_elements());
			multi::uninitialized_copy<D>(maybestd_begin(std::forward<A>(a)), maybestd_end(std::forward<A>(a)), array::begin());
		}
		return *this;
	}
	array& operator=(array const& other){return operator=<array const&>(other);}
	void swap(array& other) noexcept{
		using std::swap;
		swap(this->base_, other.base_);
		swap(static_cast<Alloc&>(*this), static_cast<Alloc&>(other));
		swap(
			static_cast<typename array_ref<T, D, typename std::allocator_traits<Alloc>::pointer>::layout_t&>(*this), 
			static_cast<typename array_ref<T, D, typename std::allocator_traits<Alloc>::pointer>::layout_t&>(other)
		);
	}
	friend void swap(array& a, array& b){a.swap(b);}
	array& operator=(array&& other){clear(); swap(other); return *this;}
	template<class It>
	array& assign(It first, It last){
		using std::next;
		using std::all_of;
		if(distance(first, last) == array::size() and multi::extensions(*first) == multi::extensions(*array::begin())){
			array_ref<T, D, typename std::allocator_traits<Alloc>::pointer>::assign(first, last);
		}else{
			clear();
			this->layout_t<D>::operator=(layout_t<D>{std::tuple_cat(std::make_tuple(index_extension{array::extension().front(), array::extension().front() + distance(first, last)}), multi::extensions(*first))});
			using std::next;
			using std::all_of;
			if(first!=last) assert( all_of(next(first), last, [x=multi::extensions(*first)](auto& e){return extensions(e)==x;}) );
			this->base_ = allocate();
			multi::uninitialized_copy<D>(first, last, array::begin());
		}
		return *this;
	}
	array& operator=(std::initializer_list<value_type> il){return assign(begin(il), end(il));}
	void reextent(typename array::extensions_type const& e, typename array::element const& v = {}){
		array tmp(e, v, static_cast<Alloc const&>(*this));
		tmp.intersection_assign_(*this);
		swap(tmp);
	}
	allocator_type get_allocator() const{return static_cast<allocator_type const&>(*this);}

	using element_const_ptr = typename std::pointer_traits<typename array::element_ptr>::template rebind<typename array::element const>;
	using const_reference = std::conditional_t<
		array::dimensionality!=1, 
		basic_array<typename array::element, array::dimensionality-1, element_const_ptr>, 
	//	typename std::iterator_traits<element_ptr>::reference const	//	
		typename pointer_traits<typename array::element_ptr>::element_type const&
	>;
	using const_iterator = typename std::conditional<
		array::dimensionality != 1,
		basic_array_ptr<const_reference, typename layout_t<D>::sub_t>,
		typename basic_array<typename array::element, dimensionality_type{1}, typename array::element_ptr>::template basic_iterator<element_const_ptr, const_reference>
	>::type;

	reference       operator[](index i)      {return ref::operator[](i);}
	typename array::const_reference operator[](index i) const{return ref::operator[](i);}

	typename array::element_ptr       data()      {return ref::data();}
	typename array::element_const_ptr data() const{return ref::data();}
	friend typename array::element_ptr       data(array&       s){return s.data();}
	friend typename array::element_const_ptr data(array const& s){return s.data();}

	typename array::element_ptr       origin()      {return ref::origin();}
	typename array::element_const_ptr origin() const{return ref::origin();}
	friend typename array::element_ptr       origin(array&       s){return s.origin();}
	friend typename array::element_const_ptr origin(array const& s){return s.origin();}

//	using const_reverse_iterator = basic_reverse_iterator<const_iterator>;

	typename array::iterator begin(){return ref::begin();}
	typename array::iterator end()  {return ref::end();}
//	typename array::iterator begin() &&{return ref::begin();}
//	typename array::iterator end()   &&{return ref::end();}

	typename array::iterator begin() const{return ref::begin();}
	typename array::iterator end()   const{return ref::end();}
#if 0
	friend auto begin(array&       s){return s.begin();}
	friend auto begin(array const& s){return s.begin();}
	friend auto end  (array&       s){return s.end();}
	friend auto end  (array const& s){return s.end();}

	auto cbegin() const{return begin();}
	auto cend()   const{return end();}
	friend auto cbegin(array const& s){return s.cbegin();}
	friend auto cend  (array const& s){return s.cend();}

	typename array::reverse_iterator rbegin(){return ref::rbegin();}
	typename array::reverse_iterator rend() {return ref::rend();}
	typename array::const_reverse_iterator rbegin() const{return ref::rbegin();}
	typename array::const_reverse_iterator rend() const{return ref::rend();}
#endif

#if 0
	friend auto begin(array const& s){return s.begin();}
	friend auto end(array const& s){return s.end();}
	friend auto begin(array& s){return s.begin();}
	friend auto end(array& s){return s.end();}
	typename array::const_reference front() const{assert(not this->empty()); return *this->cbegin();}
	typename array::const_reference back() const{assert(not this->empty()); return *(this->cbegin() + (this->size() - 1));}
	typename array::reference front(){return array_ref<T, D, typename std::allocator_traits<Alloc>::pointer>::front();}
	typename array::reference back(){return array_ref<T, D, typename std::allocator_traits<Alloc>::pointer>::back();}
#endif
	void clear() noexcept{
		destroy(); 
		deallocate();
		layout_t<D>::operator=({});
	}
	friend void clear(array& self) noexcept{self.clear();}
	~array() noexcept{clear();}	
private:
	void destroy(){destroy_n(this->data(), this->num_elements());}
	template<typename It>
	auto uninitialized_copy(It it){
		using std::uninitialized_copy_n;
		return uninitialized_copy_n(it, this->num_elements(), this->data());
	}
	auto uninitialized_default_construct(){
	//	using std::uninitialized_default_construct_n;
		return uninitialized_default_construct_n(this->base_, this->num_elements());
	}
	auto uninitialized_value_construct(){
	//	using std::uninitialized_value_construct_n;
		return uninitialized_value_construct_n(this->base_, this->num_elements());
	}
	auto uninitialized_fill(typename array::element const& el){
		using std::uninitialized_fill_n;
		return uninitialized_fill_n(this->base_, this->num_elements(), el);
	}
	typename array::element_ptr allocate(typename array::size_type n){
		return alloc_traits::allocate(*this, n);
	}
	auto allocate(){return allocate(this->num_elements());}
	void deallocate(){
		alloc_traits::deallocate(*this, this->base_, this->num_elements());
		this->base_ = nullptr;
	}
};

//template<class T, dimensionality_type D, class Alloc>
//array<T, D, Alloc>::array(size_type, array<T, D-1> const& arr, allocator_type const& a){}
//template<class T, class Alloc>
//array<T, 1u, Alloc>::array(size_type, T const& e, Alloc const& a){}

#if __cpp_deduction_guides
#define IL std::initializer_list
// clang cannot recognize templated-using, so don't replace IL<T>->IL1<T>
template<class T, class A=std::allocator<T>> array(IL<T>                , A={})->array<T,1,A>; 
template<class T, class A=std::allocator<T>> array(IL<IL<T>>            , A={})->array<T,2,A>;
template<class T, class A=std::allocator<T>> array(IL<IL<IL<T>>>        , A={})->array<T,3,A>; 
template<class T, class A=std::allocator<T>> array(IL<IL<IL<IL<T>>>>    , A={})->array<T,4,A>; 
template<class T, class A=std::allocator<T>> array(IL<IL<IL<IL<IL<T>>>>>, A={})->array<T,5,A>; 
#undef IL
#endif

}}

namespace std{
	template<class T, boost::multi::dimensionality_type N, class... Ts> 
	struct rank<boost::multi::array<T, N, Ts...>> 
	: public std::integral_constant<
		boost::multi::dimensionality_type, 
		boost::multi::array<T, N, Ts...>::rank
	>{};
}

#if _TEST_BOOST_MULTI_ARRAY

#include<cassert>
#include<numeric> // iota
#include<iostream>
#include<algorithm>
#include<vector>

#include <random>
#include <boost/timer/timer.hpp>
#include<boost/multi_array.hpp>
using std::cout;
namespace multi = boost::multi;

#if 1
template<class Matrix, class Vector>
void solve(Matrix& m, Vector& y){
//	using std::size; // assert(size(m) == std::ptrdiff_t(size(y)));
	std::ptrdiff_t msize = size(m); 
	for(auto r = 0; r != msize; ++r){ //	auto mr = m[r]; //  auto const mrr = mr[r];// assert( mrr != 0 ); // m[r][r] = 1;
		auto mr = m[r];
		auto mrr = mr[r];
		for(auto c = r + 1; c != msize; ++c) mr[c] /= mrr;
		auto yr = (y[r] /= mrr);
		for(auto r2 = r + 1; r2 != msize; ++r2){ //	auto mr2 = m[r2]; //	auto const mr2r = mr2[r]; // m[r2][r] = 0;
			auto mr2 = m[r2];
			auto const& mr2r = mr2[r];
			auto const& mr = m[r];
			for(auto c = r + 1; c != msize; ++c) mr2[c] -= mr2r*mr[c];
			y[r2] -= mr2r*yr;
		}
	}
	for(auto r = msize - 1; r > 0; --r){ //	auto mtr = m.rotated(1)[r];
		auto const& yr = y[r];
		for(auto r2 = r-1; r2 >=0; --r2)
			y[r2] -= yr*m[r2][r];
	}
}
#endif

void f(boost::multi::array<double, 4> const& A){
	A[1][2];
	auto&& a = A[1][2]; (void)a; // careful, a is a reference here, don't use auto, 
	auto const& b = A[1][2]; (void)b; // use auto const& if possible
//	A[1][2][3][4] = 5; // fail, element is read-only
}

template<class C>
void set_99(C&& c){
	for(auto j : c.extension(0))
		for(auto k : c.extension(1))
				c[j][k] = 99.;
}

namespace multi = boost::multi;

template<class T> void fun(T const& t){
	std::cout << typeid(t).name() << std::endl;
}

template<class T> struct extension{};
template<class T> void gun(extension<T>){
	std::cout << typeid(T).name() << std::endl;
}

typedef double a1010[10][10];

struct A{
	double const* p;
	A(std::initializer_list<double> il){ p = &*(il.begin() + 1); };
};

double f(){return 5.;}
int main(){

#if __cpp_deduction_guides
{
	multi::array<double, 1> A1 = {1.,2.,3.}; 
					assert(A1.dimensionality==1 and A1.num_elements()==3);
	multi::array<double, 2> A2 {
		 {1.,2.,3.},
		 {4.,5.,6.}
	};
	*A2.begin()->begin() = 99;
	assert(A2[0][0] == 99 );
//					assert(A2.dimensionality==2 and A2.num_elements()==2*3);
//	multi::array A3 = {
//		{{ 1.2,  0.}, { 2.4, 1.}},
//		{{11.2,  3.}, {34.4, 4.}},
//		{{15.2, 99.}, {32.4, 2.}}
//	};
//					assert(A3.dimensionality==3 and A3.num_elements()==3*2*2);
}
#endif


#if 0
A a{1.,2.,3.};
assert( *a.p == 2. );

A aa{1., f(), f()};

assert( *aa.p == 5. );

#if __cpp_deduction_guides
{
multi::array A1 = {1.,2.,3.}; 
                     assert(A1.dimensionality==1 and A1.num_elements()==3);
multi::array A2 {
	 {1.,2.,3.},
	 {4.,5.,6.}
};                   assert(A2.dimensionality==2 and A2.num_elements()==2*3);
multi::array A3 = {
    {{ 1.2,  0.}, { 2.4, 1.}},
    {{11.2,  3.}, {34.4, 4.}},
    {{15.2, 99.}, {32.4, 2.}}
};                   assert(A3.dimensionality==3 and A3.num_elements()==3*2*2);
}
#endif


{
	multi::array<double, 1> const A = {1.2,3.4,5.6}; assert( size(A) == 3 );
	assert( A[2] == 5.6 );
}
{
	multi::array<double, 2> const A = {
		multi::array<double, 1>{ 1.2,  2.4, 3.6},
		multi::array<double, 1>{11.2, 34.4, 5.6},
		multi::array<double, 1>{15.2, 32.4, 5.6}
	};
//	multi::array<double, 2> Acopy{cbegin(A), cend(A)};
//	assert( A == Acopy );
//	assert( size(A) == 3 and size(A[0]) == 3 );
//	assert( A[1][1] == 34.4 );
}
{
	multi::array<double, 2> const A = {
		{ 1.2, 2.4 },
		{11.2,34.4},
		{15.2,32.4}
	};
	assert( size(A) == 3 and size(A[0]) == 2 );
	assert( A[1][1] == 34.4 );
}
{
	multi::array<double, 3> const A = {
		{{ 1.2,  0.}, { 2.4, 1.}},
		{{11.2,  3.}, {34.4, 11.}},
		{{15.2, 99.}, {32.4, 2.}}
	};
	assert( A[1][1][0] == 34.4 );
	assert( A[1][1][1] == 11.   );
}
#if 0
{
	multi::array<double, 2> A = {
		multi::array<double, 1>{1.2,3.4,5.6},
		multi::array<double, 1>{1.2,3.4,5.6},
		multi::array<double, 1>{1.2,3.4,5.6}
	};
	std::initializer_list<multi::array<double, 1>> il = {
		multi::array<double, 1>{1.2,3.4,5.6},
		multi::array<double, 1>{1.2,3.4,5.6},
		multi::array<double, 1>{1.2,3.4,5.6}
	};
//	multi::array<double, 1> AA1({boost::multi::index_extension{3}}); 
//	assert( *(il.begin()->begin() + 1) == 3.4 );
	assert( size(A) == 3 );
//	std::initializer_list<double> il = {1.2,3.4,5.6};
//	boost::multi::detail::uninitialized_copy_from_il<2>::call(il.begin(), il.end(), A.begin());
//	cout<< A[0][0] <<std::endl;
//	assert( A[1][1] == 3.4 );
//	auto it = A.begin();
//	[]{}(it, *it, (*it).begin(), it->begin());
//	auto it2 = it->begin();
//	*it2 = 2.;
//	auto it2 = (*it).begin();
//	*(A.begin()->begin()) = 444.;
}
#if 0
	{
		multi::array<double, 1> A1({{3}}, 1.1); 
		multi::array<double, 2> A2({3,4}, 1.1);
		multi::array<double, 1> AA1 = {1.2, 3.4, 5.6};
		std::initializer_list<double> il = {1.2,3.4,5.6};
		boost::multi::detail::uninitialized_copy_from_il<1>::call(il.begin(), il.end(), AA1.begin());
		assert( AA1[1] == 3.4 );
		multi::array<double, 2> AA = {
			{0.0, 0.1, 0.2, 0.3},
			{1.0, 1.1, 1.2, 1.3},
			{2.0, 2.1, 2.2, 2.3}
		};
		std::initializer_list<multi::array<double, 1>> il2 = {
			multi::array<double, 1>{0.0, 0.1, 0.2, 0.3},
			multi::array<double, 1>{1.0, 1.1, 1.2, 1.3},
			multi::array<double, 1>{2.0, 2.1, 2.2, 2.3}
		};
		cout<< *((il2.begin()+1)->begin()) <<std::endl; assert(0);

		assert(*((il2.begin()+1)->begin())==1.0);
		boost::multi::detail::uninitialized_copy_from_il<2>::call(il.begin(), il.end(), AA.begin());
		assert( AA[1][1] == 1.1 );
		multi::array<double, 3> AAA = {
			{{0.00,0.01},
			 {0.10,0.11}},
			{{1.00,1.01},
			 {1.10,1.11}}
		};
		assert( size(AAA) == 2 );
		cout<< AAA[1][0][1] <<std::endl;
		assert( AAA[1][0][1] == 1.01 );
//		multi::array<double, 3> BBB(2,
//			{{0.00,0.01},
//			 {0.10,0.11}}
//		);
//		assert( BBB[0][1][0] == 0.10 );
//		assert( BBB[1][1][0] == 0.10 );

		assert( A1[1] == 1.1);
		assert( AA[1][2] == 1.2 );	
//		multi::array<double, 3> A3(10, {20, {0.0, 0.1, 0.2}});
//		assert( A3.size(0) == 10 and A3.size(1) == 20 and A3.size(2) == 3 );
//		assert( A3[8][11][1] == 0.1 );
//		assert(A1.size() == 10 and A1[0]==1.1);
//		assert(A2.size() == 10 and A2[0].size() == 20 and A2[0][0]==1.1);
//		assert(A3.size() == 10 and A3[0].size() == 20 and A3[0][0].size()==30 and A2[0][0]==1.1);
	}
	{
		multi::array<double, 2> A({10, 20});
		assert( A.num_elements() == 200 );
	//	auto const& B = 
		//	A[indices[range(2,8)][range(4,16)]];                // Boost.MultiArray only syntax
		//	A.sliced(2,8).rotated(1).sliced(4, 16).rotated(-1); // advanced interface
	//		A({2,8}, {4, 16});                                  // simple interface
	//	typedef boost::multi_array_types::index_range range;
		auto const& B = 
		//	A[boost::indices[range(2,8)][range(4,16)]];                // Boost.MultiArray only syntax
		//	A.sliced(2,8).rotated(1).sliced(4, 16).rotated(-1); // advanced interface
			A({2,8}, {4, 16}); // simple interface
		auto const& C = 
		//	A[indices[range()][range(4,16)]];      // Boost.MultiArray only syntax
		//	A.rotated(1).sliced(4, 16).rotated(-1) // advanced interface
			A( A.extension(0), {4, 16} ) // simple interface
		;	
		assert(A.size(0) == 10 and A.size(1) == 20);
		assert(B.size(0) == 6 and B.size(1) == 12);
		assert(C.size(0) == 10 and B.size(1) == 12);
		assert(&B[0][0] == &A[2][4] and &B[5][11] == &A[7][15]);
		assert(&C[0][0] == &A[0][4] and &C[9][11] == &A[9][15]);
	}
	{
		multi::array<double, 3> A({30,40,50});

		for(auto i : A.extension(0))
			for(auto j : A.extension(1))
				for(auto k : A.extension(2))
					A[i][j][k] = 5.1;

		set_99(A[1]);
		for(auto j : A.extension(1))
			for(auto k : A.extension(2))
				A[1][j][k] = 99.;
	}	
	assert(( boost::multi::list_extensions({1.,2.,3.})[0] == boost::multi::index_range{0, 3} ));
	assert(( boost::multi::list_extensions(
		{{1.,2.,3.}, {1.,2.,3.}}
	)[0] == boost::multi::index_range{0, 2} ));
	assert(( boost::multi::list_extensions(
		{{1.,2.,3.}, {1.,2.,3.}}
	)[1] == boost::multi::index_range{0, 3} ));

//	assert(( boost::multi::list_extensions(
//		{{1.,2.,3.}, {1.,2.,3.}}
//	)[1] == boost::multi::index_range{0, 3} ));

//	assert((boost::multi::extensions<1>({1.,2.})[0] == boost::multi::index_range{0,2}));
//	extensions<double>( {{1.,2.},{3.,4.}} );


//	assert( MA0[1][2] == 12. );

	multi::array<double, 2> MA0({
		multi::index_extension{0,3},
		multi::index_extension{0,3}
	});
	assert(size(MA0) == 3 and size(MA0[0])==3);
	MA0 = {
		{0.1, 01., 02.},
		{10., 11., 12.},
		{20., 21., 22.}
	};
	
	multi::array<double, 1> VV0({{{0, 3}}}); assert(size(VV0) == 3); 
	VV0 = {1.,2.,3.};
	multi::array<double, 1> VV1({{0, 3}}); assert(size(VV1) == 2); 
	VV1 = {0.,3.};
	multi::array<double, 1> VV2({0, 3}); assert(size(VV2) == 2); 
	VV2 = {0.,3.};
	multi::array<double, 1> VV3 = {0, 3}; assert(size(VV3) == 2); 
	VV2 = {0.,3.};

	return 0;


	for(auto i : MA0.extension(0)){
		for(auto j : MA0.extension(1)) 
			cout << MA0[i][j] << ' ';
		cout << '\n';
	}
	MA0.reextent({10,10});
	
	for(auto i : MA0.extension(0)){
		for(auto j : MA0.extension(1)) 
			cout << MA0[i][j] << ' ';
		cout << '\n';
	}

	return 0;
/*	for(auto i : MA0.extension(0)){
		for(auto j : MA0.extension(1))
			cout << MA0[i][j] << '\t';
		cout <<"\t|"<< VV0[i] <<'\n';
	}
	cout << "--\n";*/

/*	double MA00[3][3] = {
		{0.1, 01., 02.},
		{10., 11., 12.},
		{20., 21., 22.}
	};*/
	multi::array<double, 2> MA00({{0,3},{0,3}}); MA00 = {
		{0.1, 01., 02.},
		{10., 11., 12.},
		{20., 21., 22.}
	};
	std::vector<double> VV00 = {1.,2.,3.};
	solve(MA00, VV00);

	for(int i = 0; i != int(VV00.size()); ++i){//VV0.extension(0)){
	//	for(auto j : MA0.extension(1))
	//		cout << MA0[i][j] << '\t';
		cout <<"\t|"<< VV00[i] << std::endl;
	}
	for(auto n = 1; n < 10000; n = n*2)
	{
		multi::array<double, 2> MA( {multi::index_extension{0, n}, multi::index_extension{0, n}} );
		assert( MA.size() == n and MA[0].size() == n );
//		boost::multi_array<double, 2> MA(boost::extents[3000][3000]);
		std::vector<double> V(n);
    	std::mt19937 rng;
		std::generate_n(MA.data(), MA.num_elements(), [&]{return rng();});
		std::generate_n(V.data(), V.size(), [&]{return rng();});
		{
			auto MA2 = MA;
			auto V2 = V;
			solve(MA2, V2);
			boost::timer::cpu_timer timer;
	//		boost::timer::auto_cpu_timer t;
			solve(MA, V);
			cout << n <<'\t' << timer.elapsed().user << std::endl;
		}
//		cout << "some " << V[13] << std::endl;
	}

	return 0;

	multi::array<double, 2> MA({{0,2}, {0, 3}}); assert( MA.size()==2 and MA[0].size()== 3 );
	MA[1][1] = 11.;
	assert( MA[1][1] == 11.);
	multi::array<double, 2> MA2({{0, 4}, {0, 5}}); assert(MA2.size()==4 and MA2[0].size()==5);
	using std::swap;
	swap(MA, MA2);
	assert( MA.size() == 4 );
	assert( MA2.size() == 2 );
	cout << MA2[1][1] << std::endl;
	assert( MA2[1][1] == 11. );
	multi::array<double, 2> MA3 = MA2;//({2, 3});
//	MA3 = MA2;
	cout << MA3[1][1] << std::endl;
	assert(MA3[1][1] == 11.);
#if 0
	multi::array<double, 2> MAswap({4, 5});
	multi::array<double, 1> MA1({3});
	using std::swap;
	swap(MA, MAswap);
	assert(MA[2][2] == 0.);
	MA[1][3] = 7.1;
	assert(MA[1][3] == 7.1);
	cout << MA.stride() << '\n';	
	cout << MA.strides()[0] << '\n';
	cout << MA.strides()[1] << '\n';
	cout << "s = " << MA.size() << std::endl;
	assert( MA.size() == 4 );
	assert( MA.size(0) == 4 );
	assert( MA.size(1) == 5 );

	double d2D[4][5] {
		{ 0,  1,  2,  3,  4}, 
		{ 5,  6,  7,  8,  9}, 
		{10, 11, 12, 13, 14}, 
		{15, 16, 17, 18, 19}
	};
	multi::array_ref<double, 2> d2D_ref{&d2D[0][0], {4, 5}};
	d2D_ref[1][1] = 8.1;
	multi::array<double, 2> MA2;
	MA2 = d2D_ref;
	d2D_ref[1][1] = 8.2;
	assert(MA2.extensions() == d2D_ref.extensions());
	cout << MA2[1][1] << std::endl;
	assert(MA2[1][1] == 8.2);
#endif
#endif
#endif
#endif
}
#endif
#endif
//	array(std::initializer_list<typename array::value_type> il, allocator_type const& a = {}) : 
//		array_ref<T, D, typename array::element_ptr>{
//			nullptr,
//			layout_t<D>{std::tuple_cat(std::make_tuple(index_extension{index(il.size())}), extensions<D-1>(*il.begin()))}
//		},
//		allocator_{a}
//	{
//		if(il.size()) assert( std::all_of(il.begin()+1, il.end(), [x=extensions<D-1>(*il.begin())](auto& e){return extensions<D-1>(e)==x;}) );
//		this->base_ = alloc_traits::allocate(allocator_, array::num_elements());
//		detail::uninitialized_copy_from_il<D>::call(il.begin(), il.end(), array::begin()); // don't use begin(il) because begin() is a member
//	}

