#ifdef COMPILATION_INSTRUCTIONS
(echo "#include\""$0"\"" > $0x.cpp) && clang++ `#-DNDEBUG` -std=c++14 -Wall -Wextra `#-Wfatal-errors` -lboost_timer -I${HOME}/prj -D_TEST_BOOST_MULTI_ARRAY $0x.cpp -o $0x.x && time $0x.x $@ && rm -f $0x.x $0x.cpp; exit
#endif
#ifndef BOOST_MULTI_ARRAY_HPP
#define BOOST_MULTI_ARRAY_HPP

//#include "../multi/array_ref.hpp"
#include "./array_ref.hpp"

#include<algorithm>
#include<array>
//#include<iostream> // cerr
#include<numeric>
#include<vector>

namespace boost{
namespace multi{

template<class T>
inline std::array<index_extension, 1> extensions(std::vector<T> const& v){
	return {index_extension{0, long(v.size())}};
}

template<class TT> auto list_extensions(std::initializer_list<TT> const& il){
	return std::array<index_extension, 1>{index_extension{0, size_type(il.size())}};
}
template<class TT> auto list_extensions(std::initializer_list<std::initializer_list<TT>> il){
	return std::array<index_extension, 2>{
		index_extension{0, size_type(il.size())}, index_extension{0, size_type(il.begin()->size())} 
	};
}
template<class TT> auto list_extensions(std::initializer_list<std::initializer_list<std::initializer_list<TT>>> il){
	return std::array<index_extension, 3>{
		index_extension{0, size_type(il.size())}, index_extension{0, size_type(il.begin()->size())}, index_extension{0, size_type(il.begin()->begin()->size())} 
	};
}

namespace detail{

template<class T> void destroy_at(T* p){p->~T();}

template<class ForwardIt, class Size>
ForwardIt destroy_n(ForwardIt first, Size n){
	for(; n > 0; (void)++first, --n)
		destroy_at(std::addressof(*first));
	return first;
}

template<dimensionality_type N> struct uninitialized_copy_from_initializer_list_aux;

template<dimensionality_type N, class InputIt, class ForwardIt>
ForwardIt uninitialized_copy_from_initializer_list(InputIt first, InputIt last, ForwardIt dest){
	return uninitialized_copy_from_initializer_list_aux<N>::call(first, last, dest);
}

template<dimensionality_type N>
struct uninitialized_copy_from_initializer_list_aux{
	template<class InputIt, class ForwardIt>
	static auto call(InputIt first, InputIt last, ForwardIt dest){
		while(first != last){
			uninitialized_copy_from_initializer_list<N-1>(
				first->begin(), first->end(), dest->begin()
			);
			++first;
			++dest;
		}
		return dest;		
	}
};

template<>
struct uninitialized_copy_from_initializer_list_aux<1u>{
	template<class InputIt, class ForwardIt>
	static auto call(InputIt first, InputIt last, ForwardIt dest){
		while(first != last){
			construct_from_initializer_list(std::addressof(*dest), *first);
		//	using T = typename std::iterator_traits<ForwardIt>::value_type;
		//	::new(static_cast<void*>(std::addressof(*dest))) T(*first);
		    ++first;
		    ++dest;
		}
		return dest;
	}
};

template<dimensionality_type N> struct uninitialized_copy_aux;

template<dimensionality_type N, class InputIt, class ForwardIt>
ForwardIt uninitialized_copy(InputIt first, InputIt last, ForwardIt dest){
	return uninitialized_copy_aux<N>::call(first, last, dest);
}

template<dimensionality_type N>
struct uninitialized_copy_aux{
	template<class InputIt, class ForwardIt>
	static auto call(InputIt first, InputIt last, ForwardIt dest){
		while(first != last){
			uninitialized_copy<N-1>(
				first->begin(), first->end(), dest->begin()
			);
			++first;
			++dest;
		}
		return dest;		
	}
};

template<>
struct uninitialized_copy_aux<1u>{
	template<class InputIt, class ForwardIt>
	static auto call(InputIt first, InputIt last, ForwardIt dest){
		using std::uninitialized_copy;
		return uninitialized_copy(first, last, dest);
	}
};

}

template<class T>
auto extensions(T const& t)
->decltype(t.extensions()){
	return t.extensions();}
inline std::array<index_extension, 0> 
extensions(...){return {};};

template<dimensionality_type> struct extension_aux;

template<dimensionality_type D, class T>
auto extensions(T const& t);

template<dimensionality_type D>
struct extensions_aux{
	template<class T>
	static auto call(T const& t){
		return std::tuple_cat(std::make_tuple(t.extension()), extensions<D-1>(t));
	}
};

template<> struct extensions_aux<0>{
	template<class T> static auto call(T const& ){return std::make_tuple();}
};

template<dimensionality_type D, class T>
auto extensions(T const& t){
	return extensions_aux<D>::call(t);
//	if constexpr(D != 0)
//		return std::tuple_cat(std::make_tuple(t.extension()), extensions<D-1>(t));
//	else
//		return std::make_tuple();
}

template<class T1> struct extensions_t_aux;

template<class T1, class T2> auto extensions_t(T2 const& t2){
	return extensions_t_aux<T1>::call(t2);
}
/*template<class T1, class T2> auto extensions(T2 const& t2){
	if constexpr(std::is_same<T1, T2>{})
		return std::make_tuple();
	else
		return std::tuple_cat(std::make_tuple(t2.extension()), extensions<T1>(*begin(t2)));
}*/

template<class T1> struct extension_t_aux{
	static auto call(T1 const&){return std::make_tuple();}
	template<class T2>
	static auto call(T2 const& t2){return std::tuple_cat(std::make_tuple(t2.extension()), extensions_t<T1>(*begin(t2)));}
};

template<class T, typename = decltype(std::declval<T const&>().layout())>
std::true_type has_layout_member_aux(T const&);
std::false_type has_layout_member_aux(...);

template<class T>
struct has_layout_member : decltype(has_layout_member_aux(std::declval<T const&>())){};

template<class T, typename = std::enable_if_t<has_layout_member<T const&>{}> >
auto layout(T const& t)
->decltype(t.layout()){
	return t.layout();}

template<class T, typename = std::enable_if_t<not has_layout_member<T const&>{}> >
layout_t<0> layout(T const&){return {};};

template<class T, dimensionality_type D, class Alloc>
struct array : array_ref<T, D, typename std::allocator_traits<Alloc>::pointer>{
	using Allocator = Alloc;
	using allocator_type = Alloc;
	using add_lvalue_reference = array_ref<T, D, typename std::allocator_traits<Alloc>::pointer&>;
	using value_type = typename std::conditional<
		D == 1,
		typename array::element,
		array<typename array::element, D-1, allocator_type>
	>::type;
private:
	allocator_type allocator_;
	using alloc_traits = std::allocator_traits<allocator_type>;
public:
	array(array const& other) : array(extensions(other), other.get_allocator()){
		using std::copy_n;
		copy_n(other.data(), other.num_elements(), this->data()); // TODO use uninitialized_fill_n
	}
	array(array const& o, allocator_type const& a) : array(extensions(o), a){
		using std::copy_n;
		copy_n(o.data(), o.num_elements(), this->data()); // TODO use uninitialized_fill_n
	}
	array(array&& other) : 
		array_ref<T, D, typename std::allocator_traits<Alloc>::pointer>{
			std::exchange(other.data_, nullptr), 
			std::exchange(static_cast<layout_t<D>&>(other), {}) // extensions(other)
		},
		allocator_{other.allocator_}
	{}
	array(array&& other, allocator_type const& a) : 
		array_ref<T, D, typename std::allocator_traits<Alloc>::pointer>{
			std::exchange(other.data_, nullptr), 
			std::exchange(static_cast<layout_t<D>&>(other), {}) // extensions(other)
		},
		allocator_{a}
	{}
	template<class Array, typename = decltype(data(std::declval<Array const&>()), extension(std::declval<Array const&>()))>
	array(Array const& o, allocator_type const& a = {}) : array{extensions(o), a}{
		using std::copy_n;
		copy_n(data(o), num_elements(o), this->data());		
	}
	template<class TT, dimensionality_type DD, class... Args>
	array(basic_array<TT, DD, Args...> const& o) 
	: //array{extensions(o), 
		array_ref<T, D, typename array::element_ptr>{nullptr, extensions(o)}, 
		allocator_{pointer_traits<typename array::element_ptr>::allocator_of(o.base())}
	{
std::cout<<" here 0 " <<std::endl;
		this->data_ = alloc_traits::allocate(allocator_, array::num_elements());
std::cout<<" here 1 " <<std::endl;
		array_ref<T, D, typename std::allocator_traits<Alloc>::pointer>::operator=(o);
std::cout<<" here 2 " <<std::endl;
		detail::uninitialized_copy<D>(o.begin(), o.end(), this->begin());
std::cout<<" here 3 " <<std::endl;
	}
	explicit array(allocator_type const& alloc) : array(typename array::extensions_type{}, alloc)
	{}
	array() : 
		array_ref<T, D, typename array::element_ptr>{nullptr, {}}, 
		allocator_{} // TODO, make allocator a base?
	{}
	array(typename array::extensions_type const& e, allocator_type const& alloc = {}) : 
		array_ref<T, D, typename array::element_ptr>{nullptr, e}, 
		allocator_{alloc} // TODO, make allocator a base
	{
		this->data_ = alloc_traits::allocate(allocator_, array::num_elements());
		uninitialized_construct();
	}
	array(typename array::extensions_type e, typename array::element const& el) : 
		array_ref<T, D, typename array::element_ptr>(nullptr, e), allocator_{}
	{
		this->data_ = alloc_traits::allocate(allocator_, array::num_elements());
		uninitialized_construct(el);
	}

	array(typename array::extensions_type e, typename array::element const& el, allocator_type const& a) : 
		array_ref<T, D, typename array::element_ptr>(nullptr, e), allocator_{a}
	{
		this->data_ = alloc_traits::allocate(allocator_, array::num_elements());
		uninitialized_construct(el);
	}
	array(std::initializer_list<typename array_ref<T, D, typename std::allocator_traits<Alloc>::pointer>::value_type> il, allocator_type const& a = {})
	: 
		array_ref<T, D, typename array::element_ptr>{
			nullptr, 
			layout_t<D>(std::tuple_cat(std::make_tuple(index_extension{0, index(il.size())}), extensions<D-1>(*il.begin())))
		},
		allocator_{a}
	{
		this->data_ = alloc_traits::allocate(allocator_, array::num_elements());
		detail::uninitialized_copy_from_initializer_list<D>(il.begin(), il.end(), this->begin());
	}
	allocator_type get_allocator() const{return allocator_;}
//	friend allocator_type const& allocator(array const& s){return s.allocator_;}
	array& operator=(array const& other){
		if(other.extensions() == this->extensions() and other.get_allocator() == get_allocator()){
			array_ref<T, D, typename std::allocator_traits<Allocator>::pointer>::operator=(
				static_cast<const_array_ref<T, D, typename std::allocator_traits<Allocator>::pointer> const&>(other)
			);
		}else{
			array tmp(extensions(other), other.get_allocator());//other.extensions(), allocator_);
			tmp.array_ref<T, D, typename std::allocator_traits<Allocator>::pointer>::operator=(
				static_cast<const_array_ref<T, D, typename std::allocator_traits<Allocator>::pointer> const&>(other)
			);
			swap(tmp);
		}
		return *this;
	}
	array& operator=(array&& other){
		swap(other);
		other.clear();
		return *this;
	}
	template<class Array>
	array& operator=(Array const& a){
		array tmp(extensions(a));
		tmp.array_ref<T, D, typename std::allocator_traits<Allocator>::pointer>::operator=(a);
		swap(tmp);
		return *this;
	}
	void operator=(typename array::initializer_list il){
	//	reextent(list_extensions<typename array::element>(il));
		this->recursive_assign_(il.begin(), il.end());
	}
	void reextent(typename array::extensions_type const& e){
		array tmp(e, allocator_);
		tmp.intersection_assign_(*this);
		swap(tmp);
	}
	friend void reextent(array& self, typename array::extensions_type const& e){self.reextent(e);}
	void swap(array& other) noexcept{
		using std::swap;
		swap(this->data_, other.data_);
		swap(this->allocator_, other.allocator_);
		swap(
			static_cast<typename array_ref<T, D, typename std::allocator_traits<Allocator>::pointer>::layout_t&>(*this), 
			static_cast<typename array_ref<T, D, typename std::allocator_traits<Allocator>::pointer>::layout_t&>(other)
		);
	}
	friend void swap(array& self, array& other){self.swap(other);}
	~array(){clear();}
	void clear(){
		destroy();
		allocator_.deallocate(this->data(), this->num_elements());
		this->data_ = nullptr;
		layout_t<D>::operator=({});
	}
	friend void clear(array& self){self.clear();}
	typename array::reference operator[](index i){
		return array_ref<T, D, typename std::allocator_traits<Allocator>::pointer>::operator[](i);
	}
	typename array::const_reference operator[](index i) const{
		return array_ref<T, D, typename std::allocator_traits<Allocator>::pointer>::operator[](i);
	}
//	typename array::element_ptr data(){return array_ref<T, D, typename std::allocator_traits<Allocator>::pointer>::data();}
//	typename array::element_const_ptr data() const{return array_ref<T, D, typename std::allocator_traits<Allocator>::pointer>::data();}

	typename array::element_ptr origin(){return array_ref<T, D, typename std::allocator_traits<Allocator>::pointer>::origin();}
	typename array::element_const_ptr origin() const{return array_ref<T, D, typename std::allocator_traits<Allocator>::pointer>::origin();}
private:
	void destroy(){
		using boost::multi::detail::destroy_n;	//C++11 workaround for C++17's `using std::destroy_n;`
		destroy_n(this->data(), this->num_elements());
	}
	template<class... Args>
	auto uninitialized_construct(Args&&... args){
		using std::uninitialized_fill_n;
		return uninitialized_fill_n(
			this->data_, this->num_elements(), 
			typename array::element(std::forward<Args>(args)...)
		);
	}
};

//template<class T, dimensionality_type D, class Alloc>
//array<T, D, Alloc>::array(size_type, array<T, D-1> const& arr, allocator_type const& a){}
//template<class T, class Alloc>
//array<T, 1u, Alloc>::array(size_type, T const& e, Alloc const& a){}

}}

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

int main(){
	{
		multi::array<double, 1> A1({{3}}, 1.1); 
		multi::array<double, 2> A2({{3,4}}, 1.1);
		multi::array<double, 1> AA1 = {1.2, 3.4, 5.6};
		multi::array<double, 2> AA = {
			{0.0, 0.1, 0.2, 0.3},
			{1.0, 1.1, 1.2, 1.3},
			{2.0, 2.1, 2.2, 2.3}
		};
		multi::array<double, 3> AAA = {
			{{0.00,0.01},
			 {0.10,0.11}},
			{{1.00,1.01},
			 {1.10,1.11}}
		};
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
}
#endif
#endif

