#ifdef COMPILATION_INSTRUCTIONS
$CXX -O3 -std=c++14 -Wall -Wextra -Wpedantic `#-Wfatal-errors` $0 -o $0.x && time $0.x $@ && rm -f $0.x; exit
#endif

#include<iostream>
#include<cassert>

namespace fancy{

template<class T = void> struct ptr{ // minimal fancy ptr
	static double value;
	using difference_type = std::ptrdiff_t;
	using value_type = T;
	using element_type = T;
	using pointer = T*;
	using reference = T&;
	using iterator_category = std::random_access_iterator_tag;
	reference operator*() const{return value;}
	ptr operator+(difference_type) const{return *this;}
	ptr& operator+=(difference_type){return *this;}
	friend difference_type operator-(ptr const&, ptr const&){return 0;}
	bool operator==(ptr const&) const{return true;}
	explicit operator T*() const{return &value;}
};
#if 1
template<> struct ptr<void>{ // minimal fancy ptr
	static double value;
	using difference_type = std::ptrdiff_t;
	using value_type = void;
	using element_type = void;
	using pointer = void*;
	using reference = void;
	using iterator_category = std::random_access_iterator_tag;
	double& operator*() const{return value;}
	ptr operator+(difference_type) const{return *this;}
	ptr& operator+=(difference_type){return *this;}
	friend difference_type operator-(ptr const&, ptr const&){return 0;}
	bool operator==(ptr const&) const{return true;}
};
double ptr<void>::value = 42.;
#endif
template<> double ptr<double>::value = 42.;

// all these are optional, depending on the level of specialization needed
template<class Ptr, class T, class Size>
ptr<T> copy_n(Ptr, Size, ptr<T> d){ // custom copy_n, Boost.Multi uses copy_n
	std::cerr << "called Pointer-based copy_n(Ptr, n, fancy::ptr)" << std::endl; 
	return d;
}
template<class Ptr, class T, class Size>
Ptr copy_n(ptr<T>, Size, Ptr d){ // custom copy_n, Boost.Multi uses copy_n
	std::cerr << "called Pointer-based copy_n(fancy::ptr, n, Ptr)" << std::endl; 
	return d;
}
template<class T1, class T2, class Size>
ptr<T2> copy_n(ptr<T1>, Size, ptr<T2> d){ // custom copy_n, Boost.Multi uses copy_n
	std::cerr << "called Pointer-based copy_n(fancy::ptr, n, fancy::ptr)" << std::endl; 
	return d;
}

}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// multi-fancy glue, where should this be? In boost/multi/adaptors/MyFancyApaptor.hpp if anything, or in user code if it is very special
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#include "../array.hpp"

namespace boost{
namespace multi{

template<class It, class T>  // custom copy 1D (aka strided copy)
void copy(It first, It last, multi::array_iterator<T, 1, fancy::ptr<T> > dest){
	assert( stride(first) == stride(last) );
	std::cerr << "1D copy(it1D, it1D, it1D) with strides " << stride(first) << " " << stride(dest) << std::endl;
}

template<class It, class T> // custom copy 2D (aka double strided copy)
void copy(It first, It last, multi::array_iterator<T, 2, fancy::ptr<T> > dest){
	assert( stride(first) == stride(last) );
	std::cerr << "2D copy(It, It, it2D) with strides " << stride(first) << " " << stride(dest) << std::endl;
}

}}

////////////////////////////////////////////////////////////////////////////////
// user code
////////////////////////////////////////////////////////////////////////////////

using std::cout; using std::cerr;
namespace multi = boost::multi;

int main(){

	using ptr = fancy::ptr<double>;
	ptr p1;// = new double[25];
	multi::array_ref<double, 2, ptr> A(p1, multi::iextensions<2>{5, 5});
	
	auto AA = multi::static_array_cast<double, double*>(A); (void)AA;
	
	multi::array_ref<double, 1, ptr> A1(p1, multi::iextensions<1>{25});
	auto AA1 = multi::static_array_cast<double, double*>(A1); (void)AA1;

	ptr p2;// = new double[25]; 
	multi::array_ref<double, 2, ptr> B(p2, multi::iextensions<2>{5, 5});

	A = B;                   // called fancy::copy_n;
	rotated(A)[1] = B[1];    // called multi::copy 1D multi::iterator<...>
	rotated(A) = rotated(B); // called multi::copy 2D multi::iterator<...>

}

