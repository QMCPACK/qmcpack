#ifdef COMPILATION// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4-*-
$CXX $0 -o $0x&&$0x&&rm $0x;exit
#endif
// © Alfredo A. Correa 2018-2020

#include<iostream>
#include<cassert>

namespace fancy{

template<class T> struct ref;

template<class T = void> class ptr{
	static double value;
public:
	using difference_type = std::ptrdiff_t;
	using value_type = std::decay_t<T>;
	using pointer = T*;
	using reference = ref<T>;
	using iterator_category = std::random_access_iterator_tag;
	ptr() = default;//ptr(ptr const=default; ptr& operator=(ptr const&)=default;
	ptr(std::nullptr_t){}
	template<class Other> ptr(ptr<Other> const&){}
	reference operator*() const{return reference{&value};}
	ptr operator+(difference_type) const{return *this;}
	ptr& operator+=(difference_type){return *this;}
	ptr& operator++(){return operator+=(1);}
	friend difference_type operator-(ptr const&, ptr const&){return 0;}
	bool operator==(ptr const&) const{return true;}
	bool operator!=(ptr const&) const{return false;}
//	explicit operator T*() const{return &value;}
	ptr const& operator->() const{return *this;}
	friend ptr to_address(ptr const& p){return p;}
	explicit operator bool(){return false;}
//	operator double*() const{return &value;}
	friend std::allocator<value_type> get_allocator(ptr const&){return std::allocator<value_type>{};}
};
template<> double ptr<double>::value = 42.;
template<> double ptr<double const>::value = 42.;

template<class T> struct ref;

template<class T> struct ref{
protected:
	T* p_;
	ref(T* p) : p_{p}{}
	friend class ptr<T>;
	friend struct ref<T const>;
public:
	ref(ref<std::remove_const_t<T>> const& other) : p_{other.p_}{}
	bool operator==(ref const&) const{return true;}
	bool operator!=(ref const&) const{return false;}
	using decay_t = std::decay_t<T>;
};

template<class T> struct allocator{
	using pointer = ptr<T>;
	using value_type = T;
	auto allocate(std::size_t){return pointer{};}
	void deallocate(pointer, std::size_t){}
	std::true_type operator==(allocator const&){return {};}
	allocator(){}
	template<class T2> allocator(allocator<T2> const&){}
	template<class... Args>
	void construct(pointer, Args&&...){}
	void destroy(pointer){}
};

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
// multi-fancy glue, where should this be? 
// In boost/multi/adaptors/MyFancyApaptor.hpp if anything, or in user code if it is very specialized
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#include "../array.hpp"

namespace boost{
namespace multi{

template<class It, class T>  // custom copy 1D (aka strided copy)
fancy::ptr<T> copy(It first, It last, fancy::ptr<T> dest){ (void)last;
	return copy(first, last, multi::array_iterator<T, 1, fancy::ptr<T>>{dest});
//	std::cerr << "1D copy(it1D, it1D, it1D) with strides " << stride(first) << " " << stride(dest) << std::endl;
//	return dest;
}

template<class It, class T>  // custom copy 1D (aka strided copy)
auto copy(It first, It last, multi::array_iterator<T, 1, fancy::ptr<T>> dest){ (void)last;
	std::cerr << "1D copy(it1D, it1D, it1D) with strides " << stride(first) << " " << stride(dest) << std::endl;
	return dest;
}

template<class It, class T> // custom copy 2D (aka double strided copy)
auto copy(It first, It last, multi::array_iterator<T, 2, fancy::ptr<T>> dest){ (void)last; (void)first;
	std::cerr<<"2D copy(It, It, it2D) with strides 1"<< first.stride() <<" "<< dest.stride() <<std::endl;
	return dest;
}

//template<class Alloc, class It, class T> // custom copy 2D (aka double strided copy)
//auto uninitialized_copy(Alloc&, It first, It last, multi::array_iterator<T, 2, fancy::ptr<T>> const& dest){
//	std::cerr << "2D uninitialized_copy(...) calls raw copy 2D" << std::endl;
//	return copy(first, last, dest);
//}

}}

////////////////////////////////////////////////////////////////////////////////
// user code
////////////////////////////////////////////////////////////////////////////////

using std::cout; using std::cerr;
namespace multi = boost::multi;

int main(){
#if 0
	multi::array<double, 2, fancy::allocator<double>> A( 
//		multi::index_extensions<2>
		{5, 5}
	); // default element ctor
	assert( size(A) == 5 );
	assert( A[1][1] == A[2][2] );

	multi::array<double, 2, fancy::allocator<double>> B(A); // copy ctor
	assert( A[1][1] == B[1][1] );

	multi::array<double, 2, fancy::allocator<double>> C({5, 5}, 42.); // element value ctor

	multi::array<double, 2, fancy::allocator<double>> D; assert( D.num_elements() == 0 );

	D = C; assert( D == C);

	D = std::move(B); assert( D.size() == 5 and B.size() == 0 );

//	C = A; // calls custom copy
//	multi::array_ref<double, 2, fancy::ptr<double>> AA(A.data(), {5, 5});
//	C = A({0, 5}, {0,5});

	multi::array<double, 2, fancy::allocator<double> > CC;// = A({0, 5}, {0,5});
//	CC = A({0, 5}, {0,5});

//	using std::copy_n;
//	copy_n(A.data(), 25, A.data());
	static_assert( multi::array_iterator<double, 2, double*>::rank{} == 2, "!");
#endif
}


