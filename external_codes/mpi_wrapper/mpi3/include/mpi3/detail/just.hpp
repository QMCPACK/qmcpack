// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;autowrap:nil;-*-
// Copyright 2018-2022 Alfredo A. Correa

#ifndef BOOST_JUST_HPP
#define BOOST_JUST_HPP

#include<array>
#include<utility>  // std::forward

namespace boost {

// NOLINTNEXTLINE(cppcoreguidelines-macro-usage) :
#define BOOST_JUST_REPRODUCE_OPERATOR(OperatoR) \
	template<class T...> \
	auto OperatorR(T&&... t)

template<class T>
struct wrapper /*: T*/ {
	using type = T;
	T value;  // NOLINT(misc-non-private-member-variables-in-classes) : wrapper

	template<class... Args, typename = decltype(T(std::forward<Args>(std::declval<Args>()...)...))>
	explicit wrapper(Args&&... args) : value(std::forward<Args>(args)...) {}

	template<class Args, class = decltype(std::declval<T&>() = std::declval<Args&&>())>
	wrapper& operator=(Args&& args) {
		value = std::forward<Args>(args);
		return *this;
	}

	explicit operator T      &()       {return value;}
	explicit operator T const&() const {return value;}

	template<class Arg>
	auto operator[](Arg&& arg) const
	->decltype(value[std::forward<Arg>(arg)]) {
		return value[std::forward<Arg>(arg)]; }

	auto operator*()->decltype(*std::declval<T>()) {return *value;} 
};

template<class T>
using reference = wrapper<T&>;

template<class T>
struct just :
	std::conditional<
		std::is_class<T>::value,
			T,
			typename std::conditional<std::is_array<T>::value,
				std::array<typename std::remove_extent<T>::type, std::extent<T>::value>,
				wrapper<T>
			>::type
	> /*no ::type here*/ {
};

template<class T>
using just_t = typename just<T>::type;

#undef BOOST_JUST_REPRODUCE_OPERATOR

//template<class T>
//typename just<T>::type& _(T&& t){
//	return reinterpret_cast<typename just<T>::type&>(std::forward<T>(t));
//}

//template<class T>
//typename just<T>::type& wrap(T&& t){
//	return reinterpret_cast<typename just<T>::type&>(std::forward<T>(t));
//}


/*
template<class T>
struct just<T&>{// : boost::reference_wrapper<T>{
	T* t_ptr;
	just(T& t) : t_ptr(&t){}
	typedef just<T&> type;
	just<T&>& operator=(T const& t){*t_ptr = t;}
//	typedef std::reference_wrapper<T> type;
};*/

/*
template<>
struct just<double>{
	double impl_;
	typedef just<double> type;
	just() : impl_(){}
	just(double const& d) : impl_(d){}
	operator double const&() const{return impl_;}
	double& operator+=(double const& d){return impl_+=d;}
};
template<>
struct just<bool>{
	bool impl_;
	typedef just<bool> type;
	just(bool const& d) : impl_(d){}
	operator bool const&() const{return impl_;}
//	double& operator+=(double const& d){return impl_+=d;}
};
*/

}  // end namespace boost

#endif

// NOLINTNEXTLINE(cppcoreguidelines-macro-usage) : TODO(correaa) : check if this is necessary
#define BOOST_INHERIT_UNARY_CONSTRUCTOR(MyclasS, MybaseclasS) \
	template<typename A> MyclasS(A& arg) : MybaseclasS(arg) {}

// #ifdef _TEST_BOOST_JUST

// //struct number : boost::just<double&>::type{
// //  BOOST_INHERIT_UNARY_CONSTRUCTOR(number, boost::just<double&>::type)
// //};

// #include<iostream>
// #include<vector>
// #include<cassert>

// template<class T>
// class A : boost::just<T>::type {};

// int main() {
// 	A<int> a;

// 	A<int[8]> b;
// 	assert( std::is_class<int[8]>::value == false );
// 	assert( std::is_array<int[8]>::value == true );

// 	{
// 		double d=5;
// 		boost::wrapper<double> n(d);
// 		n+=4;
// 		std::cout<< n <<std::endl;
// 		std::cout<< d <<std::endl;
// 	}
// 	{
// 		double d = 5.;
// 		boost::wrapper<double&> n(d);
// 		double aa = 6.;
// 		std::cout<< n <<std::endl;
// 		n = aa;
// 		n+= 5.;
// 		assert(&n == &d);
// 		std::cout<< n <<std::endl;
// 		std::cout<< d <<std::endl;
// 	}

// 	{
// 		double d = 5.;
// 		std::vector<boost::reference<double>> v3;
// 		v3.push_back(d);
// 		v3.push_back(d);
// 		v3[0] = 4.;
// 		std::cout << v3[0] << " " << v3[1] << std::endl;
// 	}
// }
// #endif
