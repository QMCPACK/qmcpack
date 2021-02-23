#ifdef compile_instructions
(echo "#include\""$0"\"" > $0x.cpp) && echo $0x.cpp && time clang++ -std=c++14 -Wfatal-errors -D_TEST_BOOST_JUST -lboost_system $0x.cpp -o $0x.cpp.x -lstdc++fs && ./$0x.cpp.x $@ && (xdg-open ${0%.*}.pdf 2>/dev/null &) && rm -f $0x.cpp $0x.cpp.x ; exit;
#endif

#ifndef BOOST_JUST_HPP
#define BOOST_JUST_HPP
//#include<boost/ref.hpp>
#include<utility> // std::forward
#include<array>

namespace boost{

#define BOOST_JUST_REPRODUCE_OPERATOR(OperatoR) \
	template<class T...> \
	auto OperatorR(T&&... t)

template<class T>
struct wrapper /*: T*/{
	typedef T type;
	T value;
	template<class... Args, typename = decltype(T(std::forward<Args>(std::declval<Args>()...)...))>
	wrapper(Args&&... args) : value(std::forward<Args>(args)...){}
	
	
/*	template<class Args>
	auto operator+=(Args&& args) 
	-> decltype(
		//	just<decltype(value+=(std::forward<Args>(args)))>(
				value+=(std::forward<Args>(args))
		//	)
	){
		return 
		//	just<decltype(value+=(std::forward<Args>(args)))>(
				value+=(std::forward<Args>(args))
		//	)
		;
	}*/
	template<class Args>
	auto operator=(Args&& args) -> decltype(value = std::forward<Args>(args)){
		return value = std::forward<Args>(args);
	}
	operator T&(){return value;}
	operator T const&() const{return value;}
	decltype(auto) operator&() const{return &value;}
//	operator T&() const{return *t_;}
//	operator T const&() const{return *t_;}
//	auto operator&() -> decltype(&value){return &value;}
//	auto operator&() const -> decltype(&value){return &value;}
	template<class Arg>
	auto operator[](Arg&& arg) const
	->decltype(value[std::forward<Arg>(arg)]){
		return value[std::forward<Arg>(arg)];
	}

	auto operator*()->decltype(*std::declval<T>()){return *value;} 
//	auto operator*() const->decltype(*value){return *value;}
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
	> /*no ::type here*/{
};



template<class T>
using just_t = typename just<T>::type;

#undef BOOST_JUST_REPRODUCE_OPERATOR

//template<class T>
//just<T> _(T const& t){
//	return just<T>(t);
//}

template<class T>
typename just<T>::type& _(T&& t){
	return reinterpret_cast<typename just<T>::type&>(std::forward<T>(t));
}

template<class T>
typename just<T>::type& wrap(T&& t){
	return reinterpret_cast<typename just<T>::type&>(std::forward<T>(t));
}


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

}


#endif

#define BOOST_INHERIT_UNARY_CONSTRUCTOR(MyclasS, MybaseclasS) \
	template<typename A> MyclasS(A& arg) : MybaseclasS(arg) {}

#ifdef _TEST_BOOST_JUST

//struct number : boost::just<double&>::type{
//	BOOST_INHERIT_UNARY_CONSTRUCTOR(number, boost::just<double&>::type)
//};

#include<iostream>
#include<vector>
#include<cassert>

template<class T>
class A : boost::just<T>::type{};

int main(){
	A<int> a;

	A<int[8]> b;
	assert( std::is_class<int[8]>::value == false );
	assert( std::is_array<int[8]>::value == true );

	{
	
		double d=5;
		boost::wrapper<double> n(d);
		n+=4;
		std::cout << n << std::endl;
		std::cout << d << std::endl;
	}
	{
		double d = 5.;
		boost::wrapper<double&> n(d);
		double a = 6.;
		std::cout << n << std::endl;
		n = a;
		n+= 5.;
		assert(&n == &d);
		std::cout << n << std::endl;
		std::cout << d << std::endl;
	}

	{
		double d = 5.;
		std::vector<double> v(10);
		std::vector<boost::reference<double>> v3;
		v3.push_back(d);
		v3.push_back(d);
		v3[0] = 4.;
		std::cout << v3[0] << " " << v3[1] << std::endl;
	}

	{
		double a = 5.;
		double& b = a;
		assert( &b == &a );
	//	boost::just<double&> c = a;
	//	c = 6.;		
	}

}
#endif

