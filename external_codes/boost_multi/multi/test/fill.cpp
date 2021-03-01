#ifdef COMPILATION// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;-*-
$CXX $CXXFLAGS $0 -lm -o $0x -lboost_unit_test_framework&&$0x&&rm $0x;exit
#endif

#define BOOST_TEST_MODULE "C++ Unit Tests for Multi fill"
#define BOOST_TEST_DYN_LINK
#include<boost/test/unit_test.hpp>

#include<random>

#include "../array.hpp"

#include<cstddef>
#include<limits>
#include<type_traits> // enable_if_t
#include<iostream>  //cout

//#include <range/v3/all.hpp>

#include <boost/iterator/transform_iterator.hpp>

// from Howard Hinnart hash
std::size_t fnv1a(void const* key, std::size_t len, std::size_t h = 14695981039346656037u) noexcept{
	unsigned char const* p = static_cast<unsigned char const*>(key);
	unsigned char const* const e = p + len;
	for(; p < e; ++p) h = (h ^ *p) * 1099511628211u; // prime
	return h;
}

class fnv1a_t{
	std::size_t h = 14695981039346656037u; // offset
public:
	using result_type = std::size_t;
	static constexpr result_type min(){return std::numeric_limits<std::size_t>::min();}
	static constexpr result_type max(){return std::numeric_limits<std::size_t>::max();}
	void operator()(void const* key, std::size_t len) noexcept{h = fnv1a(key, len, h);}
	template<class T, std::enable_if_t<std::is_fundamental<T>{}, int> = 0>
	fnv1a_t& operator()(T const& t) noexcept{operator()(&t, sizeof(t)); return *this;}
	result_type operator()() && noexcept{return h;}
	result_type operator()() & noexcept{return h;}
	explicit operator result_type() && noexcept{return h;}
	explicit operator result_type() & noexcept{return h;}
};

std::random_device r;

BOOST_AUTO_TEST_CASE(fill_1d){
	namespace multi = boost::multi;
	{
		multi::array<double, 1> d1D(10);
		static_assert( std::is_same<std::iterator_traits<decltype(begin(d1D))>::value_type, double>{}, "!");

		using std::copy;
		copy(begin(extension(d1D)), end(extension(d1D)), begin(d1D));
		BOOST_REQUIRE( d1D[0] == 0. );
		BOOST_REQUIRE( d1D[1] == 1. );
		BOOST_REQUIRE( d1D[9] == 9. );
		
//		auto f = [](auto e){return e + 1;};
//		copy(
//			boost::make_transform_iterator(begin(extension(d1D)), f), 
//			boost::make_transform_iterator(end  (extension(d1D)), f), 
//			begin(d1D)
//		);
//		BOOST_REQUIRE( d1D[0] == 1. );
//		BOOST_REQUIRE( d1D[1] == 2. );
//		BOOST_REQUIRE( d1D[9] == 10. );

//		ranges::copy(extension(d1D), begin(d1D));
//		BOOST_REQUIRE( d1D[0] == 0. );
//		BOOST_REQUIRE( d1D[1] == 1. );
//		BOOST_REQUIRE( d1D[9] == 9. );

//		ranges::copy(extension(d1D) | ranges::views::transform([](auto e){return e+1;}), begin(d1D));
//		BOOST_REQUIRE( d1D[0] == 1. );
//		BOOST_REQUIRE( d1D[1] == 2. );
//		BOOST_REQUIRE( d1D[9] == 10. );

		d1D.assign(extension(d1D));
		BOOST_REQUIRE( d1D[0] == 0. );
		BOOST_REQUIRE( d1D[1] == 1. );
		BOOST_REQUIRE( d1D[9] == 9. );

//		d1D.assign(extension(d1D)|ranges::views::transform([](auto e){return e+1;}));
//		BOOST_REQUIRE( d1D[0] == 1. );
//		BOOST_REQUIRE( d1D[1] == 2. );
//		BOOST_REQUIRE( d1D[9] == 10. );

/*		fnv1a_t h;
		using urd = std::uniform_real_distribution<>;
		auto crypt = [h=h(111)](auto i){return urd{}(fnv1a_t{h}(i));};
		copy(
			boost::make_transform_iterator(begin(extension(d1D)), crypt), 
			boost::make_transform_iterator(end  (extension(d1D)), crypt), 
			begin(d1D)
		);
		for(auto&& e : d1D) std::cout<< e <<std::endl;
		std::cout<<"***"<<std::endl;

		using urd = std::uniform_real_distribution<>;
		auto crypt2 = [h=h(222)](auto i){return fnv1a_t{h}(i)();};
*/
//		copy(
//			boost::make_transform_iterator(begin(extension(d1D)), crypt2), 
//			boost::make_transform_iterator(end  (extension(d1D)), crypt2), 
//			begin(d1D)
//		);
//		for(auto&& e : d1D) std::cout<< e <<std::endl;
//		std::cout<<"///"<<std::endl;
	}
	{
//		auto const& xp1 = multi::index_extension(10) | ranges::views::transform([](auto e){return e+1;});
//		multi::array<double, 1> d1D(xp1.begin(), xp1.end());
//		BOOST_REQUIRE( d1D[0] == 1. );
//		BOOST_REQUIRE( d1D[1] == 2. );
//		BOOST_REQUIRE( d1D[9] == 10. );
	}
	{
//		multi::array<double, 1> d1D{multi::index_extension(10) | ranges::views::transform([](auto e){return e+1;})};
//		BOOST_REQUIRE( d1D[0] == 1. );
//		BOOST_REQUIRE( d1D[1] == 2. );
//		BOOST_REQUIRE( d1D[9] == 10. );
	}
//	std::mt19937 e;
//	fnv1a_t h;
//	using urd = std::uniform_real_distribution<>;
	{
//		multi::array<double, 1> d1D{multi::index_extension(10) | ranges::views::transform([h=h(e())](auto i){return urd{-10, 10}(fnv1a_t{h}(i));})};
//
//		BOOST_REQUIRE( d1D[1] >= -10. and d1D[1] <   10. ); 
//		for(auto i : extension(d1D)) std::cout << d1D[i] << std::endl;
	}
	std::cout << "====\n";
	{
//		multi::array<double, 1> d1D(10);
//		d1D.assign(extension(d1D)|ranges::views::transform([h=h(e())](auto i){return urd{-10, 10}(fnv1a_t{h}(i));}));
//
//		BOOST_REQUIRE( d1D[1] >= -10. and d1D[1] < 10. );
//		for(auto&& e : d1D) std::cout<< e <<std::endl;
	}
	{
		multi::array<double, 1> d1D(begin(multi::index_extension(10)), end(multi::index_extension(10)));
		BOOST_REQUIRE( size(d1D) == 10 );
		BOOST_REQUIRE( d1D[0] == 0. );
		BOOST_REQUIRE( d1D[1] == 1. );
		BOOST_REQUIRE( d1D[9] == 9. );
	}
	{
		multi::array<double, 1> d1D(10);
		BOOST_REQUIRE( size(d1D) == 10 );

		d1D.assign(begin(extension(d1D)), end(extension(d1D)));
		BOOST_REQUIRE( d1D[0] == 0. );
		BOOST_REQUIRE( d1D[1] == 1. );
		BOOST_REQUIRE( d1D[9] == 9. );		
	}
	{
		multi::array<double, 1> d1D(10);
		d1D.assign(extension(d1D));
		BOOST_REQUIRE( d1D[0] == 0. );
		BOOST_REQUIRE( d1D[1] == 1. );
		BOOST_REQUIRE( d1D[9] == 9. );
	}
}

BOOST_AUTO_TEST_CASE(fill_member){
	namespace multi = boost::multi;
	multi::array<double, 1> d1D = {1., 2., 3., 4.};
	d1D.fill(42.);

	multi::array<double, 2> d2D = {
		{150., 16., 17., 18., 19.}, 
		{  5.,  5.,  5.,  5.,  5.}, 
		{100., 11., 12., 13., 14.}, 
		{ 50.,  6.,  7.,  8.,  9.}  
	};
	BOOST_REQUIRE( d2D.elements().size() == d2D.num_elements() );
	BOOST_REQUIRE( d2D.elements().base() == d2D.base() );
	BOOST_REQUIRE( d2D.elements()[3] == 18. );
	BOOST_REQUIRE( &*d2D.elements().begin() == d2D.data_elements() );
	BOOST_REQUIRE( &*d2D.elements().end() == d2D.data_elements() + d2D.num_elements() );
//	std::fill( d2D.elements().begin(), d2D.elements().end() , 99. );
//	multi::adl_fill_n( d2D.elements().begin(), d2D.elements().size(), 99. );
	d2D.elements().fill(99.);

	BOOST_REQUIRE( d2D[1][1] == 99. );
}

BOOST_AUTO_TEST_CASE(fill){
	namespace multi = boost::multi;

	multi::array<double, 2> d2D = {
		{150., 16., 17., 18., 19.}, 
		{  5.,  5.,  5.,  5.,  5.}, 
		{100., 11., 12., 13., 14.}, 
		{ 50.,  6.,  7.,  8.,  9.}  
	};
	using std::all_of;
	BOOST_REQUIRE( all_of(begin(d2D[1]), end(d2D[1]),[](auto& e){return e==5.;}) );

	using std::fill;
	fill(d2D[1].begin(), d2D[1].end(), 8.);

	BOOST_REQUIRE( all_of(begin(d2D[1]), end(d2D[1]), [](auto& e){return e==8.;}) );

	fill(begin(rotated(d2D)[1]), end(rotated(d2D)[1]), 8.);
	BOOST_REQUIRE( all_of(begin(rotated(d2D)[1]), end(rotated(d2D)[1]), [](auto&& e){return e==8.;}) );

	fill(begin((d2D<<1)[1]), end((d2D<<1)[1]), 8.);
	BOOST_REQUIRE( all_of(begin((d2D<<1)[1]), end((d2D<<1)[1]), [](auto&& e){return e==8.;}) );

	std::mt19937 g{r()};
	auto rand = [d=std::normal_distribution<>{}, g = std::mt19937{r()}]() mutable{return d(g);};
	multi::array<double, 2> r2D({5, 5});
	std::for_each(begin(r2D), end(r2D), [&](auto&& e){std::generate(begin(std::move(e)), end(std::move(e)), rand);});

}

