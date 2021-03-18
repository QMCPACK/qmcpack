#ifdef COMPILATION// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;-*-
$CXX $CXXFLAGS $0 -o $0.$X -lboost_unit_test_framework&&$0.$X $@&&rm $0.$X;exit
#endif
// © Alfredo A. Correa 2019-2020

#define BOOST_TEST_MODULE "C++ Unit Tests for Multi array reference"
#define BOOST_TEST_DYN_LINK
#include<boost/test/unit_test.hpp>

#include "../array.hpp"
#include "../array_ref.hpp"

#include<boost/multi_array.hpp>

namespace multi = boost::multi;

BOOST_AUTO_TEST_CASE(array_ref_from_carray){

	double a[4][5] {
		{ 0,  1,  2,  3,  4}, 
		{ 5,  6,  7,  8,  9}, 
		{10, 11, 12, 13, 14}, 
		{15, 16, 17, 18, 19}
	};

	multi::array_ptr<double, 2> map = &a; // double (*)[4][5]
	BOOST_REQUIRE( &map->operator[](1)[1] == &a[1][1] );
	BOOST_REQUIRE( (&a)[0][1][1] == 6. );
	
	multi::array_ref<double, 2>&& mar = *map;

	BOOST_REQUIRE( &mar[1][1] == &a[1][1] );

	mar[1][1] = 9.;
	BOOST_REQUIRE( &mar[1][1] == &a[1][1] );

	double const(&a_const)[4][5] = a;
	BOOST_REQUIRE( &a_const[1][1] == &a[1][1] );

	BOOST_REQUIRE( mar(2, {1, 3}).dimensionality == 1 );
	BOOST_REQUIRE( dimensionality(mar(2, {1, 3})) == 1 );

	BOOST_REQUIRE( size(mar(2, {1, 3})) == 2 );
	BOOST_REQUIRE( &mar(2, {1, 3})[1] == &a[2][2] );

}

BOOST_AUTO_TEST_CASE(array_ref_1D_reindexed){
	std::string (&&a)[5] = {"a", "b", "c", "d", "e"};

	multi::Array<std::string(&)[1]> mar = *multi::Array<std::string(*)[1]>(&a);
	BOOST_REQUIRE( &mar[1] == &a[1] );
	BOOST_REQUIRE( sizes(mar.reindexed(1)) == sizes(mar) );
	auto diff = &(mar.reindexed(1)[1]) - &mar[0];
	BOOST_REQUIRE( diff == 0 );
	
	BOOST_REQUIRE( &mar.blocked(2, 4)[2] == &mar[2] );
	for(auto i : extension(mar.stenciled({2, 4})))
		BOOST_REQUIRE( &mar.stenciled({2, 4})[i] == &mar[i] );
		
		
	multi::array<std::string, 1> arr({{2, 7}}, std::string{"xx"});
	BOOST_REQUIRE( size(arr) == 5 );
	BOOST_REQUIRE( extension(arr) == multi::iextension(2, 7) );
	arr[2] = "a";
	arr[3] = "b";
	arr[4] = "c";
	arr[5] = "d";
	arr[6] = "e";
	BOOST_REQUIRE( std::equal(arr.begin(), arr.end(), mar.begin()) );
	
	auto arrB = multi::array<std::string, 1>({"a", "b", "c", "d", "e"}).reindex(2);
	BOOST_REQUIRE( size(arrB) == 5 );
	BOOST_REQUIRE( arrB[2] == "a" );
	BOOST_REQUIRE( arrB[6] == "e" );
}

BOOST_AUTO_TEST_CASE(array_ref_reindexed){

	double (&&a)[4][5] = {
		{ 0,  1,  2,  3,  4}, 
		{ 5,  6,  7,  8,  9}, 
		{10, 11, 12, 13, 14}, 
		{15, 16, 17, 18, 19}
	};

	multi::Array<double(&)[2]> mar = *multi::Array<double(*)[2]>(&a);

	BOOST_REQUIRE( &mar[1][1] == &a[1][1] );
	BOOST_REQUIRE( size(mar   .reindexed(1)) == size(mar) );
	BOOST_REQUIRE( size(mar[0].reindexed(1)) == size(mar[0]) );

	BOOST_REQUIRE( sizes(mar.reindexed(1)) == sizes(mar) );
	BOOST_REQUIRE( &mar.reindexed(1)[1][0] == &mar[0][0] );

	BOOST_REQUIRE( sizes(mar[0].reindexed(1)) == sizes(mar[0]) );
	BOOST_REQUIRE( mar[0].reindexed(1).extension().start () == mar[0].extension().start () + 1 );
	BOOST_REQUIRE( mar[0].reindexed(1).extension().finish() == mar[0].extension().finish() + 1 );

	auto diff = &mar[0].reindexed(1)[1] - &mar[0][0];
	BOOST_REQUIRE( diff == 0 );
	
//	BOOST_REQUIRE( &(((mar<<1).reindexed(2)>>1).reindexed(1))[1][2] == &mar[0][0] );
	BOOST_REQUIRE( &mar.reindexed(1, 2)[1][2] == &mar[0][0] );
	
	BOOST_REQUIRE( &mar.reindexed(1)({1, 5})[1][0] == &mar[0][0] );
	
	BOOST_REQUIRE( sizes(mar.stenciled({2, 4})) == std::make_tuple(2, 5) );
	BOOST_REQUIRE( &mar.stenciled({2, 4})[2][0] == &mar[2][0] );
	BOOST_REQUIRE( &mar.stenciled({2, 4}, {1, 3})[2][1] == &mar[2][1] );

//	BOOST_REQUIRE( &mar[0][0] == mar.origin() ); // origin changed meaning in on 2020/Dec/16
//	BOOST_REQUIRE( mar.base() == mar.origin() );

//	BOOST_REQUIRE( mar.stenciled({2, 4}).origin() == mar.origin() ); // origin changed meaning in on 2020/Dec/16
	BOOST_REQUIRE( mar.stenciled({2, 4}).base()   != mar.base()   );

	BOOST_REQUIRE( &mar.stenciled({2, 4})[2][0] == mar.stenciled({2, 4}).base() );

	{
		multi::array<std::string, 2> arrB = 
			{{"a", "b", "c", "d", "e"},
			 {"f", "g", "h", "f", "g"},
			 {"h", "i", "j", "k", "l"}}
		;
		arrB.reindex(2);
		BOOST_REQUIRE( size(arrB) == 3 );
		BOOST_REQUIRE( arrB[2][0] == "a" );
	}
	{
		multi::array<std::string, 2> arrB = 
			{{"a", "b", "c", "d", "e"},
			 {"f", "g", "h", "f", "g"},
			 {"h", "i", "j", "k", "l"}}
		;
		arrB.reindex(2, 1);
		BOOST_REQUIRE( size(arrB) == 3 );
		BOOST_REQUIRE( arrB[2][1] == "a" );
	}
	{
		multi::array<std::string, 2> arrB = (multi::array<std::string, 2>
			{{"a", "b", "c", "d", "e"},
			 {"f", "g", "h", "f", "g"},
			 {"h", "i", "j", "k", "l"}})//.reindex(2, 1);
		;
		BOOST_REQUIRE( arrB.reindex(2).extension() == multi::iextension(2, 5) );
		auto x = arrB.reindexed(2).extensions();

		multi::array<std::string, 2> arrC(x);
		BOOST_REQUIRE( size(arrC) == 3 );
		BOOST_REQUIRE( size(arrC) == size(arrB) );

		BOOST_REQUIRE( arrC.extension().start()  == 2 );
		BOOST_REQUIRE( arrC.extension().finish() == 5 );
		
		std::cout<< arrC.extension().start()<<"..."<< arrC.extension().finish() <<std::endl;
	//	BOOST_REQUIRE( arrC[2][1] == "a" );
	}	
	

}

BOOST_AUTO_TEST_CASE(array_ref_with_stencil){

	double a[4][5] = {
		{ 0,  1,  2,  3,  4}, 
		{ 5,  6,  7,  8,  9}, 
		{10, 11, 12, 13, 14}, 
		{15, 16, 17, 18, 19}
	};
	auto const& mar = *multi::Array<double(*)[2]>(&a);

	multi::Array<double[2]> s = {
		{ 0, +1,  0},
		{+1, -4, +1},
		{ 0, +1,  0}
	};
	auto const& st = s.reindexed(-1, -1);
	
//	BOOST_REQUIRE( &st[ 0][ 0] == st.origin() ); // origin changed meaning in on 2020/Dec/16
	BOOST_REQUIRE( &st[-1][-1] == st.base() );
	
	multi::array<double, 2> gy(extensions(mar), 0.);

	{
		auto xs = extensions(mar);
		for(auto i = std::get<0>(xs).start()+1; i != std::get<0>(xs).finish()-1; ++i)
			for(auto j = std::get<1>(xs).start()+1; j != std::get<1>(xs).finish()-1; ++j){
				auto xt = extensions(st);
				for(auto b: std::get<0>(xt))
					for(auto d: std::get<1>(xt))
						gy[i][j] += st[b][d]*a[i + b][j + d];
			}
	}

}

BOOST_AUTO_TEST_CASE(array_ref_1D_from_vector){
	std::vector<double> v = {1, 2, 3};
	multi::array_ref<double, 1> aref({{1, 3}}, v.data());
	BOOST_REQUIRE( aref.extension() == multi::iextension(1, 3) );
	BOOST_REQUIRE( &aref[1] == &v[0] );
}

BOOST_AUTO_TEST_CASE(array_ref_2D_from_vector){
	std::vector<double> v = {1, 2, 3, 4, 5, 6};
	multi::array_ref<double, 2> aref({2, 3}, v.data());
	BOOST_REQUIRE( &aref[1][0] == &v[3] );
}

BOOST_AUTO_TEST_CASE(array_ref_2D_from_vector_with_offset){
	std::vector<double> v = {
		1, 2, 3, 
		4, 5, 6
	};
	multi::array_ref<double, 2> aref({multi::iextension(1, 3), multi::iextension(1, 4)}, v.data());

	auto x = aref.extensions();
	BOOST_REQUIRE( std::get<0>(x) == multi::iextension(1, 3) );
	BOOST_REQUIRE( std::get<1>(x).start()  == 1 );
	BOOST_REQUIRE( std::get<1>(x).finish() == 4 );
	BOOST_REQUIRE( std::get<1>(x) == multi::iextension(1, 4) );
	BOOST_REQUIRE( x == decltype(x)(multi::iextension(1, 3), multi::iextension(1, 4)) );

	BOOST_REQUIRE( &aref[1][1] == &v[0] );
}

BOOST_AUTO_TEST_CASE(array_2D_with_offset){
	multi::array<double, 2> a({multi::iextension(1, 3), multi::iextension(2, 5)}, 1.2);

	BOOST_REQUIRE( a.extension().start() == 1 );
	BOOST_REQUIRE( a.extension().finish() == 3 );
}

BOOST_AUTO_TEST_CASE(array_ref_1D){

	std::array<std::string, 5> a = {"a", "b", "c", "d", "e"};

	multi::array_ref<std::string, 1>&& mar = *multi::array_ptr<std::string, 1>{&a};
//	multi::Array<std::string(&)[1]> mar = *multi::Array<std::string(*)[1]>(&a);

	BOOST_REQUIRE(  extension(mar).first() == 0 );
	BOOST_REQUIRE(  extension(mar).last()  == 5 );

	auto&& mar1 = mar.reindexed(1);

	BOOST_REQUIRE( extension(mar1).size() == extension(mar).size() );

	BOOST_REQUIRE( mar1.extension() == extension(mar1) );
	BOOST_REQUIRE(  extension(mar1).first() == 1 );
	BOOST_REQUIRE(  mar1.extension().first() == 1 );
	BOOST_REQUIRE(  mar1.extension().last()  == 6 );
	BOOST_REQUIRE( *extension(mar1).begin() == 1 );

	BOOST_REQUIRE( size(mar1) == size(mar) );
	BOOST_REQUIRE( mar1.layout().extension().start() == 1 );
//	BOOST_REQUIRE( extension(mar1).start() == 1 );
	BOOST_REQUIRE( &mar1[1] == &a[0] );
	BOOST_REQUIRE( mar1.base() == &a[0] );

}

//BOOST_AUTO_TEST_CASE(array_from_iterator_range){
//	std::vector<double> v(20);
//	multi::array_ref<double, 2, std::vector<double>::iterator> R(begin(v), {4, 5});
//	assert( P.base() == biit );
//	auto at = [](auto&& pp, std::ptrdiff_t idx) -> multi::array_ref<double, 2, decltype(biit)>::reference{
//		return {pp.layout().sub, pp.base() + pp.layout().operator()(idx)};
//	};
//	at(P, 0);
//	P[0];
//	P[0][1];
//		P[0];
//		P[0];
//	for(int i = 0; i != 2; ++i){
//		for(int j = 0; j != 2; ++j)
//			std::cout << P[i][j] << ',';
//		std::cout << std::endl;
//	}
//}

#if 0

BOOST_AUTO_TEST_CASE(iterator_1d){
	{

		multi::biiterator<std::decay_t<decltype(std::move(A).begin())>> biit{std::move(A).begin(), 0, size(*std::move(A).begin())};
		for(int i = 0; i!=num_elements(A); ++i, ++biit)
			cout << i << "->" << *biit << ", ";
		cout << std::endl;
	
//		auto biit2 = biit + 2;
		multi::array_ref<double, 2, decltype(biit)> P(biit, {2, 10});

	//	cout << std::endl;	
	//	P[0][0];
		for(int i = 0; i != 2; ++i){
			for(int j = 0; j != 10; ++j)
				std::cout << "j = "<< j << "->" << P[i][j] << " ,";
			std::cout << std::endl;
		}
#if 0
		std::vector<double> v(20);
		multi::array_ref<double, 2, std::vector<double>::iterator> R(begin(v), extensions(A));
		assert( P.base() == biit );
	//	auto at = [](auto&& pp, std::ptrdiff_t idx) -> multi::array_ref<double, 2, decltype(biit)>::reference{
	//		return {pp.layout().sub, pp.base() + pp.layout().operator()(idx)};
	//	};
	//	at(P, 0);
		P[0];
	//	P[0][1];
//		P[0];
//		P[0];
	//	for(int i = 0; i != 2; ++i){
	//		for(int j = 0; j != 2; ++j)
	//			std::cout << P[i][j] << ',';
	//		std::cout << std::endl;
	//	}
#endif
	}
//	return 0;
	double const d2D[4][5] {
		{ 0,  1,  2,  3,  4}, 
		{ 5,  6,  7,  8,  9}, 
		{10, 11, 12, 13, 14}, 
		{15, 16, 17, 18, 19}
	};
	multi::array_ref<double, 2, double const*> d2D_cref(&d2D[0][0], {4, 5});

	decltype(d2D_cref)::value_type a_row = d2D_cref[2];
	decltype(d2D_cref[2])::decay_type b_row = d2D_cref[2];
	auto c_row = d2D_cref[2].decay();
	auto d_row = decay(d2D_cref[2]);
	auto const& e_row = d2D_cref[2];
	assert( e_row[3] == d2D_cref[2][3] );
	std::vector<double> d2Dv(20); iota(d2Dv.begin(), d2Dv.end(), 0);
	multi::array_ref<double, 2, double const*> d2Dv_cref(d2Dv.data(), {4, 5});
	assert( d2Dv.size() == std::size_t(num_elements(d2Dv_cref)) );
	assert( d2Dv_cref[1][1] == 6);
	assert( std::move(d2D_cref).data_elements() == data_elements(std::move(d2D_cref)) );
	assert( data_elements(std::move(d2D_cref)) == &d2D[0][0] );
	assert( d2D_cref.num_elements() == num_elements(d2D_cref) );
	assert( num_elements(d2D_cref) == 4*5 );
	assert( d2D_cref.size() == size(d2D_cref) );
	assert( d2D_cref.size() == 4 );
	assert( std::get<0>(sizes(d2D_cref)) == size(d2D_cref) );
	assert( std::get<0>(sizes(d2D_cref)) == 4 );
	assert( std::get<1>(sizes(d2D_cref)) == 5 );
	assert( d2D_cref.shape() == sizes(d2D_cref) );


	assert( stride(d2D_cref) == 5 );
	assert( std::get<1>(strides(d2D_cref)) == 1 );
	assert( strides(d2D_cref) == d2D_cref.strides() );	
	assert( std::get<0>(extensions(d2D_cref)) == extension(d2D_cref) );
	assert( &data_elements(std::move(d2D_cref))[2] == &d2D_cref[0][2] );
	assert( &data_elements(std::move(d2D_cref))[6] == &d2D_cref[1][1] );

	assert( std::move(d2D_cref).begin() == begin(std::move(d2D_cref)) );
	assert( std::move(d2D_cref).end() == end(std::move(d2D_cref)) );
	assert( begin(std::move(d2D_cref)) != end(std::move(d2D_cref)) );
	assert( begin(std::move(d2D_cref)) + size(d2D_cref) == end(std::move(d2D_cref)) );
//	assert( end(d2D_cref) - begin(d2D_cref) == size(d2D_cref) );
	using std::distance;
	assert( distance(begin(std::move(d2D_cref)), end(std::move(d2D_cref))) == size(d2D_cref));
	assert( size(*begin(std::move(d2D_cref))) == 5 );
//	assert( distance((std::move(d2D_cref).begin())->begin(), ((std::move(d2D_cref)).begin())->end()) == ((std::move(d2D_cref)).begin())->size() );
	assert( distance(begin(*begin(std::move(std::move(d2D_cref)))), end(*begin(std::move(d2D_cref)))) == size(*begin(std::move(d2D_cref))) );

	assert( size(std::move(d2D_cref)[0]) == 5 );
	assert( num_elements(std::move(d2D_cref)[0]) == 5 );

	for(auto i=0; i != d2D_cref.size(0) ||!endl(cout); ++i)
		for(auto j=0; j != d2D_cref.size(1) ||!endl(cout); ++j)
			cout<< d2D_cref[i][j] <<' ';
			
	{
		auto ext = extensions(d2D_cref);
		using std::get;
		for(auto i : get<0>(ext)){
			for(auto j : get<1>(ext)) cout<< d2D_cref[i][j] <<' ';
			cout<<'\n';
		}
		cout<<'\n';
	}

	auto const& d2D_cref_sliced = d2D_cref.sliced(1, 3);
	assert( size(d2D_cref_sliced) == 2 );
	auto const& d2D_cref_sliced_strided = d2D_cref_sliced.strided(2);
	assert( size(d2D_cref_sliced_strided) == 1 );

	for(auto i = 0; i != std::get<0>(sizes(d2D_cref_sliced)) ||!endl(cout); ++i)
		for(auto j = 0; j != std::get<1>(sizes(d2D_cref_sliced)) ||!endl(cout); ++j)
			cout<< d2D_cref_sliced[i][j] <<' ';
			
	assert(std::move(d2D_cref).begin() == std::move(d2D_cref).begin(0));
	assert(std::move(d2D_cref).begin() != std::move(d2D_cref).begin(1));
	for(auto it1 = std::move(d2D_cref).begin(1); it1 != std::move(d2D_cref).end(1)||!endl(cout); ++it1)
		for(auto it2 = std::move(*it1).begin()   ; it2 != std::move(*it1).end()   ||!endl(cout); ++it2)
			cout << *it2 << ' ';

	multi::array_ref<double, 2, double const*> d2D_crefref(std::move(d2D_cref).data_elements(), extensions(d2D_cref));
	multi::array_cref<double, 2> d2D_crefcref(data_elements(std::move(d2D_cref)), extensions(d2D_cref));

	using std::for_each;
	using std::begin;
	using std::end;
	for_each(
		begin(std::move(d2D_cref)), end(std::move(d2D_cref)), 
		[](auto&& row){
			for_each(begin(std::move(row)), end(std::move(row)), [](auto&& e){cout<<' '<< e;})("\n");
		}
	)("\n");
	
//	for(auto&& row: std::move(d2D_cref)){for(auto&& el: row) cout<< el <<' '; cout<<'\n';}
	for(auto itrow = std::move(d2D_cref).begin(); itrow != std::move(d2D_cref).end(); ++itrow){
		for(auto itel = (*itrow).begin(); itel != (*itrow).end(); ++itel) cout<< *itel <<' ';
		cout<<'\n';
	}
//	for(decltype(d2D_cref)::reference r : d2D_cref){
//		for(decltype(d2D_cref)::element const& e: r) cout << e <<' '; 
//		cout <<'\n';
//	}
	for(auto it=begin(std::move(d2D_cref)); it!=end(std::move(d2D_cref)); ++it){
//		for(decltype(d2D_cref)::element const& e: *it) cout << e <<' '; 
		for(auto it2=begin(*it); it2!=end(*it); ++it2) cout<< *it2 <<' ';
		cout<<'\n';
	}
	
	for(auto it1=begin(std::move(d2D_cref)); it1 != end(std::move(d2D_cref)) ||!endl(cout); ++it1)
		for(auto it2=(*it1).begin(); it2!=(*it1).end() ||!endl(cout); ++it2)
			cout<< *it2 <<' ';

//	for(auto it1=(&std::move(d2D_cref))->begin(); it1 != (&std::move(d2D_cref))->end() ||!endl(cout); ++it1)
//		for(auto it2=it1->begin()   ; it2 != it1->end()    ||!endl(cout); ++it2)
//			cout<< *it2 <<' ';

	auto print = [](auto&& arr){
		for(auto it1 = std::forward<decltype(arr)>(arr).begin(); it1 != std::forward<decltype(arr)>(arr).end()||!endl(cout); ++it1)
			for(auto it2 = (*it1).begin()   ; it2 != (*it1).end()   ||!endl(cout); ++it2)
				cout << *it2 << ' ';
	};
	print(std::move(d2D_cref));
	print(std::move(d2D_cref).range({0, 2}));
	print(std::move(d2D_cref).rotated(1).range({0, 2}).rotated(1));


	using std::is_sorted;
	assert( is_sorted(begin(std::move(d2D_cref)), end(std::move(d2D_cref))) ); 

	double const d2D_prime[4][5] {
		{ 0,  1,  2,  3,  4}, 
		{ 5,  6,  7,  8,  9}, 
		{10, 11, 12, 13, 14}, 
		{15, 16, 17, 18, 19}
	};

	multi::array_cref<double, 2> d2D_prime_cref(&d2D_prime[0][0], {4, 5});
//	assert( std::move(d2D_cref) == std::move(d2D_prime_cref) ); // deep comparison
	assert( std::move(d2D_cref)[1][2] == std::move(d2D_prime_cref)[1][2] );
	assert( &std::move(d2D_cref)[1][2] != &std::move(d2D_prime_cref)[1][2] );
	assert( not(std::move(d2D_cref) != std::move(d2D_prime_cref)) );
	assert( not(std::move(d2D_cref) < std::move(d2D_cref)) );
	assert( not(std::move(d2D_cref) > std::move(d2D_cref)) );
//	assert( std::move(d2D_cref) <= std::move(d2D_cref) );
//	assert( std::move(d2D_cref) >= std::move(d2D_cref) );

	assert(( d2D_prime_cref[std::array<int, 2>{2, 3}] == 13 ));

	double const d2D_null[4][5] {
		{ 0,  0,  0,  0,  0}, 
		{ 0,  0,  0,  0,  0}, 
		{ 0,  0,  0,  0,  0}, 
		{ 0,  0,  0,  0,  0}
	};
#if defined(__INTEL_COMPILER) or (defined(__GNUC__) && (__GNUC__<6))
	auto&& d2D_null_cref = multi::make_array_ref<2>(&d2D_null[0][0], multi::iextensions<2>{4, 5});
#else
	auto&& d2D_null_cref = multi::make_array_ref(&d2D_null[0][0], {4, 5});
#endif
	(void)d2D_null_cref;
	

	using std::min;
//	assert( &min(d2D_null_cref, d2D_cref) == &d2D_null_cref );
	using std::max;
//	assert( &max(d2D_null_cref, d2D_cref) == &d2D_cref );

	using std::find;
	auto f = find(begin(std::move(d2D_cref)), end(std::move(d2D_cref)), std::move(d2D_cref)[2]);
	assert( f != end(std::move(d2D_cref)) );
	cout<<"distance = " << std::distance(begin(std::move(d2D_cref)), f) <<std::endl;
	assert( f == begin(std::move(d2D_cref)) + 2 );
	assert( *f == *(begin(std::move(d2D_cref)) + 2) );
	assert( *f == std::move(d2D_cref)[2] );
	assert( & (*f)[3] == &std::move(d2D_cref)[2][3] );

	
	using std::find_if;
	auto fif1 = find_if(begin(std::move(d2D_cref)), end(std::move(d2D_cref)), [](auto&& e){return e[3] == 8.111;});
	assert( fif1 == end(std::move(d2D_cref)) );

	using std::find_if;
	auto fif2 = find_if(begin(std::move(d2D_cref)), end(std::move(d2D_cref)), [](auto&& row){return std::move(row)[3] == 8.;});
	assert( fif2 != end(std::move(d2D_cref)) );
	assert( fif2->operator[](4) == 9. );


	using std::count;
//	assert( count(begin(d2D_cref), end(d2D_cref), d2D_prime_cref[3]) == 1 );


	using std::min_element;
	using std::max_element;

	assert( min_element(begin(std::move(d2D_cref)), end(std::move(d2D_cref))) == begin(std::move(d2D_cref)) );
	assert( max_element(begin(std::move(d2D_cref)), end(std::move(d2D_cref))) == begin(std::move(d2D_cref)) + size(d2D_cref) - 1 );

	using std::minmax_element;
	assert( minmax_element(begin(std::move(d2D_cref)), end(std::move(d2D_cref))).first == min_element(begin(std::move(d2D_cref)), end(std::move(d2D_cref))) );
	assert( minmax_element(begin(std::move(d2D_cref)), end(std::move(d2D_cref))).first == min_element(begin(std::move(d2D_cref)), end(std::move(d2D_cref))) );

	{
		decltype(d2D_cref)::iterator it; 
		decltype(d2D_cref)::iterator it2{};
		decltype(d2D_cref)::iterator it3{0};
		decltype(d2D_cref)::iterator it4{nullptr};
		decltype(d2D_cref)::iterator it5(0);
		assert( it == it2 and it2 == it3 and it3 == it4 and it4 == it5 );
	//	assert(not it); // there are not null iterator

	//	return 0;
	//	assert( std::addressof(it->operator[](0)) == nullptr);
		it = begin(std::move(d2D_cref));
		assert(it == begin(std::move(d2D_cref)));
		it = decltype(it){};

		std::vector<double>::iterator vit;
		std::list<double>::iterator lit{nullptr};
		assert( std::addressof(*vit) == nullptr );
	}

	auto NX = 2, NY = 2, NZ = 2;
	std::vector<double> v(NX*NY*NZ);
	iota(begin(v), end(v), 0.);

	multi::array_cref<double, 3> v3D_cref(v.data(), {NX, NY, NZ});
	assert( num_elements(v3D_cref) == multi::size_type(v.size()) );

	{
	auto x = extensions(v3D_cref); using std::get;
	for(auto i : get<0>(x)) for(auto j : get<1>(x)) for(auto k : get<2>(x))
		cout<< i <<' '<< j <<' '<< k <<'|'<< v3D_cref[i][j][k] <<'\n';
	}

	cout<< v3D_cref[1][1][1] <<"\n\n";
//	cout<< *std::move(v3D_cref).begin()->begin()->begin()   <<'\n\n';
//	cout<< std::move(v3D_cref).begin()->begin()->begin()[0] <<'\n\n';
	{
		std::vector<double> v(100, 1.);
		std::vector<double> w(100, 2.);
		auto f = [](std::vector<double>& a){return multi::array_ref<double,1>(a.data(), {100});}; //	auto V = f(v); //	V = f(w);
		f(v) = f(w);
		assert( v[10] == 2. );

		std::vector<double> x(100, 3.);
		auto g = [](std::vector<double>& a){return multi::array_ref<double,2>(a.data(), {10, 10});}; //	auto V = f(v); //	V = f(w);
		g(v) = g(x);
		assert( v[10] == 3. );
	}
	{
		std::vector<double> a(100, 1.), b(100, 2.);
		multi::array_ref<double, 1> A1(a.data(), {100}), B1(b.data(), {100});
		std::move(A1) = B1;
		assert( a[10] == b[10] );
		assert( A1[10] == B1[10] and &A1[10] != &B1[10] );

		iota(begin(a), end(a), 0.); iota(begin(b), end(b), 10.);
		multi::array_ref<double, 2> A2(a.data(), {10, 10}), B2(b.data(), {10, 10});
		std::move(A2) = std::move(B2);
		assert( a[10] == b[10] );
		assert( A2[4][5] == B2[4][5] and &A2[4][5] != &B2[4][5] );
	}

}
#endif

