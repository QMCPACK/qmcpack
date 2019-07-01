#ifdef COMPILATION_INSTRUCTIONS
$CXX -O3 -std=c++17 -Wall -Wextra -Wpedantic `#-Wfatal-errors` $0 -o $0.x && $0.x $@ && rm -f $0.x; exit
#endif

#include "../array_ref.hpp"
#include "../array.hpp"

#include<algorithm>
#include<cassert>
#include<iostream>
#include<cmath>
#include<vector>
#include<list>
#include<numeric> //iota

using std::cout; using std::cerr;
namespace multi = boost::multi;

double f(){return 5.;}

int main(){
//	using multi::size;
	{
		double a[4][5] {
			{ 0,  1,  2,  3,  4}, 
			{ 5,  6,  7,  8,  9}, 
			{10, 11, 12, 13, 14}, 
			{15, 16, 17, 18, 19}
		};
		double b[4][5];
		multi::array_ref<double, 2, double*> A(&a[0][0], multi::iextensions<2>{4, 5});
		multi::array_ref<double, 2, double*> B(&b[0][0], multi::iextensions<2>{4, 5});
		rotated(A) = rotated(B);
	}

	double const d2D[4][5] {
		{ 0,  1,  2,  3,  4}, 
		{ 5,  6,  7,  8,  9}, 
		{10, 11, 12, 13, 14}, 
		{15, 16, 17, 18, 19}
	};
	multi::array_ref<double, 2, double const*> d2D_cref(&d2D[0][0], multi::iextensions<2>{4, 5});

	decltype(d2D_cref)::value_type a_row = d2D_cref[2];
	decltype(d2D_cref[2])::decay_type b_row = d2D_cref[2];
	auto c_row = d2D_cref[2].decay();
	auto d_row = decay(d2D_cref[2]);
	auto const& e_row = d2D_cref[2];
	assert( e_row[3] == d2D_cref[2][3] );

	std::vector<double> d2Dv(20); iota(d2Dv.begin(), d2Dv.end(), 0);
	multi::array_ref<double, 2, double const*> d2Dv_cref(d2Dv.data(), {4, 5});
	assert( d2Dv.size() == std::size_t(num_elements(d2Dv_cref)) );
	assert(d2Dv_cref[1][1] == 6);

	assert( d2D_cref.data_elements() == data_elements(d2D_cref) );
	assert( data_elements(d2D_cref) == &d2D[0][0] );
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
	assert( &data_elements(d2D_cref)[2] == &d2D_cref[0][2] );
	assert( &data_elements(d2D_cref)[6] == &d2D_cref[1][1] );

	assert( d2D_cref.begin() == begin(d2D_cref) );
	assert( d2D_cref.end() == end(d2D_cref) );
	assert( begin(d2D_cref) != end(d2D_cref) );
	assert( begin(d2D_cref) + size(d2D_cref) == end(d2D_cref) );
//	assert( end(d2D_cref) - begin(d2D_cref) == size(d2D_cref) );
	using std::distance;
	assert( distance(begin(d2D_cref), end(d2D_cref)) == size(d2D_cref));
	assert( size(*begin(d2D_cref)) == 5 );
	assert( distance(begin(d2D_cref)->begin(), begin(d2D_cref)->end()) == begin(d2D_cref)->size() );
	assert( distance(begin(*begin(d2D_cref)), end(*begin(d2D_cref))) == size(*begin(d2D_cref)) );
	assert( size(d2D_cref[0]) == 5 );
	assert( num_elements(d2D_cref[0]) == 5 );

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
			
	assert(d2D_cref.begin() == d2D_cref.begin(0));
	assert(d2D_cref.begin() != d2D_cref.begin(1));
	for(auto it1 = d2D_cref.begin(1); it1 != d2D_cref.end(1)||!endl(cout); ++it1)
		for(auto it2 = it1->begin()   ; it2 != it1->end()   ||!endl(cout); ++it2)
			cout << *it2 << ' ';

	multi::array_ref<double, 2, double const*> d2D_crefref(data_elements(d2D_cref), extensions(d2D_cref));
	multi::array_cref<double, 2> d2D_crefcref(data_elements(d2D_cref), extensions(d2D_cref));

	using std::for_each;
	using std::begin;
	using std::end;
	for_each(
		begin(d2D_cref), end(d2D_cref), 
		[](auto&& row){
			for_each(begin(row), end(row), [](auto&& e){cout<<' '<< e;})("\n");
		}
	)("\n");
	
	for(auto& row: d2D_cref){for(auto& el: row) cout<< el <<' '; cout<<'\n';}

//	for(decltype(d2D_cref)::reference r : d2D_cref){
//		for(decltype(d2D_cref)::element const& e: r) cout << e <<' '; 
//		cout <<'\n';
//	}
	for(auto it = begin(d2D_cref); it != end(d2D_cref); ++it){
		for(decltype(d2D_cref)::element const& e: *it) cout << e <<' '; 
		cout <<'\n';
	}

	
	for(auto it1=begin(d2D_cref); it1 != end(d2D_cref) ||!endl(cout); ++it1)
		for(auto it2=it1->begin()   ; it2 != it1->end()    ||!endl(cout); ++it2)
			cout<< *it2 <<' ';

	for(auto it1=(&d2D_cref)->begin(); it1 != (&d2D_cref)->end() ||!endl(cout); ++it1)
		for(auto it2=it1->begin()   ; it2 != it1->end()    ||!endl(cout); ++it2)
			cout<< *it2 <<' ';

	auto print = [](auto&& arr){
		for(auto it1 = arr.begin(); it1 != arr.end()||!endl(cout); ++it1)
			for(auto it2 = it1->begin()   ; it2 != it1->end()   ||!endl(cout); ++it2)
				cout << *it2 << ' ';
	};
	print(d2D_cref);
	print(d2D_cref.range({0, 2}));
	print(d2D_cref.rotated(1).range({0, 2}).rotated(1));


	using std::is_sorted;
	assert( is_sorted(begin(d2D_cref), end(d2D_cref)) ); 

//	d2D_cref[3][1] = 3.; // cannot assign to const value

	double const d2D_prime[4][5] {
		{ 0,  1,  2,  3,  4}, 
		{ 5,  6,  7,  8,  9}, 
		{10, 11, 12, 13, 14}, 
		{15, 16, 17, 18, 19}
	};

	multi::array_cref<double, 2> d2D_prime_cref(&d2D_prime[0][0], {4, 5});
	assert( d2D_cref == d2D_prime_cref ); // deep comparison
	assert( d2D_cref[1][2] == d2D_prime_cref[1][2] );
	assert( &d2D_cref[1][2] != &d2D_prime_cref[1][2] );
	assert( not(d2D_cref != d2D_prime_cref) );
	assert( not(d2D_cref < d2D_cref) );
	assert( not(d2D_cref > d2D_cref) );
	assert( d2D_cref <= d2D_cref );
	assert( d2D_cref >= d2D_cref );

	assert(( d2D_prime_cref[std::array<int, 2>{2, 3}] == 13 ));

	double const d2D_null[4][5] {
		{ 0,  0,  0,  0,  0}, 
		{ 0,  0,  0,  0,  0}, 
		{ 0,  0,  0,  0,  0}, 
		{ 0,  0,  0,  0,  0}
	};
#if defined(__INTEL_COMPILER) or (defined(__GNUC__) && (__GNUC__<6))
	auto d2D_null_cref = multi::make_array_ref<2>(&d2D_null[0][0], multi::iextensions<2>{4, 5});
#else
	auto d2D_null_cref = multi::make_array_ref(&d2D_null[0][0], {4, 5});
#endif
	using std::min;
	assert( &min(d2D_null_cref, d2D_cref) == &d2D_null_cref );
	using std::max;
	assert( &max(d2D_null_cref, d2D_cref) == &d2D_cref );

	using std::find;
	auto f = find(begin(d2D_cref), end(d2D_cref), d2D_cref[2]);
	assert( f == begin(d2D_cref) + 2 );
	assert( *f == *(begin(d2D_cref) + 2) );
	assert( *f == d2D_cref[2] );
	assert( & (*f)[3] == &d2D_cref[2][3] );
	
	using std::find_if;
	auto fif1 = find_if(begin(d2D_cref), end(d2D_cref), [](auto&& e){return e[3] == 8.111;});
	assert( fif1 == end(d2D_cref) );

	using std::find_if;
	auto fif2 = find_if(begin(d2D_cref), end(d2D_cref), [](auto&& row){return row[3] == 8.;});
	assert( fif2 != end(d2D_cref) );
	assert( fif2->operator[](4) == 9. );

	using std::count;
	assert( count(begin(d2D_cref), end(d2D_cref), d2D_prime_cref[3]) == 1 );

	using std::min_element;
	using std::max_element;

	assert( min_element(begin(d2D_cref), end(d2D_cref)) == begin(d2D_cref) );
	assert( max_element(begin(d2D_cref), end(d2D_cref)) == begin(d2D_cref) + size(d2D_cref) - 1 );

	using std::minmax_element;
	assert( minmax_element(begin(d2D_cref), end(d2D_cref)).first == min_element(begin(d2D_cref), end(d2D_cref)) );
	assert( minmax_element(begin(d2D_cref), end(d2D_cref)).first == min_element(begin(d2D_cref), end(d2D_cref)) );

	{
		decltype(d2D_cref)::iterator it; 
		decltype(d2D_cref)::iterator it2{};
		decltype(d2D_cref)::iterator it3{0};
		decltype(d2D_cref)::iterator it4{nullptr};
		decltype(d2D_cref)::iterator it5(0);
		assert( it == it2 and it2 == it3 and it3 == it4 and it4 == it5 );
	//	assert(not it); // there are not null iterator
		assert( std::addressof(it->operator[](0)) == nullptr);
		it = begin(d2D_cref);
		assert(it == begin(d2D_cref));
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

	cout << v3D_cref[1][1][1] << "\n\n";
	cout << *v3D_cref.begin()->begin()->begin() << "\n\n";
	cout << v3D_cref.begin()->begin()->begin()[0] << "\n\n";
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
		A1 = B1;
		assert( a[10] == b[10] );
		assert( A1[10] == B1[10] and &A1[10] != &B1[10] );

		iota(begin(a), end(a), 0.); iota(begin(b), end(b), 10.);
		multi::array_ref<double, 2> A2(a.data(), {10, 10}), B2(b.data(), {10, 10});
		A2 = B2;
		assert( a[10] == b[10] );
		assert( A2[4][5] == B2[4][5] and &A2[4][5] != &B2[4][5] );
	}
}

