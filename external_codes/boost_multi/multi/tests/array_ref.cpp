#ifdef COMPILATION_INSTRUCTIONS
c++ -O3 `#-DNDEBUG` -std=c++14 -Wall -Wextra -Wfatal-errors -I$HOME/prj $0 -o $0.x && $0.x $@ && rm -f $0.x; exit
#endif

#include "../array_ref.hpp"
#include "../array.hpp"

#include<algorithm> // for sort
#include<iostream> // for print
#include<vector>
#include<cmath>

namespace multi = boost::multi;
using std::cout; using std::cerr;

#include<numeric>
#include<vector>

auto f(){
	std::vector<std::vector<double>> v(10, std::vector<double>(3));
	iota(v[5].begin(), v[5].end(), 1);
	return v;
}

#include<string>

template<class M>
void some_assign21(M&& m, double d){m[2][1] = d;}
template<class M>
void some_assign1(M&& m, double d){m[1] = d;}

void some_assignref21(boost::multi::array_ref<double, 2>& m, double d){m[2][1] = d;}


template<class T>
struct vect{
	T* data_;
	T& operator[](std::size_t idx){return *(data_ + idx);}
//	T const& operator[](std::size_t idx) const{assert(0); return *(data_ + idx);}
//	T&& operator[](std::size_t idx)&&{return std::move(*(data_ + idx));}
	vect() : data_(new T[10]){}
	~vect(){delete[] data_;}
};

template<class M1, class M2>
void assign4(M1&& m1, M2 const& m2){m1 = m2;}

template<class T>
struct ptr{
	T* impl_;
	T& operator*() const{return *impl_;}
	auto operator+(std::ptrdiff_t n) const{return ptr{impl_ + n};}
	T& operator[](std::ptrdiff_t n) const{return impl_[n];}
};

int f(int a){return a + 5;}

//#include <boost/timer/timer.hpp>

int main(){
	{
		std::unique_ptr<double[]> ubuf{new double[10000]};
		multi::array_ref<double, 2> const mbuf(ubuf.get(), {100, 100});

		double* Pbuf = new double[10000]; std::fill_n(Pbuf, 10000, 3.);

		multi::array_ref<double, 1, double*&> mbuf2{Pbuf, {{{10000}}}};

		assert(mbuf2.size() == 10000);

		assert( mbuf2[10] == 3. );

		Pbuf = new double[10000]; 
		std::fill_n(Pbuf, 10000, 4.);
		cout << mbuf2[10] << std::endl;
		assert( mbuf2[10] == 4. );

	}
	return 0;
#if 1
//	std::unique_ptr<double[]> 
	auto ubuf{new double[400000000]};
//	std::fill_n(ubuf.get(), 400000000, 99.);
//	std::iota(ubuf.get(), 
	multi::array_ref<double, 2> mbuf(ubuf, {20000,20000});
	std::iota(mbuf.data(), mbuf.data() + mbuf.num_elements(), 0.);
	{
	//	std::unique_ptr<double[]> 
		auto ubuf2{new double[400000000]};
		multi::array_ref<double, 2> mbuf2(ubuf2, {20000,20000});
	//	std::iota(mbuf2.data(), mbuf2.data() + mbuf2.num_elements(), 0.);
		std::iota(mbuf2.data(), mbuf2.data() + mbuf2.num_elements(), 122223.);
	//	std::copy_n(mbuf.data(), mbuf.num_elements(), mbuf2.data());
	//	{
	//		boost::timer::auto_cpu_timer t;
		//	std::copy_n(mbuf.data(), mbuf.num_elements(), mbuf2.data());
	//		mbuf2 = mbuf;
	//	}
		assert(mbuf2[123][456] == mbuf[123][456]);
	//	assert( mbuf2[123][456] == mbuf[123][456] );
		multi::array<double, 2> Mbuf = mbuf2;
	}
	std::vector<double> v(100);
//	auto v2b{new double[100]};
//	multi::array_ref<double, 1> v2(v2b, {100});
//	v2 = v;
	multi::array<double, 1> v2({multi::index_extension{0, (int)v.size()}});
	copy(begin(v), end(v), begin(v2));
	return 0;

	double* buffer = new double[100];

	multi::array_ref<double, 2, ptr<double> > CC(ptr<double>{buffer}, {10, 10} );
	CC[2];
	CC[1][1];
	CC[1][1] = 9;
	assert(CC[1][1] == 9);


	multi::array<double, 2> C3( {4, 3} );
	multi::array<double, 2> C6( {4, 3} );

	assert(C3.size(0)==4 and C3.size(1)==3);

	C6[1][1] = 99.;
	assign4(C3({0,4},{0,3}), C6({0,4},{0,3}));
	assert(C3[1][1] == 99.);
	C3[1][1] = 88.;
	assert(C6[1][1]==99.);

	C3({0,4},{0,3}) = C6({0,4},{0,3});
	C3({0,4}, 1) = C6({0,4}, 2);
	C3({0,4}, 0) = C6({0,4}, 1);

	return 0;
#if 0

	some_assign21(C3, 4.1); assert(C3[2][1]==4.1);
	some_assign21(C3({0,4},{0,3}), 4.2); assert(C3[2][1]==4.2);
	some_assign1(C3[2], 4.); assert(C3[2][1]==4.);
	some_assign1(C3({0,4},2), 4.2); assert(C3[1][2]==4.2);

	multi::array<double, 2> const C4( {4, 3} ); 
	// C4[2][1] = 5.; error C4 is const

	multi::array<double, 2> C5( {4, 3}, 1. ); 
	some_assignref21(C5, 3.);

	using boost::multi::index_extension;
	multi::array<double, 2> A(multi::array<double, 2>::extensions_type{{{0, 4}, {0, 3}}});
//	assert(A.size(0) == 4 and A.size(1) == 3);
	
//	multi::array<double, 2> A(multi::array<double, 2>::extensions_type{{{0, 4}, {5, 8}}}); assert(A.shape()[0] == 4 and A.shape()[1] == 3);
//	multi::array<double, 2> B(multi::layout_t<2>::extensions_type{{{0, 4}, {5, 8}}});
//	assert(B.size(0) == 4 and B.size(1) == 3);

	multi::array<double, 2> C({ {{0, 4}, {5, 8}} }); // C({0,4})
// 	assert(C.size(0) == 4 and C.size(1) == 3);

	multi::array<double, 2> C0( A.extensions() ); assert(C0.size(0)==A.size(0) and C0.size(1)==A.size(1));
	multi::array<double, 2> C1( {4, A.extension(1)} ); assert(C1.size(0)==4 and C1.size(1)==3);
	multi::array<double, 2> C2( {4, {0, 3}} ); assert(C2.size(0)==4 and C2.size(1)==3);

//	multi::array<double, 2> D = {{0., 4.}, {5., 8.}}; 
//	assert(B[1][1] == 8 and B.shape()[0] == 2 and B.shape()[1] == 2);
	return 0;
//	std::pointer_traits<double const*>::rebind<double const> p = 0;//::rebind<double const>::type p;
//	(void)p;

	double d2D[4][5] = {
		{150, 16, 17, 18, 19},
		{ 30,  1,  2,  3,  4}, 
		{100, 11, 12, 13, 14}, 
		{ 50,  6,  7,  8,  9} 
	};
	multi::array_ref<double, 2> d2D_ref{&d2D[0][0], {4, 5}};

	swap(d2D_ref[0], d2D_ref[3]);

	cout << "--\n";
	for(auto i : d2D_ref.extension(0)){
		for(auto j : d2D_ref.extension(1)) 
			cout << d2D_ref[i][j] << ' ';
		cout << '\n';
	}
	cout << '\n';

	multi::array_ref<double, 2>::const_reference crow1 = d2D_ref[1];
	assert( crow1[3] == 3 );
	multi::array_ref<double, 2>::iterator it = begin(d2D_ref);
	multi::array_ref<double, 2>::const_iterator cit = begin(d2D_ref);
	assert(it == cit);

	multi::array_ref<double, 2> d2D_ref2{&d2D[0][0], {4, 5}};
	d2D_ref = d2D_ref2;

	for(auto& r: d2D_ref) for(auto& e: r) e = -e;
	
	for(auto i : d2D_ref.extension(0))
		for(auto j : d2D_ref.extension(1))
			d2D_ref[i][j] = -d2D_ref[i][j];

	for(auto i : d2D_ref.extension(0)){
		for(auto j : d2D_ref.extension(1))
			cout << d2D_ref[i][j] <<' ';
		cout <<'\n';
	}

	using std::stable_sort;
	stable_sort( d2D_ref.begin(0), d2D_ref.end(0) );

	cout << "--\n";
	for(auto i : d2D_ref.extension(0)){
		for(auto j : d2D_ref.extension(1))
			cout << d2D_ref[i][j] << ' ';
		cout << '\n';
	}
	stable_sort( d2D_ref.begin(1), d2D_ref.end(1) );

	cout << "--\n";	
	for(auto i : d2D_ref.extensions()[0]){
		for(auto j : d2D_ref.extensions()[1])
			cout << d2D_ref[i][j] << ' ';
		cout << '\n';
	}
	swap(*begin(d2D_ref), *(begin(d2D_ref) + 1));
//	swap(*d2D_ref.begin(), *(d2D_ref.begin() + 1));

	cout << "--\n";	
	for(auto i : d2D_ref.extension(0)){
		for(auto j : d2D_ref.extension(1))
			cout << d2D_ref[i][j] << ' ';
		cout << '\n';
	}
	std::reverse(d2D_ref.begin(1), d2D_ref.end(1));
	cout << "--\n";	
	for(auto i : d2D_ref.extension(0)){
		for(auto j : d2D_ref.extension(1))
			cout << d2D_ref[i][j] << ' ';
		cout << '\n';
	}
	{
		cout << "--\n";
		auto d2D_take = d2D_ref({1, 3}, {1, 2});
		assert(d2D_take.dimensionality == 2);
		for(auto i : d2D_take.extension(0)){
			for(auto j : d2D_take.extension(1))
				cout << d2D_take[i][j] << ' ';
			cout << '\n';
		}
	}
	{
		cout << "--\n";
		auto d1D_take = d2D_ref(2, {1, 4});
		assert(d1D_take.dimensionality == 1);
		for(auto i : d1D_take.extension(0))
			cout << d1D_take[i] << ' '; 
		cout << '\n';
		assert( d1D_take.shape()[0] == 3 );
		assert( d1D_take.sizes()[0] == 3 );
		assert( d1D_take.size(0) == 3 );
		assert( d1D_take.size() == 3 );
	}
	{
		cout << "--\n";
		auto d1D_take = d2D_ref({1, 4}, 3);
		assert(d1D_take.dimensionality == 1);
		for(auto i : d1D_take.extension(0))
			cout << d1D_take[i] << ' ';
	}
	{
		cout << "--\n";
	//	auto d1D_take = d2D_ref(2, 3);
	}

	return 0;

	cout << "--\n";
	auto const& d2D_ref_sv = d2D_ref.range({1, 3});
	for(auto i : d2D_ref_sv.extension(0)){
		for(auto j : d2D_ref_sv.extension(1))
			cout << d2D_ref_sv[i][j] << ' ';
		cout << '\n';
	}
	assert( d2D_ref.range({1,3}).origin() == &d2D_ref_sv[0][0] );
	assert( d2D_ref_sv.origin() == &d2D_ref_sv[0][0] );
	assert( *d2D_ref_sv.origin() == d2D_ref_sv[0][0] );
	assert( d2D_ref.range({1,3}).rotated(1).range({2, 5}).rotated(-1).origin() == &d2D_ref[1][2] );

	return 0;
	multi::array_ref<double, 1> d1D_ref{&d2D[0][0], {5}};

	assert( d2D_ref.cdata() == cdata(d2D_ref) );
	assert( d2D_ref.data() == data(d2D_ref) );
	assert( data(d2D_ref) == &d2D[0][0] );
	*data(d2D_ref) = 1;
	assert(d2D[0][0] == 1);
	for(auto i : d2D_ref.extension(0))
		for(auto j : d2D_ref.extension(1)) 
			d2D_ref[i][j] = 10.*i + j;

	for(auto i : d2D_ref.extension(0)){
		for(auto j : d2D_ref.extension(1))
			cout << d2D_ref[i][j] << ' ';
		cout << '\n';
	}
	cout << '\n';
	multi::array_ref<double, 2> d2D_refref{data(d2D_ref), extensions(d2D_ref)};
	multi::array_ref<double, 2> d2D_refref2{d2D_ref}; (void)d2D_refref2;

	for(auto it1 = begin(d2D_ref); it1 != end(d2D_ref) ||!endl(cout); ++it1)
		for(auto it2 = it1->begin()   ; it2 != it1->end()    ||!endl(cout); ++it2)
			*it2 = 99.;

	assert(d2D_ref[2][2] == 99.);

	double n = 0;
	for(auto it1 = d2D_ref.begin(1); it1 != d2D_ref.end(1) ||!endl(cout); ++it1)
		for(auto it2 = it1->begin()   ; it2 != it1->end()    ||!endl(cout); ++it2)
			*it2 = n++;

	for(auto i : d2D_ref.extension(0)){
		for(auto j : d2D_ref.extension(1)){
			cout << d2D_ref[i][j] << ' ';
		}
		cout << '\n';
	}

	*(d2D_ref.begin()->begin())=88.;
	assert( d2D_ref[0][0] == 88. );

	for(auto i : d2D_ref.extension(0)){
		for(auto j : d2D_ref.extension(1)){
			cout << d2D_ref[i][j] << ' ';
		}
		cout << '\n';
	}

	assert( d2D_ref[0] > d2D_ref[1] );
	using std::swap;
	swap( d2D_ref[2], d2D_ref[3] );

//	using std::stable_sort;
	std::stable_sort( begin(d2D_ref), end(d2D_ref) );
//	swap( d2D_ref[0], d2D_ref[1] );
	
	for(auto i : d2D_ref.extension(0)){
		for(auto j : d2D_ref.extension(1))
			cout << d2D_ref[i][j] << ' ';
		cout << '\n';
	}

/*	using std::for_each;
    using namespace std::string_literals; //""s
	for_each(begin(d2D_ref), end(d2D_ref), [](auto&& row){
		for_each(begin(row), end(row), [](auto&& element){
			element = 99;
		});
	});
*/
#if 0

	for(auto it1 = begin(d2D_cref); it1 != end(d2D_cref) ||!endl(cout); ++it1)
		for(auto it2 = it1->begin()   ; it2 != it1->end()    ||!endl(cout); ++it2)
			cout << *it2 << ' ';

//	d2D_cref[3][1] = 3.; // error const ref not assignable

	double const d2D_prime[4][5] {
		{ 0,  1,  2,  3,  4}, 
		{ 5,  6,  7,  8,  9}, 
		{10, 11, 12, 13, 14}, 
		{15, 16, 17, 18, 19}
	};

	multi::array_cref<double, 2> d2D_prime_cref{&d2D_prime[0][0], {4, 5}};
//	multi::array_cref<double, 2> d2D_prime_cref{&d2D_prime[0][0], extensions(d2D_cref)};
	assert( d2D_cref == d2D_prime_cref ); // deep comparison
//	assert( d2D_cref == d2D_prime );
	assert( not(d2D_cref != d2D_prime_cref) );
	assert( not(d2D_cref < d2D_cref) );
	assert( not(d2D_cref > d2D_cref) );
	assert( d2D_cref <= d2D_cref );
	assert( d2D_cref >= d2D_cref );

	double const d2D_null[4][5] {
		{ 0,  0,  0,  0,  0}, 
		{ 0,  0,  0,  0,  0}, 
		{ 0,  0,  0,  0,  0}, 
		{ 0,  0,  0,  0,  0}
	};
	multi::array_cref<double, 2> d2D_null_cref{&d2D_null[0][0], {4, 5}};

	using std::min;
	assert( &min(d2D_null_cref, d2D_cref) == &d2D_null_cref );
	using std::max;
	assert( &max(d2D_null_cref, d2D_cref) == &d2D_cref );
	
	using std::find;
	auto f = find(begin(d2D_cref), end(d2D_cref), d2D_cref[2]);
	assert( f != end(d2D_cref) );
	assert( (*f)[3] == d2D_cref[2][3] );

	using std::find_if;
	auto fif1 = find_if(begin(d2D_cref), end(d2D_cref), [](auto&& e){return e[3] == 8.111;});
	assert( fif1 == end(d2D_cref) );

	using std::find_if;
	auto fif2 = find_if(begin(d2D_cref), end(d2D_cref), [](auto&& e){return e[3] == 8.;});
	assert( fif2 != end(d2D_cref) );
	assert( fif2->operator[](4) == 9. );

	using std::count;
	assert( count(begin(d2D_cref), end(d2D_cref), d2D_prime_cref[3]) == 1 );
	assert( count(begin(d2D_cref), end(d2D_cref), d2D_prime[3]     ) == 1 );

	using std::min_element;
	using std::max_element;

	assert( min_element(begin(d2D_cref), end(d2D_cref)) == begin(d2D_cref) );
	assert( max_element(begin(d2D_cref), end(d2D_cref)) == begin(d2D_cref) + size(d2D_cref) - 1 );

	using std::minmax_element;
	assert( minmax_element(begin(d2D_cref), end(d2D_cref)).first == min_element(begin(d2D_cref), end(d2D_cref)) );
	assert( minmax_element(begin(d2D_cref), end(d2D_cref)).first == min_element(begin(d2D_cref), end(d2D_cref)) );
	decltype(d2D_cref)::const_iterator it; // it{} == it{0} == it{nullptr} = it(0);
//	assert(not it); // there are not null iterators
	assert( std::addressof(it->operator[](0)) == nullptr);
	it = begin(d2D_cref);
	assert(it == begin(d2D_cref));
	it = decltype(it){};

	std::vector<double>::iterator vit;
	std::list<double>::iterator lit{nullptr};
	assert( std::addressof(*vit) == nullptr );

	std::ptrdiff_t NX = 2;
	std::ptrdiff_t NY = 2;
	std::ptrdiff_t NZ = 2;
	std::vector<double> v(NX*NY*NZ);
	iota(begin(v), end(v), 0.);

	multi::array_cref<double, 3> v3D_cref{v.data(), {NX, NY, NZ}};

	assert( v3D_cref.num_elements() == multi::size_type(v.size()) );
	for(auto i : v3D_cref.extension(0))
		for(auto j : v3D_cref.extension(1))
			for(auto k : v3D_cref.extension(2))
				cout << i << ' ' << j << ' ' << k << ' ' 
					<< v3D_cref[i][j][k] << '\n';

	cout << v3D_cref[9][9][9] << "\n\n";

	assert(d2D_cref.begin() == d2D_cref.begin(0));
	assert(d2D_cref.begin() != d2D_cref.begin(1));
	for(auto it1 = d2D_cref.begin(1); it1 != d2D_cref.end(1)||!endl(cout); ++it1)
		for(auto it2 = it1->begin()   ; it2 != it1->end()   ||!endl(cout); ++it2)
			cout << *it2 << ' ';
#endif
#endif
#endif
}

