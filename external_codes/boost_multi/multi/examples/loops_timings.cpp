#ifdef COMPILATION// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;-*-
$CXX $0 -o $0x -ltbb&&$0x&&rm $0x;exit
#endif
//  © Alfredo A. Correa 2018-2020

#include "../array.hpp"

#include<numeric>
#include<algorithm> // generate
#include<chrono>
#include<complex>
#include<iostream>
#include<random>
#include<vector>
#if __cplusplus >= 201703L
#include<execution>
#endif
#include<map>
#include<numeric>

using watch = std::chrono::high_resolution_clock;
using ns = std::chrono::nanoseconds;

namespace multi = boost::multi;
using std::cout;
using std::vector;
using complex = std::complex<double>;

const int nbasis = 1000000;
const int nstates = 32;

struct tic{
	std::string s;
	decltype(watch::now()) t;
	tic(std::string s) : s{s}, t{watch::now()}{}
	~tic(){std::cout<<s<<ns{watch::now() - t}.count()/1e9 <<"sec"<<std::endl;}
};

int main(){

std::map<decltype(watch::now() - watch::now()), std::string> timings;

std::random_device dev;
std::mt19937 eng{dev()};
std::uniform_real_distribution<double> dist{-1., 1.};
auto gen = [&dist, &eng](){return complex(dist(eng), dist(eng));};
{
	multi::array<complex, 2> Asb({nstates, nbasis}, complex{});
	for_each(begin(Asb), end(Asb), [&](auto&& e){generate(begin(std::move(e)), end(std::move(e)), gen);});
	vector<double> d(nbasis);
	{
		auto tic = watch::now();
		fill(begin(d), end(d), 0.);
		for(std::size_t s = 0; s != nstates; ++s)
			for(std::size_t b = 0; b != nbasis; ++b) 
				d[b] += norm(Asb[s][b]);
		timings[watch::now() - tic] = "array[state][basis] storage, raw loop state, raw loop basis";
		cout<< d[d.size()/3] <<std::endl;
	}
#if 0
	{
		auto tic = watch::now();
		fill(begin(d), end(d), 0.);
#pragma omp parallel for 
		for(std::size_t s = 0; s != nstates; ++s)
			for(std::size_t b = 0; b != nbasis; ++b) 
				d[s] += norm(Asb[s][b]);
		timings[watch::now() - tic] = "array[state][basis] storage, raw omp loop state, raw loop basis";
		cout<< d[d.size()/3] <<std::endl;
	}
#endif
	{
		auto tic = watch::now();
		transform(begin(Asb), end(Asb), begin(d),
			[](auto const& a){
				double d = 0; 
				for(std::size_t b = 0; b != nbasis; ++b) d += norm(a[b]);
				return d;
			}
		);
		timings[watch::now() - tic] = "array[state][basis] storage, transform states, raw loop basis";
		cout<< d[d.size()/3] <<std::endl;
	}
	{
		auto tic = watch::now();
		transform(begin(Asb), end(Asb), begin(d),
			[](auto&& a){
				return std::accumulate(begin(std::move(a)), end(std::move(a)), 0., [](auto&& acc, auto&& e){return acc + norm(e);});
			}
		);
		timings[watch::now() - tic] = "array[state][basis] storage, transform states, accumulate basis";
		cout<< d[d.size()/3] <<std::endl;
	}
#if __cpp_lib_execution >= 201603
	{
		auto tic = watch::now();
		transform(begin(Asb), end(Asb), begin(d),
			[](auto&& a){
				return std::reduce(std::execution::seq, begin(std::move(a)), end(std::move(a)), 0., [](auto&& acc, auto&& e){return acc + norm(e);});
			}
		);
		timings[watch::now() - tic] = "array[state][basis] storage, transform states, sequential reduce basis"; 
		cout<< d[d.size()/3] <<std::endl;
	}
#endif
#if __cpp_lib_execution >= 201603
	{
		auto tic = watch::now();
		transform(std::execution::par, begin(Asb), end(Asb), begin(d),
			[](auto&& a){
				return std::accumulate(begin(std::move(a)), end(std::move(a)), 0., [](auto&& acc, auto&& e){return acc + norm(e);});
			}
		);
		timings[watch::now() - tic] = "array[state][basis] storage, par transform states, accumulate basis"; 
		cout<< d[d.size()/3] <<std::endl;
	}
#endif
#if __cpp_lib_execution >= 201603
	{
		auto tic = watch::now();
		transform(std::execution::par, begin(Asb), end(Asb), begin(d),
			[](auto&& state){
				auto const norm = [](auto const& e){return std::norm(e);};
				return transform_reduce(std::execution::par, 
					begin(std::move(state)), end(std::move(state)), 0., std::plus{}, norm
				);
			}
		);
		timings[watch::now() - tic] = "array[state][basis] storage, par transform states, par transform reduce basis"; 
		cout<< d[d.size()/3] <<std::endl;
	}
#endif
#if __cpp_lib_execution >= 201603
	{
		auto tic = watch::now();
		transform(begin(Asb), end(Asb), begin(d),
			[](auto&& a){
				return transform_reduce(std::execution::par, begin(std::move(a)), end(std::move(a)), 0., std::plus<>{}, [](auto&& e){return std::norm(e);});
			}
		);
		timings[watch::now() - tic] = "array[state][basis] storage, seq transform states, par transform reduce basis"; 
		cout<< d[d.size()/3] <<std::endl;
	}
#endif
#if __cpp_lib_execution >= 201603
	{
		auto tic = watch::now();
		transform(std::execution::seq, begin(Asb), end(Asb), begin(d),
			[](auto&& a){
				return transform_reduce(std::execution::seq, begin(std::move(a)), end(std::move(a)), 0., std::plus<>{}, [](auto&& e){return std::norm(e);});
			}
		);
		timings[watch::now() - tic] = "array[state][basis] storage, seq transform states, seq par transform reduce basis"; 
		cout<< d[d.size()/3] <<std::endl;
	}
#endif
}

cout<<"------------"<<std::endl;
{
	multi::array<complex, 2> Abs({nbasis, nstates}, complex{});
	for_each(begin(Abs), end(Abs), [&](auto&& e){generate(begin(std::move(e)), end(std::move(e)), gen);});
	vector<double> d(nstates);
	{
		auto tic = watch::now();
		fill(begin(d), end(d), 0.);
		for(std::size_t s = 0; s != nstates; ++s){
			for(std::size_t b = 0; b != nbasis; ++b) d[s] += norm(Abs[b][s]);
		}
		timings[watch::now() - tic] = "array[basis][state] storage, raw loop states, raw loop basis"; 
		cout<< d[d.size()/3] <<std::endl;
	}
	{
		auto tic = watch::now();
		fill(begin(d), end(d), 0.);
		for(std::size_t b = 0; b != nbasis; ++b){
			for(std::size_t s = 0; s != nstates; ++s) d[s] += norm(Abs[b][s]);
		}
		timings[watch::now() - tic] = "array[basis][state] storage, raw loop basis, raw loop states"; 
		cout<< d[d.size()/3] <<std::endl;
	}
#if __cpp_lib_execution >= 201603
	{
		auto tic = watch::now();
		fill(std::execution::par, begin(d), end(d), 0.);
		for(std::size_t b = 0; b != nbasis; ++b){
			for(std::size_t s = 0; s != nstates; ++s) d[s] += norm(Abs[b][s]);
		}
		timings[watch::now() - tic] = "array[basis][state] storage, par fill, raw loop basis, raw loop states"; 
		cout<< d[d.size()/3] <<std::endl;
	}
#endif
	{
		auto tic = watch::now();
		fill(begin(d), end(d), 0.);
		for(std::size_t b = 0; b != nbasis; ++b)
			transform(begin(d), end(d), begin(Abs[b]), begin(d), [](auto&& x, auto&& y){return x + norm(y);});
		timings[watch::now() - tic] = "array[basis][state] storage, raw loop basis, transform states"; 
		cout<< d[d.size()/3] <<std::endl;
	}
	{
/*		auto tic = watch::now();
		fill(begin(d), end(d), 0.);
		for_each(begin(Abs), end(Abs), [&](auto&& e){
			transform(std::execution::par_unseq, begin(d), end(d), begin(std::move(e)), begin(d), [](auto&& x, auto&& y){return x + norm(y);});
		});
		timings[watch::now() - tic] = "array[basis][state] storage, for_each, parunseq transform"; 
		cout<< d[d.size()/3]<<std::endl;
*/	}
#if __cpp_lib_execution >= 201603
	{
		auto tic = watch::now();
		fill(begin(d), end(d), 0.);
		for_each(begin(Abs), end(Abs), [&](auto&& e){
			transform(std::execution::par, begin(d), end(d), begin(std::move(e)), begin(d), [](auto&& x, auto&& y){return x + norm(y);});
		});
		timings[watch::now() - tic] = "array[basis][state] storage, fill, for_each, par transform"; 
		cout<<"\t\t\t"<< d[d.size()/3]<<std::endl;
	}
	{
		auto tic = watch::now();
		fill(begin(d), end(d), 0.);
		for_each(std::execution::par, begin(Abs), end(Abs), [&](auto&& e){
			transform(std::execution::par, begin(d), end(d), begin(std::move(e)), begin(d), [](auto&& x, auto&& y){return x + norm(y);});
		});
		timings[watch::now() - tic] = "array[basis][state] storage, fill, for_each, par transform"; 
		cout<<"\t\t\t"<< d[d.size()/3]<<std::endl;
	}
#endif
}
for(auto&& p : timings) cout<< p.first.count()/1e9 <<"\t...."<< p.second  <<std::endl;
return 0;
{
	cout<<"v[basis][state] storage\n";
	vector<vector<complex>> vbs(nbasis, vector<complex>(nstates));	// v[basis][state]
	for_each(begin(vbs), end(vbs), [&](auto& e){generate(begin(e), end(e), gen);});
	vector<double> d(nstates);
	{
	 	cout<<"\traw loops state/basis\n";
		auto tic = watch::now();
		fill(begin(d), end(d), 0.);
		for(std::size_t s = 0; s != nstates; ++s){
			for(std::size_t b = 0; b != nbasis; ++b) d[s] += norm(vbs[b][s]);
		}
		cout<<"\t\t"<< ns{watch::now()-tic}.count()/1e9 <<" seconds\t"<<d[d.size()/2]<<std::endl;
	}
	{
	 	cout<<"\traw loops basis/state\n";
		auto tic = watch::now();
		fill(begin(d), end(d), 0.);
		for(std::size_t b = 0; b != nbasis; ++b){
			for(std::size_t s = 0; s != nstates; ++s) d[s] += norm(vbs[b][s]);
		}
		cout<<"\t\t"<< ns{watch::now()-tic}.count()/1e9 <<" sec\t"<<d[d.size()/2]<<std::endl;
	}
	{
	 	cout<<"\traw loop basis, transform states\n";
		auto tic = watch::now();
		fill(begin(d), end(d), 0.);
		for(std::size_t b = 0; b != nbasis; ++b)
			transform(begin(d), end(d), begin(vbs[b]), begin(d), [](auto const& x, auto const& y){return x + norm(y);});
		cout<<"\t\t"<< ns{watch::now()-tic}.count()/1e9 <<" sec\t"<<d[d.size()/2]<<std::endl;
	}
#if __cpp_lib_execution >= 201603
	{
	 	cout<<"\traw loop basis, parallel transform states\n";
		auto tic = watch::now();
		fill(begin(d), end(d), 0.);
		for(std::size_t b = 0; b != nbasis; ++b) // cannot parallelize this loop
			transform(std::execution::par, begin(d), end(d), begin(vbs[b]), begin(d), [](auto const& x, auto const& y){return x + norm(y);});
		cout<<"\t\t"<< ns{watch::now()-tic}.count()/1e9 <<" sec  \t"<<d[d.size()/2]<<std::endl;
	}
#endif
}
{
	cout<<"v[state][basis] storage\n";
	vector<vector<complex>> vsb(nstates, vector<complex>(nbasis));	// v[state][basis]
	for_each(begin(vsb), end(vsb), [&](auto& e){generate(begin(e), end(e), gen);});
	vector<double> d(nstates);
	{
	 	cout<<"\traw loops state/basis\n";
		auto tic = watch::now();
		fill(begin(d), end(d), 0.);
		for(std::size_t s = 0; s != nstates; ++s){
			for(std::size_t b = 0; b != nbasis; ++b) d[s] += norm(vsb[s][b]);
		}
		auto toc = watch::now();
		cout<<"\t\t"<< ns{toc-tic}.count()/1e9 <<" sec\t"<<d[d.size()/2]<<std::endl;
	}
	{
	 	cout<<"\traw loops basis/state\n";
		auto tic = watch::now();
		fill(begin(d), end(d), 0.);
		for(std::size_t b = 0; b != nbasis; ++b){
			for(std::size_t s = 0; s != nstates; ++s) d[s] += norm(vsb[s][b]);
		}
		auto toc = watch::now();
		cout<<"\t\t"<< ns{toc-tic}.count()/1e9 <<" sec\t"<<d[d.size()/2]<<std::endl;
	}
	{
	 	cout<<"\traw loop basis transform state\n";
		auto tic = watch::now();
		fill(begin(d), end(d), 0.);
		for(std::size_t b = 0; b != nbasis; ++b)
			transform(begin(d), end(d), begin(vsb), begin(d), [&](auto const& x, auto const& y){return x + norm(y[b]);});
		auto toc = watch::now();
		cout<<"\t\t"<< ns{toc-tic}.count()/1e9 <<" sec\t"<<d[d.size()/2]<<std::endl;
	}
#if __cpp_lib_execution >= 201603
	{
	 	cout<<"\traw loop basis transform state\n";
		auto tic = watch::now();
		fill(begin(d), end(d), 0.);
		for(std::size_t b = 0; b != nbasis; ++b)
			transform(std::execution::par, begin(d), end(d), begin(vsb), begin(d), [&](auto const& x, auto const& y){return x + norm(y[b]);});
		auto toc = watch::now();
		cout<<"\t\t"<< ns{toc-tic}.count()/1e9 <<" sec\t"<<d[d.size()/2]<<std::endl;
	}
#endif
	{
	 	cout<<"\traw loop basis transform state\n";
		auto tic = watch::now();
		fill(begin(d), end(d), 0.);
		for(std::size_t b = 0; b != nbasis; ++b)
			transform(begin(d), end(d), begin(vsb), begin(d), [&](auto const& x, auto const& y){return x + norm(y[b]);});
		auto toc = watch::now();
		cout<<"\t\t"<< ns{toc-tic}.count()/1e9 <<" sec\t"<<d[d.size()/2]<<std::endl;
	}
#if 0
	{
	 	cout<<"\traw loops basis/state\n";
		auto tic = watch::now();
		fill(begin(d), end(d), 0.);
		for(std::size_t b = 0; b != nbasis; ++b){
			for(std::size_t s = 0; s != nstates; ++s) d[s] += norm(vbs[b][s]);
		}
		cout<<"\t\t"<< ns{watch::now()-tic}.count()/1e9 <<" seconds\t"<<d[d.size()/2]<<std::endl;
	}
	{
	 	cout<<"\traw loop basis, transform states\n";
		auto tic = watch::now();
		fill(begin(d), end(d), 0.);
		for(std::size_t b = 0; b != nbasis; ++b)
			transform(begin(d), end(d), begin(vbs[b]), begin(d), [](auto const& a, auto const& b){return a + norm(b);});
		cout<<"\t\t"<< ns{watch::now()-tic}.count()/1e9 <<" seconds\t"<<d[d.size()/2]<<std::endl;
	}
	while(0){
	 	cout<<"\traw loop basis, parallel transform states\n";
		auto tic = watch::now();
		fill(begin(d), end(d), 0.);
		for(std::size_t b = 0; b != nbasis; ++b)
			transform(std::execution::par, begin(d), end(d), begin(vbs[b]), begin(d), [](auto const& a, auto const& b){return a + norm(b);});
		cout<<"\t\t"<< ns{watch::now()-tic}.count()/1e9 <<" seconds  \t"<<d[d.size()/2]<<std::endl;
	}
#endif
}

#if 0	


	{
		multi::array<double, 1> arr(100, 99.); assert(size(arr) == 100);
		assert( begin(arr) < end(arr) );

	}
	{
		multi::array<double, 2> arr({100, 100}, 99.); assert(size(arr) == 100);
		assert( cbegin(arr) < cend(arr) );
	}
	{
		std::vector<double> v(10000);
		multi::array_ref<double, 2> A(v.data(), {100, 100}); assert(size(A) == 100);
		begin(A)[4][3] = 2.; // ok 
		using multi::static_array_cast;
	//	auto const& A_const = static_array_cast<double const>(A);
	//	begin(A_const)[4][3] = 2.; // error, read only
	}
	{
		std::vector<double> dd(10000);
		multi::array_ref<double, 2, std::vector<double>::iterator> arr(begin(dd), {100, 100}); assert(size(arr) == 100);
		begin(arr)[4][3] = 2.;
	//	assert( cbegin(arr)/2 );
	//	assert( cbegin(arr) < cend(arr) );
	}
	return 0;

	multi::array<double, 3>::reverse_iterator rit;
	assert(( rit.base() == multi::array<double, 3>::reverse_iterator{}.base() ));
	assert(( multi::array<double, 3>::reverse_iterator{}.base() == multi::array<double, 3>::reverse_iterator{}.base() ));
	assert(( multi::array<double, 3>::reverse_iterator{} == multi::array<double, 3>::reverse_iterator{} ));
	assert(( multi::array<double, 3>::reverse_iterator{} == multi::array<double, 3>::reverse_iterator{} ));

	multi::array<double, 3> A =
	#if defined(__INTEL_COMPILER)
		(double[3][2][2])
	#endif
		{
			{{ 1.2,  1.1}, { 2.4, 1.}},
			{{11.2,  3.0}, {34.4, 4.}},
			{{ 1.2,  1.1}, { 2.4, 1.}}
		}
	;
	assert( begin(A) < end(A) );
	assert( cbegin(A) < cend(A) );
	assert( begin(A[0]) < end(A[0]) );

	multi::array<double, 1>::const_iterator i;
	assert( begin(A[0]) < end(A[0]) );

	assert( size(A) == 3 and size(A[0]) == 2 and size(A[0][0]) == 2 and A[0][0][1] == 1.1 );
	assert(( multi::array<double, 3>::reverse_iterator{A.begin()} == rend(A) ));

	assert( begin(A) < end(A) );
	assert( cbegin(A) < cend(A) );
//	assert( crbegin(A) < crend(A) );
//	assert( crend(A) > crbegin(A) );
	assert( end(A) - begin(A) == size(A) );
	assert( rend(A) - rbegin(A) == size(A) );

	assert( size(*begin(A)) == 2 );
	assert( size(begin(A)[1]) == 2 );
	assert( &(A[1][1].begin()[0]) == &A[1][1][0] );
	assert( &A[0][1][0] == &A[0][1][0] );
	assert( &((*A.begin())[1][0]) == &A[0][1][0] );
	assert( &((*A.begin()).operator[](1)[0]) == &A[0][1][0] );
	assert( &(A.begin()->operator[](1)[0]) == &A[0][1][0] );
	assert( &(A.begin()->operator[](1).begin()[0]) == &A[0][1][0] );
	assert( &((A.begin()+1)->operator[](1).begin()[0]) == &A[1][1][0] );
	assert( &((begin(A)+1)->operator[](1).begin()[0]) == &A[1][1][0] );
	assert( &((cbegin(A)+1)->operator[](1).begin()[0]) == &A[1][1][0] );

	multi::array<double, 3>::iterator it; assert(( it == multi::array<double, 3>::iterator{} ));
	--it;
	it = begin(A);                                    assert( it == begin(A) );
	multi::array<double, 3>::iterator it2 = begin(A); assert(it == it2);
	it = end(A);                                      assert(it != it2);
	assert(it > it2);
	multi::array<double, 3>::iterator it3{it};        assert( it3 == it );
	multi::array<double, 3>::const_iterator cit;
	cit = it3;                                        assert( cit == it3 );
	assert((begin(A) == multi::array<double, 3>::iterator{rend(A)}));
	{
		std::vector<double> vv = {1.,2.,3.};
		auto it = vv.begin();
		auto rit = vv.rend();
		assert(std::vector<double>::reverse_iterator{it} == rit);
	}
#endif

}

