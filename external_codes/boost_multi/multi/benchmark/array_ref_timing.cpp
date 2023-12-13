#ifdef COMPILATION_INSTRUCTIONS
clang++ -O3 -Ofast -std=c++14 -DNDEBUG -Wall -Wextra -Wpedantic -Werror $0 -o $0.x -lboost_timer && time $0.x $@ && rm -f $0.x; exit
#endif

#include "../array_ref.hpp"

#include <boost/timer/timer.hpp>

#include<iostream> 
#include<numeric>
#include<vector>
#include<numeric>

using std::cout; using std::cerr;
namespace multi = boost::multi;

int main(){

	assert(0); // check that NDEBUG is on
{
	std::ptrdiff_t NX = 20000;
	std::ptrdiff_t NY = 20000;

	std::vector<double> data(NX*NY); 
	std::iota(begin(data), end(data), 0.);
	multi::array_cref<double, 2> data2D_cref({NX, NY}, data.data());

	cout << "size " << data.size()*sizeof(double)/1e6 << "MB\n";

	iota(begin(data), end(data), 1.2); iota(begin(data), end(data), 10.1); data[1234] = 929.1;
	double sum_raw;
	{
		double sum = 0.;
		boost::timer::auto_cpu_timer t{std::cerr, 3, "sum: raw %t seconds\n"};
		for(auto const& e : data) sum += e;
		sum_raw = sum;
	}
	iota(begin(data), end(data), 1.2); iota(begin(data), end(data), 202.2); data[1234] = 399.1;
	double sum_2D;
	{
		double sum = 0.;
		boost::timer::auto_cpu_timer t{std::cerr, 3, "sum: 2D %t seconds\n"};
		auto ext = extensions(data2D_cref);
		for(auto i : std::get<0>(ext)){
		//	auto const& data2D_crefi = data2D_cref[i]; // not necessary in any clang or gcc
		//	for(auto j : std::get<1>(ext)) sum += data2D_crefi[j];
			for(auto j : std::get<1>(ext)) sum += data2D_cref(i, j);
		}
		sum_2D = sum;
	}
	iota(begin(data), end(data), 1.2); iota(begin(data), end(data), 2.21); data[1234] = 3299.1;
	double sum_2D_acc;
	{
		double sum = 0.;
		boost::timer::auto_cpu_timer t{cerr, 3, "sum: 2D acc %t seconds\n"};
		sum = std::accumulate(
			begin(data2D_cref), end(data2D_cref), 0.,
			[](auto&& a, auto&& b){return a + std::accumulate(begin(b), end(b), 0.);}
		);
		sum_2D_acc = sum;
	}
	iota(begin(data), end(data), 1.2); iota(begin(data), end(data), 2.15); data[1234] = 3419.1;
	double sum_2Dwrong_acc;
	{
		boost::timer::auto_cpu_timer t{cerr, 3, "sum: 2Dwrong acc %t seconds\n"};
		sum_2Dwrong_acc = std::accumulate(
			begin(rotated(data2D_cref)), end(rotated(data2D_cref)), 0.,
			[](auto&& a, auto&& b){return a + std::accumulate(begin(b), end(b), 0.);}
		);
	}
	iota(begin(data), end(data), 1.2); iota(begin(data), end(data), 11.2); data[1234] = 199.1;
	double sum_2Dwrong;
	{
		double sum = 0.;
		boost::timer::auto_cpu_timer t{std::cerr, 3, "sum: 2Dwrong sum %t seconds\n"};
		auto ext = extensions(data2D_cref);
		for(auto j : std::get<1>(ext)){
		//	auto const& data2D_crefXj = data2D_cref.rotated(1)[j]; // not necessary in clang or gcc
			for(auto i : std::get<0>(ext)){
			//	sum += data2D_crefXj[i];
				sum += data2D_cref[i][j];
			}
		}
		sum_2Dwrong = sum;
	}
	iota(begin(data), end(data), 1.2); iota(begin(data), end(data), 10.112); data[1234] = 99.1;
	double sum_raw2;
	{
		double sum = 0.;
		boost::timer::auto_cpu_timer t{cerr, 3, "sum: raw %t seconds\n"};
		for(auto const& e : data) sum += e;
		sum_raw2 = sum;
	}
	iota(begin(data), end(data), 1.2); iota(begin(data), end(data), 33.112); data[1234] = 1199.1;
	double sum_raw_acc;
	{
		boost::timer::auto_cpu_timer t{cerr, 3, "sum: raw acc %t seconds\n"};
		sum_raw_acc = std::accumulate(data.data(), data.data() + data.size(), 0.);
	}
	cout<< sum_2D + sum_2D_acc + sum_2Dwrong_acc + sum_2Dwrong + sum_raw + sum_raw2 + sum_raw_acc <<std::endl;
}
cout<<'\n';	
{	
	std::ptrdiff_t NX = 700, NY = 700, NZ = 700;
	std::vector<double> v(NX*NY*NZ);
	cout<<"3D data "<< v.size()*sizeof(double)/1e6 <<"MB\n";
	iota(begin(v), end(v), 0.1);

	multi::array_cref<double, 3> v3D_cref({NX, NY, NZ}, v.data());
	assert( num_elements(v3D_cref) == std::ptrdiff_t(v.size()) );
	{
		double sum = 0.;
		boost::timer::auto_cpu_timer t{std::cerr, 3, "sum: 3D raw %t seconds\n"};
		for(auto const& e : v) sum += e;
		cout << sum << '\n';
	}
	iota(begin(v), end(v), 1.2);
	{
		double sum = 0.;
		boost::timer::auto_cpu_timer t{std::cerr, 3, "sum: 3D indexed %t seconds\n"};
		auto ext = extensions(v3D_cref);
		for(auto i : std::get<0>(ext)){
			auto const& v3D_crefi = v3D_cref[i];
			for(auto j : std::get<1>(ext)){
				auto const& v3D_crefij = v3D_crefi[j];
				for(auto k : std::get<2>(ext))
					sum += v3D_crefij[k];
			}
		}
		cout << sum << '\n';
	}
	iota(begin(v), end(v), 4444.5);
	{
		double sum = 0.;
		boost::timer::auto_cpu_timer t{std::cerr, 3, "sum: 3Dwrong indexed %t seconds\n"};
		for(auto k : v3D_cref.extension(2))
			for(auto j : v3D_cref.extension(1))
				for(auto i : v3D_cref.extension(0))
					sum += v3D_cref[i][j][k];
		cout << sum << '\n';
	}
}
{	
	std::ptrdiff_t NX = 150, NY = 150, NZ = 150, NA = 150;
	std::vector<double> v(NX*NY*NZ*NA);
	multi::array_cref<double, 4> v4D_cref({NX, NY, NZ, NA}, v.data());
	assert( v4D_cref.num_elements() == std::ptrdiff_t(v.size()) );
	cout<<"4D data "<< num_elements(v4D_cref)*sizeof(double)/1e6 <<"MB\n";
	iota(begin(v), end(v), 0.1);
	{
		double sum = 0.;
		boost::timer::auto_cpu_timer t{std::cerr, 3, "sum: 4D raw %t seconds\n"};
		for(auto const& e : v) sum += e;
		cout<< sum <<'\n';
	}
	iota(begin(v), end(v), 1222.1);
	{
		double sum = 0.;
		boost::timer::auto_cpu_timer t{std::cerr, 3, "sum: 4D indexed %t seconds\n"};
		auto ext = extensions(v4D_cref);
		for(auto i : std::get<0>(ext)){
		//	auto const& v4D_crefi = v4D_cref[i]; // not necessary in clang or gcc
			for(auto j : std::get<1>(ext)){
			//	auto const& v4D_crefij = v4D_crefi[j]; // not necessary in clang or gcc
				for(auto k : std::get<2>(ext)){
				//	auto const& v4D_crefijk = v4D_crefij[k]; // not necessary in clang or gcc
					for(auto l : std::get<3>(ext))
						sum += v4D_cref[i][j][k][l];
				}
			}
		}
		cout<< std::to_string(sum)[0] <<'\n';
	}
}
	return 0;
}
