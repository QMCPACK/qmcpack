#ifdef COMPILATION_INSTRUCTIONS
clang++ -O3 -Ofast -std=c++14 -DNDEBUG -Wall -Wextra -Wpedantic -Werror $0 -o $0.x && time $0.x $@ && rm -f $0.x; exit
#endif

#include "../array_ref.hpp"

#include<chrono>
#include<iostream>
#include<numeric>
#include<string>
#include<utility>  // move
#include<vector>

using std::cout; using std::cerr;
namespace multi = boost::multi;

namespace {
// minimal RAII wall-clock timer (replaces boost::timer::auto_cpu_timer)
class auto_timer {
	std::string                           label_;
	std::chrono::steady_clock::time_point start_ = std::chrono::steady_clock::now();

 public:
	explicit auto_timer(std::string label = {}) : label_{std::move(label)} {}
	auto_timer(auto_timer const&)                    = delete;
	auto_timer(auto_timer&&)                         = delete;
	auto operator=(auto_timer const&) -> auto_timer& = delete;
	auto operator=(auto_timer&&) -> auto_timer&      = delete;
	~auto_timer() {
		std::cerr << label_ << std::chrono::duration<double>(std::chrono::steady_clock::now() - start_).count() << " s (wall)\n";
	}
};
}  // namespace

int main(){

	assert(0); // check that NDEBUG is off
{
	std::ptrdiff_t NX = 20000;
	std::ptrdiff_t NY = 20000;

	std::vector<double> data(NX*NY); 
	std::iota(begin(data), end(data), 0.);
	multi::array_cref<double, 2> data2D_cref({NX, NY}, data.data());

	cout << "size " << static_cast<double>(data.size()*sizeof(double))/1.0e6 << "MB\n";

	iota(begin(data), end(data), 1.2); iota(begin(data), end(data), 10.1); data[1234] = 929.1;
	double sum_raw;
	{
		double sum = 0.0;
		auto_timer t{"sum: raw: "};
		for(auto const& e : data) sum += e;
		sum_raw = sum;
	}
	iota(begin(data), end(data), 1.2); iota(begin(data), end(data), 202.2); data[1234] = 399.1;
	double sum_2D;
	{
		double sum = 0.0;
		auto_timer t{"sum: 2D: "};
		auto ext = extensions(data2D_cref);
		for(auto i : std::get<0>(ext)){
			for(auto j : std::get<1>(ext)) sum += data2D_cref(i, j);
		}
		sum_2D = sum;
	}
	iota(begin(data), end(data), 1.2); iota(begin(data), end(data), 2.21); data[1234] = 3299.1;
	double sum_2D_acc;
	{
		double sum = 0.0;
		auto_timer t{"sum: 2D acc: "};
		sum = std::accumulate(
			begin(data2D_cref), end(data2D_cref), 0.0,
			[](auto const& a, auto const& b){return a + std::accumulate(begin(b), end(b), 0.0);}
		);
		sum_2D_acc = sum;
	}
	iota(begin(data), end(data), 1.2); iota(begin(data), end(data), 2.15); data[1234] = 3419.1;
	double sum_2Dwrong_acc;
	{
		auto_timer t{"sum: 2Dwrong acc: "};
		sum_2Dwrong_acc = std::accumulate(
			begin(rotated(data2D_cref)), end(rotated(data2D_cref)), 0.0,
			[](auto const& a, auto const& b){return a + std::accumulate(begin(b), end(b), 0.0);}
		);
	}
	iota(begin(data), end(data), 1.2); iota(begin(data), end(data), 11.2); data[1234] = 199.1;
	double sum_2Dwrong;
	{
		double sum = 0.0;
		auto_timer t{"sum: 2Dwrong sum: "};
		auto const [is, js] = extensions(data2D_cref);
		for(auto j : js){
			for(auto i : is){
				sum += data2D_cref[i][j];
			}
		}
		sum_2Dwrong = sum;
	}
	iota(begin(data), end(data), 1.2); iota(begin(data), end(data), 10.112); data[1234] = 99.1;
	double sum_raw2;
	{
		double sum = 0.0;
		auto_timer t{"sum: raw: "};
		for(auto const& e : data) sum += e;
		sum_raw2 = sum;
	}
	iota(begin(data), end(data), 1.2); iota(begin(data), end(data), 33.112); data[1234] = 1199.1;
	double sum_raw_acc;
	{
		auto_timer t{"sum: raw acc: "};
		sum_raw_acc = std::accumulate(data.data(), data.data() + data.size(), 0.);
	}
	cout<< sum_2D + sum_2D_acc + sum_2Dwrong_acc + sum_2Dwrong + sum_raw + sum_raw2 + sum_raw_acc <<std::endl;
}
cout<<'\n';
{
	std::ptrdiff_t NX = 700;
	std::ptrdiff_t NY = 700;
	std::ptrdiff_t NZ = 700;

	std::vector<double> v(NX*NY*NZ);
	cout<<"3D data "<< static_cast<double>(v.size()*sizeof(double))/1.0e6 <<"MB\n";
	iota(begin(v), end(v), 0.1);

	multi::array_cref<double, 3> v3D_cref({NX, NY, NZ}, v.data());
	assert( num_elements(v3D_cref) == std::ptrdiff_t(v.size()) );
	{
		double sum = 0.0;
		auto_timer t{"sum: 3D raw: "};
		for(auto const& e : v) sum += e;
		cout << sum << '\n';
	}
	iota(begin(v), end(v), 1.2);
	{
		double sum = 0.0;
		auto_timer t{"sum: 3D indexed: "};
		auto const [is, js, ks] = extensions(v3D_cref);
		for(auto i : is) {
			auto const& v3D_crefi = v3D_cref[i];
			for(auto j : js) {
				auto const& v3D_crefij = v3D_crefi[j];
				for(auto k : ks) {
					sum += v3D_crefij[k];
				}
			}
		}
		cout << sum << '\n';
	}
	iota(begin(v), end(v), 4444.5);
	{
		double sum = 0.0;
		auto_timer t{"sum: 3Dwrong indexed: "};
		for(auto k : v3D_cref.extension(2))  // TODO(correaa) this doesn't work anymore
			for(auto j : v3D_cref.extension(1))
				for(auto i : v3D_cref.extension(0))
					sum += v3D_cref[i][j][k];
		cout << sum << '\n';
	}
}
{
	std::ptrdiff_t NX = 150;
	std::ptrdiff_t NY = 150;
	std::ptrdiff_t NZ = 150;
	std::ptrdiff_t NA = 150;

	std::vector<double> v(NX*NY*NZ*NA);
	multi::array_cref<double, 4> v4D_cref({NX, NY, NZ, NA}, v.data());
	assert( v4D_cref.num_elements() == std::ptrdiff_t(v.size()) );
	cout<<"4D data "<< num_elements(v4D_cref)*sizeof(double)/1e6 <<"MB\n";
	iota(begin(v), end(v), 0.1);
	{
		double sum = 0.0;
		auto_timer t{"sum: 4D raw: "};
		for(auto const& e : v) sum += e;
		cout<< sum <<'\n';
	}
	iota(begin(v), end(v), 1222.1);
	{
		double sum = 0.0;
		auto_timer t{"sum: 4D indexed: "};
		auto const [is, js, ks, ls] = extensions(v4D_cref);
		for(auto i : is) {
		//	auto const& v4D_crefi = v4D_cref[i]; // not necessary in clang or gcc
			for(auto j : js) {
			//	auto const& v4D_crefij = v4D_crefi[j]; // not necessary in clang or gcc
				for(auto k : ks) {
				//	auto const& v4D_crefijk = v4D_crefij[k]; // not necessary in clang or gcc
					for(auto l : ls) {
						sum += v4D_cref[i][j][k][l];
					}
				}
			}
		}
		cout<< std::to_string(sum)[0] <<'\n';
	}
}
}
