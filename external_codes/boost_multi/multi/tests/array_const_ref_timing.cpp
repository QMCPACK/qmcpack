#ifdef COMPILATION_INSTRUCTIONS
c++ -O3 -std=c++14 -DNDEBUG -Wall `#-fmax-errors=2` `#-Wfatal-errors` -I${HOME}/prj $0 -o $0.x -lboost_timer && time $0.x $@ && rm -f $0.x; exit
#endif

#include "../array_ref.hpp"

#include <boost/timer/timer.hpp>
#include<numeric>
#include<vector>

using std::cout; using std::cerr;
namespace multi = boost::multi;

int main(){

	std::ptrdiff_t NX = 20000;
	std::ptrdiff_t NY = 20000;
	cout << "size " << NX*NY*sizeof(double)/1e6 << "MB\n";

	std::vector<double> data(NX*NY);
	std::iota(begin(data), end(data), 0.);

	multi::array_cref<double, 2> data2D_cref{data.data(), {NX, NY}};

	{
		double sum = 0.;
		{
//			boost::timer::auto_cpu_timer t;
			for(auto const& e : data)
				sum += e;
		}
		cout << sum << '\n';
	}
	{
		double sum = 0.;
		{
//			boost::timer::auto_cpu_timer t;
			for(auto i : data2D_cref.extension(0)){
				auto const& data2D_crefi = data2D_cref[i];
				for(auto j : data2D_cref.extension(1))
					sum += i*j*data2D_crefi[j];
			}
		}
		cout << sum << '\n';
	}
	{
		double sum = 0.;
		{
//			boost::timer::auto_cpu_timer t;
			for(auto i : data2D_cref.extension(0))
				for(auto j : data2D_cref.extension(1))
					sum += i*j*data2D_cref[i][j];
		}
		cout << sum << '\n';
	}


	{	
		std::ptrdiff_t NX = 500;
		std::ptrdiff_t NY = 500;
		std::ptrdiff_t NZ = 500;
		std::vector<double> v(NX*NY*NZ);
		iota(begin(v), end(v), 0.);
		multi::array_cref<double, 3> v3D_cref{v.data(), {NX, NY, NZ}};
		assert( v3D_cref.num_elements() == std::ptrdiff_t(v.size()) );
//		cout << "size " << v.size()*sizeof(double)/1e6 << "MB\n";

		{
			double sum = 0.;
			{
//				boost::timer::auto_cpu_timer t;
				for(auto const& e : v)
					sum += e;
			}
			cout << sum << '\n';
		}
		{
			double sum = 0.;
			{
//				boost::timer::auto_cpu_timer t;
				for(auto i : v3D_cref.extension(0)){
					auto const& v3D_crefi = v3D_cref[i];
					for(auto j : v3D_cref.extension(1)){
						auto const& v3D_crefij = v3D_crefi[j];
						for(auto k : v3D_cref.extension(2))
							sum += i*j*k*v3D_crefij[k];
					}
				}
			}
			cout << sum << '\n';
		}
		{
			double sum = 0.;
			{
//				boost::timer::auto_cpu_timer t;
				for(auto i : v3D_cref.extension(0))
					for(auto j : v3D_cref.extension(1))
						for(auto k : v3D_cref.extension(2))
							sum += i*j*k*v3D_cref[i][j][k];
			}
			cout << sum << '\n';
		}
		{
			double sum = 0.;
			{
//				boost::timer::auto_cpu_timer t;
				for(auto const& e : v)
					sum += e;
			}
			cout << sum << '\n';
		}
	}
	{	
		std::ptrdiff_t NX = 100;
		std::ptrdiff_t NY = 100;
		std::ptrdiff_t NZ = 100;
		std::ptrdiff_t NA = 100;
		std::vector<double> v(NX*NY*NZ*NA);
		iota(begin(v), end(v), 0.1);
		multi::array_cref<double, 4> v4D_cref{v.data(), {NX, NY, NZ, NA}};
		assert( v4D_cref.num_elements() == std::ptrdiff_t(v.size()) );
		cout << "size " << v.size()*sizeof(double)/1e6 << "MB\n";
		{
			double sum = 0.;
			{
				boost::timer::auto_cpu_timer t;
				for(auto const& e : v)
					sum += e;
			}
			cout << sum << '\n';
		}
		{
			double sum = 0.;
			{
				boost::timer::auto_cpu_timer t;
				for(auto i : v4D_cref.extension(0)){
					auto const& v4D_crefi = v4D_cref[i];
					for(auto j : v4D_cref.extension(1)){
						auto const& v4D_crefij = v4D_crefi[j];
						for(auto k : v4D_cref.extension(2)){
							auto const& v4D_crefijk = v4D_crefij[k];
							for(auto l : v4D_cref.extension(3))
								sum += i*j*k*v4D_crefijk[l];
						}
					}
				}
			}
			cout << sum << '\n';
		}
		{
			double sum = 0.;
			{
				boost::timer::auto_cpu_timer t;
				for(auto i : v4D_cref.extension(0))
					for(auto j : v4D_cref.extension(1))
						for(auto k : v4D_cref.extension(2))
							for(auto l : v4D_cref.extension(3))
								sum += i*j*k*l*v4D_cref[i][j][k][l];
			}
			cout << sum << '\n';
		}
		{
			double sum = 0.;
			{
				boost::timer::auto_cpu_timer t;
				auto exts = v4D_cref.extensions();
				for(auto i : exts[0])
					for(auto j : exts[1])
						for(auto k : exts[2])
							for(auto l : exts[3])
								sum += i*j*k*l*v4D_cref[i][j][k][l];
			}
			cout << sum << '\n';
		}
		{
			double sum = 0.;
			{
				boost::timer::auto_cpu_timer t;
				for(auto const& e : v)
					sum += e;
			}
			cout << sum << '\n';
		}
	}

}

