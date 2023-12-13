#if COMPILATION_INSTRUCTIONS
OMPI_CXX=$CXX mpicxx $CXXFLAGS -I$HOME/prj/alf/boost $0 -o $0x -lboost_serialization&&mpirun -n 2 $0x&&rm $0x;exit
#endif

#include<complex>

#include "mpi3/main_environment.hpp"
#include "mpi3/process.hpp"

#include "../array.hpp"

namespace mpi3 = boost::mpi3;
namespace multi = boost::multi;

void test_1D(mpi3::communicator& comm){

	switch(comm.rank()){
		case 0:{
			multi::array<double, 1> v(100);
			std::iota(v.begin(), v.end(), 0.);
			assert( v.strided(2).size() == 50 and v.strided(2)[9] == 18 );
			comm.send(v.strided(2).begin(), v.strided(2).end(), 1);
			return;
		}
		case 1:{
			multi::array<double, 1> v(50);
			comm.receive(v.begin(), v.end(), 0);
			assert( v[9] == 18 );
			return;
		}
	}
	assert(0);

}

void test_2D(mpi3::communicator& comm){

	auto const v = []{
		multi::array<double, 2> v({4, 5});
		std::iota(begin(v.elements()), end(v.elements()), 0.);
		return v;
	}();

	auto&& vpart = v({1, 4}, {1, 3});
	switch(comm.rank()){
		case 0:
			comm.send(begin(vpart), end(vpart), 1);
			return;
		case 1:
			multi::array<double, 2> w(extensions(vpart), 99.);
			comm.receive(begin(w), end(w), 0); 
			assert( w == vpart );
			return;
	}
	assert(0);

}

void test_2D_complex(mpi3::communicator& comm){

	using complex = std::complex<double>;
	
	multi::array<complex, 2> v({4, 5}); using std::get;
	if(auto x = v.extensions())
		for(auto i: get<0>(x))
			for(auto j: get<1>(x))
				v[i][j] = complex(i, j);

	switch(comm.rank()){
		case 0:
			comm.send(begin(v), end(v), 1); break;
		case 1:
			multi::array<complex, 2> w(extensions(v), 99.);
			comm.receive(begin(w), end(w), 0); 
			assert( w[2][3] == std::complex<double>(2., 3.) );
			break;
	}

}

void test_3D(mpi3::communicator& comm){

	auto const v = []{
		multi::array<double, 3> v({4, 5, 7});
		std::iota(begin(v.elements()), end(v.elements()), 0.);
		return v;
	}();

	auto&& vpart = v({1, 4}, {1, 3}, {3, 6});
	switch(comm.rank()){
		case 0:
			comm.send(begin(vpart), end(vpart), 1);
			return;
		case 1:
			multi::array<double, 3> w(extensions(vpart), 99.);
			comm.receive(begin(w), end(w), 0); assert( w == vpart );
			return;
	}
	assert(0);

}

void test_2D_strides(mpi3::communicator& comm){

	multi::array<double, 2> v({4, 5});
	std::iota(v.elements().begin(), v.elements().end(), 0.);

	std::cout << std::endl;
	switch(comm.rank()){
		case 0:{
			comm.send_n(v({1, 3}, {2, 4}).begin(), v({1, 3}, {2, 4}).size(), 1);
			return;
		}
		case 1:{
			multi::array<double, 2> w({4, 5}, 99.);
			comm.receive(w({1, 3}, {2, 4}).begin(), w({1, 3}, {2, 4}).end(), 0);
			assert( w({1, 3}, {2, 4}) == v({1, 3}, {2, 4}) );
			return;
		}
	}
	assert(0);

}

/*
void test_3D(mpi3::communicator& comm){

	switch(comm.rank()){
		case 0:{
			multi::array<double, 3> v({4, 5, 6});
			std::iota(v.elements().begin(), v.elements().end(), 0.);
			comm.send_n(v.begin(), v.size(), 1);
			return;
		}
		case 1:{
			multi::array<double, 3> v({4, 5, 6});
			comm.receive_n(v.begin(), v.size(), 0);
			for(auto const& e : v.elements()) std::cout<< e <<',';
			std::cout<<std::endl;
		//	for(auto&& row: v){for(auto&& e: row) std::cout<< e <<','; std::cout<<std::endl;}
			return;
		}
	}
	assert(0);

}
*/

void test_vector_nonpod(mpi3::communicator& comm){

	switch(comm.rank()){
		case 0:{
			std::vector<std::string> v(10);
			v[2] = "hola";
			comm.send_n(v.begin(), v.size(), 1);
			return;
		}
		case 1:{
			std::vector<std::string> v(10);
			comm.receive_n(v.begin(), v.size(), 0);
			assert( v[2] == "hola" );
			return;
		}
	}
	assert(0);

}

#if 0
void test_1D_nonpod(mpi3::communicator& comm){

	switch(comm.rank()){
		case 0:{
			multi::array<std::string, 1> v(10);
			v[2] = "hola";
			comm.send_n(v.begin(), v.size(), 1);
			return;
		}
		case 1:{
			multi::array<std::string, 1> v(10);
			comm.receive_n(v.begin(), v.size(), 0);
			assert( v[2] == "hola" );
			return;
		}
	}
	assert(0);

}
#endif

int mpi3::main(int, char*[], mpi3::environment& env){

	auto world = env.world();

	test_1D(world);
	test_2D_strides(world);
	test_2D(world);
	test_2D_complex(world);
	test_3D(world);
	
	{
		auto self = env.get_self_instance();
		auto const v = []{
			multi::array<double, 2> v({4, 5});
			std::iota(v.elements().begin(), v.elements().end(), 0.);
			return v;
		}();
		multi::array<double, 2> w({4, 2}, 0.);

		self.gather_n(v({0, 4}, {2, 4}).begin(), v({0, 4}, {2, 4}).size(), w.begin());

		using std::cout; using std::endl;
		std::cout<<"v()="<<std::endl;
		for(auto&& r: v){for(auto&& e: r) cout<< e <<'\t'; cout<<endl;}

		auto vv = v({0, 4}, {2, 4}).decay();
		std::cout<<"vv()="<<std::endl;
		for(auto&& r: vv){for(auto&& e: r) cout<< e <<'\t'; cout<<endl;}

		std::cout<<"w="<<std::endl;
		for(auto&& r : w){for(auto&& e: r) cout<< e <<'\t'; cout<<endl;}
		std::cout<<std::endl;
	}
	{
		auto self = env.get_self_instance();
		auto const v = []{
			multi::array<double, 2> v({4, 5});
			std::iota(v.elements().begin(), v.elements().end(), 0.);
			return v;
		}();
		multi::array<double, 2> w({5, 4}, 0.);
		self.gather_n(v.rotated().begin(), v.rotated().size(), w.begin());
		assert( v.rotated() == w );
	}

	multi::array<double, 1> v(100);
	{
		mpi3::type t(v.begin());
		assert( t.size() == sizeof(double) );
		assert( t.extent() == sizeof(double) );
	}
	{
		assert( v.strided(2).size() == 50 );
		mpi3::type t(v.strided(2).begin());
		assert( t.size() == sizeof(double) );
		assert( t.extent() == sizeof(double)*2 );
	}
	
	{
		auto self = env.get_self_instance();
		multi::array<double, 1> v(10);
		std::iota(v.elements().begin(), v.elements().end(), 0.);
		multi::array<double, 1> w(10);
		self.gather_n(v.begin(), v.size(), w.begin());
		assert( w == v );
	}
	{
		auto self = env.get_self_instance();
		auto const v = []{
			multi::array<double, 1> v(60);
			std::iota(v.elements().begin(), v.elements().end(), 0.);
			return v;
		}();
		multi::array<double, 1> w(30);
		assert( v.strided(2).size() == 30 );

		self.gather_n(v.strided(2).begin(), v.strided(2).size(), w.data_elements());
		assert( v.strided(2) == w );
	}
	{
		auto self = env.get_self_instance();
		auto const v = []{
			multi::array<double, 1> v(30);
			std::iota(v.elements().begin(), v.elements().end(), 0.);
			return v;
		}();
		multi::array<double, 1> w(10);
		self.gather_n(v({10, 30}).strided(2).begin(), v({10, 30}).strided(2).size(), w.data_elements());
		assert( v({10, 30}).strided(2) == w );
	}
	{
		auto self = env.get_self_instance();
		auto const v = []{
			multi::array<double, 2> v({4, 5});
			std::iota(v.elements().begin(), v.elements().end(), 0.);
			return v;
		}();
		multi::array<double, 2> w({2, 5}, 0.);

		self.gather_n(v.strided(2).begin(), v.strided(2).size(), w.begin());

		assert( v.strided(2) == w );		
	}

	{
		auto self = env.get_self_instance();
		auto const v = []{
			multi::array<double, 3> v({6, 4, 5});
			std::iota(v.elements().begin(), v.elements().end(), 0.);
			return v;
		}();
		multi::array<double, 3> w({3, 4, 5}, 0.);

		for(auto const& e : v.elements()) std::cout<< e <<',';
		std::cout << std::endl;

		self.gather_n(v.strided(2).begin(), v.strided(2).size(), w.begin());

		assert( v.strided(2) == w );		
	}
	return 0;
}

