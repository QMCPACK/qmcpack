#if COMPILATION_INSTRUCTIONS
mpic++ -O3 -std=c++17 -Wall -Wextra `#-Wfatal-errors` $0 -o $0x.x && time mpirun -n 1 $0x.x $@ && rm -f $0x.x; exit
#endif

#include "../../mpi3/environment.hpp"
#include "../../mpi3/communicator.hpp"

#include<any>

namespace mpi3 = boost::mpi3;
using std::cout;

mpi3::environment e;
mpi3::communicator::keyval<int> k;
mpi3::communicator::keyval<std::map<std::string, mpi3::any>> k_map;

int main(){

	auto world = e.get_world_instance();

	world.set_attribute(k, 10);
	world.get_attribute(k) = 20;
	int const& v = world.get_attribute(k);
	assert(v == 20);

	auto w2 = world;
	assert(w2.get_attribute(k) == 20);

	world.delete_attribute(k);
	assert(not world.has_attribute(k));
	assert(w2.has_attribute(k));

//	assert(not world.has_attribute(k_map));
//	world.attribute(k_map);
//	assert(world.has_attribute(k_map));
//	world.get_attribute(k_map)["leader"] = std::string{"red"};
//	assert(mpi3::any_cast<std::string>(world.get_attribute(k_map)["leader"]) == "red");
//	world.attribute(mpi3::environment::named_attributes_key());
	world.attribute("lop");
//	world.attribute("color") = std::string("red");
//	std::string s = mpi3::any_cast<std::string>(world.attribute("color"));
//	assert( s == "red" );

//	assert(mpi3::any_cast<bool>(world.named_attributes()["leader"]));

	std::map<std::string, std::any> m = {{"node", 555}};
	m["color"] = 11;
//	int n = 
	std::any_cast<int&>(m["color"]) = 99;
	int n = std::any_cast<int>(m["color"]);
	assert( n == 99 );

	assert( std::any_cast<int>(m["node"]) == 555 );

	return 0;
}

