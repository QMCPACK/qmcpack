#if COMPILATION_INSTRUCTIONS
(echo "#include\""$0"\"" > $0x.cpp) && mpic++ -O3 -std=c++17 -Wfatal-errors -D_TEST_BOOST_MPI3_CARTESIAN_COMMUNICATOR $0x.cpp -o $0x.x && time mpirun -n 12 $0x.x $@ && rm -f $0x.x $0x.cpp; exit
#endif
#ifndef BOOST_MPI3_CARTESIAN_COMMUNICATOR_HPP
#define BOOST_MPI3_CARTESIAN_COMMUNICATOR_HPP

#include "../mpi3/communicator.hpp"
#include "../mpi3/process.hpp"

namespace boost{
namespace mpi3{

struct cartesian_communicator : communicator{
	private:
	cartesian_communicator() : communicator(){}
	public:
	template<class Shape, class Period>
	cartesian_communicator(communicator& comm_old, Shape const& s, Period const& p){
		assert(s.size() == p.size());
		int status = MPI_Cart_create(&comm_old, s.size(), s.data(), p.data(), false, &impl_);
		if(status != MPI_SUCCESS) throw std::runtime_error("cannot create cart comm ");
		assert(impl_ != MPI_COMM_NULL);
		// there is an bug in mpich, in which the the remaining dim are none then the communicator is not well defined.
	}
	template<class Shape>
	cartesian_communicator(communicator& comm_old, Shape const& s) : cartesian_communicator(comm_old, s, std::vector<int>(s.size(), 0)){}

	cartesian_communicator(communicator& comm_old, std::initializer_list<int> shape) 
		: cartesian_communicator(comm_old, std::vector<int>(shape)){}
	cartesian_communicator(communicator& comm_old, std::initializer_list<int> shape, std::initializer_list<int> period) 
		: cartesian_communicator(comm_old, std::vector<int>(shape), std::vector<int>(period)){}

	int dimension() const{
		int ret;
		MPI_Cartdim_get(impl_, &ret);
		return ret;
	}
	std::vector<int> coordinates(){
		std::vector<int> ret(dimension());
		MPI_Cart_coords(impl_, rank(), dimension(), ret.data());
		return ret;
	}
	template<class Coord>
	auto operator()(Coord const& coord){
		int rank = -1;
		MPI_Cart_rank(impl_, coord.data(), &rank);
		return (*this)[rank];
	//	return operator[](rank);
	}
	// int MPI_Cart_map not implemented
	template<class RemainDim>
	cartesian_communicator
	sub(RemainDim const& remain_dims) const{
		assert(remain_dims.size() == dimension());
		cartesian_communicator ret;
		MPI_Cart_sub(impl_, remain_dims.data(), &ret.impl_);
		return ret;
	}
//	cartesian_communicator
	cartesian_communicator sub() const{
		std::vector<int> remain(1, dimension());
		return sub(remain);
	}
};

}}

#ifdef _TEST_BOOST_MPI3_CARTESIAN_COMMUNICATOR

#include<iostream>
#include "../mpi3/main.hpp"
#include "../mpi3/version.hpp"
#include "../mpi3/ostream.hpp"

using std::cout;
namespace mpi3 = boost::mpi3;
int mpi3::main(int, char*[], boost::mpi3::communicator world){

	if(world.size() != 12) throw std::runtime_error("run with 12 procs!");

	mpi3::cartesian_communicator comm(world, {4, 3}, {1, 0});
	cout <<"I am rank "<< comm.rank() <<" and have coordinates "<< comm.coordinates()[0] <<", "<< comm.coordinates()[1] <<"\n";
	
	return 0;
}

#endif
#endif


