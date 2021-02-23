#if COMPILATION// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4-*-
OMPI_CXX=$CXXX OMPI_CXXFLAGS=$CXXFLAGS mpic++  $0 -o $0x&&mpirun -n 6 --oversubscribe $0x;exit
#endif
// © Alfredo A. Correa 2018-2020

#ifndef BOOST_MPI3_CARTESIAN_COMMUNICATOR_HPP
#define BOOST_MPI3_CARTESIAN_COMMUNICATOR_HPP

#include "../mpi3/communicator.hpp"
#include "../mpi3/process.hpp"

#include "../mpi3/detail/call.hpp"

namespace boost{
namespace mpi3{

using dimensionality_type = int;
static constexpr dimensionality_type dynamic_extent = -1;

template<dimensionality_type D = dynamic_extent> struct cartesian_communicator;

template<>
struct cartesian_communicator<dynamic_extent> : communicator{
	private:
	cartesian_communicator() : communicator(){}
	public:
	template<class Shape, class Period>
	cartesian_communicator(communicator& comm_old, Shape const& s, Period const& p){
		assert(s.size() == p.size());
		MPI_(Cart_create)(comm_old.get(), s.size(), s.data(), p.data(), false, &impl_);
	//	assert(impl_ != MPI_COMM_NULL); // null communicator is a valid outcome
		// TODO try with mpich, WAS: there is an bug in mpich, in which if the remaining dim are none then the communicator is not well defined.
	}
	template<class Shape>
	cartesian_communicator(communicator& comm_old, Shape const& s) : cartesian_communicator(comm_old, s, std::vector<int>(s.size(), true)){}
	
	cartesian_communicator(communicator& comm_old, std::initializer_list<int> shape) 
		: cartesian_communicator(comm_old, std::vector<int>(shape)){}
	cartesian_communicator(communicator& comm_old, std::initializer_list<int> shape, std::initializer_list<int> period) 
		: cartesian_communicator(comm_old, std::vector<int>(shape), std::vector<int>(period)){}

	[[deprecated("use dimensionality() instead of dimension")]] 
	int dimension() const{int ret; MPI_Cartdim_get(impl_, &ret); return ret;}

	int dimensionality() const{int ret; MPI_(Cartdim_get)(impl_, &ret); return ret;}
	std::vector<int> coordinates() const{
		std::vector<int> ret(dimensionality());
		MPI_(Cart_coords)(impl_, rank(), dimensionality(), ret.data());
		return ret;
	}
	auto topology() const{
		auto maxdims = dimensionality();
		struct topology_t{
			std::vector<int> dimensions;
			std::vector<int> periods;
			std::vector<int> coordinates;
			topology_t(std::size_t n) : dimensions(n), periods(n), coordinates(n){}
		} ret(maxdims);
		MPI_(Cart_get)(impl_, maxdims, ret.dimensions.data(), ret.periods.data(), ret.coordinates.data());
		assert( ret.coordinates == coordinates() );
		return ret;
	}
	std::vector<int> shape() const{return topology().dimensions;}
	std::vector<bool> periods() const{auto ps = topology().periods; return {ps.begin(), ps.end()};}
	auto num_elements() const{return size();}

	template<class Coord>
	auto operator()(Coord const& coord){
		int rank = -1;
		MPI_(Cart_rank)(impl_, coord.data(), &rank);
		return (*this)[rank];
	//	return operator[](rank);
	}
	// int MPI_Cart_map not implemented
	cartesian_communicator sub(std::vector<int> const& remain_dims) const{
		assert( static_cast<dimensionality_type>(remain_dims.size()) == dimensionality() );
		cartesian_communicator ret; MPI_(Cart_sub)(impl_, remain_dims.data(), &ret.impl_); return ret;
	}
	template<class RemainDim = std::initializer_list<bool>>
	cartesian_communicator sub(RemainDim const& remain_dims) const{
		return sub(std::vector<int>(remain_dims.begin(), remain_dims.end()));
	}
	cartesian_communicator sub() const{
		assert( dimensionality()>1 );
		std::vector<int> remain(dimensionality(), true); remain[0] = false;
		return sub(remain);
	}
};

enum fill_t{fill = 0};

template<dimensionality_type D>
struct cartesian_communicator : cartesian_communicator<>{
	static std::array<int, D> division(int nnodes, std::array<int, D> suggest = {}){
		return MPI_(Dims_create)(nnodes, D, suggest.data()), suggest;
	}
	constexpr static dimensionality_type dimensionality = D;
	cartesian_communicator(communicator& other, std::array<int, D> dims) try: 
		cartesian_communicator<>(other, division(other.size(), dims))
	{}catch(std::runtime_error& e){
		std::ostringstream ss;
		std::copy(dims.begin(), dims.end(), std::ostream_iterator<int>{ss, " "});
		throw std::runtime_error{"cannot create cartesian communicator with constrains "+ss.str()+" from communicator of size "+std::to_string(other.size())+" because "+e.what()};
	}
	auto topology() const{
		struct topology_t{
			std::array<int, dimensionality> dimensions, periods, coordinates;
		} ret;
		MPI_(Cart_get)(
			impl_, dimensionality, 
			ret.dimensions.data(), ret.periods.data(), ret.coordinates.data()
		);
		return ret;
	}
	auto dimensions() const{return topology().dimensions;}
	cartesian_communicator<D-1> sub() const{
		static_assert( D != 1 , "!");
		auto comm_sub = cartesian_communicator<>::sub();
		return static_cast<cartesian_communicator<D-1>&>(comm_sub);
//		return cartesian_communicator<D-1>(comm_sub, comm_sub.shape());
	}
	cartesian_communicator<1> axis(int d) const{
		std::vector<int> remains(D, false); remains[d] = true;
		auto comm_sub = cartesian_communicator<>::sub(remains);
		return static_cast<cartesian_communicator<1>&>(comm_sub);
//		return cartesian_communicator<1>(comm_sub, {comm_sub.shape()[d]});				
	}
};

#ifdef __cpp_deduction_guides
template<class T> cartesian_communicator(T, std::initializer_list<int>, std::initializer_list<bool>)
	->cartesian_communicator<dynamic_extent>;
template<class T> cartesian_communicator(T, std::initializer_list<int>, std::initializer_list<int>)
	->cartesian_communicator<dynamic_extent>;
template<class T> cartesian_communicator(T, std::initializer_list<int>)
	->cartesian_communicator<dynamic_extent>;
template<class... As> cartesian_communicator(As...)
	->cartesian_communicator<dynamic_extent>;
#endif

}}

#if not __INCLUDE_LEVEL__ // def _TEST_BOOST_MPI3_CARTESIAN_COMMUNICATOR

#include<iostream>

#include "../mpi3/main.hpp"
#include "../mpi3/version.hpp"
#include "../mpi3/ostream.hpp"

using std::cout;
using std::cerr;
namespace mpi3 = boost::mpi3;

int mpi3::main(int, char*[], boost::mpi3::communicator world){
{
	auto div = mpi3::cartesian_communicator<2>::division(6);
	assert( div[0]*div[1] == 6 );
	assert( div[0] == 3 );
	assert( div[1] == 2 );
}
{
	auto div = mpi3::cartesian_communicator<2>::division(6, {});
	assert( div[0]*div[1] == 6 );
	assert( div[0] == 3 );
	assert( div[1] == 2 );
}
{
	auto div = mpi3::cartesian_communicator<2>::division(6, {mpi3::fill});
	assert( div[0]*div[1] == 6 );
	assert( div[0] == 3 );
	assert( div[1] == 2 );
}
{
	auto div = mpi3::cartesian_communicator<2>::division(6, {mpi3::fill, mpi3::fill});
	assert( div[0]*div[1] == 6 );
	assert( div[0] == 3 );
	assert( div[1] == 2 );
}
{
	assert(world.size() == 6);
	auto div = mpi3::cartesian_communicator<2>::division(6, {2});
	assert( div[0]*div[1] == 6 );
	assert( div[0] == 2 );
	assert( div[1] == 3 );
}
{
	auto div = mpi3::cartesian_communicator<2>::division(6, {2, mpi3::fill});
	assert( div[0]*div[1] == 6 );
	assert( div[0] == 2 );
	assert( div[1] == 3 );
}
{
	auto div = mpi3::cartesian_communicator<2>::division(6, {mpi3::fill, 3});
	assert( div[0]*div[1] == 6 );
	assert( div[0] == 2 );
	assert( div[1] == 3 );
}
{
	auto div = mpi3::cartesian_communicator<2>::division(7);
	assert( div[0]*div[1] == 7 );
	assert( div[0] == 7 );
	assert( div[1] == 1 );
}
{
	auto div = mpi3::cartesian_communicator<2>::division(7);
	assert( div[0]*div[1] == 7 );
	assert( div[0] == 7 );
	assert( div[1] == 1 );
}
{
	auto div = mpi3::cartesian_communicator<2>::division(7, {mpi3::fill, mpi3::fill});
	assert( div[0]*div[1] == 7 );
	assert( div[0] == 7 );
	assert( div[1] == 1 );
}
{
	assert(world.size() == 6);
	mpi3::cartesian_communicator<2> cart_comm(world, {2, 3});
	assert( cart_comm );
	assert( cart_comm.dimensions()[0] == 2 );
	assert( cart_comm.dimensions()[1] == 3 );
	auto row = cart_comm.axis(0);
	auto col = cart_comm.axis(1);
	assert( row.size() == 2 );
	assert( col.size() == 3 );
}
{
	assert(world.size() == 6);
	if(mpi3::cartesian_communicator<2> cart_comm{world, {2, 2}}){
		auto row = cart_comm.axis(0);
		auto col = cart_comm.axis(1);
	}
}
try{
	assert(world.size() == 6);
	mpi3::cartesian_communicator<2> cart_comm(world, {4});
	assert(cart_comm.dimensions()[0] == 2);
	assert(cart_comm.dimensions()[1] == 3);
}catch(...){}
{
	mpi3::cartesian_communicator<2> cart_comm(world, {2, mpi3::fill});
	assert(cart_comm.dimensions()[0] == 2);
	assert(cart_comm.dimensions()[1] == 3);
}
{
	mpi3::cartesian_communicator<2> cart_comm(world, {mpi3::fill, 2});
	assert(cart_comm.dimensions()[0] == 3);
	assert(cart_comm.dimensions()[1] == 2);
}
{
	return 0;

	mpi3::cartesian_communicator<> comm(world, {4, 3}, {true, false});
	assert( comm.dimensionality() == 2 );
	cerr <<"= I am rank "<< comm.rank() <<" and have coordinates "<< comm.coordinates()[0] <<", "<< comm.coordinates()[1] <<"\n";
	auto comm_sub = comm.sub();
	assert( comm_sub.dimensionality() == 1 );
}
	if(world.root()) cerr<<"---"<<std::endl;
#ifdef __cpp_deduction_guides
{
	assert(world.size() == 12);

	mpi3::cartesian_communicator comm(world, {4, 3}, {true, false});
	assert( comm.dimensionality() == 2 );
	cerr <<"- I am rank "<< comm.rank() <<" and have coordinates "<< comm.coordinates()[0] <<", "<< comm.coordinates()[1] <<"\n";
}
#endif
	if(world.root()) cerr<<"---"<<std::endl;
#ifdef __cpp_deduction_guides
{
	assert(world.size() == 12);

	mpi3::cartesian_communicator comm(world, {4, 3}, {true, false});
	assert( comm.dimensionality() == 2 );
	cerr <<"I am rank "<< comm.rank() <<" and have coordinates "<< comm.coordinates()[0] <<", "<< comm.coordinates()[1] <<"\n";
}
#endif
{
	assert(world.size() == 12);

	mpi3::cartesian_communicator<3> comm(world, {});
	static_assert( mpi3::cartesian_communicator<3>::dimensionality == 3, "!");
	assert( comm.cartesian_communicator<>::dimensionality() == 3 );
	assert( comm.num_elements() == world.size() );
	assert( comm.shape()[0] == 3 );
	assert( comm.shape()[1] == 2 );
	assert( comm.shape()[2] == 2 );
	cerr<<"+ I am rank "<< comm.rank() <<" and have coordinates "<< comm.coordinates()[0] <<", "<< comm.coordinates()[1] <<", "<< comm.coordinates()[2] <<'\n';

	auto comm_sub = comm.sub();
	static_assert( comm_sub.dimensionality == 2 , "!" );
	std::cout << "numelements " << comm_sub.num_elements() << std::endl;
	assert( comm_sub.num_elements() == 4 );

	assert( comm_sub.shape()[0] == 2 );
	assert( comm_sub.shape()[1] == 2 );
	{
		auto comm_sub0 = comm.axis(0);
		assert( comm_sub0.shape()[0] == 3 );
		assert( comm_sub0.size() == 3 );
	}
	{
		auto comm_sub1 = comm.axis(1);
		assert( comm_sub1.shape()[0] == 2 );
		assert( comm_sub1.size() == 2 );
	}
	{
		auto comm_sub2 = comm.axis(2);
		assert( comm_sub2.shape()[0] == 2 );
		assert( comm_sub2.size() == 2 );
	}

}
	return 0;
}

#endif
#endif

