// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;autowrap:nil;-*-
// Copyright 2022 Alfredo A. Correa

#ifndef MPI3_NCCL_UNIVERSAL_COMMUNICATOR_HPP_
#define MPI3_NCCL_UNIVERSAL_COMMUNICATOR_HPP_

#include "../../mpi3/communicator.hpp"
#include "../../mpi3/nccl/communicator.hpp"

#define BMPI3_DECLRET(ExpR) ->decltype(ExpR) {return ExpR;}  // NOLINT(cppcoreguidelines-macro-usage) saves a lot of typing
#define BMPI3_JUSTRET(ExpR)                  {return ExpR;}  // NOLINT(cppcoreguidelines-macro-usage) saves a lot of typing

namespace boost {
namespace mpi3 {
namespace nccl {

template<std::size_t I> struct priority : std::conditional_t<I==0, std::true_type, priority<I-1>> {};

struct universal_communicator : mpi3::communicator, private mpi3::nccl::communicator {
	template<class... As>
	universal_communicator(As&&... as)
	: mpi3::communicator{std::forward<As>(as)...}
	, mpi3::nccl::communicator{static_cast<mpi3::communicator&>(*this)} {}

	using mpi3::communicator::rank;
	using mpi3::communicator::size;

 private:
	template<class P>
	using comm_system = typename std::conditional<
		std::is_same_v<typename ::thrust::iterator_system<P>::type, thrust::cuda_cub::tag>,
		mpi3::nccl::communicator,
		mpi3::      communicator
	>::type;

 public:
	template<class... As>          auto all_reduce_n(priority<0>, As... as         ) BMPI3_JUSTRET(mpi3::communicator::all_reduce_n(               as...       ))
	template<class P, class... Rs> auto all_reduce_n(priority<1>, P first, Rs... rs) BMPI3_DECLRET(comm_system<P>    ::all_reduce_n(               first, rs...))
	template<class... As         > auto all_reduce_n(             As... as         ) {return this->all_reduce_n(priority<1>{}, as...       );}

	template<class... As>          auto broadcast_n(priority<0>, As... as         ) BMPI3_JUSTRET(mpi3::communicator::broadcast_n(               as...       ))
	template<class P, class... Rs> auto broadcast_n(priority<1>, P first, Rs... rs) BMPI3_DECLRET(comm_system<P>    ::broadcast_n(               first, rs...))
	template<class... As         > auto broadcast_n(             As... as         ) {return this->                    broadcast_n(priority<1>{}, as...       );}

	template<class... As>          auto send_n(priority<0>, As... as         ) BMPI3_JUSTRET(mpi3::communicator::send_n(               as...       ))
	template<class P, class... Rs> auto send_n(priority<1>, P first, Rs... rs) BMPI3_DECLRET(comm_system<P>    ::send_n(               first, rs...))
	template<class... As         > auto send_n(             As... as         ) {return this->                    send_n(priority<1>{}, as...       );}

	template<class... As>          auto receive_n(priority<0>, As... as         ) BMPI3_JUSTRET(mpi3::communicator::receive_n(               as...       ))
	template<class P, class... Rs> auto receive_n(priority<1>, P first, Rs... rs) BMPI3_DECLRET(comm_system<P>    ::receive_n(               first, rs...))
	template<class... As         > auto receive_n(             As... as         ) {return this->                    receive_n(priority<1>{}, as...       );}
};

}  // end namespace nccl
}  // end namespace mpi3
}  // end namespace boost
#endif  // MPI3_NCCL_UNIVERSAL_COMMUNICATOR_HPP_
