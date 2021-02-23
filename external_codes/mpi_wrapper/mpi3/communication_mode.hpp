//  (C) Copyright Alfredo A. Correa 2018.
#ifndef MPI3_DETAIL_COMMUNICATION_MODE
#define MPI3_DETAIL_COMMUNICATION_MODE

#define OMPI_SKIP_MPICXX 1//https://github.com/open-mpi/ompi/issues/5157
#include <mpi.h> // if you get an error here probably you need to compile with and MPI compiler wrapper

#include<utility> // forward

namespace boost{
namespace mpi3{

struct blocking_mode{};
struct nonblocking_mode{};

struct standard_communication_mode{
	template<class... Args>
	int Send(Args&&... args) const{return MPI_Send(std::forward<Args>(args)...);}
	template<class... Args>
	int Recv(Args&&... args) const{return MPI_Recv(std::forward<Args>(args)...);}
	template<class... Args>
	int ISend(Args&&... args) const{return MPI_Isend(std::forward<Args>(args)...);}
	template<class... Args>
	int IRecv(Args&&... args) const{return MPI_Irecv(std::forward<Args>(args)...);}
};
struct buffered_communication_mode{
	template<class... Args>
	int Send(Args&&... args) const{return MPI_Bsend(std::forward<Args>(args)...);}
	template<class... Args>
	int Recv(Args&&... args) const{return MPI_Brecv(std::forward<Args>(args)...);}
	template<class... Args>
	int ISend(Args&&... args) const{return MPI_Ibsend(std::forward<Args>(args)...);}
	template<class... Args>
	int IRecv(Args&&... args) const{return MPI_Ibrecv(std::forward<Args>(args)...);}
};
struct synchronous_communication_mode{
	template<class... Args>
	int Send(Args&&... args) const{return MPI_Ssend(std::forward<Args>(args)...);}
	template<class... Args>
	int Recv(Args&&... args) const{return MPI_Srecv(std::forward<Args>(args)...);}
	template<class... Args>
	int ISend(Args&&... args) const{return MPI_Issend(std::forward<Args>(args)...);}
	template<class... Args>
	int IRecv(Args&&... args) const{return MPI_Isrecv(std::forward<Args>(args)...);}
};
struct ready_communication_mode{
	template<class... Args>
	int Send(Args&&... args) const{return MPI_Rsend(std::forward<Args>(args)...);}
	template<class... Args>
	int Recv(Args&&... args) const{return MPI_Rrecv(std::forward<Args>(args)...);}
	template<class... Args>
	int ISend(Args&&... args) const{return MPI_Irsend(std::forward<Args>(args)...);}
	template<class... Args>
	int IRecv(Args&&... args) const{return MPI_Irrecv(std::forward<Args>(args)...);}
};

struct gather_mode{
	template<class... Args>
	int operator()(Args&&... args) const{return MPI_Gather(std::forward<Args>(args)...);}
};
struct igather_mode{
	template<class... Args>
	int operator()(Args&&... args) const{return MPI_Igather(std::forward<Args>(args)...);}
};
struct all_gather_mode{
	template<class... Args>
	int operator()(Args&&... args) const{return MPI_Allgather(std::forward<Args>(args)...);}
};

struct reduce_mode{
	template<class... Args>
	int operator()(Args&&... args) const{return MPI_Reduce(std::forward<Args>(args)...);}
};
struct ireduce_mode{
	template<class... Args>
	int operator()(Args&&... args) const{return MPI_Ireduce(std::forward<Args>(args)...);}
};
struct all_reduce_mode{
	template<class... Args>
	int operator()(Args&&... args) const{return MPI_Allreduce(std::forward<Args>(args)...);}
};

}}
#endif


