// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;autowrap:nil;-*-
// Copyright 2022- Alfredo A. Correa

#ifndef MPI3_NCCL_COMMUNICATOR_HPP_
#define MPI3_NCCL_COMMUNICATOR_HPP_

#include "../../mpi3/communicator.hpp"
#include "../../mpi3/nccl/detail/basic_datatype.hpp"
#include "../../mpi3/nccl/detail/basic_reduction.hpp"

#include <thrust/system/cuda/memory.h>
#include <thrust/complex.h>

#include <functional>  // for plus
#include <iostream>

#include <nccl.h>

namespace boost {
namespace mpi3 {
namespace nccl {

namespace detail {

template<class T> auto datatype(T const&) -> decltype(basic_datatype_t<T>::value) {return basic_datatype<T>;}
template<class T> auto reduction(T const&) -> decltype(basic_reduction<T>) {return basic_reduction<T>;}

}

class communicator {
	static auto get_unique_id() -> ncclUniqueId {
		ncclUniqueId nccl_id;
		auto r = ncclGetUniqueId(&nccl_id);
		switch(r) {
			case ncclSuccess: break;
			case ncclUnhandledCudaError: assert(0);
			case ncclSystemError: assert(0);
			case ncclInternalError: assert(0);
			case ncclInvalidArgument: assert(0);
			case ncclInvalidUsage: assert(0 && "likely \"Duplicate GPU detected\", for example if rank 0 and rank 1 both on CUDA device 1000");
			case ncclRemoteError: assert(0);
			case ncclNumResults: assert(0);
		}
		return nccl_id;
	}

 public:
	explicit communicator(mpi3::communicator& mpi) : impl_{nullptr} {
		if(mpi.empty()) {return;}
		ncclUniqueId nccl_id = mpi.root()?get_unique_id():ncclUniqueId{};
		mpi.broadcast_n(reinterpret_cast<char*>(&nccl_id), sizeof(ncclUniqueId));
		// TODO(correaa) may need mpi.barrier(); here
		{
			ncclResult_t r = ncclCommInitRank(&impl_, mpi.size(), nccl_id, mpi.rank());
			switch(r) {
				case ncclSuccess: break;
				case ncclUnhandledCudaError: assert(0);
				case ncclSystemError: assert(0);
				case ncclInternalError: assert(0);
				case ncclInvalidArgument: assert(0);
				case ncclInvalidUsage: assert(0 && "likely \"Duplicate GPU detected\", for example if rank 0 and rank 1 both on CUDA device 1000");
				case ncclRemoteError: assert(0);
				case ncclNumResults: assert(0);
			}
		}
	}

	communicator(communicator const&) = delete;
//  [[deprecated("experimental")]]
	communicator(communicator& other) : impl_{nullptr} {
		if(other.empty()) {return;}
		ncclUniqueId nccl_id = other.root()?get_unique_id():ncclUniqueId{};
		{
			ncclResult_t r = ncclBcast(&nccl_id, sizeof(ncclUniqueId), ncclChar, 0, other.impl_, NULL);
			switch(r) {
				case ncclSuccess: break;
				case ncclUnhandledCudaError: assert(0);
				case ncclSystemError: assert(0);
				case ncclInternalError: assert(0);
				case ncclInvalidArgument: assert(0);
				case ncclInvalidUsage: assert(0);
				case ncclRemoteError: assert(0);
				case ncclNumResults: assert(0);
			}
		}
		cudaStreamSynchronize(NULL);
		// TODO(correaa) may need mpi.barrier(); here
		{
			ncclResult_t r = ncclCommInitRank(&impl_, other.count(), nccl_id, other.rank());
			switch(r) {
				case ncclSuccess: break;
				case ncclUnhandledCudaError: assert(0);
				case ncclSystemError: assert(0);
				case ncclInternalError: assert(0);
				case ncclInvalidArgument: assert(0);
				case ncclInvalidUsage: assert(0 && "likely \"Duplicate GPU detected\", for example if rank 0 and rank 1 both on CUDA device 1000");
				case ncclRemoteError: assert(0);
				case ncclNumResults: assert(0);
			}
		}
	}
	communicator(communicator&& other) : impl_{std::exchange(other.impl_, nullptr)} {}
	// moved from communicators are left in a partially formed, since there is no assigmment it cannot be used

	auto duplicate() {return communicator{*this};}

	template<
		class Op = std::plus<>,
		class P1, class Size, class P2,
		typename = decltype(
			*thrust::raw_pointer_cast(P2{}) = Op{}(*thrust::raw_pointer_cast(P1{}), *thrust::raw_pointer_cast(P1{})),
			detail::datatype(*raw_pointer_cast(P1{}))
		)
	>
	auto all_reduce_n(P1 first, Size count, P2 dest, Op op = {}) {
		ncclResult_t r = ncclAllReduce(
			thrust::raw_pointer_cast(first), thrust::raw_pointer_cast(dest), count, 
			detail::datatype(*raw_pointer_cast(first)),
			detail::reduction(op), impl_, NULL
		);
		switch(r) {
			case ncclSuccess: break;
			case ncclUnhandledCudaError: assert(0);
			case ncclSystemError: assert(0);
			case ncclInternalError: assert(0);
			case ncclInvalidArgument: assert(0);
			case ncclInvalidUsage: assert(0);
			case ncclRemoteError: assert(0);
			case ncclNumResults: assert(0);
		}
		return dest + count;
	}

	template<class P, class Size, typename = decltype(detail::datatype(*raw_pointer_cast(P{})))>
	auto send_n(P first, Size n, int peer) {
		// ncclGroupStart();
		ncclResult_t r = ncclSend(thrust::raw_pointer_cast(first), n, detail::datatype(*raw_pointer_cast(first)), peer, impl_, NULL);
		switch(r) {
			case ncclSuccess: break;
			case ncclUnhandledCudaError: assert(0);
			case ncclSystemError: assert(0);
			case ncclInternalError: assert(0);
			case ncclInvalidArgument: assert(0);
			case ncclInvalidUsage: assert(0);
			case ncclRemoteError: assert(0);
			case ncclNumResults: assert(0);
		}
		// ncclGroupEnd();
		// cudaStreamSynchronize(NULL);
		return first + n;
	}
	template<class P, class Size, typename = decltype(detail::datatype(*raw_pointer_cast(P{})))>
	auto receive_n(P first, Size n, int peer) {
		// ncclGroupStart();
		ncclResult_t r = ncclRecv(thrust::raw_pointer_cast(first), n, detail::datatype(*raw_pointer_cast(first)), peer, impl_, NULL);
		switch(r) {
			case ncclSuccess: break;
			case ncclUnhandledCudaError: assert(0);
			case ncclSystemError: assert(0);
			case ncclInternalError: assert(0);
			case ncclInvalidArgument: assert(0);
			case ncclInvalidUsage: assert(0);
			case ncclRemoteError: assert(0);
			case ncclNumResults: assert(0);
		}
		// ncclGroupEnd();
		// cudaStreamSynchronize(NULL);
		return first + n;
	}
	template<
		class P, //typename = decltype(detail::datatype(*raw_pointer_cast(P{}))),
		class Size
	>
	P broadcast_n(P first, Size n, int root = 0) {
		// ncclGroupStart();
		using thrust::raw_pointer_cast;
		ncclResult_t r = ncclBcast(raw_pointer_cast(first), n, detail::datatype(*raw_pointer_cast(first)), root, impl_, NULL);
		switch(r) {
			case ncclSuccess: break;
			case ncclUnhandledCudaError: assert(0);
			case ncclSystemError: assert(0);
			case ncclInternalError: assert(0);
			case ncclInvalidArgument: assert(0);
			case ncclInvalidUsage: assert(0);
			case ncclRemoteError: assert(0);
			case ncclNumResults: assert(0);
		}
		// ncclGroupEnd();
		// cudaStreamSynchronize(NULL);
		return first + n;
	}

 private:
	template<class Complex, class Value> static constexpr bool is_numeric_complex_of =
		std::is_same_v<decltype(Complex{Complex{}.real(), Complex{}.imag()}), Complex>
		and sizeof(Complex) == sizeof(Complex{}.real()) + sizeof(Complex{}.imag())
	;

 public:
	template<class Complex> static constexpr bool is_numeric_complex = is_numeric_complex_of<Complex, typename Complex::value_type>;
	template<
		class P, class Size,
		class PT = std::pointer_traits<P>, class Complex = typename PT::element_type, 
			class Value = typename std::conditional<
				std::is_const_v<Complex>,
				std::add_const_t<typename Complex::value_type>,
			                     typename Complex::value_type
			>::type,
		std::enable_if_t<is_numeric_complex_of<std::decay_t<Complex>, std::decay_t<Value>>, int> =0
	>
	auto send_n(P first, Size n, int peer) {
		send_n(thrust::reinterpret_pointer_cast<typename PT::template rebind<Value>>(first), n*2, peer);
		return first + n;
	}
	template<
		class P, class Size,
		class PT = std::pointer_traits<P>, class Complex = typename PT::element_type,
			class Value = typename std::conditional<
				std::is_const_v<Complex>,
				std::add_const_t<typename Complex::value_type>,
			                     typename Complex::value_type
			>::type,
		std::enable_if_t<is_numeric_complex_of<std::decay_t<Complex>, std::decay_t<Value>>, int> =0
	>
	auto receive_n(P first, Size n, int peer) {
		receive_n(thrust::reinterpret_pointer_cast<typename PT::template rebind<Value>>(first), n*2, peer);
		return first + n;
	}

//	template<
//		class P1, class Size, class P2
////		, typename = decltype(
////			*thrust::raw_pointer_cast(P2{}) = std::plus<>{}(*thrust::raw_pointer_cast(P1{}), *thrust::raw_pointer_cast(P1{})),
////			detail::datatype(*raw_pointer_cast(P1{}))
////		)
//		, class PT1 = std::pointer_traits<P1>
////		, class Complex1 = typename PT1::element_type,
////			class Value1 = typename std::conditional<
////				std::is_const_v<Complex1>,
////				std::add_const_t<typename Complex1::value_type>,
////			                     typename Complex1::value_type
////			>::type
////		, std::enable_if_t<is_numeric_complex_of<std::decay_t<Complex1>, std::decay_t<Value1>>, int> =0,
//		, class PT2 = std::pointer_traits<P2>
////		, class Complex2 = typename PT2::element_type,
////			class Value2 = typename std::conditional<
////				std::is_const_v<Complex2>,
////				std::add_const_t<typename Complex2::value_type>,
////			                     typename Complex2::value_type
////			>::type,
////		std::enable_if_t<is_numeric_complex_of<std::decay_t<Complex2>, std::decay_t<Value2>>, int> =0
//	>
//	auto all_reduce_n(P1 first, Size count, P2 dest, std::plus<> op = {}) {
//		using Value1 = double;
//		using Value2 = double;
//		all_reduce_n(
//			reinterpret_pointer_cast<typename PT1::template rebind<Value1>>(first), count*2,
//			reinterpret_pointer_cast<typename PT2::template rebind<Value2>>(dest ), op
//		);
//		return dest + count;
//	}
	template<class Size>
	auto all_reduce_n(thrust::cuda::pointer<thrust::complex<double>> first, Size n, thrust::cuda::pointer<thrust::complex<double>> dest) {
		using thrust::reinterpret_pointer_cast;
		all_reduce_n(
			reinterpret_pointer_cast<thrust::cuda::pointer<double>>(first), n*2,
			reinterpret_pointer_cast<thrust::cuda::pointer<double>>(dest ), std::plus<>{}
		);
	}
	template<class Size>
	auto all_reduce_n(thrust::cuda::universal_pointer<thrust::complex<double>> first, Size n, thrust::cuda::universal_pointer<thrust::complex<double>> dest) {
		using thrust::reinterpret_pointer_cast;
		all_reduce_n(
			reinterpret_pointer_cast<thrust::cuda::universal_pointer<double>>(first), n*2,
			reinterpret_pointer_cast<thrust::cuda::universal_pointer<double>>(dest ), std::plus<>{}
		);
	}

	~communicator() {
		if(impl_) {
			ncclCommDestroy(impl_);  // call ncclCommFinalize internally if necessary
		}
	}

	int rank() const {  // aka user_rank()
		int ret;
		ncclResult_t r = ncclCommUserRank(impl_, &ret);
		switch(r) {
			case ncclSuccess: break;
			case ncclUnhandledCudaError: assert(0);
			case ncclSystemError: assert(0);
			case ncclInternalError: assert(0);
			case ncclInvalidArgument: assert(0);
			case ncclInvalidUsage: assert(0);
			case ncclRemoteError: assert(0);
			case ncclNumResults: assert(0);
		}
		return ret;
	}
	int count() const {
		int ret;
		ncclResult_t r = ncclCommCount(impl_, &ret);
		switch(r) {
			case ncclSuccess: break;
			case ncclUnhandledCudaError: assert(0);
			case ncclSystemError: assert(0);
			case ncclInternalError: assert(0);
			case ncclInvalidArgument: assert(0);
			case ncclInvalidUsage: assert(0);
			case ncclRemoteError: assert(0);
			case ncclNumResults: assert(0);
		}
		return ret;
	}
	[[deprecated("in NCCL nomenclature `.size` is called `.count`")]] int size() const {return count();}
	[[nodiscard]] bool    empty() const {return not count();}
	[[nodiscard]] bool is_empty() const {return not count();}

	[[nodiscard]] bool root() const {return not rank();}

	[[deprecated("using comm handle, try implementating")]] ncclComm_t operator&() {return impl_;}

 private:
	ncclComm_t impl_;
};

struct group {
	static void start() {
		ncclResult_t r = ncclGroupStart();
		switch(r) {
			case ncclSuccess: break;
			case ncclUnhandledCudaError: assert(0);
			case ncclSystemError: assert(0);
			case ncclInternalError: assert(0);
			case ncclInvalidArgument: assert(0);
			case ncclInvalidUsage: assert(0);
			case ncclRemoteError: assert(0);
			case ncclNumResults: assert(0);
		}
	}
	static void end() {
		ncclResult_t r = ncclGroupEnd();
		switch(r) {
			case ncclSuccess: break;
			case ncclUnhandledCudaError: assert(0);
			case ncclSystemError: assert(0);
			case ncclInternalError: assert(0);
			case ncclInvalidArgument: assert(0);
			case ncclInvalidUsage: assert(0);
			case ncclRemoteError: assert(0);
			case ncclNumResults: assert(0);
		}
	}
};

class group_block {
	group g_;

 public:
	[[nodiscard]] group_block(group g = {}) : g_{std::move(g)} {g_.start();}
	~group_block() {g_.end();}
};

}  // end namespace nccl
}  // end namespace mpi3
}  // end namespace boost
#endif  // MPI3_NCCL_COMMUNICATOR_HPP_
