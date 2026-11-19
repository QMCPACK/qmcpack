#include <boost/multi/adaptors/thrust.hpp>

#include <assert.h>
#include <cuda_runtime.h>
#include <cutensor.h>
#include <stdio.h>
#include <stdlib.h>

#include <iostream>
#include <unordered_map>
#include <vector>

// Handle cuTENSOR errors
#define HANDLE_ERROR(x)                                         \
	{                                                           \
		const auto err = x;                                     \
		if(err != CUTENSOR_STATUS_SUCCESS) {                    \
			printf("Error: %s\n", cutensorGetErrorString(err)); \
			exit(-1);                                           \
		}                                                       \
	};

#define HANDLE_CUDA_ERROR(x)                                \
	{                                                       \
		const auto err = x;                                 \
		if(err != cudaSuccess) {                            \
			printf("Error: %s\n", cudaGetErrorString(err)); \
			exit(-1);                                       \
		}                                                   \
	};

namespace boost::multi::cutensor {

template<class T> struct data_type;
template<> struct data_type<float> : std::integral_constant<cutensorDataType_t, CUTENSOR_R_32F> {};

template<class T> struct compute_type;
template<> struct compute_type<float> {
	operator cutensorComputeDescriptor_t() const { return CUTENSOR_COMPUTE_DESC_32F; }
};

class context {
	cutensorHandle_t handle_;

 public:
	context() {
		HANDLE_ERROR(cutensorCreate(&handle_));
	}
	context(context const&) = delete;
	~context() { HANDLE_ERROR(cutensorDestroy(handle_)); }

	auto operator&() { return handle_; }
};

auto default_context() -> context& {
	thread_local context ctxt;
	return ctxt;
}

class stream {
	cudaStream_t stream_;

 public:
	stream() { HANDLE_CUDA_ERROR(cudaStreamCreate(&stream_)); }
	stream(context const&) = delete;
	~stream() { HANDLE_CUDA_ERROR(cudaStreamDestroy(stream_)); }

	auto operator&() { return stream_; }
};

auto default_stream() -> stream& {
	thread_local stream instance;
	return instance;
}

template<class T = void, ::boost::multi::dimensionality_t D = 0> class descriptor;

template<>
class descriptor<void> {
	cutensorTensorDescriptor_t desc_;

 public:
	template<class StridedLayout>
	descriptor(StridedLayout const& lyt, cutensorDataType_t type, uint32_t kAlignment, context& ctxt) {
		auto extents = std::apply([](auto... es) { return std::array{static_cast<int64_t>(es)...}; }, lyt.sizes());
		cutensorCreateTensorDescriptor(&ctxt, &desc_, lyt.dimensionality, extents.data(), NULL, /*stride*/
									   type, kAlignment);
	}

	template<class StridedLayout>
	descriptor(StridedLayout const& lyt, cutensorDataType_t type, uint32_t kAlignment)
	: descriptor(lyt, type, kAlignment, default_context()) {}

	descriptor(descriptor const&) = delete;
	~descriptor() { HANDLE_ERROR(cutensorDestroyTensorDescriptor(desc_)); }

	static auto default_alignment() { return 128; }

	auto operator&() { return desc_; }
};

template<class T, ::boost::multi::dimensionality_t D>
class descriptor : public descriptor<void> {
 public:
	using type = T;

	constexpr static dimensionality_t dimensionality = D;

	template<class StridedLayout>
	descriptor(StridedLayout const& lyt, uint32_t kAlignment)
	: descriptor<void>(lyt, data_type<T>::value, kAlignment) {}

	template<class StridedLayout>
	descriptor(StridedLayout const& lyt) : descriptor(lyt, default_alignment()) {}
};

struct operation {
	cutensorOperationDescriptor_t desc_;

 public:
	~operation() {
		HANDLE_ERROR(cutensorDestroyOperationDescriptor(desc_));
	}
	auto operator&() { return desc_; }

	auto get_scalar_type(context& ctxt) const {
		cutensorDataType_t scalarType;
		HANDLE_ERROR(cutensorOperationDescriptorGetAttribute(&ctxt, desc_, CUTENSOR_OPERATION_DESCRIPTOR_SCALAR_TYPE, (void*)&scalarType, sizeof(scalarType)));
		return scalarType;
	}

	auto get_scalar_type() const { return get_scalar_type(default_context()); }

	static auto contraction(descriptor<>& dA, std::vector<int> const& modeA, descriptor<>& dB, std::vector<int> const& modeB, descriptor<>& dC, std::vector<int> const& modeC, cutensorComputeDescriptor_t descCompute, context& ctxt) {
		operation ret;
		HANDLE_ERROR(cutensorCreateContraction(&ctxt, &ret.desc_, &dA, modeA.data(), /* unary operator A*/ CUTENSOR_OP_IDENTITY, &dB, modeB.data(), /* unary operator B*/ CUTENSOR_OP_IDENTITY, &dC, modeC.data(), /* unary operator C*/ CUTENSOR_OP_IDENTITY, &dC, modeC.data(), descCompute));
		return ret;
	}

	static auto contraction(descriptor<>& dA, std::vector<int> const& modeA, descriptor<>& dB, std::vector<int> const& modeB, descriptor<>& dC, std::vector<int> const& modeC, cutensorComputeDescriptor_t descCompute) {
		return contraction(dA, modeA, dB, modeB, dC, modeC, descCompute, default_context());
	}

	template<class ComputeType>
	static auto contraction(descriptor<>& dA, std::vector<int> const& modeA, descriptor<>& dB, std::vector<int> const& modeB, descriptor<>& dC, std::vector<int> const& modeC) {
		return contraction(dA, modeA, dB, modeB, dC, modeC, multi::cutensor::compute_type<ComputeType>{});
	}
};

struct plan_preference {
	cutensorPlanPreference_t planPref_;
	plan_preference(cutensorAlgo_t const algo, context& ctxt) {
		HANDLE_ERROR(cutensorCreatePlanPreference(
			&ctxt,
			&planPref_,
			algo,
			CUTENSOR_JIT_MODE_NONE
		));
	}
	plan_preference(cutensorAlgo_t const algo = CUTENSOR_ALGO_DEFAULT)
	: plan_preference(algo, default_context()) {}

	plan_preference(plan_preference const&) = delete;
	~plan_preference() {
		HANDLE_ERROR(cutensorDestroyPlanPreference(planPref_));
	}
	auto operator&() { return planPref_; }
};

template<class Allocator = ::thrust::cuda::allocator<std::byte>> class plan;

template<>
class plan<void> {
	cutensorPlan_t plan_;

 public:
	static auto estimate_workspace_size(operation& op, plan_preference const& pp, cutensorWorksizePreference_t workspacePref, context& ctxt) {
		uint64_t ret;
		HANDLE_ERROR(cutensorEstimateWorkspaceSize(&ctxt, &op, &const_cast<plan_preference&>(pp), workspacePref, &ret));
		return ret;
	}

	static auto estimate_workspace_size(operation& op, plan_preference const& pp, cutensorWorksizePreference_t workspacePref) {
		return estimate_workspace_size(op, pp, workspacePref, default_context());
	}

	static auto estimate_workspace_size(operation& op, plan_preference const& pp) {
		return estimate_workspace_size(op, pp, CUTENSOR_WORKSPACE_DEFAULT);
	}

	static auto estimate_workspace_size(operation& op) {
		return estimate_workspace_size(op, plan_preference{});
	}

	plan(multi::cutensor::operation& op, plan_preference const& pp, uint64_t workspaceSizeEstimate, context& ctxt) {
		HANDLE_ERROR(cutensorCreatePlan(&ctxt, &plan_, &op, &const_cast<plan_preference&>(pp), workspaceSizeEstimate));
	}

	plan(multi::cutensor::operation& op, plan_preference const& pp, cutensorWorksizePreference_t workspacePref, context& ctxt)
	: plan(op, pp, estimate_workspace_size(op, pp, workspacePref, ctxt), ctxt) {}

	plan(multi::cutensor::operation& op, plan_preference const& pp, cutensorWorksizePreference_t workspacePref)
	: plan(op, pp, workspacePref, default_context()) {}

	plan(multi::cutensor::operation& op, plan_preference const& pp)
	: plan(op, pp, CUTENSOR_WORKSPACE_DEFAULT) {}

	plan(multi::cutensor::operation& op)
	: plan(op, plan_preference{}) {}

	plan(plan const&) = delete;
	~plan() {
		HANDLE_ERROR(cutensorDestroyPlan(plan_));
	}

	auto operator&() { return plan_; }

	auto get_required_workspace(context& ctxt) const {
		uint64_t actualWorkspaceSize;
		HANDLE_ERROR(cutensorPlanGetAttribute(&ctxt, plan_, CUTENSOR_PLAN_REQUIRED_WORKSPACE, &actualWorkspaceSize, sizeof(actualWorkspaceSize)));
		return actualWorkspaceSize;
	}

	auto get_required_workspace() const { return get_required_workspace(default_context()); }

	template<class ComputeType>
	void contract(ComputeType const& alpha, void* A_d, void* B_d, ComputeType const& beta, void* C_d, void* work, uint64_t actualWorkspaceSize, cudaStream_t stream, context& ctxt) const {

		assert(uintptr_t(A_d) % multi::cutensor::descriptor<>::default_alignment() == 0);
		assert(uintptr_t(B_d) % multi::cutensor::descriptor<>::default_alignment() == 0);
		assert(uintptr_t(C_d) % multi::cutensor::descriptor<>::default_alignment() == 0);

		HANDLE_ERROR(cutensorContract(&ctxt, plan_, (void*)&alpha, A_d, B_d, (void*)&beta, C_d, C_d, (void*)::thrust::raw_pointer_cast(work), actualWorkspaceSize, stream));
	}

	template<class ComputeType>
	void contract(ComputeType const& alpha, void* A_d, void* B_d, ComputeType const& beta, void* C_d, void* work, uint64_t actualWorkspaceSize, cudaStream_t stream) const {
		return contract(alpha, A_d, B_d, beta, C_d, work, actualWorkspaceSize, stream, default_context());
	}

	template<class ComputeType>
	void contract(ComputeType const& alpha, void* A_d, void* B_d, ComputeType const& beta, void* C_d, void* work, uint64_t actualWorkspaceSize) const {
		contract(alpha, A_d, B_d, beta, C_d, work, actualWorkspaceSize, &default_stream());
	}

	template<class ComputeType, class Allocator = ::thrust::cuda::allocator<std::byte>>
	void contract(ComputeType const& alpha, void* A_d, void* B_d, ComputeType const& beta, void* C_d, Allocator alloc = {}) const {
		uint64_t actualWorkspaceSize = this->get_required_workspace();

		auto work = actualWorkspaceSize ? alloc.allocate(actualWorkspaceSize) : nullptr;
		assert(work == nullptr || uintptr_t(raw_pointer_cast(work)) % multi::cutensor::descriptor<>::default_alignment() == 0);  // workspace must be aligned to 128 byte-boundary

		contract(alpha, A_d, B_d, beta, C_d, (void*)raw_pointer_cast(work), actualWorkspaceSize);
		if(work) {
			alloc.deallocate(work, actualWorkspaceSize);
		}
	}
};

template<class Allocator>
class plan : plan<void> {
	Allocator alloc_;

	typename std::allocator_traits<Allocator>::size_type work_size_;
	typename std::allocator_traits<Allocator>::pointer   work_;

 public:
	plan(operation& op, Allocator const& alloc = {})
	: plan<void>(op), alloc_{alloc}, work_size_{plan<void>::get_required_workspace()}, work_{work_size_ ? alloc_.allocate(work_size_) : nullptr} {}

	template<class ComputeType>
	void contract(ComputeType const& alpha, void* A_d, void* B_d, ComputeType const& beta, void* C_d) {  // no const! (workspace is mutable)
		plan<void>::contract(alpha, A_d, B_d, beta, C_d, (void*)raw_pointer_cast(work_), work_size_);
	}

	~plan() {
		if(work_) {
			alloc_.deallocate(work_, work_size_);
		}
	}
};

}  // namespace boost::multi::cutensor

namespace multi = boost::multi;

int main() {
	using floatTypeCompute = float;

	// Create vector of modes
	std::vector<int>
		modeA{'m', 'h', 'k', 'n'},
		modeB{'u', 'k', 'v', 'h'},
		modeC{'m', 'u', 'n', 'v'};

	// Extents
	std::unordered_map<int, int64_t> extent = {
		{'m', 96},
		{'n', 96},
		{'u', 96},
		{'v', 64},
		{'h', 64},
		{'k', 64},
	};

	multi::thrust::cuda::array<float, 4>
		A_dev = +([](auto...) { return (static_cast<float>(rand()) / static_cast<float>(RAND_MAX) - 0.5) * 100; } ^ multi::extents_t<4>{extent['m'], extent['h'], extent['k'], extent['n']}),
		B_dev = +([](auto...) { return (static_cast<float>(rand()) / static_cast<float>(RAND_MAX) - 0.5) * 100; } ^ multi::extents_t<4>{extent['u'], extent['k'], extent['v'], extent['h']}),
		C_dev = +([](auto...) { return (static_cast<float>(rand()) / static_cast<float>(RAND_MAX) - 0.5) * 100; } ^ multi::extents_t<4>{extent['m'], extent['u'], extent['n'], extent['v']})
	;

	multi::cutensor::descriptor<decltype(A_dev)::element, decltype(A_dev)::dimensionality> descA(A_dev.layout());
	multi::cutensor::descriptor<decltype(B_dev)::element, decltype(B_dev)::dimensionality> descB(B_dev.layout());
	multi::cutensor::descriptor<decltype(C_dev)::element, decltype(C_dev)::dimensionality> descC(C_dev.layout());

	using floatTypeCompute = float;

	multi::cutensor::operation opc = multi::cutensor::operation::contraction<floatTypeCompute>(
		descA, modeA,
		descB, modeB,
		descC, modeC
	);

	auto alpha = static_cast<floatTypeCompute>(1.1f);
	auto beta  = static_cast<floatTypeCompute>(0.0f);

	assert(opc.get_scalar_type() == multi::cutensor::data_type<decltype(alpha)>::value);
	assert(opc.get_scalar_type() == multi::cutensor::data_type<decltype(beta)>::value);

	multi::cutensor::plan<>(opc).contract(
		alpha, (void*)raw_pointer_cast(A_dev.base()), (void*)raw_pointer_cast(B_dev.base()),
		beta, (void*)raw_pointer_cast(C_dev.base())
	);
}
