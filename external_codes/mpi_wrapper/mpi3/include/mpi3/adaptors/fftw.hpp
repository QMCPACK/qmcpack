// Copyright 2019-2024 Alfredo A. Correa

#include <fftw3-mpi.h>

#include <multi/adaptors/fftw.hpp>

#include <multi/array_ref.hpp>

#include <mpi3/allocator.hpp>
#include <mpi3/communicator.hpp>

namespace boost::mpi3::fftw {

struct environment {
	environment() { fftw_mpi_init(); }

	environment(environment const&) = delete;
	environment(environment&&)      = delete;

	environment& operator=(environment const&) = delete;
	environment& operator=(environment&&)      = delete;

	~environment() { fftw_mpi_cleanup(); }

	template<class... Args>
	static auto local_size_2d(Args... args) { return fftw_mpi_local_size_2d(args...); }

	template<class... Args>
	static auto local_size_many(Args... args) { return fftw_mpi_local_size_many(args...); }
};

template<class T>
using default_allocator =
	// std::allocator<T>
	boost::multi::fftw::allocator<T>
	// boost::mpi3::allocator<T>
	;

struct local_2d {
 private:
	std::ptrdiff_t             n0_     = -1;
	std::ptrdiff_t             start0_ = -1;
	std::ptrdiff_t             count_;
	boost::mpi3::communicator* handle_;
	multi::extensions_t<2>     exts_;

 public:
	local_2d(multi::extensions_t<2> const& exts, boost::mpi3::communicator& comm)
	: count_{environment::local_size_2d(std::get<0>(exts).size(), std::get<1>(exts).size(), &comm, &n0_, &start0_)}, handle_{&comm}, exts_{exts} {}

	auto count() const { return count_; }
	auto extension() const { return multi::extension_t{start0_, start0_ + n0_}; }
	auto comm() const -> communicator& { return *handle_; }
	auto global_extensions() const { return exts_; }
    auto static block() {return FFTW_MPI_DEFAULT_BLOCK;}
};

struct local_2d_many {
 private:
	std::ptrdiff_t n0_     = -1;
	std::ptrdiff_t start0_ = -1;
	std::ptrdiff_t count_;

	boost::mpi3::communicator* handle_;
	multi::extensions_t<2>     exts_;
	std::ptrdiff_t             block_;

 public:
	local_2d_many(multi::extensions_t<2> const& exts, boost::mpi3::communicator& comm) : handle_{&comm}, exts_{exts}, block_{FFTW_MPI_DEFAULT_BLOCK} {
		count_ = environment::local_size_many(2, std::array<std::ptrdiff_t, 2>{std::get<0>(exts).size(), std::get<1>(exts).size()}.data(), /*howmany*/ 1, FFTW_MPI_DEFAULT_BLOCK /*std::get<0>(ext).size()*/, &comm, &n0_, &start0_);
	}
	local_2d_many(multi::extensions_t<2> const& exts, boost::mpi3::communicator& comm, std::ptrdiff_t block) : handle_{&comm}, exts_{exts}, block_{block} {
		count_ = environment::local_size_many(2, std::array<std::ptrdiff_t, 2>{std::get<0>(exts).size(), std::get<1>(exts).size()}.data(), /*howmany*/ 1, block, &comm, &n0_, &start0_);
	}

	auto count() const { return count_; }
	auto extension() const { return multi::extension_t{start0_, start0_ + n0_}; }
	auto comm() const -> communicator& { return *handle_; }
	auto global_extensions() const { return exts_; }
	auto block() const {return block_;}
};

template<
	class T,
	boost::multi::dimensionality_type D,
	class LocalLayout = local_2d,
	class Alloc       = default_allocator<T>>
class array;

template<
	class T,
	boost::multi::dimensionality_type D,
	class Alloc = default_allocator<T>>
class unbalanced_array;

// namespace bmpi3 = boost::mpi3;

template<class T, class LocalLayout, class Alloc>
class array<T, multi::dimensionality_type{2}, LocalLayout, Alloc> {
	LocalLayout local_layout_;

	Alloc alloc_;

	boost::multi::array_ptr<T, 2, typename std::allocator_traits<Alloc>::pointer> local_ptr_;

 public:
	using element_type = T;


	array(multi::extensions_t<2> exts, boost::mpi3::communicator& comm, Alloc alloc = Alloc{})
	: alloc_{alloc},
	  local_layout_(exts, comm),
	  local_ptr_{alloc_.allocate(local_layout_.count()), multi::extensions_t<2>(local_layout_.extension(), std::get<1>(exts))} {}

	array(multi::extensions_t<2> exts, element_type const& e, boost::mpi3::communicator& comm, Alloc alloc = Alloc{})
	: alloc_{alloc},
	  local_layout_(exts, comm),
	  local_ptr_{alloc_.allocate(local_layout_.count()), multi::extensions_t<2>(local_layout_.extension(), std::get<1>(exts))} {
		std::uninitialized_fill_n(local_ptr_.base(), local_ptr_->num_elements(), e);
	}

	array(array const&) = delete;
	array(array&&) = delete;

	auto operator=(array const&) -> array& = delete;
	auto operator=(array&&) -> array& = delete;

	boost::multi::array_ref<T, 2>  local_cutout() & { return *local_ptr_; }
	boost::multi::array_cref<T, 2> local_cutout() const& { return *local_ptr_; }

	auto local_layout() const { return local_layout_; }

	ptrdiff_t local_count() const& { return local_layout_.count(); }

	// auto extensions() const& { return multi::extensions_t<2>{local_layout_.n0_, std::get<1>(local_cutout().extensions())}; }

	// ptrdiff_t num_elements() const& { return multi::layout_t<2>(extensions()).num_elements(); }

	template<class Array>
	static auto from_scatter(Array const& snd) -> array {
		array ret(snd.extensions());
		ret.scatter(snd);
		return ret;
	}

	// template<class Array>
	// void scatter(Array const& snd) & {
	//  auto& comm = reinterpret_cast<boost::mpi3::communicator&>(handle_);

	//  auto const sendcounts = comm |= static_cast<int>(local_cutout().num_elements());
	//  auto const displs     = comm |= static_cast<int>(snd[local_cutout().extension().front()].base() - snd.base());

	//  MPI_Scatterv(
	//    snd.base(), sendcounts.data(), displs.data(), MPI_DOUBLE_COMPLEX,
	//    local_cutout().base(), local_cutout().num_elements(), MPI_DOUBLE_COMPLEX,
	//    0, &comm
	//  );
	// }

	// auto communicator() const -> boost::mpi3::communicator& {
	//  return const_cast<boost::mpi3::communicator&>(reinterpret_cast<boost::mpi3::communicator const&>(handle_));
	// }

	// template<class Array>
	// void all_gather(Array&& rcv) const& {
	//  assert(rcv.extensions() == extensions());

	//  auto& comm = const_cast<boost::mpi3::communicator&>(reinterpret_cast<boost::mpi3::communicator const&>(handle_));

	//  auto const recvcounts = comm |= static_cast<int>(local_cutout().num_elements());
	//  auto const displs     = comm |= static_cast<int>(rcv[local_cutout().extension().front()].base() - rcv.base());

	//  MPI_Allgatherv(
	//    local_cutout().base(), local_cutout().num_elements(), MPI_DOUBLE_COMPLEX,
	//    rcv.base(),
	//    recvcounts.data(), displs.data(), MPI_DOUBLE_COMPLEX,
	//    handle_
	//  );
	// }

	// template<class Alloc2 = Alloc>
	// explicit operator multi::array<T, 2, Alloc2>() const& {
	//  multi::array<T, 2, Alloc2> ret(extensions());
	//  all_gather(ret);
	//  return ret;
	// }

	auto extensions() const { return local_layout_.global_extensions(); }

	template<class OtherLayout, class OtherAlloc>
	array& operator=(array<T, 2, OtherLayout, OtherAlloc> const& other) {
		int P = -1;
		MPI_Comm_size(&local_layout_.comm(), &P);
		fftw_plan plan = fftw_mpi_plan_many_transpose(
			6, 6, /*howmany*/ 2 /*2 for complex*/, 1, 1,
			const_cast<double*>(reinterpret_cast<double const*>(other.local_cutout().base())),  // NOLINT(cppcoreguidelines-pro-type-const-cast,cppcoreguidelines-pro-type-reinterpret-cast)
			reinterpret_cast<double*>(this->local_cutout().base()),  // NOLINT(cppcoreguidelines-pro-type-reinterpret-cast)
			&local_layout_.comm(),
			FFTW_ESTIMATE
		);
		fftw_execute(plan);
		fftw_destroy_plan(plan);
		return *this;
	}

	// array& operator=(multi::array<T, 2> const& other) & {
	//  if(other.extensions() == extensions())
	//    local_cutout() = other.stenciled(std::get<0>(local_cutout().extensions()), std::get<1>(local_cutout().extensions()));
	//  else {
	//    array tmp{other};
	//    std::swap(*this, tmp);
	//  }
	//  return *this;
	// }
	// bool operator==(multi::array<T, 2> const& other) const&{
	//  if(other.extensions() != extensions()) return false;
	//  return comm_&=(local_cutout() == other.stenciled(std::get<0>(local_cutout().extensions()), std::get<1>(local_cutout().extensions())));
	// }
	// friend bool operator==(multi::array<T, 2> const& other, array const& self){
	//  return self.operator==(other);
	// }
	// bool operator==(array<T, 2> const& other) const&{assert(comm_==other.comm_);
	//  return comm_&=(local_cutout() == other.local_cutout());
	// }
	// array& operator=(array const& other)&{
	//  if(other.extensions() == this->extensions() and other.comm_ == other.comm_)
	//      local_cutout() = other.local_cutout();
	//  else assert(0);
	//  return *this;
	// }
	~array() { alloc_.deallocate(local_cutout().base(), local_layout_.count()); }

	template<class Array>
	void all_gather(Array&& rcv) const& {
		assert(rcv.extensions() == extensions());

		auto const recvcounts = local_layout_.comm() |= static_cast<int>(local_cutout().num_elements());
		auto const displs     = local_layout_.comm() |= static_cast<int>(rcv[local_cutout().extension().front()].base() - rcv.base());

		MPI_Allgatherv(
			local_cutout().base(), local_cutout().num_elements(), MPI_DOUBLE_COMPLEX,
			rcv.base(),
			recvcounts.data(), displs.data(), MPI_DOUBLE_COMPLEX,
			&local_layout_.comm()
		);
	}

	template<class Alloc2 = Alloc>
	explicit operator multi::array<T, 2, Alloc2>() const& {
		multi::array<T, 2, Alloc2> ret(extensions());
		all_gather(ret);
		return ret;
	}
};

template<class Array>
auto scatter(Array const& arr) {
	return array<typename Array::element_type, Array::dimensionality>::from_scatter(arr);
}

template<class MPIArrayIn, class MPIArrayOut>
auto dft_forward(MPIArrayIn const& A, MPIArrayOut& B) -> MPIArrayOut& {
	assert( &A.local_layout().comm() == &B.local_layout().comm() );

	// fftw_plan p = fftw_mpi_plan_dft_2d(
	//  std::get<0>(A.extensions()).size(), std::get<1>(A.extensions()).size(),
	//  (fftw_complex*)A.local_cutout().base(), (fftw_complex*)B.local_cutout().base(),
	//  &A.local_layout().comm(),
	//  FFTW_FORWARD, FFTW_ESTIMATE
	// );

	fftw_plan p = fftw_mpi_plan_many_dft(
		2, std::array<std::ptrdiff_t, 2>{std::get<0>(A.extensions()).size(), std::get<1>(A.extensions()).size()}.data(),
		1,
        A.local_layout().block(), B.local_layout().block(),  // FFTW_MPI_DEFAULT_BLOCK, FFTW_MPI_DEFAULT_BLOCK,
		const_cast<fftw_complex*>(reinterpret_cast<fftw_complex const*>(A.local_cutout().base())),  // NOLINT(cppcoreguidelines-pro-type-reinterpret-cast,cppcoreguidelines-pro-type-const-cast)
		                          reinterpret_cast<fftw_complex      *>(B.local_cutout().base()) ,  // NOLINT(cppcoreguidelines-pro-type-reinterpret-cast)
		&A.local_layout().comm(),
		FFTW_FORWARD, FFTW_ESTIMATE
	);
	fftw_execute(p);
	fftw_destroy_plan(p);

	// // B = BT;

	// fftw_plan t = fftw_mpi_plan_transpose(
	//  std::get<0>(A.extensions()).size(), std::get<1>(A.extensions()).size(),
	//  (double*)(fftw_complex*)BT.local_cutout().data_elements(), (double*)(fftw_complex*)B.local_cutout().data_elements(),
	//  &A.communicator(), FFTW_ESTIMATE);
	// fftw_execute(t);
	// fftw_destroy_plan(t);

	return B;
}

#if 0
template<class T, class Alloc>
class array<T, multi::dimensionality_type{2}, Alloc> {
	boost::mpi3::communicator* commhandle_;

	Alloc alloc_;

	class local_2d_type {
		std::ptrdiff_t n0_     = -1;
		std::ptrdiff_t start0_ = -1;
		std::ptrdiff_t count_;

	 public:
		local_2d_type(multi::extensions_t<2> const& ext, boost::mpi3::communicator& comm)
		: count_{environment::local_size_2d(std::get<0>(ext).size(), std::get<1>(ext).size(), &comm, &n0_, &start0_)} {}

		auto count() const { return count_; }
		auto extension() const { return multi::extension_t{start0_, start0_ + n0_}; }
	} local_;

 public:
	using element_type = T;
	using element_ptr  = typename std::allocator_traits<Alloc>::pointer;

 private:
	boost::multi::array_ptr<element_type, 2, element_ptr> local_ptr_;

 public:
    auto local_count() const {return local_.count();}

	array(multi::extensions_t<2> ext, element_type const& e, boost::mpi3::communicator& comm, Alloc alloc = Alloc{})
	: commhandle_{&comm},
	  alloc_{alloc},
	  local_{ext, comm},
	  local_ptr_{alloc_.allocate(local_.count()), multi::extensions_t<2>(local_.extension(), std::get<1>(ext))} {
		std::uninitialized_fill(local_ptr_->elements().begin(), local_ptr_->elements().end(), e);  // TODO(correaa) use adl_uninit_fill or uninitialized_fill member
		// std::uninitialized_fill_n(local_ptr_->base(), local_ptr_->num_elements(), e);
	}

    array(array const&);

	boost::multi::array_ref<T, 2>  local() & { return *local_ptr_; }
	boost::multi::array_cref<T, 2> local() const& { return *local_ptr_; }

	boost::multi::array_cref<T, 2> clocal() const { return local(); }

	// template<class Array>
	// void scatter(Array const& snd) & {
	//     auto& comm = reinterpret_cast<boost::mpi3::communicator&>(handle_);

	//     auto const sendcounts = comm |= static_cast<int>(local_cutout().num_elements());
	//     auto const displs     = comm |= static_cast<int>(snd[local_cutout().extension().front()].base() - snd.base());

	//     MPI_Scatterv(
	//         snd.base(), sendcounts.data(), displs.data(), MPI_DOUBLE_COMPLEX,
	//         local_cutout().base(), local_cutout().num_elements(), MPI_DOUBLE_COMPLEX,
	//         0, &comm
	//     );
	// }

	// auto communicator() const -> boost::mpi3::communicator& {
	//     return const_cast<boost::mpi3::communicator&>(reinterpret_cast<boost::mpi3::communicator const&>(handle_));
	// }

	// template<class Array>
	// void all_gather(Array&& rcv) const& {
	//     assert(rcv.extensions() == extensions());

	//     auto& comm = const_cast<boost::mpi3::communicator&>(reinterpret_cast<boost::mpi3::communicator const&>(handle_));

	//     auto const recvcounts = comm |= static_cast<int>(local_cutout().num_elements());
	//     auto const displs     = comm |= static_cast<int>(rcv[local_cutout().extension().front()].base() - rcv.base());

	//     MPI_Allgatherv(
	//         local_cutout().base(), local_cutout().num_elements(), MPI_DOUBLE_COMPLEX,
	//         rcv.base(),
	//         recvcounts.data(), displs.data(), MPI_DOUBLE_COMPLEX,
	//         handle_
	//     );
	// }

	// template<class Alloc2 = Alloc>
	// explicit operator multi::array<T, 2, Alloc2>() const& {
	//     multi::array<T, 2, Alloc2> ret(extensions());
	//     all_gather(ret);
	//     return ret;
	// }

	// array& operator=(multi::array<T, 2> const& other) & {
	//     if(other.extensions() == extensions())
	//         local_cutout() = other.stenciled(std::get<0>(local_cutout().extensions()), std::get<1>(local_cutout().extensions()));
	//     else {
	//         array tmp{other};
	//         std::swap(*this, tmp);
	//     }
	//     return *this;
	// }
	// // bool operator==(multi::array<T, 2> const& other) const&{
	// //  if(other.extensions() != extensions()) return false;
	// //  return comm_&=(local_cutout() == other.stenciled(std::get<0>(local_cutout().extensions()), std::get<1>(local_cutout().extensions())));
	// // }
	// // friend bool operator==(multi::array<T, 2> const& other, array const& self){
	// //  return self.operator==(other);
	// // }
	// // bool operator==(array<T, 2> const& other) const&{assert(comm_==other.comm_);
	// //  return comm_&=(local_cutout() == other.local_cutout());
	// // }
	// // array& operator=(array const& other)&{
	// //  if(other.extensions() == this->extensions() and other.comm_ == other.comm_)
	// //      local_cutout() = other.local_cutout();
	// //  else assert(0);
	// //  return *this;
	// // }
	~array() { alloc_.deallocate(local().base(), local_.count()); }
};
#endif

template<class T, class Alloc>
class unbalanced_array<T, multi::dimensionality_type{2}, Alloc> {
	boost::mpi3::communicator* commhandle_;

	Alloc alloc_;

	class local_2d_type {
		std::ptrdiff_t n0_     = -1;
		std::ptrdiff_t start0_ = -1;
		std::ptrdiff_t count_;

	 public:
		local_2d_type(multi::extensions_t<2> const& ext, boost::mpi3::communicator& comm)
		: count_{environment::local_size_2d(std::get<0>(ext).size(), std::get<1>(ext).size(), &comm, &n0_, &start0_)} {}

		auto count() const { return count_; }
		auto extension() const { return multi::extension_t{start0_, start0_ + n0_}; }
	} local_;

 public:
	using element_type = T;
	using element_ptr  = typename std::allocator_traits<Alloc>::pointer;

 private:
	boost::multi::array_ptr<element_type, 2, element_ptr> local_ptr_;

 public:
	unbalanced_array(multi::extensions_t<2> ext, element_type const& e, boost::mpi3::communicator& comm, Alloc alloc = Alloc{})
	: commhandle_{&comm},
	  alloc_{alloc},
	  local_{ext, *commhandle_},
	  local_ptr_{alloc_.allocate(local_.count()), multi::extensions_t<2>(local_.extension(), std::get<1>(ext))} {
		std::uninitialized_fill(local_ptr_->elements().begin(), local_ptr_->elements().end(), e);  // TODO(correaa) use adl_uninit_fill or uninitialized_fill member
		// std::uninitialized_fill_n(local_ptr_->base(), local_ptr_->num_elements(), e);
	}

	boost::multi::array_ref<T, 2>  local() & { return *local_ptr_; }
	boost::multi::array_cref<T, 2> local() const& { return *local_ptr_; }

	boost::multi::array_cref<T, 2> clocal() const { return local(); }
};

}  // namespace boost::mpi3::fftw
