// Copyright 2018-2023 Alfredo A. Correa

#ifndef BOOST_MPI3_ENVIRONMENT_HPP
#define BOOST_MPI3_ENVIRONMENT_HPP

#pragma once

#include "./communicator.hpp"
#include "./core.hpp"
#include "./wall_clock.hpp"

#include "./detail/call.hpp"
#include "./version.hpp"

#include <mpi.h>

#include <string>

namespace boost {
namespace mpi3 {

enum class thread_level : int {
	single     = MPI_THREAD_SINGLE,
	funneled   = MPI_THREAD_FUNNELED,
	serialized = MPI_THREAD_SERIALIZED,
#if not defined(EXAMPI)
	multiple   = MPI_THREAD_MULTIPLE
#endif
};

using thread = thread_level;

inline void finalize() {
	assert(initialized());
	assert(not finalized());

	if(int const count = std::uncaught_exceptions()) {
		std::cerr << "finalizing MPI environment with " << count << " uncaught exceptions";
	}

	std::set_terminate(&std::abort);
	int const s = MPI_Finalize();  // TODO(correaa) modernize call?
	if(s != MPI_SUCCESS) {
		std::terminate();
	}  //{throw std::runtime_error{"cannot finalize"};}
}
inline void myterminate() {
	std::cerr << "myterminate handler called" << '\n';
	finalize();
	std::abort();
}

inline void initialize(int& argc, char**& argv) {
	assert(not initialized());  // avoid double initialization`
#if not defined(EXAMPI)
	assert(not finalized());
#endif

	if(mpi3::version() != mpi3::Version()) {
		std::cerr << "WARNING: MPI version inconsistency\n";
		std::cerr << "Compile version (" << mpi3::version() <<") and linked version (" << mpi3::Version() << ") do not agree. Likely a link error.";
	}

#ifndef OMPI_MAJOR_VERSION
	if(char const* ompi_size_cstr = std::getenv("OMPI_COMM_WORLD_SIZE")) {  // NOLINT(concurrency-mt-unsafe)
		if(char const* ompi_rank_cstr = std::getenv("OMPI_COMM_WORLD_RANK")) {  // NOLINT(concurrency-mt-unsafe)
			if(std::string{ompi_rank_cstr} == "0") {
				std::cerr << "WARNING: MPI environment inconsistency?\n";
				std::cerr << "running program " << *argv << " compiled with " << boost::mpi3::library_version_short() << " but allocated with OMPI (np = " << ompi_size_cstr << "). Try running with `mpirun.mpich` or load MPICH module.\n\n";
			}
		} else {
			std::cerr << "WARNING: MPI environment inconsistency?\n";
			std::cerr << "running program " << *argv << " compiled with " << boost::mpi3::library_version_short() << " but allocated with OMPI (np = " << ompi_size_cstr << "). Try running with `mpirun.mpich` or load MPICH module.\n\n";
		}
		using namespace std::chrono_literals;
		std::this_thread::sleep_for(1s);
	}
#endif
#ifndef MPICH_VERSION
	if(char const* pmi_size_cstr = std::getenv("PMI_SIZE")) {  // NOLINT(concurrency-mt-unsafe)
		if(char const* pmi_rank_cstr = std::getenv("PMI_RANK")) {  // NOLINT(concurrency-mt-unsafe)
			if(std::string{pmi_rank_cstr} == "0") {
				std::cerr << "WARNING: MPI environment inconsistency?\n";
				std::cerr << "running program " << *argv << " compiled with " << boost::mpi3::library_version_short() << " but allocated with PMI (Hydra or MPICH) (np = " << pmi_size_cstr << "). Try running with `mpirun.openmpi` or load OpenMPI module.\n\n";
			}
		} else {
			std::cerr << "WARNING: MPI environment inconsistency?\n";
			std::cerr << "running program " << *argv << " compiled with " << boost::mpi3::library_version_short() << " but allocated with PMI (Hydra or MPICH) (np = " << pmi_size_cstr << "). Try running with `mpirun.openmpi` or load OpenMPI module.\n\n";
		}
		using namespace std::chrono_literals;
		std::this_thread::sleep_for(1s);
	}
#endif

	MPI_(Init)(&argc, &argv);  // TODO(correaa) modernize call?

	int nprocs = -1;
	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

	int rank = -1;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	if(char const* ompi_size_cstr = std::getenv("OMPI_COMM_WORLD_SIZE")) {  // NOLINT(concurrency-mt-unsafe)
		char const* ompi_rank_cstr = std::getenv("OMPI_COMM_WORLD_RANK");  // NOLINT(concurrency-mt-unsafe)
		if(std::to_string(nprocs) != ompi_size_cstr and std::string{ompi_rank_cstr} == "0" and rank == 0) {
			std::cerr << "WARNING: MPI size inconsistency?\n";
			std::cerr << "running program " << *argv << " in " << std::to_string(nprocs) << " processes but allocated with " << ompi_size_cstr << " processes \n\n";
			MPI_Barrier(MPI_COMM_WORLD);
			using namespace std::chrono_literals;
			std::this_thread::sleep_for(1s);
		}
	}
	if(char const* pmi_size_cstr = std::getenv("PMI_SIZE")) {  // NOLINT(concurrency-mt-unsafe)
		char const* pmi_rank_cstr = std::getenv("PMI_RANK");  // NOLINT(concurrency-mt-unsafe)
		if(std::to_string(nprocs) != pmi_size_cstr and std::string{pmi_rank_cstr} == "0" and rank == 0) {
			std::cerr << "WARNING: MPI size inconsistency?\n";
			std::cerr << "running program " << *argv << " in " << std::to_string(nprocs) << " processes but allocated with " << pmi_size_cstr << " processes \n\n";
			MPI_Barrier(MPI_COMM_WORLD);
			using namespace std::chrono_literals;
			std::this_thread::sleep_for(1s);
		}
	}
}

inline thread_level initialize(int& argc, char**& argv, thread_level required) {
	int provided;  // NOLINT(cppcoreguidelines-init-variables) : delayed initialization
	MPI_(Init_thread)(&argc, &argv, static_cast<int>(required), &provided);
	return static_cast<thread_level>(provided);
}

inline thread_level initialize_thread(thread_level required) {
	int provided;  // NOLINT(cppcoreguidelines-init-variables) delayed initialization
	MPI_(Init_thread)(nullptr, nullptr, static_cast<int>(required), &provided);
	return static_cast<thread_level>(provided);
}

inline thread_level initialize(thread_level required = thread_level::single) {
	return initialize_thread(required);
}

// inline void throw_error_fn(MPI_Comm* comm, int* errorcode, ...) {
//  char name[MPI_MAX_OBJECT_NAME];
//  int rlen;  // NOLINT(cppcoreguidelines-init-variables) delayed initialization
//  int status = MPI_Comm_get_name(*comm, name, &rlen);
//  if(status != MPI_SUCCESS) {throw std::runtime_error{"cannot get name"};}
//  std::string sname(name, rlen);
//  throw std::runtime_error{"error code "+ std::to_string(*errorcode) +" from comm "+ sname};
// }

inline thread_level initialize_thread(
	int& argc, char**& argv, thread_level required
) {
	int ret;  // NOLINT(cppcoreguidelines-init-variables) delayed initialization
	MPI_(Init_thread)
	(&argc, &argv, static_cast<int>(required), &ret);
	return static_cast<thread_level>(ret);
}

#if not defined(EXAMPI)
inline thread_level thread_support() {
	return static_cast<thread_level>(MPI_(Query_thread)());
}
#endif

inline bool cuda_support() {
	#if defined(CUDA_AWARE_SUPPORT) and CUDA_AWARE_SUPPORT
		static bool const ret = Query_cuda_support();
	#else
		static bool const ret = false;
	#endif
	return ret;
}

inline bool is_thread_main() {
	int flag;  // NOLINT(cppcoreguidelines-init-variables) : delayed initialization
	MPI_(Is_thread_main)
	(&flag);
	return flag != 0;
}

inline std::string processor_name() { return detail::call<&MPI_Get_processor_name>(); }

inline std::string get_processor_name() { return detail::call<&MPI_Get_processor_name>(); }

class environment {
 public:
	environment() {
#if not defined(EXAMPI)
		initialize_thread(thread_level::multiple);
		named_attributes_key_f() = std::make_unique<communicator::keyval<std::map<std::string, mpi3::any>>>();
#else
		initialize_thread(thread_level::serialized);
#endif
	}
	explicit environment(thread_level required) {
		initialize_thread(required);
	#if not defined(EXAMPI)
		named_attributes_key_f() = std::make_unique<communicator::keyval<std::map<std::string, mpi3::any>>>();
	#endif
	}
	explicit environment(int& argc, char**& argv) {  // cppcheck-suppress [constParameter, constParameterReference] ; bug in cppcheck 2.3 and 2.9 or it can't see through the MPI C-API
		initialize(argc, argv);  // initialize(argc, argv); // TODO have an environment_mt/st version?
	#if not defined(EXAMPI)
		named_attributes_key_f() = std::make_unique<communicator::keyval<std::map<std::string, mpi3::any>>>();
	#endif
	}
	explicit environment(int& argc, char**& argv, thread_level required) {  // cppcheck-suppress [constParameter, constParameterReference] ; bug in cppcheck 2.3 and 2.9 or it can't see through the MPI C-API
		initialize(argc, argv, required);
	#if not defined(EXAMPI)
		named_attributes_key_f() = std::make_unique<communicator::keyval<std::map<std::string, mpi3::any>>>();
	#endif
	}

	environment(environment const&) = delete;
	environment(environment&&)      = delete;

	environment& operator=(environment const&) = delete;
	environment& operator=(environment&&)      = delete;

	~environment() noexcept {  // NOLINT(bugprone-exception-escape) finalizes throws as an instance of UB
	#if not defined(EXAMPI)
		named_attributes_key_f().reset();
	#endif
		finalize();  // cppcheck-suppress throwInNoexceptFunction ; finalizes throws as an instance of UB
	}

#if not defined(EXAMPI)
	inline static thread_level thread_support() { return mpi3::thread_support(); }
#endif
	//  static /*inline*/ communicator::keyval<int> const* color_key_p;
	//  static communicator::keyval<int> const& color_key(){return *color_key_p;}
	//  static /*inline*/ communicator::keyval<std::map<std::string, mpi3::any>> const* named_attributes_key_p;

#if not defined(EXAMPI)
	static std::unique_ptr<communicator::keyval<std::map<std::string, mpi3::any>> const>& named_attributes_key_f() {
		static std::unique_ptr<communicator::keyval<std::map<std::string, mpi3::any>> const> named_attributes_key_p;
		return named_attributes_key_p;
	}
	static communicator::keyval<std::map<std::string, mpi3::any>> const& named_attributes_key() {
		//  static communicator::keyval<std::map<std::string, mpi3::any>> const named_attributes_key_p;
		//  return named_attributes_key_p;
		return *named_attributes_key_f();
	}
#endif
	static bool cuda_support() { return mpi3::cuda_support(); }  // cppcheck-suppress knownConditionTrueFalse ; might be known at compile time

	static inline void initialize() { mpi3::initialize(); }
	static inline void initialize(int argc, char** argv) { mpi3::initialize(argc, argv); }

	static inline thread_level initialize(thread_level required) { return mpi3::initialize_thread(required); }
	static inline thread_level initialize(int argc, char** argv, thread_level req) { return mpi3::initialize_thread(argc, argv, req); }

	static inline void finalize() { mpi3::finalize(); }

	static inline bool is_initialized() { return mpi3::initialized(); }
	static inline bool is_finalized() { return mpi3::finalized(); }

	using wall_clock = mpi3::wall_clock;
	explicit operator bool() const { return initialized(); }

	static bool is_thread_main() { return mpi3::is_thread_main(); }

	static inline communicator& get_self_instance() {
		assert(initialized());
		static communicator instance = [] {
			//  MPI_Comm_create_errhandler(&throw_error_fn, &throw_error_);
			//  MPI_Comm_set_errhandler(MPI_COMM_WORLD, throw_error_);
			MPI_Comm_set_errhandler(MPI_COMM_SELF, MPI_ERRORS_RETURN);
			//  MPI_Comm_set_errhandler(MPI_COMM_WORLD, MPI_ERRORS_ARE_FATAL);
			return communicator{MPI_COMM_SELF};
		}();
		return instance;
	}

	static communicator self() {  // returns a copy!
		MPI_Comm_set_errhandler(MPI_COMM_SELF, MPI_ERRORS_RETURN);
		return communicator{MPI_COMM_SELF};
	}
	static inline communicator& get_world_instance() {
		assert(initialized());
		static communicator instance = [] {
			//  MPI_Comm_create_errhandler(&throw_error_fn, &throw_error_);
			//  MPI_Comm_set_errhandler(MPI_COMM_WORLD, throw_error_);
			MPI_Comm_set_errhandler(MPI_COMM_WORLD, MPI_ERRORS_RETURN);
			//  MPI_Comm_set_errhandler(MPI_COMM_WORLD, MPI_ERRORS_ARE_FATAL);
			return communicator{MPI_COMM_WORLD};
		}();
		return instance;
	}

	[[nodiscard]] communicator world() {  // NOLINT(readability-convert-member-functions-to-static) to force instance
		communicator ret{get_world_instance()};
		ret.set_name("world");
		return ret;
	}

	static std::string processor_name() { return get_processor_name(); }
	static std::string get_processor_name() { return mpi3::get_processor_name(); }

	inline static auto wall_time() { return mpi3::wall_time(); }
	inline static auto wall_tick() { return mpi3::wall_tick(); }

	template<class Duration = wall_clock::duration>
	static auto wall_sleep_for(Duration d) { return mpi3::wall_sleep_for(d); }
};

#if not defined(EXAMPI)
inline mpi3::any& communicator::attribute(std::string const& s) {
	return attribute(environment::named_attributes_key())[s];
}
#endif

}  // end namespace mpi3
}  // end namespace boost
#endif
