//  -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;autowrap:nil;-*-
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

// #include <iostream>  // for std::clog
#include <string>

namespace boost {
namespace mpi3 {

enum class thread_level : int {
	single     = MPI_THREAD_SINGLE,
	funneled   = MPI_THREAD_FUNNELED,
	serialized = MPI_THREAD_SERIALIZED,
	multiple   = MPI_THREAD_MULTIPLE
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
	assert(not initialized());  // double initialize
	assert(not finalized());

	if(mpi3::version() != mpi3::Version()) {
		std::cerr << "WARNING: MPI version inconsistency\n";
		std::cerr << "Compile version (" << mpi3::version() <<") and linked version (" << mpi3::Version() << ") do not agree. Likely a link error.";
	}

	if([[maybe_unused]] char const* ompi_size_cstr = std::getenv("OMPI_COMM_WORLD_SIZE")) {  // NOLINT(concurrency-mt-unsafe)
#ifndef OMPI_MAJOR_VERSION
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
#endif
	}
	if([[maybe_unused]] char const* pmi_size_cstr = std::getenv("PMI_SIZE")) {  // NOLINT(concurrency-mt-unsafe)
#ifndef MPICH_VERSION
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
#endif
	}

	MPI_(Init)(&argc, &argv);  // TODO(correaa) modernize call?

	int nprocs = -1;
	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

	int rank = -1;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	if([[maybe_unused]] char const* ompi_size_cstr = std::getenv("OMPI_COMM_WORLD_SIZE")) {  // NOLINT(concurrency-mt-unsafe)
		char const* ompi_rank_cstr = std::getenv("OMPI_COMM_WORLD_RANK");  // NOLINT(concurrency-mt-unsafe)
		if(std::to_string(nprocs) != ompi_size_cstr and std::string{ompi_rank_cstr} == "0" and rank == 0) {
			std::cerr << "WARNING: MPI size inconsistency?\n";
			std::cerr << "running program " << *argv << " in " << std::to_string(nprocs) << " processes but allocated with " << ompi_size_cstr << " processes \n\n";
			MPI_Barrier(MPI_COMM_WORLD);
			using namespace std::chrono_literals;
			std::this_thread::sleep_for(1s);
		}
	}
	if([[maybe_unused]] char const* pmi_size_cstr = std::getenv("PMI_SIZE")) {  // NOLINT(concurrency-mt-unsafe)
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
//	char name[MPI_MAX_OBJECT_NAME];
//	int rlen;  // NOLINT(cppcoreguidelines-init-variables) delayed initialization
//	int status = MPI_Comm_get_name(*comm, name, &rlen);
//	if(status != MPI_SUCCESS) {throw std::runtime_error{"cannot get name"};}
//	std::string sname(name, rlen);
//	throw std::runtime_error{"error code "+ std::to_string(*errorcode) +" from comm "+ sname};
// }

inline thread_level initialize_thread(
	int& argc, char**& argv, thread_level required
) {
	int ret;  // NOLINT(cppcoreguidelines-init-variables) delayed initialization
	MPI_(Init_thread)
	(&argc, &argv, static_cast<int>(required), &ret);
	return static_cast<thread_level>(ret);
}

inline thread_level thread_support() {
	int r;  // NOLINT(cppcoreguidelines-init-variables) : delayed initialization
	MPI_(Query_thread)
	(&r);
	return static_cast<thread_level>(r);
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
		initialize_thread(thread_level::multiple);
		// std::clog << "ctor() environment" << std::endl;
		named_attributes_key_f() = std::make_unique<communicator::keyval<std::map<std::string, mpi3::any>>>();
	}
	explicit environment(thread_level required) {
		// std::clog << "ctor(thread_level) environment" << std::endl;
		initialize_thread(required);
		named_attributes_key_f() = std::make_unique<communicator::keyval<std::map<std::string, mpi3::any>>>();
	}
	explicit environment(int& argc, char**& argv) {  // cppcheck-suppress constParameter ; bug in cppcheck 2.3 or it can't see through the MPI C-API
		initialize(argc, argv);  // initialize(argc, argv); // TODO have an environment_mt/st version?
		named_attributes_key_f() = std::make_unique<communicator::keyval<std::map<std::string, mpi3::any>>>();
	}
	explicit environment(int& argc, char**& argv, thread_level required) {  // cppcheck-suppress constParameter ; bug in cppcheck 2.3
		initialize(argc, argv, required);
		named_attributes_key_f() = std::make_unique<communicator::keyval<std::map<std::string, mpi3::any>>>();
	}

	environment(environment const&) = delete;
	environment(environment&&)      = delete;

	environment& operator=(environment const&) = delete;
	environment& operator=(environment&&)      = delete;

	~environment() noexcept {  // NOLINT(bugprone-exception-escape) finalizes throws as an instance of UB
		// std::clog << "dtor environment" << std::endl;
		named_attributes_key_f().reset();
		finalize();  // cppcheck-suppress throwInNoexceptFunction ; finalizes throws as an instance of UB
	}

	inline static thread_level thread_support() { return mpi3::thread_support(); }
	//	static /*inline*/ communicator::keyval<int> const* color_key_p;
	//	static communicator::keyval<int> const& color_key(){return *color_key_p;}
	//	static /*inline*/ communicator::keyval<std::map<std::string, mpi3::any>> const* named_attributes_key_p;
	static std::unique_ptr<communicator::keyval<std::map<std::string, mpi3::any>> const>& named_attributes_key_f() {
		static std::unique_ptr<communicator::keyval<std::map<std::string, mpi3::any>> const> named_attributes_key_p;
		return named_attributes_key_p;
	}
	static communicator::keyval<std::map<std::string, mpi3::any>> const& named_attributes_key() {
		//	static communicator::keyval<std::map<std::string, mpi3::any>> const named_attributes_key_p;
		//	return named_attributes_key_p;
		return *named_attributes_key_f();
	}

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
			//	MPI_Comm_create_errhandler(&throw_error_fn, &throw_error_);
			//	MPI_Comm_set_errhandler(MPI_COMM_WORLD, throw_error_);
			MPI_Comm_set_errhandler(MPI_COMM_SELF, MPI_ERRORS_RETURN);
			//	MPI_Comm_set_errhandler(MPI_COMM_WORLD, MPI_ERRORS_ARE_FATAL);
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
			//	MPI_Comm_create_errhandler(&throw_error_fn, &throw_error_);
			//	MPI_Comm_set_errhandler(MPI_COMM_WORLD, throw_error_);
			MPI_Comm_set_errhandler(MPI_COMM_WORLD, MPI_ERRORS_RETURN);
			//	MPI_Comm_set_errhandler(MPI_COMM_WORLD, MPI_ERRORS_ARE_FATAL);
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

inline mpi3::any& communicator::attribute(std::string const& s) {
	return attribute(environment::named_attributes_key())[s];
}

}  // end namespace mpi3
}  // end namespace boost

//#if not __INCLUDE_LEVEL__ // _TEST_BOOST_MPI3_ENVIRONMENT

//#include<thread> // this_tread::sleep_for
//#include<chrono>

// namespace mpi3 = boost::mpi3;
// using std::cout;
// using namespace std::chrono_literals; // 2s

// int main(){//int argc, char* argv[]){
//	mpi3::environment::initialize(mpi3::thread_level::multiple);//argc, argv); // same as MPI_Init(...);
//	assert( mpi3::environment::thread_support() == mpi3::thread_level::multiple );
//	assert(mpi3::environment::is_initialized());
//	{
//		mpi3::communicator world = mpi3::environment::get_world_instance(); // a copy
//		auto then = mpi3::environment::wall_time();
//		std::this_thread::sleep_for(2s);
//	//	cout<< (mpi3::environment::wall_time() - then).count() <<" seconds\n";
//	}
//	mpi3::environment::finalize(); // same as MPI_Finalize()
//	assert(mpi3::environment::is_finalized());
//// or better:
////	mpi3::environment env(argc, argv);
////	auto world = env.world();
//}
//#endif
#endif
