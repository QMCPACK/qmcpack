#if COMPILATION//-*- indent-tabs-mode:t;c-basic-offset:4;tab-width:4; -*-
$CXXX `mpicxx -showme:compile|sed 's/-pthread/ /g'` -std=c++14 $0 -o $0x `mpicxx -showme:link|sed 's/-pthread/ /g'`&&mpirun -n 4 $0x&&rm $0x;exit
#endif
// © Alfredo A. Correa 2018-2019

#ifndef BOOST_MPI3_ENVIRONMENT_HPP
#define BOOST_MPI3_ENVIRONMENT_HPP

#include "./core.hpp"
#include "./communicator.hpp"
#include "./wall_clock.hpp"
#include "./detail/call.hpp"

#include<mpi.h>

#include<string>

namespace boost{
namespace mpi3{

enum thread_level : int{
	single     = MPI_THREAD_SINGLE,
	funneled   = MPI_THREAD_FUNNELED,
	serialized = MPI_THREAD_SERIALIZED, 
	multiple   = MPI_THREAD_MULTIPLE
};

inline void finalize(){
	std::set_terminate(&std::abort);
	int s = MPI_Finalize();
	if(s != MPI_SUCCESS) throw std::runtime_error{"cannot finalize"};
}
inline void myterminate(){
	std::cerr << "myterminate handler called" << '\n';
	finalize();
	std::abort();
//	exit(1);  // forces abnormal termination
}

inline void initialize(int& argc, char**& argv){
	int s = MPI_Init(&argc, &argv);
	if(s != MPI_SUCCESS) throw std::runtime_error{"cannot initialize"};
//	std::set_terminate(&finalize);
//	std::set_terminate(myterminate);
}

inline thread_level initialize(int& argc, char**& argv, thread_level required){
	int provided;
	int s = MPI_Init_thread(&argc, &argv, static_cast<int>(required), &provided);
	if(s != MPI_SUCCESS) throw std::runtime_error{"cannot thread-initialize"};
	return static_cast<thread_level>(provided);
}

inline thread_level initialize_thread(thread_level required){
	int provided;
	int s = MPI_Init_thread(nullptr, nullptr, static_cast<int>(required), &provided);
	if(s != MPI_SUCCESS) throw std::runtime_error{"cannot thread-initialize"};
	return static_cast<thread_level>(provided);
}

inline thread_level initialize(thread_level required = thread_level::single){
	return initialize_thread(required);
}

inline void throw_error_fn(MPI_Comm* comm, int* errorcode, ...){
	char name[MPI_MAX_OBJECT_NAME];
	int rlen;
	int status = MPI_Comm_get_name(*comm, name, &rlen);
	if(status != MPI_SUCCESS) throw std::runtime_error{"cannot get name"};
	std::string sname(name, rlen);
	throw std::runtime_error{"error code "+ std::to_string(*errorcode) +" from comm "+ sname};
}

inline thread_level initialize_thread(
	int& argc, char**& argv, thread_level required
){
	int ret;
	int status = MPI_Init_thread(&argc, &argv, static_cast<int>(required), &ret);
	if(status != MPI_SUCCESS) throw std::runtime_error("cannot thread-initialize");
	return static_cast<thread_level>(ret);
}

inline thread_level thread_support(){
	int r; MPI_(Query_thread)(&r); return static_cast<thread_level>(r);
}

inline bool is_thread_main(){
	int flag = -1;
	int s = MPI_Is_thread_main(&flag);
	if(s != MPI_SUCCESS) throw std::runtime_error{"cannot determine is thread main"};
	return flag;
}

inline std::string processor_name(){return detail::call<&MPI_Get_processor_name>();}

inline std::string get_processor_name(){return detail::call<&MPI_Get_processor_name>();}

class environment{
	public:
	environment(){
		initialize_thread(thread_level::multiple);
		named_attributes_key_f() = new communicator::keyval<std::map<std::string, mpi3::any>>;
	}
	environment(thread_level required){
		initialize_thread(required);
		named_attributes_key_f() = new communicator::keyval<std::map<std::string, mpi3::any>>;
	}
	environment(int& argc, char**& argv){
		initialize(argc, argv); // initialize(argc, argv); // TODO have an environment_mt/st version?
		named_attributes_key_f() = new communicator::keyval<std::map<std::string, mpi3::any>>;
	}
	environment(int& argc, char**& argv, thread_level required){
		initialize(argc, argv, required);
		named_attributes_key_f() = new communicator::keyval<std::map<std::string, mpi3::any>>;		
	}
	environment(environment const&) = delete;
	environment& operator=(environment const&) = delete;
	~environment(){
		delete named_attributes_key_f();
		finalize();
	}
	inline static thread_level thread_support(){return mpi3::thread_support();}
//	static /*inline*/ communicator::keyval<int> const* color_key_p;
//	static communicator::keyval<int> const& color_key(){return *color_key_p;}
//	static /*inline*/ communicator::keyval<std::map<std::string, mpi3::any>> const* named_attributes_key_p;
	static communicator::keyval<std::map<std::string, mpi3::any>> const*& named_attributes_key_f(){
		static communicator::keyval<std::map<std::string, mpi3::any>> const* named_attributes_key_p;
		return named_attributes_key_p;
	}
	static communicator::keyval<std::map<std::string, mpi3::any>> const& named_attributes_key(){
	//	static communicator::keyval<std::map<std::string, mpi3::any>> const named_attributes_key_p;
	//	return named_attributes_key_p;
		return *named_attributes_key_f();
	}

	static inline void initialize(){mpi3::initialize();}
	static inline void initialize(thread_level required){mpi3::initialize_thread(required);}
	static inline void initialize(int argc, char** argv){mpi3::initialize(argc, argv);}
	static inline void initialize(int argc, char** argv, thread_level req){mpi3::initialize_thread(argc, argv, req);}

	static inline void finalize(){mpi3::finalize();}

	static inline bool is_initialized(){return mpi3::initialized();}
	static inline bool is_finalized(){return mpi3::finalized();}
	using wall_clock = mpi3::wall_clock;
	operator bool() const{return initialized();}
	bool is_thread_main() const{return mpi3::is_thread_main();}

	static inline communicator& get_self_instance(){
		assert(initialized());
		static communicator instance = []{
		//	MPI_Comm_create_errhandler(&throw_error_fn, &throw_error_);
		//	MPI_Comm_set_errhandler(MPI_COMM_WORLD, throw_error_);
			MPI_Comm_set_errhandler(MPI_COMM_SELF, MPI_ERRORS_RETURN);
		//	MPI_Comm_set_errhandler(MPI_COMM_WORLD, MPI_ERRORS_ARE_FATAL);
			return communicator{MPI_COMM_SELF};
		}();
		return instance;
	}

	static communicator self(){ // returns a copy!
		MPI_Comm_set_errhandler(MPI_COMM_SELF, MPI_ERRORS_RETURN);
		return communicator{MPI_COMM_SELF};
	}
	static inline communicator& get_world_instance(){
		assert(initialized());
		static communicator instance = []{
		//	MPI_Comm_create_errhandler(&throw_error_fn, &throw_error_);
		//	MPI_Comm_set_errhandler(MPI_COMM_WORLD, throw_error_);
			MPI_Comm_set_errhandler(MPI_COMM_WORLD, MPI_ERRORS_RETURN);
		//	MPI_Comm_set_errhandler(MPI_COMM_WORLD, MPI_ERRORS_ARE_FATAL);
			return communicator{MPI_COMM_WORLD};
		}();
		return instance;
	}
	communicator world(){// const{ // returns a copy!
		communicator ret{get_world_instance()}; ret.name("world");
		return ret;
	}
	std::string processor_name() const{return get_processor_name();}
	std::string get_processor_name() const{return mpi3::get_processor_name();}

	inline static auto wall_time(){return mpi3::wall_time();}
	inline static auto wall_tick(){return mpi3::wall_tick();}
	template<class Duration = wall_clock::duration>
	static auto wall_sleep_for(Duration d){return mpi3::wall_sleep_for(d);}
};

inline mpi3::any& communicator::attribute(std::string const& s){
	return attribute(environment::named_attributes_key())[s];
}

}}

#if not __INCLUDE_LEVEL__ // _TEST_BOOST_MPI3_ENVIRONMENT

#include<thread> // this_tread::sleep_for
#include<chrono>

namespace mpi3 = boost::mpi3;
using std::cout;
using namespace std::chrono_literals; // 2s

int main(){//int argc, char* argv[]){
	mpi3::environment::initialize(mpi3::thread_level::multiple);//argc, argv); // same as MPI_Init(...);
	assert( mpi3::environment::thread_support() == mpi3::thread_level::multiple );
	assert(mpi3::environment::is_initialized());
	{
		mpi3::communicator world = mpi3::environment::get_world_instance(); // a copy
		auto then = mpi3::environment::wall_time();
		std::this_thread::sleep_for(2s);
	//	cout<< (mpi3::environment::wall_time() - then).count() <<" seconds\n"; 
	}
	mpi3::environment::finalize(); // same as MPI_Finalize()
	assert(mpi3::environment::is_finalized());
// or better:
//	mpi3::environment env(argc, argv);
//	auto world = env.world();
}
#endif
#endif

