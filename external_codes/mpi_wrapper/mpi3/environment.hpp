#if COMPILATION_INSTRUCTIONS
(echo "#include\""$0"\"" > $0x.cpp) && mpic++ -Wall -Wextra -O3 -std=c++14 -Wfatal-errors -D_TEST_BOOST_MPI3_ENVIRONMENT $0x.cpp -o $0x.x && time mpirun -n 4 $0x.x ddd $@ && rm -f $0x.x $0x.cpp; exit
#endif
//  (C) Copyright Alfredo A. Correa 2018.
#ifndef BOOST_MPI3_ENVIRONMENT_HPP
#define BOOST_MPI3_ENVIRONMENT_HPP

#include "../mpi3/communicator.hpp"
#include "../mpi3/wall_clock.hpp"
#include "../mpi3/detail/call.hpp"

#include<mpi.h>

#include<string>

namespace boost{
namespace mpi3{

struct thread_level{
//	decltype(MPI_THREAD_SINGLE)
	int impl_;
	bool operator==(thread_level const& other) const{
		return impl_ == other.impl_;
	}
};

/*[[maybe_unused]]*/ static thread_level single{MPI_THREAD_SINGLE};
/*[[maybe_unused]]*/ static thread_level funneled{MPI_THREAD_FUNNELED};
/*[[maybe_unused]]*/ static thread_level serialized{MPI_THREAD_SERIALIZED};
/*[[maybe_unused]]*/ static thread_level multiple{MPI_THREAD_MULTIPLE};

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
inline void initialize(){
	int s = MPI_Init(nullptr, nullptr);
	if(s != MPI_SUCCESS) throw std::runtime_error{"cannot initialize"};
	std::set_terminate(&finalize);
}
inline void initialize(int& argc, char**& argv){
	int s = MPI_Init(&argc, &argv);
	if(s != MPI_SUCCESS) throw std::runtime_error{"cannot initialize"};
	std::set_terminate(&finalize);
//	std::set_terminate(myterminate);
}
inline thread_level initialize_thread(thread_level required){
	int provided;
	int s = MPI_Init_thread(nullptr, nullptr, required.impl_, &provided);
	if(s != MPI_SUCCESS) throw std::runtime_error{"cannot thread-initialize"};
	return {provided};
}
inline thread_level initialize(thread_level required){
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
	(void)single;
	(void)funneled;
	(void)serialized;
	(void)multiple;
	thread_level ret;
	int status = MPI_Init_thread(&argc, &argv, required.impl_, reinterpret_cast<int*>(&ret));
	if(status != MPI_SUCCESS) throw std::runtime_error("cannot thread-initialize");
	return ret;
}
inline thread_level initialize(int& argc, char**& argv, thread_level required){
	return initialize_thread(argc, argv, required);
}
inline bool initialized(){
	int flag = -1; 
	int s = MPI_Initialized(&flag); 
	if(s != MPI_SUCCESS) throw std::runtime_error{"cannot probe initialization"};
	return flag;
}
inline bool finalized(){
	int flag = -1; 
	int s = MPI_Finalized(&flag);
	if(s != MPI_SUCCESS) throw std::runtime_error{"cannot probe finalization"};
	return flag;
}
inline bool is_thread_main(){
	int flag = -1;
	int s = MPI_Is_thread_main(&flag);
	if(s != MPI_SUCCESS) throw std::runtime_error{"cannot determine is thread main"};
	return flag;
}
inline thread_level query_thread(){
	int ret;
	int status = MPI_Query_thread(&ret);
	if(status != MPI_SUCCESS) throw std::runtime_error("cannot query thread level");
	return {ret};
}

inline std::string processor_name(){return detail::call<&MPI_Get_processor_name>();}
inline std::string get_processor_name(){return detail::call<&MPI_Get_processor_name>();}

class environment{
	public:
	environment(){
		initialize();
	//	color_key_p = new communicator::keyval<int>;
		named_attributes_key_f() = new communicator::keyval<std::map<std::string, mpi3::any>>;
	}
	environment(thread_level required){initialize_thread(required);}
	environment(int argc, char** argv){
		initialize(argc, argv);
	//	color_key_p = new communicator::keyval<int>;
		named_attributes_key_f() = new communicator::keyval<std::map<std::string, mpi3::any>>;
	}
	environment(int argc, char** argv, thread_level required){initialize_thread(argc, argv, required);}
	environment(environment const&) = delete;
	environment& operator=(environment const&) = delete;
	~environment(){
		delete named_attributes_key_f();
	//	delete color_key_p;
		finalize();
	}
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
	thread_level query_thread() const{return mpi3::query_thread();}

//	communicator& null() const{return mpi3::communicator::null;}
	static communicator self(){ // returns a copy!
		MPI_Comm_set_errhandler(MPI_COMM_SELF, MPI_ERRORS_RETURN);
		return communicator{MPI_COMM_SELF};
	}
	static inline communicator& get_world_instance(){
		assert(initialized());
	//	static communicator instance{MPI_COMM_WORLD};
		static communicator instance = []{
		//	MPI_Comm_create_errhandler(&throw_error_fn, &throw_error_);
		//	MPI_Comm_set_errhandler(MPI_COMM_WORLD, throw_error_);
			MPI_Comm_set_errhandler(MPI_COMM_WORLD, MPI_ERRORS_RETURN);
		//	MPI_Comm_set_errhandler(MPI_COMM_WORLD, MPI_ERRORS_ARE_FATAL);
			return communicator{MPI_COMM_WORLD};
		}();
		//	MPI_Comm_create_errhandler(&throw_error_fn, &throw_error_);
		//	MPI_Comm_set_errhandler(MPI_COMM_NULL, MPI_ERRORS_RETURN);
		return instance;
	}
	communicator world() const{ // returns a copy!
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

#ifdef _TEST_BOOST_MPI3_ENVIRONMENT
#include<iostream>
#include <thread> // this_tread::sleep_for
#include<chrono>

namespace mpi3 = boost::mpi3;
using std::cout;
using namespace std::chrono_literals; // 2s

int main(int argc, char* argv[]){
	mpi3::environment::initialize(argc, argv); // same as MPI_Init(...);
	assert(mpi3::environment::is_initialized());
	{
		mpi3::communicator world = mpi3::environment::get_world_instance(); // a copy
		auto then = mpi3::environment::wall_time();
		std::this_thread::sleep_for(2s);
		cout<< (mpi3::environment::wall_time() - then).count() <<" seconds\n"; 
	}
	mpi3::environment::finalize(); // same as MPI_Finalize()
	assert(mpi3::environment::is_finalized());
// or better:
//	mpi3::environment env(argc, argv);
//	auto world = env.world();
	return 0;
}
#endif
#endif

