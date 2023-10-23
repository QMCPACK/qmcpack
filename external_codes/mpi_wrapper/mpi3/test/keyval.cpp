#if COMPILATION_INSTRUCTIONS
mpicxx -O3 -std=c++14 `#-Wfatal-errors` $0 -o $0x.x && time mpirun -np 8s $0x.x $@ && rm -f $0x.x; exit
#endif

#include "alf/boost/mpi3/environment.hpp"
#include "alf/boost/mpi3/communicator.hpp"

namespace boost{
namespace mpi3{
	struct keyval{
		template<class T>
		static int default_MPI_Comm_copy_attr_function(
			MPI_Comm oldcomm, int comm_keyval, void *extra_state, void *attribute_val_in, void *attribute_val_out, int *flag)
		{
		//	std::memcpy(attribute_val_in, attribute_val_out, (std::size_t)(*extra_state));
			attribute_val_out = new T(*reinterpret_cast<T*>(attribute_val_in));
			*flag = 1;
			return MPI_SUCCESS;
		}
		template<class T>
		static int default_MPI_Comm_delete_attr_function(MPI_Comm old, int comm_keyval, void *attribute_val, void *extra_state){
			delete reinterpret_cast<T*>(attribute_val);
			return MPI_SUCCESS;
		}
		std::vector<int> key_;
		std::vector<int> attrval_;
		keyval(int size) : key_(size), attrval_(size){}
		int size(){return key_.size();}
		~keyval(){
			for(int i = 0; i != size(); ++i){
				MPI_Comm_free_keyval(&key_[i]);
			}
		}
		template<class T>
		void set(int idx){
			MPI_Comm_create_keyval(
				default_MPI_Comm_copy_attr_function<T>,// MPI_NULL_COPY_FN, 
				default_MPI_Comm_delete_attr_function<T>,// MPI_NULL_DELETE_FN, 
				&key_[idx], (void*)0
			);
		}
		template<class T1, class T2>
		void set(){
			set<T1>(0);
			set<T2>(1);
		}
		template<class T1, class T2, class T3>
		void set(){
			set<T1>(0);
			set<T2>(1);
			set<T3>(2);
		}
	};


template<typename T>
void communicator::set_attribute(keyval const& kv, int idx, T const& val){
	int status = MPI_Comm_set_attr(impl_, kv.key_[idx], new T(val));
	if(status != MPI_SUCCESS) throw std::runtime_error("cannot set attribute");
}
void communicator::delete_attribute(keyval const& kv, int idx){
	int status = MPI_Comm_delete_attr(impl_, kv.key_[idx]);
}

template<typename T>
void communicator::get_attribute(keyval const& kv, int idx, T& val){
	int flag = 0;
	T* p = nullptr;
	int status = MPI_Comm_get_attr(impl_, kv.key_[idx], &p, &flag);
	assert(flag);
	val = *p;
}

bool communicator::has_attribute(keyval const& kv, int idx){
	int flag = 0;
	void* p = nullptr;
	int status = MPI_Comm_get_attr(impl_, kv.key_[idx], &p, &flag);
	return flag;
}

template<typename T>
T const& communicator::get_attribute_as(keyval const& kv, int idx){
	int flag = 0;
	T* p = nullptr;
	int status = MPI_Comm_get_attr(impl_, kv.key_[idx], &p, &flag);
	assert(flag);
	return *p;
}

}}


namespace mpi3 = boost::mpi3;
using std::cout;

int main(int argc, char* argv[]){

	mpi3::environment env(argc, argv);
	mpi3::communicator comm = env.world();

	mpi3::keyval kv(3);
	kv.set<std::string, double, int>();

	comm.set_attribute(kv, 0, std::string("hola")); 
	comm.set_attribute(kv, 1, 5.1); 
	comm.set_attribute(kv, 2, 8); 

	assert( comm.has_attribute(kv, 0) );
	std::string ret;
	comm.get_attribute(kv, 0, ret);
	assert(ret == "hola");
	assert( comm.get_attribute_as<std::string>(kv, 0) == "hola" );

	mpi3::communicator comm2 = comm;
	assert( comm2.has_attribute(kv, 0) );
	comm2.set_attribute(kv, 0, std::string("chau")); 
	assert( comm2.get_attribute_as<std::string>(kv, 0) == "chau" );

	assert( comm.get_attribute_as<std::string>(kv, 0) == "hola" );

	comm.delete_attribute(kv, 0);
	comm.delete_attribute(kv, 1);
	comm.delete_attribute(kv, 2);

	assert( not comm.has_attribute(kv, 0) );

	cout << "TAG_UB = " << comm.attribute_as<int>(MPI_TAG_UB) << std::endl;
	cout << "HOST = " << comm.attribute_as<int>(MPI_HOST) << std::endl;
	cout << "IO = " << comm.attribute_as<int>(MPI_IO) << std::endl;
	cout << "WTIME_IS_GLOBAL = " << comm.attribute_as<int>(MPI_WTIME_IS_GLOBAL) << std::endl;

}

