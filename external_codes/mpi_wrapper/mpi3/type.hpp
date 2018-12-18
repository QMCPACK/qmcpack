#if COMPILATION_INSTRUCTIONS
(echo "#include\""$0"\"" > $0x.cpp) && mpic++ -O3 -std=c++14 -Wfatal-errors -D_TEST_BOOST_MPI3_TYPE $0x.cpp -o $0x.x -lboost_serialization && time mpirun -n 2 $0x.x $@ && rm -f $0x.x $0x.cpp; exit
#endif
//  (C) Copyright Alfredo A. Correa 2018.

#ifndef BOOST_MPI3_TYPE_HPP
#define BOOST_MPI3_TYPE_HPP

#define OMPI_SKIP_MPICXX 1  // https://github.com/open-mpi/ompi/issues/5157
#include<mpi.h>
//#include "../mpi3/communicator.hpp"
//#include "../mpi3/detail/datatype.hpp"
#include<typeindex>
#include<map>
#include "detail/datatype.hpp"
#include<cassert>
#include<iostream>
#include<vector>

namespace boost{
namespace mpi3{

struct type{
private:
public:
	type(MPI_Datatype dt) : impl_(dt){}
	type() = default;// delete; //: impl_(MPI_DATATYPE_NULL){}
	MPI_Datatype impl_ = MPI_DATATYPE_NULL;
	type(type const& other){MPI_Type_dup(other.impl_, &impl_);}
	type& operator=(type const& other){
		type tmp(other);
		swap(tmp);
		return *this;
	}
	void swap(type& other){std::swap(impl_, other.impl_);}
	void commit() const{MPI_Type_commit(const_cast<MPI_Datatype*>(&impl_));}
	template<class T> void commit_as() const{
		commit();
		registered.emplace(typeid(T), *this);
	}
	template<class T> void commit_as(T const&){return commit_as<T>();}
	~type(){
	//	if(impl_!=MPI_DATATYPE_NULL) MPI_Type_free(&impl_);
	}
	type contiguous(int count) const{
		type ret;
		int s = MPI_Type_contiguous(count, impl_, &ret.impl_);
		if(s != MPI_SUCCESS) throw std::runtime_error("cannot build contiguous type");
		ret.set_name("(" + name() + ")[" + std::to_string(count) + "]");
		return ret;
	}
	type vector(int count, int block_length, int stride) const{
		type ret;
		MPI_Type_vector(count, block_length, stride, impl_, &ret.impl_);
		using std::to_string;
		ret.set_name("("+ name() +")["+to_string(count)+","+to_string(block_length)+":"+to_string(stride)+"]");
		return ret;
	}
	// MPI_Type_struct is deprecated
	static type struct_(std::initializer_list<type> il){
		type ret;
		std::vector<int> blocklen(il.size(), 1);
		std::vector<MPI_Aint> disp;
		disp.reserve(il.size());
	//	std::vector<MPI_Datatype> array_of_types; 
	//	array_of_types.reserve(il.size());
		MPI_Aint current_disp = 0;
		std::string new_name = "{";
		for(auto& e : il){
			disp.push_back(current_disp);
			current_disp += e.size();
			new_name+=(&e!=il.begin()?", ":"") + e.name();
	//		array_of_types.push_back(e.impl_);
		}
		MPI_Type_create_struct(
			il.size(), 
			blocklen.data(), 
			disp.data(), 
			&il.begin()->impl_,
	//		array_of_types.data(), 
			&ret.impl_
		);
		ret.name(new_name);
		return ret;
	}
//  this is giving a problem with intel compiler, commenting for the moment
//	template<class T, class... Ts>
//	static auto struct_(T&& t, Ts&&... ts) -> decltype(struct_({std::forward<T>(t), std::forward<Ts>(ts)...})){
//		return struct_({std::forward<T>(t), std::forward<Ts>(ts)...});
//	}

	type operator[](int count) const{return contiguous(count);}
	type operator()(int stride) const{return vector(1, 1, stride);}
	MPI_Aint size() const{
		int ret;
		MPI_Type_size(impl_, &ret);
		return ret;
	}
	std::string name() const{
		char name[MPI_MAX_OBJECT_NAME];
		int namelen;
		MPI_Type_get_name(impl_, name, &namelen);
		return std::string(name, namelen);
	}
	void name(std::string s){set_name(s);}
	void set_name(std::string s){MPI_Type_set_name(impl_, s.c_str());}

//	int extent() const; // MPI_Type_get_extent
	MPI_Aint extent() const{
		MPI_Aint lb, ext;
		int s = MPI_Type_get_extent(impl_, &lb, &ext);
		if(s != MPI_SUCCESS) throw std::runtime_error{"cannot extent"};
		return ext;
	}
	MPI_Aint lower_bound() const{
		MPI_Aint lb, ext;
		int s = MPI_Type_get_extent(impl_, &lb, &ext);
		if(s != MPI_SUCCESS) throw std::runtime_error{"cannot lower bound"};
		return lb;
	}
	MPI_Aint upper_bound() const{
		MPI_Aint lb, ext;
		int s = MPI_Type_get_extent(impl_, &lb, &ext);
		if(s != MPI_SUCCESS) throw std::runtime_error{"cannot lower bound"};
		return lb + ext;
	}
/*	MPI_Aint upper_bound() const{ // this function is replaced by get_extent
		MPI_Aint ub;
		int s = MPI_Type_ub(impl_, &ub); // http://www.mpi.deino.net/mpi_functions/MPI_Type_ub.html
		if(s != MPI_SUCCESS) throw std::runtime_error("cannot upper bound");
		return ub;
	}*/
	type operator,(type const& other) const{
		type ret;
		int count = 2;
		int blocklen[2] = {1,1};
		MPI_Aint disp[2] = {0, this->size()};
		MPI_Datatype array_of_types[2] = {impl_, other.impl_};
		MPI_Type_create_struct(count, blocklen, disp, array_of_types, &ret.impl_);
		std::string newname = name() + ", " + other.name();
		MPI_Type_set_name(ret.impl_, newname.c_str());
		return ret;
	}
	static std::map<std::type_index, type const&> registered;

#if 0
	enum struct code : MPI_Datatype{
		char_ = MPI_CHAR, 
		unsigned_char = MPI_UNSIGNED_CHAR, unsigned_char_ = MPI_UNSIGNED_CHAR,
		byte = MPI_BYTE, byte_ = MPI_BYTE,
	//	wchar_ = MPI_WCHAR_T,
		short_ = MPI_SHORT, 
		unsigned_short = MPI_UNSIGNED_SHORT, unsigned_short_ = MPI_UNSIGNED_SHORT,
		int_ = MPI_INT, 
		unsigned_ = MPI_UNSIGNED, unsigned_int = MPI_UNSIGNED, unsigned_int_ = MPI_UNSIGNED,
		long_ = MPI_LONG,
		unsigned_long = MPI_UNSIGNED_LONG, unsigned_long_ = MPI_UNSIGNED_LONG,
		float_ = MPI_FLOAT, 
		double_ = MPI_DOUBLE,
		long_double = MPI_LONG_DOUBLE, long_double_ = MPI_LONG_DOUBLE,
		long_long_int = MPI_LONG_LONG_INT, long_long_int_ = MPI_LONG_LONG_INT, 
		long_long = MPI_LONG_LONG, long_long_ = MPI_LONG_LONG, 
		float_int = MPI_FLOAT_INT, float_int_ = MPI_FLOAT_INT,
		double_int = MPI_DOUBLE_INT, double_int_ = MPI_DOUBLE_INT,
		long_int = MPI_LONG_INT, long_int_ = MPI_LONG_INT,
		_2int = MPI_2INT, int2 = MPI_2INT, int2_ = MPI_2INT, _2int_ = MPI_2INT,
		long_double_int = MPI_LONG_DOUBLE_INT, long_double_int_ = MPI_LONG_DOUBLE_INT,

		lb = MPI_LB, ub = MPI_UB
	};
#endif
};

static type const char_{MPI_CHAR};
static type const unsigned_char{MPI_UNSIGNED_CHAR}; static type const& unsigned_char_ = unsigned_char;
static type const short_{MPI_SHORT};
static type const unsigned_short{MPI_UNSIGNED_SHORT}; static type const& unsigned_short_ = unsigned_short;
static type const int_{MPI_INT};
static type const unsigned_int_{MPI_UNSIGNED}; static type const& unsigned_int = unsigned_int_; static type const& unsigned_ = unsigned_int_;
static type const long_{MPI_LONG};
static type const unsigned_long_{MPI_UNSIGNED_LONG}; static type const& unsigned_long = unsigned_long_;
static type const float_{MPI_FLOAT};
static type const double_{MPI_DOUBLE};
static type const long_double_{MPI_LONG_DOUBLE}; static type const& long_double = long_double_;
static type const long_long_int{MPI_LONG_DOUBLE_INT};
static type const float_int{MPI_FLOAT_INT};
static type const long_int{MPI_LONG_INT};
static type const double_int{MPI_DOUBLE_INT};
static type const short_int{MPI_SHORT_INT};
static type const int_int{MPI_2INT}; static type const& _2int = int_int;
static type const long_double_int{MPI_LONG_DOUBLE_INT};

//static std::map<std::type_index, type const&> type::registered;

template<class T>
type make_type();

template<> inline type make_type<double>(){return mpi3::double_;}
template<> inline type make_type<int>(){return mpi3::int_;}

}}

#ifdef _TEST_BOOST_MPI3_TYPE

#include "../mpi3/main.hpp"
#include "../mpi3/communicator.hpp"
#include "../mpi3/process.hpp"

#include<iostream>

namespace mpi3 = boost::mpi3;
using std::cout;

struct A{
	double d[6];
	long l[20];
};

struct B{
	double d[7];
	long l[9];
};

template<class Archive>
void serialize(Archive& ar, A& a, const unsigned int){
	ar & boost::serialization::make_array(a.d, 6);
	ar & boost::serialization::make_array(a.l, 20);
}

int mpi3::main(int argc, char* argv[], mpi3::communicator world){

	assert(world.size() >= 2);
	{
		int value = -1;
		if(world.rank() == 0){
			value = 5;
			world[1] << value; //	world.send_value(value, 1);
		}else{
			assert(value == -1);
			world[0] >> value; //world.receive(&value, &value + 1, 0);
			assert(value == 5);
		}
	}
	{
		int buffer[100]; std::fill(&buffer[0], &buffer[100], 0);
		if(world.rank() == 0){
			std::iota(&buffer[0], &buffer[100], 0);
			world[1] << buffer; // world.send_value(buffer, 1);
		}else{
			assert(buffer[11]==0);
			world[0] >> buffer; //	world.receive_value(buffer, 0);
		//	assert(buffer[11]==11);
		}
	}
	{
	//	auto Atype = (
	//		mpi3::double_[6], 
	//		mpi3::long_[20]
	//	);
	//	Atype.commit_as<A>();

		A particle;
		particle.d[2] = 0.;
		if(world.rank()==0){
			particle.d[2] = 5.1;
			world[1] << particle;
		}else{
			assert(particle.d[2]==0.);
			world[0] >> particle;
		}
	}
/*	{
		auto Btype = (
			mpi3::double_[6], 
			mpi3::long_[20]
		);
		Btype.commit_as<B>();
		B b;
		b.d[2] = 0.;
		if(world.rank()==0){
			b.d[2] = 5.1;
			world[1] << b;
		}else{
			assert(b.d[2]==0.);
			world[0] >> b;
		}
	}*/
	return 0;
#if 0
	{
		std::vector<double> v(100);
		mpi3::type d100 = mpi3::type::double_[100];
		d100.commit();
		if(world.rank()==0){
			v[5] = 6.;
		//	world.send_n(v.data(), 1, d100);
		}else{
		//	world.receive_n(v.data(), 1, d100);
			assert(v[5] == 6.);
		}
	}
#endif
	return 0;
}

#endif
#endif

