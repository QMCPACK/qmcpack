#if COMPILATION_INSTRUCTIONS /* -*- indent-tabs-mode: t -*- */
(echo "#include\""$0"\"" > $0x.cpp) && mpicxx -O3 -std=c++14 -Wfatal-errors -D_TEST_BOOST_MPI3_DETAIL_HANDLE $0x.cpp -o $0x.x && time mpirun -np 4 $0x.x $@ && rm -f $0x.cpp; exit
#endif

#ifndef BOOST_MPI3_DETAIL_HANDLE_HPP
#define BOOST_MPI3_DETAIL_HANDLE_HPP

//#include<boost/exception/to_string.hpp>

#include<cassert>
#include<stdexcept> // runtime_error
#include<string>

#include<mpi.h>

namespace boost{
namespace mpi3{
namespace detail{

template<class Self, class Impl>
struct caller{
	Impl& impl(){return static_cast<Self&>(*this).impl_;}
	Impl const& impl() const{return static_cast<Self const&>(*this).impl_;}

	template<class F, class... Args>
	void static_call(F f, Args&&... args){
		int status = f(std::forward<Args>(args)...);
		if(status != 0) throw std::runtime_error{"error "+ std::to_string(status)};
	}
	template<int(*F)(Impl, char const*, char const*)> void call(
		char const* c1, char const* c2
	){
		int status = F(impl(), c1, c2);
		if(status != 0) throw std::runtime_error{"error "+ std::to_string(status)};
	}
	template<int(*F)(Impl, char const*, char const*)> void call(
		std::string const& s1, std::string const& s2
	){
		return call<F>(s1.c_str(), s2.c_str());
	}
	template<int(*F)(Impl, int*)> int call() const{
		int ret = -1;
	//	static_call(F, impl_, &ret);
		int status = F(impl(), &ret);
		if(status != MPI_SUCCESS) throw std::runtime_error{"error " + std::to_string(status)};
		return ret;
	}
	template<int(*F)(Impl, int, char*)> std::string call(int n) const{
		char ret[MPI_MAX_INFO_KEY];
		int status = F(impl(), n, ret);
		if(status != 0) throw std::runtime_error{"error "+ std::to_string(status)};
		return ret;
	}
	template<int(*F)(Impl, char const*, int*, int*)> std::pair<int, int> call(std::string const& key) const{
		int flag;
		int valuelen;
		int status = F(impl(), key.c_str(), &valuelen, &flag);
		if(status != 0) throw std::runtime_error{"error "+ std::to_string(status)};
		return {valuelen, flag};
	}
	template<int(*F)(Impl, char const*, int, char*, int*)> std::pair<std::string, int> call(std::string const& key, int valuelen) const{
		char value[MPI_MAX_INFO_VAL];
		int flag;
		int status = F(impl(), key.c_str(), valuelen, value, &flag);
		if(status != 0) throw std::runtime_error{"error "+ std::to_string(status)};
		return {std::string(value), flag};
	}
	template<int(*F)(Impl, char const*)> void call(std::string const& key) const{
		int status = F(impl(), key.c_str());
		if(status != 0) throw std::runtime_error("error "+ std::to_string(status));
	}
};

template<
	class Self, 
	class Impl, 
	int(*CreateFunction)(Impl*), 
	int(*DupFunction)(Impl, Impl*), 
	int(*FreeFunction)(Impl*)
>
struct regular_handle : caller<regular_handle<Self, Impl, CreateFunction, DupFunction, FreeFunction>, Impl>{
	using caller<regular_handle<Self, Impl, CreateFunction, DupFunction, FreeFunction>, Impl>::call;
	using impl_t = Impl;
	impl_t impl_;
	regular_handle(){CreateFunction(&impl_);}
	regular_handle(Self const& other){
		int status = DupFunction(other.impl_, &impl_);
		if(status != MPI_SUCCESS) throw std::runtime_error("cannot copy handle");
	}
	~regular_handle(){
		assert(impl_ != MPI_INFO_NULL);
		if(impl_ != MPI_INFO_NULL) FreeFunction(&impl_);
	}
	void swap(Self& other){std::swap(impl_, other.impl_);}
	Self& operator=(Self const& other){
		Self tmp = other;
		swap(tmp);
		return static_cast<Self&>(*this);
	}
};

template<class Self, class Impl, int(*CreateFunction)(Impl*), int(*FreeFunction)(Impl*)>
struct noncopyable_handle : caller<noncopyable_handle<Self, Impl, CreateFunction, FreeFunction>, Impl>{
	using impl_t = Impl;
	impl_t impl_;
	bool predefined_ = false;
	noncopyable_handle(Impl code) : impl_(code), predefined_(true){}
	noncopyable_handle(){CreateFunction(&impl_);}
	noncopyable_handle(noncopyable_handle const&) = delete;
	~noncopyable_handle(){
		assert(impl_ != MPI_INFO_NULL);
	//	if(impl_ != MPI_INFO_NULL) 
		if(not predefined_) FreeFunction(&impl_);
	}
	void swap(noncopyable_handle& other){std::swap(impl_, other.impl_);}
	Self& operator=(noncopyable_handle const& other) = delete;
};

struct uninitialized{};

template<class Self, class Impl, int(*FreeFunction)(Impl*)>
struct nondefault_handle : caller<nondefault_handle<Self, Impl, FreeFunction>, Impl>{
	using impl_t = Impl;
	impl_t impl_;
	bool predefined_ = false;
	nondefault_handle(Impl code) : impl_(code), predefined_(true){}
	nondefault_handle() = delete;
	nondefault_handle(uninitialized){};
	nondefault_handle(nondefault_handle const&) = delete;
	~nondefault_handle(){if(not predefined_) FreeFunction(&impl_);}
	void swap(nondefault_handle& other){std::swap(impl_, other.impl_);}
	Self& operator=(nondefault_handle const& other) = delete;
};

template<class Self, class Impl>
struct persistent_handle : caller<persistent_handle<Self, Impl>, Impl>{
	using impl_t = Impl;
	impl_t impl_;
	persistent_handle() = delete;
	persistent_handle(uninitialized){};
	persistent_handle(Impl code) : impl_(code){}
	persistent_handle(persistent_handle const&) = delete;
	~persistent_handle() = default;
};

}}}

#ifdef _TEST_BOOST_MPI3_DETAIL_HANDLE
int main(){
}
#endif
#endif

