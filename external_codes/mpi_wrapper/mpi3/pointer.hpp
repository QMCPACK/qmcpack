#if COMPILATION_INSTRUCTIONS /* -*- indent-tabs-mode: t -*- */
(echo "#include\""$0"\"" > $0x.cpp) && mpic++ -O3 -std=c++14 -Wfatal-errors -D_TEST_BOOST_MPI3_POINTER $0x.cpp -o $0x.x && time mpirun -n 4 $0x.x $@ && rm -f $0x.x $0x.cpp; exit
#endif
#ifndef BOOST_MPI3_POINTER_HPP
#define BOOST_MPI3_POINTER_HPP

#include "../mpi3/window.hpp"
#include<memory>

namespace boost{
namespace mpi3{

template<> struct pointer<void>{
	std::shared_ptr<window<void>> winSP_;
};

template<class T>
struct pointer{
	std::shared_ptr<window<T>> winSP_;
	~pointer(){}
//	pointer() : winSP_(){}
//	pointer(pointer const& other) : winSP_(other.win_){}
//	pointer(pointer&& other) : win_(std::move(other.win_)){}
//	pointer(pointer&& other) : ptr_(other.ptr_), pimpl_(other.pimpl_){
//		other.pimpl_ = nullptr;
//	}
//	pointer& operator=(pointer&& other){
//		pimpl_ = other.pimpl_;
//		other.pimpl_ = nullptr;
//	}
//	pointer(pointer const& other) = delete;
//	pointer& operator=(pointer const& other) = delete;
//	~pointer(){if(pimpl_) delete pimpl_;}

	T* local_ptr() const{
		void* base;
		int flag;
		MPI_Win_get_attr(winSP_->impl_, MPI_WIN_BASE, &base, &flag);
		return static_cast<T*>(base);
	}
/*
	MPI_Aint local_size() const{
		MPI_Aint* size_p;
		int flag;
		int i = MPI_Win_get_attr(pimpl_->impl_, MPI_WIN_SIZE, &size_p, &flag);
		assert(i==0);
		assert(flag);
		return *size_p/local_disp_unit();
	}
	int const& local_disp_unit() const{
		int* disp_unit_p;
		int flag;
		int i = MPI_Win_get_attr(pimpl_->impl_, MPI_WIN_DISP_UNIT, &disp_unit_p, &flag);
		assert(i==0);
		assert(flag);
		return *disp_unit_p;
	}
*/
	void source(T& t) const{
		winSP_->lock_exclusive(0);
		winSP_->get_n(&t, sizeof(T), 0, 0*sizeof(T));
		winSP_->unlock(0);
	}
	reference<T> operator*() const{
		return {*this};
	}
};

template<class T>
struct reference{
	pointer<T> p_;
	T buffer_;
	reference const& operator=(T const& t) const{
		p_.winSP_->lock_exclusive(0);
		p_.winSP_->put_n(&t, sizeof(T), 0, 0*sizeof(T));
		p_.winSP_->unlock(0);
		return *this;
	}
	operator T(){
		p_.source(buffer_);
		return buffer_;
	}
};

#if 0
template<class T = void>
pointer<T> communicator::allocate(MPI_Aint size) const{
	pointer<T> ret;
	ret.pimpl_ = new window;
	int i = MPI_Win_allocate(
		size*sizeof(T), sizeof(T), MPI_INFO_NULL, impl_, 
		&ret.ptr_, //&static_cast<window&>(ret).impl_
		&ret.pimpl_->impl_
	);
	if(size == 0) ret.ptr_ = nullptr;
	if(i!=0) assert(0);
	return ret;
}
#endif

pointer<void> communicator::malloc(MPI_Aint size) const{
	pointer<void> ret;
	void* ignore; // ???
	ret.winSP_ = std::make_shared<window<>>();
	int i = MPI_Win_allocate(
		size, 1, MPI_INFO_NULL, impl_, 
		&ignore, //&ret.ptr_,
		&(ret.winSP_->operator&())
	);
//	if(size == 0) ret.ptr_ = nullptr; 
	if(i!=0) assert(0);
	return ret;
}

void communicator::free(pointer<void>& p) const{
//	p.win_.fence();
	MPI_Free_mem(p.winSP_->base()); //the window frees the memory
//	MPI_Win_free(&p.win_.impl_);
//	p.win_.impl_ = MPI_WIN_NULL;
}

template<class T>
struct pgas_allocator{
	communicator const& comm_;
	pgas_allocator(communicator const& comm) : comm_(comm){}
	pointer<T> allocate(std::size_t size){
		pointer<T> ret;
		void* ignore;
		int local_size = size/comm_.size() + (comm_.rank() < (size % comm_.size()))?1:0;
		ret.winSP_ = std::make_shared<window<T>>();
		int i = MPI_Win_allocate(
			local_size*sizeof(T), sizeof(T), 
			MPI_INFO_NULL, &comm_,
			&ignore, 
			&(ret.winSP_->operator&())
		);
		assert(i==0);
		return ret;
	}
	void deallocate(pointer<T>& p, std::size_t){
	//	p = pointer<T>();
	//	MPI_Free_mem(p.local_ptr());
	//	p.win_.impl_ = MPI_WIN_NULL;
	}
};


}}

#ifdef _TEST_BOOST_MPI3_POINTER
#include<iostream>

#include "../mpi3/main.hpp"
using std::cout;
using std::endl;

namespace mpi3 = boost::mpi3;

int mpi3::main(int, char*[], mpi3::communicator world){
	{
		auto p = world.malloc(world.rank()==0?100*sizeof(double):0);
		if(world.rank() == 1){
			double cinco = 5;
			p.winSP_->lock_exclusive(0);
			p.winSP_->put_n(&cinco, sizeof(double), 0, 11*sizeof(double));
			p.winSP_->unlock(0);
		}
		p.winSP_->fence();
		if(world.rank() == 0){
			cout << static_cast<double*>(p.winSP_->base())[11] << endl;
			cout << p.winSP_->size() << endl;
			cout << p.winSP_->disp_unit() << endl;
		}
		return 0;
		if(world.rank() == 1){
			double t;
			p.winSP_->lock_exclusive(0);
			p.winSP_->get_n(&t, sizeof(double), 0, 11*sizeof(double));
			p.winSP_->unlock(0);
			assert(t == 5.);
		}
		world.free(p);
	}
	return 0;
	if(1){
		boost::mpi3::pgas_allocator<double> alloc(world);
		boost::mpi3::pointer<double> p = alloc.allocate(1);
		if(world.rank() == 0) *p = 5.1;
		p.winSP_->fence();
		if(world.rank() == 0){if(*p == 5.1) cout << "ok\n"; else cout << "BAD\n";}
		if(world.rank() == 1){if(*p == 5.1) cout << "ok\n"; else cout << "BAD\n";}
		if(world.rank() == 2){if(*p == 5.1) cout << "ok\n"; else cout << "BAD\n";}
		p.winSP_->fence();
	//	alloc.deallocate(p, 20);
		MPI_Win_free(&(p.winSP_.get()->operator&()));
	}
	double r = 5.;
	cout <<"great\n";
	
	return 0;
}

#endif
#endif

