#if COMPILATION// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;-*-
$CXXX $CXXFLAGS $0 -o $0x  -lfftw3 -lfftw3_mpi&&$0x&&rm $0x;exit
#endif
// © Alfredo A. Correa 2020
// apt-get install libfftw3-mpi-dev
// compile with: mpicc simple_mpi_example.c  -Wl,-rpath=/usr/local/lib -lfftw3_mpi -lfftw3 -o simple_mpi_example */

#ifndef MULTI_ADAPTOR_FFTW_MEMORY_HPP
#define MULTI_ADAPTOR_FFTW_MEMORY_HPP

#include <fftw3.h>

#include "../../config/NODISCARD.hpp"

#include<cassert>
#include<cstddef>
#include<complex>
#include<limits>
#include<memory>
#include<type_traits>

namespace boost{
namespace multi{
namespace fftw{

template<class T>
class allocator{
public:
	using value_type    = T;

#if 1
	using pointer       = value_type*;
	using const_pointer = typename std::pointer_traits<pointer>::template
                                                     rebind<value_type const>;
	using void_pointer       = typename std::pointer_traits<pointer>::template
                                                           rebind<void>;
	using const_void_pointer = typename std::pointer_traits<pointer>::template
                                                           rebind<const void>;
	using difference_type = typename std::pointer_traits<pointer>::difference_type;
	using size_type       = std::make_unsigned_t<difference_type>;

	template <class U> struct rebind {typedef allocator<U> other;};
#endif

	allocator() noexcept {}  // not required, unless used
	template <class U> allocator(allocator<U> const&) noexcept {}

	NODISCARD("to avoid memory leak") 
	value_type* allocate(std::size_t n) const{return static_cast<value_type*>(fftw_malloc(sizeof(T)*n));}

//	value_type*  // Use pointer if pointer is not a value_type*
//	allocate(std::size_t n){return static_cast<value_type*>(::operator new (n*sizeof(value_type)));}

	void deallocate(value_type* p, std::size_t){fftw_free(p);}
//	void deallocate(value_type* p, std::size_t) noexcept  // Use pointer if pointer is not a value_type*
//	{::operator delete(p);}

	static int alignment_of(value_type* p){return fftw_alignment_of((double*)p);}

#if 1
	value_type* allocate(std::size_t n, const_void_pointer){return allocate(n);}

	template <class U, class ...Args>
	void construct(U* p, Args&& ...args){::new(p) U(std::forward<Args>(args)...);}

	template <class U> void destroy(U* p) noexcept{p->~U();}

	std::size_t max_size() const noexcept{return std::numeric_limits<size_type>::max();}

	allocator select_on_container_copy_construction() const{return *this;}

	using propagate_on_container_copy_assignment = std::false_type;
	using propagate_on_container_move_assignment = std::false_type;
	using propagate_on_container_swap            = std::false_type;
	using is_always_equal                        = std::is_empty<allocator>;
#endif
};

template <class T, class U>
bool operator==(allocator<T> const&, allocator<U> const&) noexcept{return true;}

template <class T, class U> 
bool operator!=(allocator<T> const& x, allocator<U> const& y) noexcept{
	return !(x == y);
}

#if 0
template<typename T>
struct allocator{
	using value_type      = T;
	using pointer         = value_type*;
	using size_type 	  = std::size_t;
	using difference_type =	std::ptrdiff_t;
	using propagate_on_container_move_assignment = std::true_type;
//	NODISCARD("to avoid memory leak") 
	pointer allocate(size_type n) const{return static_cast<pointer>(fftw_malloc(sizeof(T)*n));}
	void deallocate(pointer data, size_type){fftw_free(data);}
};
#endif

//template<> allocator<std::complex<double>>::pointer allocator<std::complex<double>>::allocate(size_type n){return reinterpret_cast<std::complex<double>*>(fftw_alloc_complex(n));}
//template<> allocator<             double >::pointer allocator<             double >::allocate(size_type n){return                                         fftw_alloc_real(n)    ;}

#if 0
template<>
struct allocator<std::complex<double>>{
	using value_type      = std::complex<double>;
	using pointer         = value_type*;
	using size_type 	  = std::size_t;
	using difference_type =	std::ptrdiff_t;
	using propagate_on_container_move_assignment = std::true_type;
	NODISCARD("to avoid memory leak") 
	pointer allocate(size_type n){return reinterpret_cast<std::complex<double>*>(fftw_alloc_complex(n));}
	void deallocate(pointer data, size_type){fftw_free(data);}
};

template<>
struct allocator<double>{
	using value_type      = double;
	using pointer         = value_type*;
	using size_type 	  = std::size_t;
	using difference_type =	std::ptrdiff_t;
	using propagate_on_container_move_assignment = std::true_type;
	NODISCARD("to avoid memory leak") 
	pointer allocate(size_type n){return fftw_alloc_real(n);}
	void deallocate(pointer data, size_type){fftw_free(data);}
};
#endif

}}}


#if 0//__NVCC__
namespace std{

template<class T> struct allocator_traits<boost::multi::fftw::allocator<T>> : std::allocator_traits<std::allocator<T>>{
	using base = std::allocator_traits<std::allocator<T>>;
	template<class U> using rebind_alloc = boost::multi::fftw::allocator<U>;
	template<class A>
	static auto allocate(A& a, typename base::size_type n){return a.allocate(n);}
};

}
#endif

#if 0 //def __NVCC__
namespace std{
template<class T> struct allocator_traits<boost::multi::fftw::allocator<T>> : std::allocator_traits<std::allocator<T>>{
	template<class U> using rebind_alloc = boost::multi::fftw::allocator<U>;
};
}
#endif

#if not __INCLUDE_LEVEL__

#include "../../array.hpp"

#include<vector>

namespace multi = boost::multi;

int main(){
	{
		std::vector<double, multi::fftw::allocator<double>> v(100);
		multi::array<double, 2> arr({10, 20});
	}
	{
		std::vector<std::complex<double>, multi::fftw::allocator<std::complex<double>>> v(100);
		multi::array<std::complex<double>, 2> arr({10, 20});
	}
}
#endif
#endif

