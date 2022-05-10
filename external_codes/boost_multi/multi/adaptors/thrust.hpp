// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;autowrap:nil;-*-
// © Alfredo A. Correa 2021

#pragma once

#include "../array.hpp"

#include "./thrust/cuda/managed.hpp"

#include <thrust/device_allocator.h>
#include <thrust/system/cuda/memory.h> // ::thrust::cuda::allocator

#include <thrust/system/cuda/experimental/pinned_allocator.h>
#include <thrust/host_vector.h>

namespace boost{
namespace multi{
namespace thrust{

template<class T, multi::dimensionality_type D> using device_array = multi::array<T, D, ::thrust::device_allocator<T>>;
template<class T, multi::dimensionality_type D> using host_array   = multi::array<T, D                               >;

namespace device{

template<class T, multi::dimensionality_type D> using array = device_array<T, D>;

}

namespace host{

template<class T, multi::dimensionality_type D> using array = host_array<T, D>;

}

namespace cuda{

template<class T, multi::dimensionality_type D> using array = multi::array<T, D, ::thrust::cuda::allocator<T>>;

namespace managed{

template<class T, multi::dimensionality_type D> using array = multi::array<T, D, boost::multi::thrust::cuda::managed::allocator<T>>;

}

}

}}}

namespace thrust{

//template<>
//struct iterator_system<boost::multi::basic_array<char, 2L, char *, boost::multi::layout_t<2L, boost::multi::size_type>>::elements_iterator_t<char *>>{
//	using type = thrust::iterator_system<char *>::type;
//};

//template<>
//struct iterator_system<boost::multi::basic_array<char, 2L, std::__detected_or_t<char *, std::__allocator_traits_base::__pointer, thrust::cuda_cub::allocator<char>>, boost::multi::layout_t<2L, boost::multi::size_type>>::elements_iterator_t<std::__detected_or_t<char *, std::__allocator_traits_base::__pointer, thrust::cuda_cub::allocator<char>>>>{
//	using type = thrust::iterator_system<thrust::cuda::pointer<char>>::type;
//};

//template<class MultiIterator>
//struct elements_range{
//	using difference_type = std::ptrdiff_t;

// public:
//	struct strides_functor {
//		std::ptrdiff_t z_;
//		typename MultiIterator::layout_type l_;

//		constexpr difference_type operator()(const difference_type& n) const {
//			auto const x = l_.extensions();
//			return n/x.num_elements()*z_ + std::apply(l_, x.from_linear(n%x.num_elements()));
//		}
//	};

// protected:
//	using CountingIterator = thrust::counting_iterator<difference_type>;
//	using TransformIterator = thrust::transform_iterator<strides_functor, CountingIterator>;
//	using PermutationIterator = thrust::permutation_iterator<typename MultiIterator::element_ptr, TransformIterator>;

// public:
//	elements_range(MultiIterator it, std::ptrdiff_t count)
//	: base_{it.base()}
//	, sf_{it.stride(), it->layout()}
//	, size_{count*(it->extensions().num_elements())} {}

//	using iterator = typename thrust::permutation_iterator<typename MultiIterator::element_ptr, TransformIterator>;

//	auto begin() const {return iterator{base_, TransformIterator{CountingIterator{    0}, sf_}};}
//	auto end()   const {return iterator{base_, TransformIterator{CountingIterator{size_}, sf_}};}

//	auto size() const {return size_;}

// protected:
//	typename MultiIterator::element_ptr base_;
//	strides_functor sf_;
//	std::ptrdiff_t size_;
//};

template<class T1, class Q1, class Size, class T2, class Q2>
auto copy_n(
	boost::multi::array_iterator<T1, 1, Q1*                      >   first , Size count,
	boost::multi::array_iterator<T2, 1, thrust::cuda::pointer<Q2>> d_first
)-> boost::multi::array_iterator<T2, 1, thrust::cuda::pointer<Q2>> {
	assert(first->extensions() == d_first->extensions());

	if constexpr(std::is_trivially_assignable<Q2&, Q1&>{}) {
		if(count == 0) return d_first;
		if(first.stride() == 1 and d_first.stride() == 1) {
			auto s = cudaMemcpy  (raw_pointer_cast(d_first.base()),                                                        first.base(),                                                      sizeof(T2)* static_cast<std::size_t>(count), cudaMemcpyHostToDevice); assert( s == cudaSuccess );
		} else {
			auto s = cudaMemcpy2D(raw_pointer_cast(d_first.base()), static_cast<std::size_t>(d_first.stride())*sizeof(T2), first.base(), static_cast<std::size_t>(first.stride())*sizeof(T2), sizeof(T2), static_cast<std::size_t>(count), cudaMemcpyHostToDevice); assert( s == cudaSuccess );
		}
	} else {
		auto const& source_range = boost::multi::ref(first , first + count).elements();
		thrust::host_vector<T1, thrust::cuda::experimental::pinned_allocator<T1>> buffer(source_range.begin(), source_range.end());
		::thrust::copy_n(thrust::device,
			buffer.begin(), buffer.size(),
			boost::multi::ref(d_first, d_first + count).template reinterpret_array_cast<T2, Q2*>().elements().begin()
		);
	}
	return d_first + count;
}

template<class T1, class Q1, class Size, class T2, class Q2, boost::multi::dimensionality_type D>
auto copy_n(
	boost::multi::array_iterator<T1, D, Q1*                      >   first, Size count,
	boost::multi::array_iterator<T2, D, thrust::cuda::pointer<Q2>> d_first
)-> boost::multi::array_iterator<T2, D, thrust::cuda::pointer<Q2>> {
	assert(first->extensions() == d_first->extensions());

	auto const& source_range = boost::multi::ref(first , first + count).elements();

	if constexpr(std::is_trivially_assignable<Q2&, Q1&>{}) {
		if(first->stride() == 1 and d_first->stride() == 1 and D == 2) {
			auto s = cudaMemcpy2D(raw_pointer_cast(d_first.base()), static_cast<std::size_t>(d_first.stride())*sizeof(T2), first.base(), static_cast<std::size_t>(first.stride())*sizeof(T2), static_cast<std::size_t>(first->size())*sizeof(T2), static_cast<std::size_t>(count), cudaMemcpyHostToDevice); assert( s == cudaSuccess );
		} else {
			cudaHostRegister(
				const_cast<void*>(static_cast<void const*>(raw_pointer_cast(&source_range.front()))),
			//	(void*)&source_range.front(),
				static_cast<std::size_t>(&source_range.back() - &source_range.front())*sizeof(T1), cudaHostRegisterPortable
			);  // cudaHostRegisterReadOnly not available in cuda < 11.1
			::thrust::copy_n(
				thrust::cuda::par,
				source_range.begin(), source_range.size(),
				boost::multi::ref(d_first, d_first + count).template reinterpret_array_cast<T2, Q2*>().elements().begin()
			);
			cudaHostUnregister(const_cast<void*>(static_cast<void const*>(raw_pointer_cast(&source_range.front()))));
		}
	} else {
		thrust::host_vector<T1, thrust::cuda::experimental::pinned_allocator<T1>> buffer(source_range.begin(), source_range.end());
		::thrust::copy_n(thrust::cuda::par,
			buffer.begin(), buffer.size(),
			boost::multi::ref(d_first, d_first + count).template reinterpret_array_cast<T2, Q2*>().elements().begin()
		);
	}

	return d_first + count;
}

template<class T1, class Q1, class T2, class Q2, boost::multi::dimensionality_type D>
auto copy(
	boost::multi::array_iterator<T1, D, Q1*                      >   first,
	boost::multi::array_iterator<T1, D, Q1*                      >   last ,
	boost::multi::array_iterator<T2, D, thrust::cuda::pointer<Q2>> d_first
)-> boost::multi::array_iterator<T2, D, thrust::cuda::pointer<Q2>> {
	return copy_n(first, last - first, d_first);
}

template<class T1, class Q1, class T2, class Q2, boost::multi::dimensionality_type D>
auto uninitialized_copy(
	boost::multi::array_iterator<T1, D, Q1*                      >   first,
	boost::multi::array_iterator<T1, D, Q1*                      >   last ,
	boost::multi::array_iterator<T2, D, thrust::cuda::pointer<Q2>> d_first
)-> boost::multi::array_iterator<T2, D, thrust::cuda::pointer<Q2>> {
	if constexpr(std::is_trivially_constructible<T2>{}){
		return copy(first, last, d_first);
	}
	throw std::logic_error{"uninitialized_copy for nontrivials in cuda device not implemented"};
}

template<class T1, class Q1, class Size, class T2, class Q2, boost::multi::dimensionality_type D>
auto uninitialized_copy_n(
	boost::multi::array_iterator<T1, D, Q1*                      >   first, Size count,
	boost::multi::array_iterator<T2, D, thrust::cuda::pointer<Q2>> d_first
)-> boost::multi::array_iterator<T2, D, thrust::cuda::pointer<Q2>> {
	if constexpr(std::is_trivial_v<T2> and std::is_nothrow_assignable_v<T2&, Q2&>) {
		return copy_n(first, count, d_first);
	}
	throw std::logic_error{"uninitialized_copy_n for nontrivials in cuda device not implemented"};
}

}

