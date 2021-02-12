#ifndef MULTI_ADAPTORS_THRUST_HPP
#define MULTI_ADAPTORS_THRUST_HPP

#include <thrust/iterator/counting_iterator.h>
#include <thrust/iterator/transform_iterator.h>
#include <thrust/iterator/permutation_iterator.h>

namespace boost{
namespace multi{

template<class Difference>
struct stride_functor{
	Difference stride;
	constexpr std::ptrdiff_t operator()(Difference i) const{return stride*i;} // needs --expt-relaxed-constexpr
};

template<class It>
using strided_iterator = thrust::permutation_iterator<It, thrust::transform_iterator<stride_functor<typename std::iterator_traits<It>::difference_type>, thrust::counting_iterator<typename std::iterator_traits<It>::difference_type>>>;

template<class It>
constexpr auto make_strided_iterator(It it, typename std::iterator_traits<It>::difference_type d){
	return strided_iterator<It>{it, thrust::make_transform_iterator(thrust::make_counting_iterator(0*d), stride_functor<decltype(d)>{d})};
}

template<class T1, class Q1, class Size, class T2, class Q2>
constexpr auto copy_n(
	array_iterator<T1, 1, thrust::device_ptr<Q1>> first, Size count, 
	array_iterator<T2, 1, thrust::device_ptr<Q2>> result
){
	thrust::copy_n(
		make_strided_iterator(first.base() , first.stride() ), count, 
		make_strided_iterator(result.base(), result.stride())
	);
	return result + count;
}

template<class T1, class Q1, class T2, class Q2>
auto copy(
	array_iterator<T1, 1, thrust::device_ptr<Q1>> first, array_iterator<T1, 1, thrust::device_ptr<Q1>> last, 
	array_iterator<T2, 1, thrust::device_ptr<Q2>> result
){
	return copy_n(first, last - first, result);
}

}}

#endif

