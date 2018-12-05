#if COMPILATION_INSTRUCTIONS
(echo "#include<"$0">" > $0x.cpp) && c++ -O3 -std=c++14 -Wall -D_TEST_BOOST_MPI3_DETAIL_STRIDED $0x.cpp -o $0x.x && $0x.x $@ && rm -f $0x.cpp; exit
#endif

#ifndef ALF_BOOST_MPI3_DETAIL_STRIDED_HPP
#define ALF_BOOST_MPI3_DETAIL_STRIDED_HPP

#include<boost/iterator/iterator_adaptor.hpp>

namespace boost{
namespace mpi3{
namespace detail{

template <class Iterator>
class strided : 
	public boost::iterator_adaptor<
		strided<Iterator>, 
		Iterator,
		typename std::iterator_traits<std::decay_t<Iterator>>::value_type,
		typename std::iterator_traits<std::decay_t<Iterator>>::iterator_category,
		typename std::iterator_traits<std::decay_t<Iterator>>::reference,
		typename std::iterator_traits<std::decay_t<Iterator>>::difference_type
	>
{
	using stride_type = typename strided::iterator_adaptor_::difference_type;
	stride_type stride_;
	Iterator&& base_reference() &&{return std::move(strided::iterator_adaptor_::base_reference());}
	Iterator& base_reference() &{return strided::iterator_adaptor_::base_reference();}
	Iterator const& base_reference() const &{return strided::iterator_adaptor_::base_reference();}
public:
	strided(Iterator it, stride_type stride) : strided::iterator_adaptor_(it), stride_(stride){}
	template<class Strided>
	strided(Strided&& s) : 
		strided::iterator_adaptor_(std::forward<Strided>(s).base_reference()), 
		stride_(std::forward<Strided>(s).stride_){}
	using difference_type = typename strided::iterator_adaptor_::difference_type;
	stride_type stride() const{return stride_;}
private:
	friend class boost::iterator_core_access;
	void increment(){std::advance(this->base_reference(),  stride_);}
	void decrement(){std::advance(this->base_reference(), -stride_);}
	void advance(difference_type n){std::advance(this->base_reference(), stride_*n);}
	auto distance_to(strided const& other) const{
		return std::distance(this->base_reference(), other.base_reference())/stride_;
	}
};

template<class Iterator>
strided<Iterator> make_strided(Iterator it, typename std::iterator_traits<Iterator>::difference_type stride){
    return {it, stride};
};

template<class Iterator>
strided<Iterator> stride(Iterator&& it, typename std::iterator_traits<std::decay_t<Iterator>>::difference_type stride){
	return strided<Iterator>(std::forward<Iterator>(it), stride);
}

template<class Iterator, class = decltype(std::declval<Iterator>().stride())>
std::true_type is_strided_aux(Iterator);
std::false_type is_strided_aux(...);

template<class Iterator> struct is_strided : decltype(is_strided_aux(std::declval<Iterator>())){};

}}}

#ifdef _TEST_BOOST_MPI3_DETAIL_STRIDED

#include "alf/boost/mpi3/detail/iterator.hpp"

#include<vector>
#include<iostream>
#include<numeric>
#include<cassert>

using std::cout;

int main(){

	std::vector<double> v(9); std::iota(v.begin(), v.end(), 0);
	for(auto it = v.begin(); it != v.end(); ++it) cout << *it << " ";
	cout << '\n';

	using boost::mpi3::detail::stride;

	for(auto it = stride(v.begin(), 3); it != stride(v.end(), 3); ++it)
		cout << *it << ' '
	;
	auto it = stride(v.begin(), 3);

	assert(( std::is_same<std::iterator_traits<decltype(it)>::value_type, double>{} ));
	assert(( std::is_same<std::iterator_traits<decltype(it)>::iterator_category, std::random_access_iterator_tag>{} ));
	assert(( std::is_same<decltype(it)::base_type, std::vector<double>::iterator>{} ));
	assert( &*it.base() == &*it );
	assert( &*it == boost::mpi3::detail::data(it.base()) );
	assert( it.stride() == 3 );

}

#endif
#endif


