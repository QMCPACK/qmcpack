#ifndef AFQMC_TUPLE_ITERATORS
#define AFQMC_TUPLE_ITERATORS

#include<tuple>
#include<cstdlib>
#include<functional>
#include <algorithm> 
#include<boost/iterator/iterator_facade.hpp>
#include<boost/iterator/iterator_traits.hpp>
//#include"AFQMC/Utilities/Utils.h"


namespace qmcplusplus{

template<
	class FirstIter, class SecondIter, 
	typename R = std::pair<typename std::iterator_traits<FirstIter>::value_type&, typename std::iterator_traits<SecondIter>::value_type&>,
	typename D = typename std::iterator_traits<FirstIter>::difference_type
>
class paired_iterator : 
	public boost::iterator_facade<
        paired_iterator<FirstIter, SecondIter>,
     	std::pair<
	        typename std::iterator_traits<FirstIter>::value_type,
  	    	typename std::iterator_traits<SecondIter>::value_type 
		>,
        std::random_access_iterator_tag,
        R,
	D
    >
{
public:
   FirstIter first;
   SecondIter second;
   paired_iterator(FirstIter first_, SecondIter second_) : first(first_), second(second_){}

    using difference_type = D;
    using reference = R;

private:
    friend class boost::iterator_core_access;

    void increment(){++first; ++second;}
    void decrement(){--first; --second;}
    bool equal(paired_iterator const& other) const{
	    return first == other.first and second == other.second;
	}
    reference dereference() const{return {*first, *second};} 
    void advance(difference_type n){first += n; second += n;}
    difference_type distance_to(paired_iterator other) const{
	    return other.first - first;
	}

};

template<
        class FirstIter, class SecondIter, class ThirdIter,
        typename R = std::tuple<typename std::iterator_traits<FirstIter>::value_type&, 
                                typename std::iterator_traits<SecondIter>::value_type&,
                                typename std::iterator_traits<ThirdIter>::value_type&>,
        typename D = typename std::iterator_traits<FirstIter>::difference_type
>
class tuple_iterator :
        public boost::iterator_facade<
        tuple_iterator<FirstIter, SecondIter, ThirdIter>,
        std::tuple<
                typename std::iterator_traits<FirstIter>::value_type,
                typename std::iterator_traits<SecondIter>::value_type,
                typename std::iterator_traits<ThirdIter>::value_type
                >,
        std::random_access_iterator_tag,
        R,
        D
    >
{
public:
   FirstIter first;
   SecondIter second;
   ThirdIter third;
   tuple_iterator(FirstIter first_, SecondIter second_, ThirdIter third_) :
       first(first_), second(second_),third(third_){}

    // Need default constructor for inplace_merge in some versions of the C++ std library
    tuple_iterator() {}

    using difference_type = D;
    using reference = R;

private:
    friend class boost::iterator_core_access;

    void increment(){++first; ++second; ++third;}
    void decrement(){--first; --second; --third;}
    // is this correct for what I want
    bool equal(tuple_iterator const& other) const{
            return first == other.first and second == other.second and third == other.third;
        }
    reference dereference() const
       {return std::make_tuple(std::ref(*first), std::ref(*second), std::ref(*third));}
    void advance(difference_type n){first += n; second += n; third += n;}
    difference_type distance_to(tuple_iterator other) const{
            return other.first - first;
        }

};

template <class FirstIter, class SecondIter>
paired_iterator<FirstIter, SecondIter>
make_paired_iterator(FirstIter ci, SecondIter vi) {
    return {ci, vi};
};

template <class FirstIter, class SecondIter, class ThirdIter>
tuple_iterator<FirstIter, SecondIter, ThirdIter>
make_tuple_iterator(FirstIter ci, SecondIter vi, ThirdIter ti){
    return {ci, vi, ti};
};


}

namespace std{

    template<class It1, class It2>
    void iter_swap(qmcplusplus::paired_iterator<It1, It2> const& a, qmcplusplus::paired_iterator<It1, It2> const& b){
        using std::swap;
        swap(std::get<0>(*a), std::get<0>(*b));
        swap(std::get<1>(*a), std::get<1>(*b));
    }

    template<class It1, class It2, class It3>
    void iter_swap(qmcplusplus::tuple_iterator<It1, It2, It3> const& a, qmcplusplus::tuple_iterator<It1, It2, It3> const& b){
        using std::swap;
        swap(std::get<0>(*a), std::get<0>(*b));
        swap(std::get<1>(*a), std::get<1>(*b));
        swap(std::get<2>(*a), std::get<2>(*b));
    }

}


#endif

