#ifdef COMPILATION_INSTRUCTIONS// -*-indent-tabs-mode:t;tab-width:4;c-basic-offset:4;truncate-lines:1-*-
$CXX $0 -o $0x -DBOOST_TEST_DYN_LINK -lboost_unit_test_framework&&$0x&&rm $0x;exit
#endif

dsadsadsadsa

#include<iterator> // iterator_traits

namespace boost{
namespace multi{

template<class Reference>
struct reference_traits{
	using reference = Reference;
	using rebind_const = std::add_const_t<Reference>;
};

template<class T>
struct reference_traits<T&>{
	using reference = T&;
	using rebind_const = std::add_const_t<T&>;
};

template<class Iter>
class const_iterator : Iter{
// see notes in https://en.cppreference.com/w/cpp/iterator/move_iterator
 public:
	using base_type = Iter;
	using iterator_type = Iter;
	using iterator_category = typename std::iterator_traits<Iter>::iterator_category;
	using iterator_concept = std::input_iterator_tag;
	using value_type = typename std::iterator_traits<Iter>::value_type;
	using difference_type = typename std::iterator_traits<Iter>::difference_type;
	using pointer = Iter;
	using reference = typename multi::reference_traits<typename std::iterator_traits<Iter>::reference>::rebind_const;

 public:
	constexpr const_iterator() = default;                                 // (1)
	constexpr explicit const_iterator(iterator_type x) : Iter{x} {}       // (2)
	template<class U>                                                     // (3)
	constexpr explicit const_iterator(const_iterator<U> const& o) : Iter{o.base()} {}
	template<class U> 
	constexpr const_iterator& operator=(const_iterator<U> const& o) {
		static_cast<Iter&>(*this)=o.base();
		return *this;
	}
	constexpr base_type base() const {return static_cast<Iter const&>(*this);}
// https://en.cppreference.com/w/cpp/iterator/move_iterator/operator*
	reference operator*() const {return *static_cast<Iter const&>(*this);}
	constexpr pointer operator->() const {return &*static_cast<Iter const&>(*this);}
// https://en.cppreference.com/w/cpp/iterator/move_iterator/operator_at
	constexpr reference operator[](difference_type n) const {return static_cast<Iter&>(*this)[n];}
// https://en.cppreference.com/w/cpp/iterator/move_iterator/operator_arith
	constexpr const_iterator& operator++(){return ++static_cast<Iter&>(*this), *this;}            //(1)
	constexpr const_iterator& operator--(){return --static_cast<Iter&>(*this), *this;}            //(2)
	constexpr const_iterator  operator++(int){return const_iterator{static_cast<Iter&>(*this)++};}//(3)
	constexpr const_iterator  operator--(int){return const_iterator{static_cast<Iter&>(*this)--};}//(4)
	constexpr const_iterator  operator+(difference_type n) const{          //(5)
		return const_iterator{static_cast<Iter const&>(*this)+n};
	}
	constexpr const_iterator  operator-(difference_type n) const{          //(6)
		return const_iterator{static_cast<Iter const&>(*this)-n};
	}
	constexpr const_iterator& operator+=(difference_type n){               //(7)
		return static_cast<Iter&>(*this)+=n, *this;
	}
	constexpr const_iterator& operator-=(difference_type n){               //(8)
		return static_cast<Iter&>(*this)-=n, *this;
	}
	template<class Other, std::enable_if_t<
		std::is_same<typename std::iterator_traits<Other>::value_type, typename std::iterator_traits<Iter>::value_type>{} and
		std::is_assignable<Other&, Iter>{} and 
		not std::is_assignable<typename std::iterator_traits<Other>::reference, typename std::iterator_traits<Iter>::value_type>{}, int 
	> =0>
	operator Other() const{return base();}
	using rebind_const = const_iterator;
};

template<class Iter1, class Iter2> constexpr bool operator==(const_iterator<Iter1> const& lhs, const_iterator<Iter2> const& rhs){return lhs.base()==rhs.base();} //(1)
template<class Iter1, class Iter2> constexpr bool operator!=(const_iterator<Iter1> const& lhs, const_iterator<Iter2> const& rhs){return lhs.base()!=rhs.base();} //(2)
template<class Iter1, class Iter2> constexpr bool operator< (const_iterator<Iter1> const& lhs, const_iterator<Iter2> const& rhs){return lhs.base()< rhs.base();} //(3)
template<class Iter1, class Iter2> constexpr bool operator<=(const_iterator<Iter1> const& lhs, const_iterator<Iter2> const& rhs){return lhs.base()<=rhs.base();} //(4)
template<class Iter1, class Iter2> constexpr bool operator> (const_iterator<Iter1> const& lhs, const_iterator<Iter2> const& rhs){return lhs.base()> rhs.base();} //(5)
template<class Iter1, class Iter2> constexpr bool operator>=(const_iterator<Iter1> const& lhs, const_iterator<Iter2> const& rhs){return lhs.base()>=rhs.base();} //(6)
// TODO three way comparison for C++20

template<class It>
const_iterator<It> make_const_iterator(It it){return const_iterator<It>{it};}

template <typename Iterator, typename Base = typename std::iterator_traits<Iterator>,typename Enable = void> 
struct iterator_traits : Base{
	using rebind_const = multi::const_iterator<Iterator>;
};

template<typename Iterator, typename Base> 
struct iterator_traits<Iterator, Base, typename Iterator::rebind_const> : Base{
	using rebind_const = typename Iterator::rebind_const; 
};

template<class T>
struct iterator_traits<T*> : std::iterator_traits<T*>{
	using rebind_const = T const*;
};

}}

#if not __INCLUDE_LEVEL__  // TEST BELOW

#define BOOST_TEST_MODULE test const_iterator
#ifdef BOOST_TEST_DYN_LINK
#include<boost/test/unit_test.hpp>
#else
#include<boost/test/included/unit_test.hpp>
#endif

#include<vector>

namespace multi = boost::multi;

BOOST_AUTO_TEST_CASE(sematics) {
	std::vector<int> v(5, 9);
	std::vector<int>::iterator it = v.begin();
	*it += 1;
	BOOST_REQUIRE(v[0] == 10 );

	static_assert( std::is_same<multi::iterator_traits<double*>::rebind_const, double const*>{}, "" );
	static_assert( std::is_same<multi::iterator_traits<std::vector<double>::iterator>::rebind_const, multi::const_iterator<std::vector<double>::iterator> >{}, "" );

	std::vector<int>::const_iterator cit = multi::make_const_iterator(v.begin()); (void)cit;
}

#endif

