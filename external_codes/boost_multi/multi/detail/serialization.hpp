#ifndef _MULTI_DETAIL_SERIALIZATION_HPP
#include<boost/serialization/nvp.hpp>

namespace boost{
namespace multi{

template<class Ar, typename = decltype(Ar::make_nvp((char const*)nullptr, int{}))
std::true_type  has_make_nvp_aux(Ar const&);
std::false_type has_make_nvp_aux(...);

template<class Ar, class T>
struct has_make_nvp_t : decltype(has_make_nvp_aux(std::declval<Ar>()){};

template<class Ar, typename = has_make_nvp_aux(std::declval<Ar)>
struct archive_traits{
};

struct archive_traits{
	template<class T>
	make_nvp(char const*, T&)
};

}}
#endif

