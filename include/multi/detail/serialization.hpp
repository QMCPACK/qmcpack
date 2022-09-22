// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;autowrap:nil;-*-
// Copyright 2018-2021 Alfredo A. Correa

#ifndef MULTI_DETAIL_SERIALIZATION_HPP
#define MULTI_DETAIL_SERIALIZATION_HPP

#include<algorithm>  // for std::for_each
#include<cstdint>    // for std::uint32_t

namespace boost {  // NOLINT(modernize-concat-nested-namespaces) keep c++14 compat
namespace archive {  // NOLINT(modernize-concat-nested-namespaces) keep c++14 compat
namespace detail {

template<class Ar> class common_iarchive;
template<class Ar> class common_oarchive;

}  // end namespace detail
}  // end namespace archive

namespace serialization {  // NOLINT(modernize-concat-nested-namespaces) keep c++14 compat

template<class T> class  nvp;            // dependency "in name only"
template<class T> class  array_wrapper;  // dependency "in name only"
                  struct binary_object;  // dependency "in name only", if you get an error here, it means that eventually you need to include #include<boost/serialization/binary_object.hpp>

template<typename T> struct version;

//template<class Archive, class T>//, std::enable_if_t<std::is_same<T, std::decay_t<T>>{}, int> =0>
//auto operator>>(Archive& ar, T&& t) -> Archive& {return ar>> t;}

}  // end namespace serialization
}  // end namespace boost

namespace cereal {

template<class ArchiveType, std::uint32_t Flags> struct OutputArchive;
template<class ArchiveType, std::uint32_t Flags> struct  InputArchive;

template<class T> class NameValuePair;  // dependency "in name only", if you get an error here you many need to #include <cereal/archives/xml.hpp> at some point

}  // end namespace cereal

namespace boost {  // NOLINT(modernize-concat-nested-namespaces) keep c++14 compat
namespace multi {

template<class Archive, class MA, std::enable_if_t<std::is_same_v<MA, std::decay_t<MA>> and (MA::dimensionality > -1) , int> =0>
auto operator>>(Archive& arxiv, MA&& self)  // this is for compatibility with Archive type
->decltype(arxiv>> self) {
	return arxiv>> self; }

template<class Archive, class MA, std::enable_if_t<std::is_same_v<MA, std::decay_t<MA>> and (MA::dimensionality > -1), int> =0>
auto operator& (Archive& arxiv, MA&& self)  // this is for compatibility with Archive type
->decltype(arxiv& self) {
	return arxiv& self; }

template<class Ar, class Enable = void>
struct archive_traits {
	template<class T>
	inline static auto make_nvp  (char const* /*n*/, T&& value) noexcept {return std::forward<T>(value);}
};

template<class Ar>
struct archive_traits<Ar, typename std::enable_if<std::is_base_of_v<boost::archive::detail::common_oarchive<Ar>, Ar> || std::is_base_of_v<boost::archive::detail::common_iarchive<Ar>, Ar>>::type> {
	template<class T> using nvp           = boost::serialization::nvp          <T>;
	template<class T> using array_wrapper = boost::serialization::array_wrapper<T>;
	template<class T> struct binary_object_t {using type = boost::serialization::binary_object;};
	template<class T>        inline static auto make_nvp          (char const* name, T&  value) noexcept -> const nvp<T> {return nvp<T>{name, value};}  // NOLINT(readability-const-return-type) : original boost declaration
	template<class T>        inline static auto make_nvp          (char const* name, T&& value) noexcept -> const nvp<T> {return nvp<T>{name, value};}  // NOLINT(readability-const-return-type) : original boost declaration

	template<class T>        inline static auto make_array        (               T* first, std::size_t size) noexcept -> const array_wrapper<T> {return array_wrapper<T>{first, size};}  // NOLINT(readability-const-return-type) : original boost declaration
	template<class T = void> inline static auto make_binary_object(      const void* first, std::size_t size) noexcept -> const typename binary_object_t<T>::type {return typename binary_object_t<T>::type(first, size); }  // if you get an error here you need to eventually `#include<boost/serialization/binary_object.hpp>`// NOLINT(readability-const-return-type,clang-diagnostic-ignored-qualifiers) : original boost declaration
};

#if 1
template<class Ar>
struct archive_traits<
	Ar, 
	typename std::enable_if<
		   std::is_base_of_v<cereal::OutputArchive<Ar, 0>, Ar> or std::is_base_of_v<cereal::OutputArchive<Ar, 1>, Ar>
		or std::is_base_of_v<cereal::InputArchive <Ar, 0>, Ar> or std::is_base_of_v<cereal:: InputArchive<Ar, 1>, Ar>
	>::type
> {
	using self_t = archive_traits<Ar, typename std::enable_if<
		   std::is_base_of_v<cereal::OutputArchive<Ar, 0>, Ar> or std::is_base_of_v<cereal::OutputArchive<Ar, 1>, Ar>
		or std::is_base_of_v<cereal:: InputArchive<Ar, 0>, Ar> or std::is_base_of_v<cereal:: InputArchive<Ar, 1>, Ar>
	>::type>;

	template<class T>
	inline static auto make_nvp  (char const* name, T&& value) noexcept {return cereal::NameValuePair<T&>{name, value};}  // if you get an error here you many need to #include <cereal/archives/xml.hpp> at some point
	template<class T>
	inline static auto make_nvp  (char const* name, T&  value) noexcept {return cereal::NameValuePair<T&>{name, value};}  // if you get an error here you many need to #include <cereal/archives/xml.hpp> at some point

	template<class T>
	struct array_wrapper {
		T* p_;
		std::size_t c_;
		template<class Archive>
		void serialize(Archive& arxiv, const unsigned int /*version*/) {
			for(std::size_t i = 0; i != c_; ++i) {  // NOLINT(altera-unroll-loops) TODO(correaa) consider using an algorithm
				auto& item = p_[i];  // NOLINT(cppcoreguidelines-pro-bounds-pointer-arithmetic)
				arxiv &                                        make_nvp("item"   , item   );  // "item" is the name used by Boost.Serialization XML make_array
			//	arxiv & boost::multi::archive_traits<Archive>::make_nvp("element", element);
			//	arxiv &                                cereal::make_nvp("element", element);
			//	arxiv &                                      CEREAL_NVP(           element);
			//	arxiv &                                                            element ;
			}
		}
	};

	template<class T>
	inline static auto make_array(T* ptr, std::size_t count) -> array_wrapper<T> {return array_wrapper<T>{ptr, count};}

	template<class T>
	inline static auto make_nvp  (char const* name, array_wrapper<T>&& value) noexcept {return make_nvp(name, value);}
};

//template<class Ar, class T> auto make_nvp(char const* n, T&& v) -> decltype(auto) {return archive_traits<Ar>::make_nvp(n, std::forward<T>(v));}  // NOLINT(readability-const-return-type)
#endif

}  // end namespace multi
}  // end namespace boost

namespace boost {

template<class T, std::size_t D, class As>
class multi_array;

}  // end namespace boost

namespace boost {  // NOLINT(modernize-concat-nested-namespaces) keep c++14 compat
namespace serialization {

//template<class Archive, class T, std::size_t D, class A>
//auto serialize(Archive& ar, boost::multi_array<T, D, A>& arr, unsigned int /*version*/)
//{
//	auto x = boost::multi::extensions(arr);
//	ar & multi::archive_traits<Archive>::make_nvp("extensions", x);
//	if( x != boost::multi::extensions(arr) ) {
//		arr.resize( std::array<std::size_t, 2>{} );
//		arr.resize( std::array<std::size_t, 2>{static_cast<std::size_t>(std::get<0>(x).size()), static_cast<std::size_t>(std::get<1>(x).size())} );
//	}
//	ar & multi::archive_traits<Archive>::make_nvp("data_elements", multi::archive_traits<Archive>::make_array(boost::multi::data_elements(arr), static_cast<std::size_t>(boost::multi::num_elements(arr))));
//}

}  // end namespace serialization
}  // end namespace boost

#endif
