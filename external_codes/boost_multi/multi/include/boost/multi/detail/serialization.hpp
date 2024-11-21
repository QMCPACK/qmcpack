// Copyright 2018-2024 Alfredo A. Correa
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#ifndef BOOST_MULTI_DETAIL_SERIALIZATION_HPP_
#define BOOST_MULTI_DETAIL_SERIALIZATION_HPP_
#pragma once

#include <algorithm>    // for std::for_each  // IWYU pragma: keep  // bug in iwyu 0.18
#include <cstddef>      // for size_t, byte
#include <cstdint>      // for uint32_t
#include <iterator>     // for next
#include <type_traits>  // for enable_if_t, decay_t
#include <utility>      // for forward

namespace boost::archive::detail { template <class Ar> class common_iarchive; }  // lines 24-24
namespace boost::archive::detail { template <class Ar> class common_oarchive; }  // lines 25-25

namespace boost::serialization { struct binary_object; }
namespace boost::serialization { template <class T> class array_wrapper; }
namespace boost::serialization { template <class T> class nvp; }

namespace cereal { template <class ArchiveType, std::uint32_t Flags> struct InputArchive; }
namespace cereal { template <class ArchiveType, std::uint32_t Flags> struct OutputArchive; }
namespace cereal { template <class T> class NameValuePair; }  // if you get an error here you many need to #include <cereal/archives/xml.hpp> at some point  // IWYU pragma: keep  // bug in iwyu 0.18

// namespace boost {
// namespace serialization {

// template<class Archive, class T>//, std::enable_if_t<std::is_same<T, std::decay_t<T>>{}, int> =0>
// auto operator>>(Archive& ar, T&& t) -> Archive& {return ar>> t;}

// }  // end namespace serialization
// }  // end namespace boost

namespace boost {  // NOLINT(modernize-concat-nested-namespaces) keep c++14 compat
namespace multi {

template<class Ar, class Enable = void>
struct archive_traits {
	template<class T>
	/*inline*/ static auto make_nvp(char const* /*n*/, T&& value) noexcept { return std::forward<T>(value); }  // match original boost declaration
};

template<class Archive, class MA, std::enable_if_t<std::is_same_v<MA, std::decay_t<MA>> && (MA::dimensionality > -1), int> = 0>  // NOLINT(modernize-use-constraints) TODO(correaa)
auto operator>>(Archive& arxiv, MA&& self)  // this is for compatibility with Archive type
	-> decltype(arxiv >> static_cast<MA&>(std::forward<MA>(self))) {
	return arxiv >> static_cast<MA&>(std::forward<MA>(self));
}

template<class Archive, class MA, std::enable_if_t<std::is_same_v<MA, std::decay_t<MA>> && (MA::dimensionality > -1), int> = 0>  // NOLINT(modernize-use-constraints) TODO(correaa)
auto operator<<(Archive& arxiv, MA&& self)  // this is for compatibility with Archive type
	-> decltype(arxiv << static_cast<MA&>(std::forward<MA>(self))) {
	return arxiv << static_cast<MA&>(std::forward<MA>(self));
}

template<class Archive, class MA, std::enable_if_t<std::is_same_v<MA, std::decay_t<MA>> && (MA::dimensionality > -1), int> = 0>  // NOLINT(modernize-use-constraints) TODO(correaa)
auto operator&(Archive& arxiv, MA&& self)  // this is for compatibility with Archive type
	-> decltype(arxiv & static_cast<MA&>(std::forward<MA>(self))) {
	return arxiv & static_cast<MA&>(std::forward<MA>(self));
}

template<class Ar>
struct archive_traits<Ar, typename std::enable_if_t<std::is_base_of_v<boost::archive::detail::common_oarchive<Ar>, Ar> || std::is_base_of_v<boost::archive::detail::common_iarchive<Ar>, Ar>>> {
	template<class T> using nvp           = boost::serialization::nvp<T>;
	template<class T> using array_wrapper = boost::serialization::array_wrapper<T>;
	template<class T> struct binary_object_t {
		using type = boost::serialization::binary_object;
	};

	template<class T> /*inline*/ static auto make_nvp(char const* name, T& value) noexcept -> nvp<T> const { return nvp<T>{name, value}; }  // NOLINT(readability-const-return-type) : match original boost declaration
	template<class T> /*inline*/ static auto make_nvp(char const* name, T&& value) noexcept -> nvp<T> const { return nvp<T>{name, static_cast<T&>(std::forward<T>(value))}; }  // NOLINT(readability-const-return-type) : match original boost declaration

	template<class T>        /*inline*/ static auto make_array(T* first, std::size_t size) noexcept -> array_wrapper<T> const { return array_wrapper<T>{first, size}; }  // NOLINT(readability-const-return-type) original boost declaration
	template<class T = void> /*inline*/ static auto make_binary_object(std::byte const* first, std::size_t size) noexcept -> const typename binary_object_t<T>::type { return typename binary_object_t<T>::type(first, size); }  // if you get an error here you need to eventually `#include<boost/serialization/binary_object.hpp>`// NOLINT(readability-const-return-type,clang-diagnostic-ignored-qualifiers) : original boost declaration
};

template<class Ar>
struct archive_traits<
	Ar,
	typename std::enable_if_t<
		std::is_base_of_v<cereal::OutputArchive<Ar, 0>, Ar> || std::is_base_of_v<cereal::OutputArchive<Ar, 1>, Ar> || std::is_base_of_v<cereal::InputArchive<Ar, 0>, Ar> || std::is_base_of_v<cereal::InputArchive<Ar, 1>, Ar>>> {
	using self_t = archive_traits<Ar, typename std::enable_if_t<
		                                  std::is_base_of_v<cereal::OutputArchive<Ar, 0>, Ar> || std::is_base_of_v<cereal::OutputArchive<Ar, 1>, Ar> || std::is_base_of_v<cereal::InputArchive<Ar, 0>, Ar> || std::is_base_of_v<cereal::InputArchive<Ar, 1>, Ar>>>;

	//  template<class T>
	//  inline static auto make_nvp  (char const* name, T const& value) noexcept {return cereal::NameValuePair<T const&>{name, value};}  // if you get an error here you many need to #include <cereal/archives/xml.hpp> at some point  // TODO(correaa) replace by cereal::make_nvp from cereal/cereal.hpp
	// template<class T>
	// inline static auto make_nvp  (std::string const& name, T&& value) noexcept {return cereal::NameValuePair<T>{name.c_str(), std::forward<T>(value)};}  // if you get an error here you many need to #include <cereal/archives/xml.hpp> at some point
	template<class T>
	/*inline*/ static auto make_nvp(char const* name, T&& value) noexcept { return cereal::NameValuePair<T>{name, std::forward<T>(value)}; }  // if you get an error here you many need to #include <cereal/archives/xml.hpp> at some point
	//  template<class T>
	//  inline static auto make_nvp  (char const* name, T&  value) noexcept {return cereal::NameValuePair<T&>{name,                 value};}  // if you get an error here you many need to #include <cereal/archives/xml.hpp> at some point

	template<class T>
	struct array_wrapper {
		T*          p_;
		std::size_t c_;

		template<class Archive>
		void serialize(Archive& arxiv, unsigned int const /*version*/) {
			std::for_each(  // std::for_each_n is absent in GCC 7
				p_, std::next(p_, c_),
				[&arxiv](auto& item) { arxiv& make_nvp("item", item); }
			);
			// for(std::size_t i = 0; i != c_; ++i) {  // NOLINT(altera-unroll-loops) TODO(correaa) consider using an algorithm
			//  auto& item = p_[i];  // NOLINT(cppcoreguidelines-pro-bounds-pointer-arithmetic)
			//  arxiv &                                        make_nvp("item"   , item   );  // "item" is the name used by Boost.Serialization XML make_array
			//  // arxiv & boost::multi::archive_traits<Archive>::make_nvp("element", element);
			//  // arxiv &                                cereal::make_nvp("element", element);
			//  // arxiv &                                      CEREAL_NVP(           element);
			//  // arxiv &                                                            element ;
			// }
		}
	};

	template<class T>
	/*inline*/ static auto make_array(T* ptr, std::size_t count) -> array_wrapper<T> { return array_wrapper<T>{ptr, count}; }

	template<class T>
	/*inline*/ static auto make_nvp(char const* name, array_wrapper<T>&& value) noexcept { return make_nvp(name, static_cast<array_wrapper<T>&>(std::move(value))); }
};

}  // end namespace multi
}  // end namespace boost

// namespace boost {

// template<class T, std::size_t D, class As>
// class multi_array;

// }  // end namespace boost

namespace boost {  // NOLINT(modernize-concat-nested-namespaces) keep c++14 compat
namespace serialization {

template<class T>  // , class = std::enable_if_t<std::is_same_v<T&&, T&>> >
inline auto make_nvp(char const* name, T&& value) {
	return boost::serialization::make_nvp(name, static_cast<T&>(std::forward<T>(value)));
}

}  // end namespace serialization

using boost::serialization::make_nvp;

}  // end namespace boost

#endif  // BOOST_MULTI_DETAIL_SERIALIZATION_HPP_
