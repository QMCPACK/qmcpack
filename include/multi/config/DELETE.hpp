// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4-*-
// Copyright 2018-2021 Alfredo A. Correa

#ifndef MULTI_CONFIG_DELETE_HPP
#define MULTI_CONFIG_DELETE_HPP

namespace boost::multi {

template<bool B, class T = int> struct disable_if_impl{};
template<class T> struct disable_if_impl<false, T>{using type = T;};

template<bool B = false, class T = int> using disable_if = typename disable_if_impl<B, T>::type;

} // end namespace boost::multi

#define DELETE(ConD) boost::multi::disable_if<ConD> =0  // NOLINT(cppcoreguidelines-macro-usage) TODO(correaa) remove
#endif
