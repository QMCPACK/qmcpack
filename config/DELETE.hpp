#ifdef COMPILATION// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4-*-
$CXX $0 -o $0x &&$0x&&rm $0x;exit
#endif
// © Alfredo A. Correa 2020

#ifndef MULTI_CONFIG_DELETE_HPP
#define MULTI_CONFIG_DELETE_HPP

namespace boost{
namespace multi{

template<bool B, class T = int> struct disable_if_impl{};
template<class T> struct disable_if_impl<false, T>{ typedef T type; };
template<bool B = false, class T = int> using disable_if = typename disable_if_impl<B, T>::type;

}}

#define DELETE(ConD) boost::multi::disable_if<ConD> =0

#if not __INCLUDE_LEVEL__ // _TEST_MULTI_CONFIG_NODISCARD

int main(){
}
#endif
#endif


