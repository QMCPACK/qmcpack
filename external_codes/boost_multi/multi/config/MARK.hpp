#ifdef COMPILATION// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4-*-
$CXXX $CXXFLAGS $0 -o $0x -include"boost/log/trivial.hpp" -D'MULTI_MARK_SCOPE(MsG)=BOOST_LOG_TRIVIAL(trace)<<MsG' -DBOOST_LOG_DYN_LINK -lboost_log -lboost_thread -lboost_system -lboost_log_setup -lpthread &&$0x&&rm $0x;exit
#endif
// © Alfredo A. Correa 2020-2021

//#ifndef MULTI_CONFIG_MARK_HPP
//#define MULTI_CONFIG_MARK_HPP

#ifndef MULTI_MARK_SCOPE
#ifdef CALI_CXX_MARK_SCOPE
#define MULTI_MARK_SCOPE(MsG) CALI_CXX_MARK_SCOPE(MsG)
#else
#define MULTI_MARK_SCOPE(MsG) ((void)0)
#endif
#endif

//#ifndef MULTI_MARK_SCOPE
//#define MULTI_MARK_SCOPE(MsG) ((void)0)
//#else
//#define MULTI_MARK_SCOPE(MsG) CALI_CXX_MARK_SCOPE(MsG)
//#endif

//#ifndef MULTI_MARK_FUNCTION
//#define MULTI_MARK_FUNCTION MULTI_MARK_SCOPE(__func__)
//#endif

//#ifndef MULTI_MARK_PRETTY_FUNCTION
//#define MULTI_MARK_PRETTY_FUNCTION MULTI_MARK_SCOPE(__PRETTY_FUNCTION__)
//#endif

#if not __INCLUDE_LEVEL__

// #include <boost/log/trivial.hpp> or use command line `-include file.hpp`
// for example add this to the end of compilation line to use markings for tracing
// -include"boost/log/trivial.hpp" -D'MULTI_MARK_SCOPE(MsG)=BOOST_LOG_TRIVIAL(trace)<<MsG' -DBOOST_LOG_DYN_LINK -lboost_log -lboost_thread -lboost_system -lboost_log_setup -lpthread &&$0x&&

void multi_fun(){
	MULTI_MARK_SCOPE("function multi_fun"); //	BOOST_LOG_TRIVIAL(trace)<< "[2020-09-23 16:42:55.035034] [0x00007f6060af3a00] [trace]   function multi_fun";
}

void multi_gun(){
	MULTI_MARK_FUNCTION;
}

int main(){
	multi_fun();
	multi_gun();
// output: 
//[2020-09-24 16:36:44.683215] [0x00007f54677eed40] [trace]   function multi_fun
//[2020-09-24 16:36:44.683279] [0x00007f54677eed40] [trace]   multi_gun
}

#endif
//#endif

