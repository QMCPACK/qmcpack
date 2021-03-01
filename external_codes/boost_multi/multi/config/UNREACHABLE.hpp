#ifdef COMPILATION// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4-*-
echo $X
$CXXX $CXXFLAGS $0 -o $0.$X &&$0.$X&&rm $0.$X;exit
#endif
// © Alfredo A. Correa 2020

#ifndef MULTI_CONFIG_NODISCARD_HPP
#define MULTI_CONFIG_NODISCARD_HPP

#if defined(__GNUC__) || defined(__clang__) || defined(__INTEL_COMPILER) || defined(__NVCC__)
#define MULTI_UNREACHABLE __builtin_unreachable()
#else
#define MULTI_UNREACHABLE do { std::abort(); } while(0)
#endif

#if defined(__INCLUDE_LEVEL__) and not __INCLUDE_LEVEL__

enum color{red, green, blue};

int f(enum color c){
	switch(c){
	case red  : return 1;
	case green: return 2;
	case blue : return 3; // comment case make gcc, clang, culang causes -Wswitch warning
	} MULTI_UNREACHABLE;  // commnet unreachable in gcc and nvcc causes -Wreturn-type warning 
}

int main(){
}

#endif
#endif

