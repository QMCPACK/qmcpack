//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Paul R. C. Kent, kentpr@ornl.gov, Oak Ridge National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Paul R. C. Kent, kentpr@ornl.gov, Oak Ridge National Laboratory
//////////////////////////////////////////////////////////////////////////////////////
    
    
//           http://pathintegrals.info                     //
/////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// Standard definitions and functions that are used almost everywhere         //
//                                                                            //
// Burkhard Militzer                                        Urbana 4-9-99     //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#ifndef _STANDARD_
#define _STANDARD_

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

// Define RANGE_CHECKING for testing
// #define RANGE_CHECKING

////////////////////////////////////////////////

std::string IntToString(const int i);
std::string DoubleToString(const double d);
int StringToInt(const std::string & s);
double StringToDouble(const std::string & s);
std::string UpperCase(const std::string & s);
std::string LowerCase(const std::string & s);

#ifdef __PGI // PG compiler bug
inline bool getline( std::istream & is, std::string & s) {
  s="";
  char c;
  while(is.get(c)) {
    if (c=='\n') return true;
    s += c;
  }
  return false;
}
inline bool getline(ifstream & is, std::string & s) {
  s="";
  char c;
  while(is.get(c)) {
    if (c=='\n') return true;
    s += c;
  }
  return false;
}
#endif

////////// redefine std::cout to write to a file //////////////
// #ifndef NO_COUT
// #define std::cout COUT
// extern std::ofstream COUT;
// #endif // NO_COUT

////////////////////////////////////////////////

#include <cmath>

#ifndef pi
const double pi=3.14159265358979323846;
#endif

// Is required for KCC but not on IBM using KCC3.3
// laters ones require it too.
/* #ifndef __rs6000__  */
/* #ifndef __osf__  */
#if !defined __GNUC__ || (__GNUC__ == 2 || (__GNUC__ == 3 && __GNUC_MINOR__ < 1))
//inline double std::abs(double x) {  
//   return std::abs(x);  
//} 
#endif
/* #endif  */
/* #endif  */

// #define double long double

inline double sign(double x) {
  return (x>0.0) ? 1.0 : ((x<0.0) ? -1.0 : 0.0);
}

inline int sign(int x) {
  return (x>0) ? 1 : ((x<0) ? -1 : 0);
}

inline double nint(const double x) {
  return int(x+0.5*sign(x));
}

inline double min(double x, double y) {
  return (x<=y) ? x : y;
}

inline int min(int x, int y) {
  return (x<=y) ? x : y;
}

inline double max(double x, double y) {
  return (x>=y) ? x : y;
}

inline int max(int x, int y) {
  return (x>=y) ? x : y;
}

inline double sqr(double x) {
  return (x*x);
}

inline int sqr(int x) {
  return (x*x);
}

///////////////////////////////////////////////////////////////////////////

// Write name fo the variable and its value
#define write1(i)                    {std::cout << " "#i"= " << i; }
#define write2(i,j)                  {write1(i); write1(j);}
#define write3(i,j,k)                {write2(i,j); write1(k); }
#define write4(i,j,k,l)              {write3(i,j,k); write1(l); }
#define write5(i,j,k,l,m)            {write4(i,j,k,l); write1(m); }
#define write6(i,j,k,l,m,n)          {write5(i,j,k,l,m); write1(n); }
#define write7(i,j,k,l,m,n,o)        {write6(i,j,k,l,m,n); write1(o); }
#define write8(i,j,k,l,m,n,o,p)      {write7(i,j,k,l,m,n,o); write1(p); }
#define write9(i,j,k,l,m,n,o,p,q)    {write8(i,j,k,l,m,n,o,p); write1(q); }
#define write10(i,j,k,l,m,n,o,p,q,r) {write9(i,j,k,l,m,n,o,p,q); write1(r); }

#define BMWrite(i)                     {write1(i); std::cout << std::endl;}
#define BMWrite2(i,j)                  {write2(i,j); std::cout << std::endl;}
#define BMWrite3(i,j,k)                {write3(i,j,k); std::cout << std::endl;}
#define BMWrite4(i,j,k,l)              {write4(i,j,k,l); std::cout << std::endl;}
#define BMWrite5(i,j,k,l,m)            {write5(i,j,k,l,m); std::cout << std::endl;}
#define BMWrite6(i,j,k,l,m,n)          {write6(i,j,k,l,m,n); std::cout << std::endl;}
#define BMWrite7(i,j,k,l,m,n,o)        {write7(i,j,k,l,m,n,o); std::cout << std::endl;}
#define BMWrite8(i,j,k,l,m,n,o,p)      {write8(i,j,k,l,m,n,o,p); std::cout << std::endl;}
#define BMWrite9(i,j,k,l,m,n,o,p,q)    {write9(i,j,k,l,m,n,o,p,q); std::cout << std::endl;}
#define BMWrite10(i,j,k,l,m,n,o,p,q,r) {write10(i,j,k,l,m,n,o,p,q,r); std::cout << std::endl;}

void Terminate();

inline void WriteError(std::ostringstream & ss) {
  const std::string errorString = "Error   ";  
  //  ss << ends;
  // std::cout is redirect into a file which might not yet be opened
  if (std::cout) { 
    std::cout.precision(16);
    std::cout << errorString << ss.str() << std::endl;
  }
  std::cerr.precision(16);
  std::cerr << errorString << ss.str() << std::endl;
  Terminate();
}

inline void error(char* m){
  std::ostringstream ss;
  ss << m;
  WriteError(ss);
}

template<class T> inline 
void error(char* m, const T& n){
  std::ostringstream ss;
  ss << m << " " << n;
  WriteError(ss);
}

template<class T, class U> inline
void error(char* m, const T& t, const U& u){
  std::ostringstream ss;
  ss << m << " " << t << " " << u;
  WriteError(ss);
}

template<class T, class U, class V> inline
void error(char* m, const T& t, const U& u, const V& v){
  std::ostringstream ss;
  ss << m << " " << t << " " << u << " " << v;
  WriteError(ss);
}

template<class T, class U, class V, class W> inline
void error(char* m, const T& t, const U& u, const V& v, const W& w){
  std::ostringstream ss;
  ss << m << " " << t << " " << u << " " << v << " " << w;
  WriteError(ss);
}

template<class T, class U, class V, class W, class X> inline
void error(char* m, const T& t, const U& u, const V& v, const W& w, const X& x){
  std::ostringstream ss;
  ss << m << " " << t << " " << u << " " << v << " " << w << " " << x;
  WriteError(ss);
}

template<class T, class U, class V, class W, class X, class Y> inline
void error(char* m, const T& t, const U& u, const V& v, const W& w, const X& x, const Y& y){
  std::ostringstream ss;
  ss << m << " " << t << " " << u << " " << v << " " << w << " " << x << " " << y;
  WriteError(ss);
}

template<class T, class U, class V, class W, class X, class Y, class Z> inline
void error(char* m, const T& t, const U& u, const V& v, const W& w, const X& x, const Y& y, const Z& z){
  std::ostringstream ss;
  ss << m << " " << t << " " << u << " " << v << " " << w << " " << x << " " << y << " " << z;
  WriteError(ss);
}

inline void WriteWarning(std::ostringstream & ss) {
  const std::string warningString = "WARNING   ";  
  //  ss << ends;
  std::cout << warningString << ss.str() << std::endl;
  //  std::cerr << warningString << ss.str() << std::endl;
}

inline void warning() {
  std::ostringstream ss;
  ss << "...";
  WriteWarning(ss);
}

inline void warning(char* m){
  std::ostringstream ss;
  ss << m;
  WriteWarning(ss);
}

template<class T> inline 
void warning(const char* m, const T& t){
  std::ostringstream ss;
  ss << m << " " << t;
  WriteWarning(ss);
}

template<class T, class U> inline
void warning(const char* m, const T& t, const U& u){
  std::ostringstream ss;
  ss << m << " " << t << " " << u;
  WriteWarning(ss);
}

template<class T, class U, class V> inline
void warning(const char* m, const T& t, const U& u, const V& v){
  std::ostringstream ss;
  ss << m << " " << t << " " << u << " " << v;
  WriteWarning(ss);
}

template<class T, class U, class V, class W> inline
void warning(const char* m, const T& t, const U& u, const V& v, const W& w){
  std::ostringstream ss;
  ss << m << " " << t << " " << u << " " << v << " " << w;
  WriteWarning(ss);
}

template<class T, class U, class V, class W, class X> inline
void warning(const char* m, const T& t, const U& u, const V& v, const W& w, const X& x){
  std::ostringstream ss;
  ss << m << " " << t << " " << u << " " << v << " " << w << " " << x;
  WriteWarning(ss);
}

template<class T, class U, class V, class W, class X, class Y> inline
void warning(const char* m, const T& t, const U& u, const V& v, const W& w, const X& x, const Y& y){
  std::ostringstream ss;
  ss << m << " " << t << " " << u << " " << v << " " << w << " " << x << " " << y;
  WriteWarning(ss);
}

template<class T, class U, class V, class W, class X, class Y, class Z> inline
void warning(const char* m, const T& t, const U& u, const V& v, const W& w, const X& x, const Y& y, const Z& z){
  std::ostringstream ss;
  ss << m << " " << t << " " << u << " " << v << " " << w << " " << x << " " << y << " " << z;
  WriteWarning(ss);
}

////////////////////// Functions for array bound checking //////////////////////

// Limits is an inline function that checks indices 
// ok is 0<= n < max
inline void Limits(const int n, const int max) {
#ifdef RANGE_CHECKING
  if ((n<0) || (n>=max)) {
    error("Array Index out of range ",n,max);
    std::cerr << "Array error: Index out of range:  0<= " 
	 << n << " < " << max << "\n" ;
    Terminate();
  }
#endif // RANGE_CHECKING
}

// Limits is an inline function that checks indices 
// ok is 0<= n <= max
inline void LimitsInclusive(const int n, const int max) {
#ifdef RANGE_CHECKING
  if ((n<0) || (n>max)) {
    //    error("Array Error: Index out of range ",n,max);
    std::cerr << "Array error: Upper limit for index out of range:  0<= " 
	 << n << " <= " << max << "\n" ;
    Terminate();
  }
#endif // RANGE_CHECKING
}

inline void EqualLimits(const int max1, const int max2) {
#ifdef RANGE_CHECKING
  if (max1!=max2) {
    std::cerr << "Array copy error: array sizes not equal:" 
	 << max1 << "," << max2 << std::endl;
    Terminate();
  }
#endif // RANGE_CHECKING
}

inline void BiggerLimit(const int lower, const int upper) {
#ifdef RANGE_CHECKING
  if (lower>=upper) {
    std::cerr << "Sub-array limits error: lower limit not lower " 
	 << lower << "," << upper << std::endl;
    Terminate();
  }
#endif // RANGE_CHECKING
}


#endif // _STANDARD_
