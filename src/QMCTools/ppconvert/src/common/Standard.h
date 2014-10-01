/////////////////////////////////////////////////////////////
// Copyright (C) 2003-2006 Bryan Clark and Kenneth Esler   //
//                                                         //
// This program is free software; you can redistribute it  //
// and/or modify it under the terms of the GNU General     //
// Public License as published by the Free Software        //
// Foundation; either version 2 of the License, or         //
// (at your option) any later version.  This program is    //
// distributed in the hope that it will be useful, but     //
// WITHOUT ANY WARRANTY; without even the implied warranty //
// of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. //  
// See the GNU General Public License for more details.    //
// For more information, please see the PIMC++ Home Page:  //
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
using namespace std;

// Define RANGE_CHECKING for testing
// #define RANGE_CHECKING

////////////////////////////////////////////////

string IntToString(const int i);
string DoubleToString(const double d);
int StringToInt(const string & s);
double StringToDouble(const string & s);
string UpperCase(const string & s);
string LowerCase(const string & s);

#ifdef __PGI // PG compiler bug
inline bool getline(istream & is, string & s) {
  s="";
  char c;
  while(is.get(c)) {
    if (c=='\n') return true;
    s += c;
  }
  return false;
}
inline bool getline(ifstream & is, string & s) {
  s="";
  char c;
  while(is.get(c)) {
    if (c=='\n') return true;
    s += c;
  }
  return false;
}
#endif

////////// redefine cout to write to a file //////////////
// #ifndef NO_COUT
// #define cout COUT
// extern ofstream COUT;
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
//inline double abs(double x) {  
//   return fabs(x);  
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
#define write1(i)                    {cout << " "#i"= " << i; }
#define write2(i,j)                  {write1(i); write1(j);}
#define write3(i,j,k)                {write2(i,j); write1(k); }
#define write4(i,j,k,l)              {write3(i,j,k); write1(l); }
#define write5(i,j,k,l,m)            {write4(i,j,k,l); write1(m); }
#define write6(i,j,k,l,m,n)          {write5(i,j,k,l,m); write1(n); }
#define write7(i,j,k,l,m,n,o)        {write6(i,j,k,l,m,n); write1(o); }
#define write8(i,j,k,l,m,n,o,p)      {write7(i,j,k,l,m,n,o); write1(p); }
#define write9(i,j,k,l,m,n,o,p,q)    {write8(i,j,k,l,m,n,o,p); write1(q); }
#define write10(i,j,k,l,m,n,o,p,q,r) {write9(i,j,k,l,m,n,o,p,q); write1(r); }

#define BMWrite(i)                     {write1(i); cout << endl;}
#define BMWrite2(i,j)                  {write2(i,j); cout << endl;}
#define BMWrite3(i,j,k)                {write3(i,j,k); cout << endl;}
#define BMWrite4(i,j,k,l)              {write4(i,j,k,l); cout << endl;}
#define BMWrite5(i,j,k,l,m)            {write5(i,j,k,l,m); cout << endl;}
#define BMWrite6(i,j,k,l,m,n)          {write6(i,j,k,l,m,n); cout << endl;}
#define BMWrite7(i,j,k,l,m,n,o)        {write7(i,j,k,l,m,n,o); cout << endl;}
#define BMWrite8(i,j,k,l,m,n,o,p)      {write8(i,j,k,l,m,n,o,p); cout << endl;}
#define BMWrite9(i,j,k,l,m,n,o,p,q)    {write9(i,j,k,l,m,n,o,p,q); cout << endl;}
#define BMWrite10(i,j,k,l,m,n,o,p,q,r) {write10(i,j,k,l,m,n,o,p,q,r); cout << endl;}

void Terminate();

inline void WriteError(ostringstream & ss) {
  const string errorString = "Error   ";  
  //  ss << ends;
  // cout is redirect into a file which might not yet be opened
  if (cout) { 
    cout.precision(16);
    cout << errorString << ss.str() << endl;
  }
  cerr.precision(16);
  cerr << errorString << ss.str() << endl;
  Terminate();
}

inline void error(char* m){
  ostringstream ss;
  ss << m;
  WriteError(ss);
}

template<class T> inline 
void error(char* m, const T& n){
  ostringstream ss;
  ss << m << " " << n;
  WriteError(ss);
}

template<class T, class U> inline
void error(char* m, const T& t, const U& u){
  ostringstream ss;
  ss << m << " " << t << " " << u;
  WriteError(ss);
}

template<class T, class U, class V> inline
void error(char* m, const T& t, const U& u, const V& v){
  ostringstream ss;
  ss << m << " " << t << " " << u << " " << v;
  WriteError(ss);
}

template<class T, class U, class V, class W> inline
void error(char* m, const T& t, const U& u, const V& v, const W& w){
  ostringstream ss;
  ss << m << " " << t << " " << u << " " << v << " " << w;
  WriteError(ss);
}

template<class T, class U, class V, class W, class X> inline
void error(char* m, const T& t, const U& u, const V& v, const W& w, const X& x){
  ostringstream ss;
  ss << m << " " << t << " " << u << " " << v << " " << w << " " << x;
  WriteError(ss);
}

template<class T, class U, class V, class W, class X, class Y> inline
void error(char* m, const T& t, const U& u, const V& v, const W& w, const X& x, const Y& y){
  ostringstream ss;
  ss << m << " " << t << " " << u << " " << v << " " << w << " " << x << " " << y;
  WriteError(ss);
}

template<class T, class U, class V, class W, class X, class Y, class Z> inline
void error(char* m, const T& t, const U& u, const V& v, const W& w, const X& x, const Y& y, const Z& z){
  ostringstream ss;
  ss << m << " " << t << " " << u << " " << v << " " << w << " " << x << " " << y << " " << z;
  WriteError(ss);
}

inline void WriteWarning(ostringstream & ss) {
  const string warningString = "WARNING   ";  
  //  ss << ends;
  cout << warningString << ss.str() << endl;
  //  cerr << warningString << ss.str() << endl;
}

inline void warning() {
  ostringstream ss;
  ss << "...";
  WriteWarning(ss);
}

inline void warning(char* m){
  ostringstream ss;
  ss << m;
  WriteWarning(ss);
}

template<class T> inline 
void warning(const char* m, const T& t){
  ostringstream ss;
  ss << m << " " << t;
  WriteWarning(ss);
}

template<class T, class U> inline
void warning(const char* m, const T& t, const U& u){
  ostringstream ss;
  ss << m << " " << t << " " << u;
  WriteWarning(ss);
}

template<class T, class U, class V> inline
void warning(const char* m, const T& t, const U& u, const V& v){
  ostringstream ss;
  ss << m << " " << t << " " << u << " " << v;
  WriteWarning(ss);
}

template<class T, class U, class V, class W> inline
void warning(const char* m, const T& t, const U& u, const V& v, const W& w){
  ostringstream ss;
  ss << m << " " << t << " " << u << " " << v << " " << w;
  WriteWarning(ss);
}

template<class T, class U, class V, class W, class X> inline
void warning(const char* m, const T& t, const U& u, const V& v, const W& w, const X& x){
  ostringstream ss;
  ss << m << " " << t << " " << u << " " << v << " " << w << " " << x;
  WriteWarning(ss);
}

template<class T, class U, class V, class W, class X, class Y> inline
void warning(const char* m, const T& t, const U& u, const V& v, const W& w, const X& x, const Y& y){
  ostringstream ss;
  ss << m << " " << t << " " << u << " " << v << " " << w << " " << x << " " << y;
  WriteWarning(ss);
}

template<class T, class U, class V, class W, class X, class Y, class Z> inline
void warning(const char* m, const T& t, const U& u, const V& v, const W& w, const X& x, const Y& y, const Z& z){
  ostringstream ss;
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
    cerr << "Array error: Index out of range:  0<= " 
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
    cerr << "Array error: Upper limit for index out of range:  0<= " 
	 << n << " <= " << max << "\n" ;
    Terminate();
  }
#endif // RANGE_CHECKING
}

inline void EqualLimits(const int max1, const int max2) {
#ifdef RANGE_CHECKING
  if (max1!=max2) {
    cerr << "Array copy error: array sizes not equal:" 
	 << max1 << "," << max2 << endl;
    Terminate();
  }
#endif // RANGE_CHECKING
}

inline void BiggerLimit(const int lower, const int upper) {
#ifdef RANGE_CHECKING
  if (lower>=upper) {
    cerr << "Sub-array limits error: lower limit not lower " 
	 << lower << "," << upper << endl;
    Terminate();
  }
#endif // RANGE_CHECKING
}


#endif // _STANDARD_
