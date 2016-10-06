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

#ifndef IO_VAR_BASE_H
#define IO_VAR_BASE_H

#include "blitz/array.h"
#include <iostream>


namespace IO {
using namespace blitz;
  typedef enum { DOUBLE_TYPE, INT_TYPE, STRING_TYPE, BOOL_TYPE, COMPLEX_TYPE, INVALID } IODataType;
  typedef enum { ASCII_TYPE} IOFileType;
  
  ///////////////////////////////////////////////////////////////////
  /// The following is a template trick for counting how many     ///
  /// dimension reductions we have made to a dataset by indexing  ///
  /// by integer arguements.                                      ///
  ///////////////////////////////////////////////////////////////////
  template<typename T> class SliceCheck
  { public:  static const int isSlice = 0; };

  template<> class SliceCheck<int>
  { public:  static const int isSlice = 1; };

  class IOVarBase
  {
  private:
    nilArraySection n0;
  protected:
    std::string Name;
  public:
    virtual int GetRank()            = 0;
    virtual IODataType GetType()     = 0;
    virtual IOFileType GetFileType() = 0;
    virtual void Print(std::ofstream &out) { /* do nothing */ }
    //    virtual std::string GetTypeString()   = 0;
    std::string GetName () const
    { return Name; }
    virtual int GetExtent(int dim)   = 0;

    /// Resizes the first dimension of the variable
    virtual void Resize(int n)       = 0;

    //////////////////////
    /// Read functions ///
    //////////////////////
    template<typename T> bool Read(T &val);
    template<typename T, int  LEN> bool Read(TinyVector<T,LEN> &val);
    template<typename T, int RANK> bool Read(blitz::Array<T,RANK> &val);

    template<typename T,  int RANK,    typename T0, typename T1, typename T2,
	     typename T3, typename T4, typename T5, typename T6, typename T7,
	     typename T8, typename T9, typename T10>
    bool Read(blitz::Array<T,RANK> &val, T0 s0, T1 s1, T2 s2, T3 s3, T4 s4, T5 s5, T6 s6,
	      T7 s7, T8 s8, T9 s9, T10 s10);
    
    template<typename T,  int RANK,    typename T0, typename T1, typename T2,
	     typename T3, typename T4, typename T5, typename T6, typename T7,
	     typename T8, typename T9>
    bool Read(blitz::Array<T,RANK> &val, T0 s0, T1 s1, T2 s2, T3 s3, T4 s4, T5 s5, T6 s6,
	      T7 s7, T8 s8, T9 s9) 
    { return Read(val, s0, s1, s2, s3, s4, s5, s6, s7, s8, s9, n0); }
    
    template<typename T,  int RANK,    typename T0, typename T1, typename T2,
	     typename T3, typename T4, typename T5, typename T6, typename T7,
	     typename T8>
    bool Read(blitz::Array<T,RANK> &val, T0 s0, T1 s1, T2 s2, T3 s3, T4 s4, T5 s5, T6 s6,
	      T7 s7, T8 s8) 
    { return Read(val, s0, s1, s2, s3, s4, s5, s6, s7, s8, n0, n0); }
    
    template<typename T,  int RANK,    typename T0, typename T1, typename T2,
	     typename T3, typename T4, typename T5, typename T6, typename T7>
    bool Read(blitz::Array<T,RANK> &val, T0 s0, T1 s1, T2 s2, T3 s3, T4 s4, T5 s5, T6 s6,
	      T7 s7) 
    { return Read(val, s0, s1, s2, s3, s4, s5, s6, s7, n0, n0, n0); }
    
    template<typename T,  int RANK,    typename T0, typename T1, typename T2,
	     typename T3, typename T4, typename T5, typename T6>
    bool Read(blitz::Array<T,RANK> &val, T0 s0, T1 s1, T2 s2, T3 s3, T4 s4, T5 s5, T6 s6) 
    { return Read(val, s0, s1, s2, s3, s4, s5, s6, n0, n0, n0, n0); }
    
    template<typename T,  int RANK,    typename T0, typename T1, typename T2,
	     typename T3, typename T4, typename T5>
    bool Read(blitz::Array<T,RANK> &val, T0 s0, T1 s1, T2 s2, T3 s3, T4 s4, T5 s5) 
    { return Read(val, s0, s1, s2, s3, s4, s5, n0, n0, n0, n0, n0); }
    
    template<typename T,  int RANK,    typename T0, typename T1, typename T2,
	     typename T3, typename T4>
    bool Read(blitz::Array<T,RANK> &val, T0 s0, T1 s1, T2 s2, T3 s3, T4 s4) 
    { return Read(val, s0, s1, s2, s3, s4, n0, n0, n0, n0, n0, n0); }
    
    template<typename T,  int RANK,    typename T0, typename T1, typename T2,
	     typename T3>
    bool Read(blitz::Array<T,RANK> &val, T0 s0, T1 s1, T2 s2, T3 s3) 
    { return Read(val, s0, s1, s2, s3, n0, n0, n0, n0, n0, n0, n0); }
    
    template<typename T,  int RANK,    typename T0, typename T1, typename T2>
    bool Read(blitz::Array<T,RANK> &val, T0 s0, T1 s1, T2 s2) 
    { return Read(val, s0, s1, s2, n0, n0, n0, n0, n0, n0, n0, n0); }
    
    template<typename T,  int RANK,    typename T0, typename T1>
    bool Read(blitz::Array<T,RANK> &val, T0 s0, T1 s1) 
    { return Read(val, s0, s1, n0, n0, n0, n0, n0, n0, n0, n0, n0); }
    
    template<typename T,  int RANK,    typename T0>
    bool Read(blitz::Array<T,RANK> &val, T0 s0) 
    { return Read(val, s0, n0, n0, n0, n0, n0, n0, n0, n0, n0, n0); }
    

    ///////////////////////
    /// Write functions ///
    ///////////////////////    
    template<typename T,  int RANK,    typename T0, typename T1, typename T2,
	     typename T3, typename T4, typename T5, typename T6, typename T7,
	     typename T8, typename T9, typename T10>
    bool Write(const blitz::Array<T,RANK> &val, T0 s0, T1 s1, T2 s2, T3 s3, T4 s4, T5 s5, T6 s6,
	       T7 s7, T8 s8, T9 s9, T10 s10);
    
    template<typename T,  int RANK,    typename T0, typename T1, typename T2,
	     typename T3, typename T4, typename T5, typename T6, typename T7,
	     typename T8, typename T9>
    bool Write(const blitz::Array<T,RANK> &val, T0 s0, T1 s1, T2 s2, T3 s3, T4 s4, T5 s5, T6 s6,
	       T7 s7, T8 s8, T9 s9) 
    { return Write(val, s0, s1, s2, s3, s4, s5, s6, s7, s8, s9, n0); }
    
    template<typename T,  int RANK,    typename T0, typename T1, typename T2,
	     typename T3, typename T4, typename T5, typename T6, typename T7,
	     typename T8>
    bool Write(const blitz::Array<T,RANK> &val, T0 s0, T1 s1, T2 s2, T3 s3, T4 s4, T5 s5, T6 s6,
	       T7 s7, T8 s8) 
    { return Write(val, s0, s1, s2, s3, s4, s5, s6, s7, s8, n0, n0); }
    
    template<typename T,  int RANK,    typename T0, typename T1, typename T2,
	     typename T3, typename T4, typename T5, typename T6, typename T7>
    bool Write(const blitz::Array<T,RANK> &val, T0 s0, T1 s1, T2 s2, T3 s3, T4 s4, T5 s5, T6 s6,
	       T7 s7) 
    { return Write(val, s0, s1, s2, s3, s4, s5, s6, s7, n0, n0, n0); }
    
    template<typename T,  int RANK,    typename T0, typename T1, typename T2,
	     typename T3, typename T4, typename T5, typename T6>
    bool Write(const blitz::Array<T,RANK> &val, T0 s0, T1 s1, T2 s2, T3 s3, T4 s4, T5 s5, T6 s6) 
    { return Write(val, s0, s1, s2, s3, s4, s5, s6, n0, n0, n0, n0); }
    
    template<typename T,  int RANK,    typename T0, typename T1, typename T2,
	     typename T3, typename T4, typename T5>
    bool Write(const blitz::Array<T,RANK> &val, T0 s0, T1 s1, T2 s2, T3 s3, T4 s4, T5 s5) 
    { return Write(val, s0, s1, s2, s3, s4, s5, n0, n0, n0, n0, n0); }

    template<typename T,  int RANK,    typename T0, typename T1, typename T2,
	     typename T3, typename T4>
    bool Write(const blitz::Array<T,RANK> &val, T0 s0, T1 s1, T2 s2, T3 s3, T4 s4) 
    { return Write(val, s0, s1, s2, s3, s4, n0, n0, n0, n0, n0, n0); }

    template<typename T,  int RANK,    typename T0, typename T1, typename T2,
	     typename T3>
    bool Write(const blitz::Array<T,RANK> &val, T0 s0, T1 s1, T2 s2, T3 s3) 
    { return Write(val, s0, s1, s2, s3, n0, n0, n0, n0, n0, n0, n0); }

    template<typename T,  int RANK,    typename T0, typename T1, typename T2>
    bool Write(const blitz::Array<T,RANK> &val, T0 s0, T1 s1, T2 s2) 
    { return Write(val, s0, s1, s2, n0, n0, n0, n0, n0, n0, n0, n0); }

    template<typename T,  int RANK,    typename T0, typename T1>
    bool Write(const blitz::Array<T,RANK> &val, T0 s0, T1 s1) 
    { return Write(val, s0, s1, n0, n0, n0, n0, n0, n0, n0, n0, n0); }

    template<typename T,  int RANK,    typename T0>
    bool Write(const blitz::Array<T,RANK> &val, T0 s0) 
    { return Write(val, s0, n0, n0, n0, n0, n0, n0, n0, n0, n0, n0); }
  
    template<typename T>
    bool Append (const T val);
    template<class T>
    bool Append(const blitz::Array<T,1> &val);
    template<class T>
    bool Append(const blitz::Array<T,2> &val);
    template<class T>
    bool Append(const blitz::Array<T,3> &val);
    template<class T>
    bool Append(const blitz::Array<T,4> &val);

    virtual ~IOVarBase() { /* do nothing at this level */ }
  };


  ///////////////////////////////////////////////////////////////////
  /// The following are template tricks for converting a C++ type ///
  /// into an enumerated type variable.                           ///
  ///////////////////////////////////////////////////////////////////
  template<typename T> class TypeConvert
  { public: static const IODataType Type = INVALID; };

  template<> class TypeConvert<double>
  { public: static const IODataType Type = DOUBLE_TYPE; };

  template<> class TypeConvert<int>
  { public: static const IODataType Type = INT_TYPE; };

  template<> class TypeConvert<std::string>
  { public: static const IODataType Type = STRING_TYPE; };

  template<> class TypeConvert<bool>
  { public: static const IODataType Type = BOOL_TYPE; };

  template<> class TypeConvert<std::complex<double> >
  { public: static const IODataType Type = COMPLEX_TYPE; };
  



  template<typename T> inline std::string TypeString(T &val)     {  return "unknown"; }
  template<> std::string inline TypeString(double &val)          {  return "double";  }
  template<> std::string inline TypeString(int &val)             {  return "int";     }
  template<> std::string inline TypeString( std::string &val)          {  return "string";  }
  template<> std::string inline TypeString(bool &val)            {  return "bool";    }
  template<> std::string inline TypeString(std::complex<double> &val) {  return "complex";    }

}



#endif // ifndef IO_BASE_H
