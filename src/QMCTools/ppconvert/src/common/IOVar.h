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

#ifndef IO_VAR_H
#define IO_VAR_H

#include "IOVarASCII.h"

namespace IO {

  template<typename T> bool
  IOVarBase::Read(T &val)
  {
    if (GetFileType() == ASCII_TYPE) {
      IOVarASCII<T,0>* newVar = dynamic_cast<IOVarASCII<T,0>*>(this); 
      if (newVar == NULL) {
	cerr << "Error in dynamic cast to IOVarASCII.\n";
	abort();
      }
      return newVar->VarRead(val);
    }
    else {
      std::cerr << "Error:  unknown type in IOVarBase::Read().\n";
      abort();
    }
  }

  template<typename T, int LEN> bool
  IOVarBase::Read(TinyVector<T,LEN> &val) 
  {
    if (GetFileType() == ASCII_TYPE) {
      IOVarASCII<T,1>* newVar = dynamic_cast<IOVarASCII<T,1>*>(this); 
      if (newVar == NULL) {
	cerr << "Error in dynamic cast to IOVarASCII.\n";
	abort();
      }
      Array<T,1> aVal;
      bool success = newVar->VarRead (aVal);
      if (!success)
	return false;
      else if (aVal.size() == LEN) 
	for (int i=0; i<LEN; i++)
	  val[i] = aVal(i);
      else
	return false;
      return true;
    }
    else {
      std::cerr << "Error:  unknown type in IOVarBase::Read().\n";
      abort();
    }
  }

  template<typename T, int RANK> bool
  IOVarBase::Read(blitz::Array<T,RANK> &val)
  {
    if (GetFileType() == ASCII_TYPE) {
      IOVarASCII<T,RANK>* newVar = dynamic_cast<IOVarASCII<T,RANK>*>(this); 
      if (newVar == NULL) {
	cerr << "Error in dynamic cast to IOVarASCII.\n";
	abort();
      }
      return newVar->VarRead(val);
    }
    else {
      std::cerr << "Error:  unknown type in IOVarBase::Read().\n";
      abort();
    }
  }	


  template<typename T,  int RANK, typename T0, typename T1, typename T2, 
	   typename T3, typename T4, typename T5, typename T6, typename T7, 
	   typename T8, typename T9, typename T10> bool
  IOVarBase::Read(blitz::Array<T,RANK> &val, T0 s0, T1 s1, T2 s2, T3 s3, T4 s4,
		  T5 s5, T6 s6, T7 s7, T8 s8, T9 s9, T10 s10)
  {
    static const int numSlices = 
      SliceCheck<T0>::isSlice+SliceCheck<T1>::isSlice+SliceCheck<T2>::isSlice+
      SliceCheck<T3>::isSlice+SliceCheck<T4>::isSlice+SliceCheck<T5>::isSlice+
      SliceCheck<T6>::isSlice+SliceCheck<T7>::isSlice+SliceCheck<T8>::isSlice+
      SliceCheck<T9>::isSlice+SliceCheck<T10>::isSlice;
  
    /// The rank of the array must be the rank of the IO variable minus
    /// the number of slices by integer singlet ranges.
    static const int varRank=numSlices+RANK;
  
    if (GetFileType() == ASCII_TYPE) {
      IOVarASCII<T,varRank>* newVar = dynamic_cast<IOVarASCII<T,varRank>*>(this); 
      if (newVar == NULL) {
	cerr << "Error in dynamic cast to IOVarHDF5.\n";
	abort();
      }
      return newVar->Slice(s0, s1, s2, s3, s4, s5, s6, s7, s8, s9, s10).VarRead(val);
    }
  }	


  template<typename T,  int RANK, typename T0, typename T1, typename T2, 
	   typename T3, typename T4, typename T5, typename T6, typename T7, 
	   typename T8, typename T9, typename T10> bool
  IOVarBase::Write(const blitz::Array<T,RANK> &val, T0 s0, T1 s1, T2 s2, T3 s3, T4 s4,
		   T5 s5, T6 s6, T7 s7, T8 s8, T9 s9, T10 s10)
  {
    static const int numSlices = 
      SliceCheck<T0>::isSlice+SliceCheck<T1>::isSlice+SliceCheck<T2>::isSlice+
      SliceCheck<T3>::isSlice+SliceCheck<T4>::isSlice+SliceCheck<T5>::isSlice+
      SliceCheck<T6>::isSlice+SliceCheck<T7>::isSlice+SliceCheck<T8>::isSlice+
      SliceCheck<T9>::isSlice+SliceCheck<T10>::isSlice;
    
    /// The rank of the array must be the rank of the IO variable minus
    /// the number of slices by integer singlet ranges.
    static const int varRank=numSlices+RANK;

    if (GetFileType() == ASCII_TYPE) {
      IOVarASCII<T,varRank>* newVar = dynamic_cast<IOVarASCII<T,varRank>*>(this); 
      if (newVar == NULL) {
	cerr << "Error in dynamic cast to IOVarASCII.\n";
	abort();
      }
      return newVar->Slice(s0, s1, s2, s3, s4, s5, s6, s7, s8, s9, s10).VarWrite(val);
    }

  }

  template<typename T> bool 
  IOVarBase::Append (const T val)
  {
    assert (GetRank() == 1);
    int n = GetExtent(0);
    Resize(n+1);
    Array<T,1> v(1);
    v(0) = val;
    return Write (v, Range(n,n));
  }

//   template<typename T, int RANK> bool 
//   IOVarBase::Append(Array<T,RANK> &val)
//   {
//     assert (GetRank() == (RANK+1));
//     int n = GetExtent(0);
//     for (int i=0; i<RANK; i++)
//       assert (val.extent(i) == GetExtent(i+1));

//     if (GetFileType() == HDF5_TYPE) {
//       IOVarHDF5<T,RANK+1>* newVar = dynamic_cast<IOVarHDF5<T,RANK+1>*>(this);
//       if (newVar == NULL) {
// 	cerr << "Error in dynamic cast to IOVarHDF5.\n";
// 	abort();
//       }
//       newVar->Resize(n+1);
//       return newVar->Slice(n, n0, n0, n0, n0, n0, n0, n0, n0, n0, n0).VarWrite(val);
//     }
//     else if (GetFileType() == ASCII_TYPE) {
//       IOVarASCII<T,RANK+1>* newVar = dynamic_cast<IOVarASCII<T,RANK+1>*>(this);
//       if (newVar == NULL) {
// 	cerr << "Error in dynamic cast to IOVarASCII.\n";
// 	abort();
//       }
//       newVar->Resize(n+1);
//       return newVar->Slice(n, n0, n0, n0, n0, n0, n0, n0, n0, n0, n0).VarWrite(val);
//     }
//   }
  
  template<class T> bool
  IOVarBase::Append(const blitz::Array<T,1> &val)
  {
    assert (GetRank()==2);
    int n = GetExtent(0);
    assert (val.extent(0) == GetExtent(1));
    Resize(n+1);
    return Write(val, n, Range::all());
  } 

  template<class T> bool
  IOVarBase::Append(const blitz::Array<T,2> &val)
  {
    assert (GetRank()==3);
    int n = GetExtent(0);
    assert (val.extent(0) == GetExtent(1));
    assert (val.extent(1) == GetExtent(2));
    Resize(n+1);
    return Write(val, n, Range::all(), Range::all());
  } 

  template<class T> bool
  IOVarBase::Append(const blitz::Array<T,3> &val)
  {
    assert (GetRank()==4);
    int n = GetExtent(0);
    assert (val.extent(0) == GetExtent(1));
    assert (val.extent(1) == GetExtent(2));
    assert (val.extent(2) == GetExtent(3));
    Resize(n+1);
    return Write(val, n, Range::all(), Range::all(), Range::all());
  } 

  template<class T> bool
  IOVarBase::Append(const blitz::Array<T,4> &val)
  {
    assert (GetRank()==5);
    int n = GetExtent(0);
    assert (val.extent(0) == GetExtent(1));
    assert (val.extent(1) == GetExtent(2));
    assert (val.extent(2) == GetExtent(3));
    assert (val.extent(3) == GetExtent(4));
    Resize(n+1);
    return Write(val, n, Range::all(), Range::all(), Range::all(), Range::all());
  } 


}       /// Ends namespace IO


#endif  /// Ends ifndef IO_VAR_H
