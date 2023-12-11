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

#ifndef IO_VAR_ASCII_H
#define IO_VAR_ASCII_H

#include "IOVarBase.h"
#include <fstream>

namespace IO
{
template<typename T, int RANK>
class IOVarASCII;

//   template<typename T,  typename T0, typename T1, typename T2, typename T3, typename T4,
// 	   typename T5, typename T6, typename T7, typename T8, typename T9, typename T10>
//   class ASCIISliceMaker
//   {
//   public:
//     static const int rank =      ArraySectionInfo<T0>::rank + ArraySectionInfo<T1>::rank +
//     ArraySectionInfo<T2>::rank + ArraySectionInfo<T3>::rank + ArraySectionInfo<T4>::rank +
//     ArraySectionInfo<T5>::rank + ArraySectionInfo<T6>::rank + ArraySectionInfo<T7>::rank +
//     ArraySectionInfo<T8>::rank + ArraySectionInfo<T9>::rank + ArraySectionInfo<T10>::rank;

//     using SliceType = IOVarASCII<T,rank>;
//   };

template<typename T,
         int RANK,
         typename T0,
         typename T1,
         typename T2,
         typename T3,
         typename T4,
         typename T5,
         typename T6,
         typename T7,
         typename T8,
         typename T9,
         typename T10>
class ASCIISliceMaker
{
public:
  static const int rank = RANK -
      (SliceCheck<T0>::isSlice + SliceCheck<T1>::isSlice + SliceCheck<T2>::isSlice + SliceCheck<T3>::isSlice +
       SliceCheck<T4>::isSlice + SliceCheck<T5>::isSlice + SliceCheck<T6>::isSlice + SliceCheck<T7>::isSlice +
       SliceCheck<T8>::isSlice + SliceCheck<T9>::isSlice + SliceCheck<T10>::isSlice);

  using SliceType = IOVarASCII<T, rank>;
  using T_slice   = Array<T, rank>;
};

template<typename T, int RANK>
class IOVarASCII : public IOVarBase
{
public:
  Array<T, RANK> ArrayValue;
  template<typename T0,
           typename T1,
           typename T2,
           typename T3,
           typename T4,
           typename T5,
           typename T6,
           typename T7,
           typename T8,
           typename T9,
           typename T10>
  typename ASCIISliceMaker<T, RANK, T0, T1, T2, T3, T4, T5, T6, T7, T8, T9, T10>::SliceType Slice(T0 s0,
                                                                                                  T1 s1,
                                                                                                  T2 s2,
                                                                                                  T3 s3,
                                                                                                  T4 s4,
                                                                                                  T5 s5,
                                                                                                  T6 s6,
                                                                                                  T7 s7,
                                                                                                  T8 s8,
                                                                                                  T9 s9,
                                                                                                  T10 s10);

  int GetRank() override;
  IODataType GetType() override;
  IOFileType GetFileType() override;

  void Print(std::ofstream& out) override;

  int GetExtent(int dim) override;
  void Resize(int n) override;

  bool VarRead(Array<T, RANK>& val);
  template<class TT>
  bool VarRead(Array<TT, 1>& val)
  {
    assert(0);
    return false;
  }

  template<typename T0,
           typename T1,
           typename T2,
           typename T3,
           typename T4,
           typename T5,
           typename T6,
           typename T7,
           typename T8,
           typename T9,
           typename T10>
  bool VarRead(typename SliceInfo<T, T0, T1, T2, T3, T4, T5, T6, T7, T8, T9, T10>::T_slice& val,
               T0 s0,
               T1 s1,
               T2 s2,
               T3 s3,
               T4 s4,
               T5 s5,
               T6 s6,
               T7 s7,
               T8 s8,
               T8 s9,
               T10 s10);

  bool VarWrite(const Array<T, RANK>& val);
  template<typename T0,
           typename T1,
           typename T2,
           typename T3,
           typename T4,
           typename T5,
           typename T6,
           typename T7,
           typename T8,
           typename T9,
           typename T10>
  bool VarWrite(const typename SliceInfo<T, T0, T1, T2, T3, T4, T5, T6, T7, T8, T9, T10>::T_slice& val,
                T0 s0,
                T1 s1,
                T2 s2,
                T3 s3,
                T4 s4,
                T5 s5,
                T6 s6,
                T7 s7,
                T8 s8,
                T8 s9,
                T10 s10);

  /// Assignment operator
  IOVarASCII<T, RANK>& operator=(const IOVarASCII<T, RANK>& var)
  {
    /// This references var's value.
    ArrayValue.reference(var.ArrayValue);
  }

  IOVarASCII(std::string name, const Array<T, RANK>& val)
  {
    Name = name;
    ArrayValue.resize(val.shape());
    ArrayValue = val;
  }
  /// Default constructor
  IOVarASCII(std::string name) { Name = name; }
  /// Copy constructor
  IOVarASCII(const IOVarASCII<T, RANK>& var) : ArrayValue(var.ArrayValue)
  {
    /// This references the arguments array value.
  }
};

template<>
class IOVarASCII<double, 0> : public IOVarBase
{
public:
  double Value;

  int GetRank() override;
  IODataType GetType() override;
  IOFileType GetFileType() override;

  int GetExtent(int dim) override;
  void Resize(int n) override;

  bool VarRead(double& val)
  {
    val = Value;
    return true;
  }
  bool VarWrite(double& val);
  IOVarASCII(std::string name, double val)
  {
    Name  = name;
    Value = val;
  }
  /// Default constructor
  IOVarASCII(std::string name) { Name = name; }
};

template<>
class IOVarASCII<int, 0> : public IOVarBase
{
public:
  int Value;

  int GetRank() override;
  IODataType GetType() override;
  IOFileType GetFileType() override;

  int GetExtent(int dim) override;
  void Resize(int n) override;

  bool VarRead(int& val)
  {
    val = Value;
    return true;
  }
  bool VarWrite(int& val);
  IOVarASCII(std::string name, int val)
  {
    Name  = name;
    Value = val;
  }
  /// Default constructor
  IOVarASCII(std::string name) { Name = name; }
};

template<>
class IOVarASCII<std::string, 0> : public IOVarBase
{
public:
  std::string Value;

  int GetRank() override;
  IODataType GetType() override;
  IOFileType GetFileType() override;

  int GetExtent(int dim) override;
  void Resize(int n) override;

  bool VarRead(std::string& val);
  bool VarWrite(std::string val);
  IOVarASCII(std::string name, std::string val)
  {
    Name  = name;
    Value = val;
  }
  /// Default constructor
  IOVarASCII(std::string name) { Name = name; }
};

template<>
class IOVarASCII<bool, 0> : public IOVarBase
{
public:
  bool Value;

  int GetRank() override;
  IODataType GetType() override;
  IOFileType GetFileType() override;

  int GetExtent(int dim) override;
  void Resize(int n) override;

  bool VarRead(bool& val)
  {
    val = Value;
    return true;
  }
  bool VarWrite(bool& val);
  IOVarASCII(std::string name, bool val)
  {
    Name  = name;
    Value = val;
  }
  /// Default constructor
  IOVarASCII(std::string name) { Name = name; }
};


template<>
class IOVarASCII<std::complex<double>, 0> : public IOVarBase
{
public:
  std::complex<double> Value;

  int GetRank() override;
  IODataType GetType() override;
  IOFileType GetFileType() override;

  int GetExtent(int dim) override;
  void Resize(int n) override;

  bool VarRead(std::complex<double>& val)
  {
    val = Value;
    return true;
  }
  bool VarWrite(std::complex<double>& val);
  IOVarASCII(std::string name, std::complex<double> val)
  {
    Name  = name;
    Value = val;
  }
  /// Default constructor
  IOVarASCII(std::string name) { Name = name; }
};


template<typename T, int RANK>
inline int IOVarASCII<T, RANK>::GetRank()
{
  return RANK;
}

template<typename T, int RANK>
inline IODataType IOVarASCII<T, RANK>::GetType()
{
  return TypeConvert<T>::Type;
}

template<typename T, int RANK>
inline IOFileType IOVarASCII<T, RANK>::GetFileType()
{
  return ASCII_TYPE;
}


template<typename T, int RANK>
inline int IOVarASCII<T, RANK>::GetExtent(int dim)
{
  return ArrayValue.extent(dim);
}


template<typename T, int RANK>
inline void IOVarASCII<T, RANK>::Resize(int n)
{
  //  TinyVector<int, RANK>
  auto dims         = ArrayValue.shape();
  std::get<0>(dims) = n;
  //  dims[0]                    = n;
  ArrayValue.resizeAndPreserve(dims);
}

template<typename T, int RANK>
bool IOVarASCII<T, RANK>::VarRead(Array<T, RANK>& val)
{
  val.resize(ArrayValue.shape());
  val = ArrayValue;
  return true;
}

template<typename T, int RANK>
template<typename T0,
         typename T1,
         typename T2,
         typename T3,
         typename T4,
         typename T5,
         typename T6,
         typename T7,
         typename T8,
         typename T9,
         typename T10>
inline bool IOVarASCII<T, RANK>::VarRead(
    typename SliceInfo<T, T0, T1, T2, T3, T4, T5, T6, T7, T8, T9, T10>::T_slice& val,
    T0 s0,
    T1 s1,
    T2 s2,
    T3 s3,
    T4 s4,
    T5 s5,
    T6 s6,
    T7 s7,
    T8 s8,
    T8 s9,
    T10 s10)
{
  val.resize(ArrayValue(s0, s1, s2, s3, s4, s5, s6, s7, s8, s9, s10).shape());
  val = ArrayValue(s0, s1, s2, s3, s4, s5, s6, s7, s8, s9, s10);
  return true;
}

template<typename T, int RANK>
bool IOVarASCII<T, RANK>::VarWrite(const Array<T, RANK>& val)
{
  bool mustReshape = false;
  for (int i = 0; i < RANK; i++)
    if (ArrayValue.extent(i) != val.extent(i))
      mustReshape = true;
  if (mustReshape)
    ArrayValue.resize(val.shape());
  ArrayValue = val;
  return true;
}

template<typename T, int RANK>
template<typename T0,
         typename T1,
         typename T2,
         typename T3,
         typename T4,
         typename T5,
         typename T6,
         typename T7,
         typename T8,
         typename T9,
         typename T10>
inline bool IOVarASCII<T, RANK>::VarWrite(
    const typename SliceInfo<T, T0, T1, T2, T3, T4, T5, T6, T7, T8, T9, T10>::T_slice& val,
    T0 s0,
    T1 s1,
    T2 s2,
    T3 s3,
    T4 s4,
    T5 s5,
    T6 s6,
    T7 s7,
    T8 s8,
    T8 s9,
    T10 s10)
{
  ArrayValue(s0, s1, s2, s3, s4, s5, s6, s7, s8, s9, s10) = val;
  return true;
}


template<class T, int RANK>
template<typename T0,
         typename T1,
         typename T2,
         typename T3,
         typename T4,
         typename T5,
         typename T6,
         typename T7,
         typename T8,
         typename T9,
         typename T10>
typename ASCIISliceMaker<T, RANK, T0, T1, T2, T3, T4, T5, T6, T7, T8, T9, T10>::SliceType IOVarASCII<T, RANK>::
    Slice(T0 s0, T1 s1, T2 s2, T3 s3, T4 s4, T5 s5, T6 s6, T7 s7, T8 s8, T9 s9, T10 s10)
{
  using newSliceType = typename ASCIISliceMaker<T, RANK, T0, T1, T2, T3, T4, T5, T6, T7, T8, T9, T10>::SliceType;
  newSliceType newVar(Name);
  newVar.ArrayValue.reference(ArrayValue(s0, s1, s2, s3, s4, s5, s6, s7, s8, s9, s10));
  return newVar;
}

template<class T, int RANK>
void IOVarASCII<T, RANK>::Print(std::ofstream& out)
{
  T a;
  if (GetRank() == 0)
  {
    out << TypeString(a) << " " << Name << " = " << std::endl;
  }
  else
  {
    out << "Array<" << TypeString(a) << "," << GetRank() << ">(";
    for (int i = 0; i < RANK; i++)
    {
      out << ArrayValue.extent(i);
      if (i < (RANK - 1))
        out << ",";
    }
    out << ") = [ ";
    for (int i = 0; i < ArrayValue.size(); i++)
    {
      out << *(ArrayValue.data() + i);
      if (i < (ArrayValue.size() - 1))
        out << ", ";
    }
    out << " ];" << std::endl;
  }
}

} // namespace IO

#endif
