// -*- C++ -*-
// ACL:license
// ----------------------------------------------------------------------
// This software and ancillary information (herein called "SOFTWARE")
// called PETE (Portable Expression Template Engine) is
// made available under the terms described here.  The SOFTWARE has been
// approved for release with associated LA-CC Number LA-CC-99-5.
//
// Unless otherwise indicated, this SOFTWARE has been authored by an
// employee or employees of the University of California, operator of the
// Los Alamos National Laboratory under Contract No.  W-7405-ENG-36 with
// the U.S. Department of Energy.  The U.S. Government has rights to use,
// reproduce, and distribute this SOFTWARE. The public may copy, distribute,
// prepare derivative works and publicly display this SOFTWARE without
// charge, provided that this Notice and any statement of authorship are
// reproduced on all copies.  Neither the Government nor the University
// makes any warranty, express or implied, or assumes any liability or
// responsibility for the use of this SOFTWARE.
//
// If SOFTWARE is modified to produce derivative works, such modified
// SOFTWARE should be clearly marked, so as not to confuse it with the
// version available from LANL.
//
// For more information about PETE, send e-mail to pete@acl.lanl.gov,
// or visit the PETE web page at http://www.acl.lanl.gov/pete/.
// ----------------------------------------------------------------------
// ACL:license

#ifndef PETE_PETE_TYPECOMPUTATIONS_H
#define PETE_PETE_TYPECOMPUTATIONS_H

namespace qmcplusplus
{
//-----------------------------------------------------------------------------
//
// CLASS NAME
//    UnaryReturn<T, Op>
//
// DESCRIPTION
//    This template describes the default mechanism for calculating the
//    return type to a unary expression given the argument type T and
//    the operation type Op.
//
//    There are a several sensible things one can do:
//      o make the return type the same as the argument to the
//        function/operation. For example, operator-(T) should return a T.
//      o return a type based entirely on the operation.
//        For example, operator! always returns a bool.
//      o sythesize a type based on the type of the argument and the operation.
//    The first case is most common. We therefore make it the behavior
//    for the base template. The other cases are handled by partial
//    specialization.
//
//    For example, the abs function typically returns a double when the
//    argument is a std::complex<double>. The appropriate specialization here
//    would be:
//      template<> struct PETEUnaryReturn<std::complex<double>, FnAbs> {
//        typedef double Type_t;
//      };
//
//-----------------------------------------------------------------------------

template<class T, class Op>
struct UnaryReturn
{
  typedef T Type_t;
};

//-----------------------------------------------------------------------------
//
// CLASS NAME
//    Promote<T1, T2>
//
// DESCRIPTION
//    General template and specializations to implement C++ type promotion
//    for basic types.
//
//-----------------------------------------------------------------------------

// Base template: don't do anything by default.

template<class T1, class T2>
struct Promote
{
  typedef T1 Type_t;
};

// bool

template<>
struct Promote<bool, bool>
{
  typedef bool Type_t;
};

template<>
struct Promote<bool, char>
{
  typedef char Type_t;
};

template<>
struct Promote<bool, short>
{
  typedef short Type_t;
};

template<>
struct Promote<bool, int>
{
  typedef int Type_t;
};

template<>
struct Promote<bool, long>
{
  typedef long Type_t;
};

template<>
struct Promote<bool, float>
{
  typedef float Type_t;
};

template<>
struct Promote<bool, double>
{
  typedef double Type_t;
};

// char

template<>
struct Promote<char, bool>
{
  typedef char Type_t;
};

template<>
struct Promote<char, char>
{
  typedef char Type_t;
};

template<>
struct Promote<char, short>
{
  typedef short Type_t;
};

template<>
struct Promote<char, int>
{
  typedef int Type_t;
};

template<>
struct Promote<char, long>
{
  typedef long Type_t;
};

template<>
struct Promote<char, float>
{
  typedef float Type_t;
};

template<>
struct Promote<char, double>
{
  typedef double Type_t;
};

// short

template<>
struct Promote<short, bool>
{
  typedef short Type_t;
};

template<>
struct Promote<short, char>
{
  typedef short Type_t;
};

template<>
struct Promote<short, short>
{
  typedef short Type_t;
};

template<>
struct Promote<short, int>
{
  typedef int Type_t;
};

template<>
struct Promote<short, long>
{
  typedef long Type_t;
};

template<>
struct Promote<short, float>
{
  typedef float Type_t;
};

template<>
struct Promote<short, double>
{
  typedef double Type_t;
};

// int

template<>
struct Promote<int, bool>
{
  typedef int Type_t;
};

template<>
struct Promote<int, char>
{
  typedef int Type_t;
};

template<>
struct Promote<int, short>
{
  typedef int Type_t;
};

template<>
struct Promote<int, int>
{
  typedef int Type_t;
};

template<>
struct Promote<int, long>
{
  typedef long Type_t;
};

template<>
struct Promote<int, float>
{
  typedef float Type_t;
};

template<>
struct Promote<int, double>
{
  typedef double Type_t;
};

// long

template<>
struct Promote<long, bool>
{
  typedef long Type_t;
};

template<>
struct Promote<long, char>
{
  typedef long Type_t;
};

template<>
struct Promote<long, short>
{
  typedef long Type_t;
};

template<>
struct Promote<long, int>
{
  typedef long Type_t;
};

template<>
struct Promote<long, long>
{
  typedef long Type_t;
};

template<>
struct Promote<long, float>
{
  typedef float Type_t;
};

template<>
struct Promote<long, double>
{
  typedef double Type_t;
};

// float

template<>
struct Promote<float, bool>
{
  typedef float Type_t;
};

template<>
struct Promote<float, char>
{
  typedef float Type_t;
};

template<>
struct Promote<float, short>
{
  typedef float Type_t;
};

template<>
struct Promote<float, int>
{
  typedef float Type_t;
};

template<>
struct Promote<float, long>
{
  typedef float Type_t;
};

template<>
struct Promote<float, float>
{
  typedef float Type_t;
};

template<>
struct Promote<float, double>
{
  typedef double Type_t;
};

// double

template<>
struct Promote<double, bool>
{
  typedef double Type_t;
};

template<>
struct Promote<double, char>
{
  typedef double Type_t;
};

template<>
struct Promote<double, short>
{
  typedef double Type_t;
};

template<>
struct Promote<double, int>
{
  typedef double Type_t;
};

template<>
struct Promote<double, long>
{
  typedef double Type_t;
};

template<>
struct Promote<double, float>
{
  typedef double Type_t;
};

template<>
struct Promote<double, double>
{
  typedef double Type_t;
};

//-----------------------------------------------------------------------------
//
// CLASS NAME
//    BinaryReturn<T1, T2, Op>
//
// DESCRIPTION
//    This template describes the default mechanism for calculating the
//    return type to a binary expression given the argument types T1 and
//    T2 and the operation type Op.
//
//    There are several sensible things one can do:
//      o make the return type by promoting/converting the "simpler" type
//        into the more "complex." For example, we typically want to do
//        this with addition.
//      o return the type of the left-hand operand. For example, this is
//        what happens with operator<<.
//      o return the type of the right-hand operand.
//      o return a type based entirely on the operation.
//        For example, operator!= always returns a bool.
//      o synthesize the return type based on the operation and types.
//    The first option is most common, so we make that the behavior of the
//    base template. The other cases are handled by partial specialization.
//
//    For example, the multiplication between a matrix and a vector might do a
//    matrix/vector product, thereby returning a vector. The appropriate
//    specialization here would be:
//      struct BinaryReturn<Mat<double,3>, Vec<float,3>, OpMultiply> {
//        typedef Vec<double,3> Type_t;
//      };
//    Notice how the element type is promoted.
//
//-----------------------------------------------------------------------------

template<class T1, class T2, class Op>
struct BinaryReturn
{
  typedef typename Promote<T1, T2>::Type_t Type_t;
};


//-----------------------------------------------------------------------------
//
// CLASS NAME
//    TrinaryReturn<T1, T2, T3, Op>
//
// DESCRIPTION
//    This template describes the default mechanism for calculating the
//    return type to a trinary expression given the argument types T1, T2, and
//    T3 and the operation type Op. The only trinary expression supported
//    in C++ is the ?: operation. In this case, T1 should end up being bool
//    and the result of the calculation is of type Binary_Promotion(T2,T3)
//    with the value being that associated with T2 if T1's associated value
//    turns out to be true and T3 if T1's associated value turns out to be
//    false.
//
//-----------------------------------------------------------------------------

template<class T1, class T2, class T3, class Op>
struct TrinaryReturn
{
  typedef typename BinaryReturn<T2, T3, Op>::Type_t Type_t;
};

}

#endif // PETE_PETE_TYPE_COMPUTATIONS_H

// ACL:rcsinfo
