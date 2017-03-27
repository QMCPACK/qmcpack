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

#ifndef PETE_PETE_SCALAR_H
#define PETE_PETE_SCALAR_H

///////////////////////////////////////////////////////////////////////////////
//
// WARNING: THIS FILE IS FOR INTERNAL PETE USE. DON'T INCLUDE IT YOURSELF
//
///////////////////////////////////////////////////////////////////////////////

//-----------------------------------------------------------------------------
//
// CLASS NAME
//    Scalar<T>
//
// DESCRIPTION
//    A wrapper around a scalar to be used in PETE expressions.
//
//-----------------------------------------------------------------------------

template<class T>
class Scalar
{
public:

  //---------------------------------------------------------------------------
  // Default constructor takes no action.

  inline
  Scalar() { }

  //---------------------------------------------------------------------------
  // Constructor from a single value.

  inline
  Scalar(const T &t) : scalar_m(t) { }

  template<class T1>
  inline
  explicit Scalar(const T1 &t) : scalar_m(t) { }

  //---------------------------------------------------------------------------
  // Constructor with arbitary second/third arguments, which is/are ignored.
  // Needed for compatibility with tree node constructors taking an
  // arbitrary argument.

  template<class Arg>
  inline
  Scalar(const Scalar<T> &s, const Arg &)
    : scalar_m(s.scalar_m) { }

  template<class Arg1, class Arg2>
  inline
  Scalar(const Scalar<T> &s, const Arg1 &, const Arg2 &)
    : scalar_m(s.scalar_m) { }

  //---------------------------------------------------------------------------
  // Copy constructor

  inline
  Scalar(const Scalar<T> &s) : scalar_m(s.scalar_m) { }

  //---------------------------------------------------------------------------
  // Return value.

  inline
  const T &value() const
  {
    return scalar_m;
  }

  //---------------------------------------------------------------------------
  // Assignment operators.

  inline
  Scalar<T> &operator=(const Scalar<T> &rhs)
  {
    scalar_m = rhs.scalar_m;
    return *this;
  }

  inline
  Scalar<T> &operator=(const T &rhs)
  {
    scalar_m = rhs;
    return *this;
  }

private:

  //---------------------------------------------------------------------------
  // The scalar value is stored here.

  T scalar_m;
};


#endif // PETE_PETE_SCALAR_H

// ACL:rcsinfo
