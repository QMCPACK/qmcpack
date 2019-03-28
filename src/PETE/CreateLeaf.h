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

#ifndef PETE_PETE_CREATELEAF_H
#define PETE_PETE_CREATELEAF_H

namespace qmcplusplus
{
//-----------------------------------------------------------------------------
// Expression<T> - a class that wraps the contents of an expression
// (so that we don't need to define operators for all the tree types)
// Expression<T> provides the following interface:
// Expression<T>::Expression_t defines the expression type (T)
// const Expression_t& expression() returns a const ref to the
// expression object.
//-----------------------------------------------------------------------------

#ifdef PETE_USER_DEFINED_EXPRESSION

template<class T> class Expression;

#else

template<class T>
class Expression
{
public:
  // Type of the expression.

  typedef T Expression_t;

  // Construct from an expression.

  Expression(const T& expr) : expr_m(expr)
  { }

  // Accessor that returns the expression.

  const Expression_t& expression() const
  {
    return expr_m;
  }

private:
  // Store the expression by value since it is a temporary produced
  // by operator functions.

  T expr_m;
};

#endif // !PETE_USER_DEFINED_EXPRESSION

//-----------------------------------------------------------------------------
// CreateLeaf<T>
//
// The class CreateLeaf is used to tell PETE what to stick in expression trees
// for various objects that can appear in expressions.  Users MUST specialize
// CreateLeaf<T> for any container objects that they intend to use in
// expressions.  The typedef CreateLeaf<T>::Leaf_t is the type of object that
// appears in the expression tree (it could be T itself, or you might chose
// to just store iterators in the expression tree in which case you might
// define it to be T::const_iterator or something).  CreateLeaf also needs
// to provide a function make(const T&) which converts an object of type T
// into the object of type CreateLeaf<T>::Leaf_t which goes into the tree.
//-----------------------------------------------------------------------------

// The general case is assumed to be a scalar, since users need to specialize
// CreateLeaf for their container classes.

template<class T>
struct CreateLeaf
{
  typedef Scalar<T> Leaf_t;

  inline static
  Leaf_t make(const T &a)
  {
    return Scalar<T>(a);
  }
};

#ifndef PETE_USER_DEFINED_EXPRESSION

// For Expressions, we strip the Expression<> wrapper since it is intended
// to wrap the whole expression. (Expression<Scalar<>>+Expression<BinaryNode<>>
// becomes Expression<BinaryNode<OpAdd,Scalar<>,BinaryNode<>>>)

template<class T>
struct CreateLeaf<Expression<T> >
{
  typedef typename Expression<T>::Expression_t Leaf_t;

  inline static
  const Leaf_t &make(const Expression<T> &a)
  {
    return a.expression();
  }
};

#endif // !PETE_USER_DEFINED_EXPRESSION

//-----------------------------------------------------------------------------
// MakeReturn<T>
//
// MakeReturn is used to wrap expression objects (UnaryNode, BinaryNode etc.)
// inside an Expression<T> object.  Usually this indirection is unnecessary,
// but the indirection allows users to specify their own approach to storing
// trees.  By specializing MakeReturn<UnaryNode<>>, MakeReturn<BinaryNode<>>,
// etc. you could cause the expression trees to be stored in another format.
// For example, POOMA stores expressions inside Arrays, so the result of
// Array+Array is another Array.
//-----------------------------------------------------------------------------

template<class T>
struct MakeReturn
{
  typedef Expression<T> Expression_t;
  inline static
  Expression_t make(const T &a)
  {
    return Expression_t(a);
  }
};

}

#endif // PETE_PETE_CREATELEAF_H

// ACL:rcsinfo
