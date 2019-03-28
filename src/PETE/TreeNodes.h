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

#ifndef PETE_PETE_TREENODES_H
#define PETE_PETE_TREENODES_H

namespace qmcplusplus
{
///////////////////////////////////////////////////////////////////////////////
//
// WARNING: THIS FILE IS FOR INTERNAL PETE USE. DON'T INCLUDE IT YOURSELF
//
///////////////////////////////////////////////////////////////////////////////

//-----------------------------------------------------------------------------
//
// CLASS NAME
//   Reference<T>
//
// DESCRIPTION
//   Reference is a special kind of node that contains a reference to an object
//   of type T.  It can be converted to a (const T &), and other tree objects
//   will perform this conversion before returning their elements.
//
//-----------------------------------------------------------------------------

template<class T>
struct Reference
{
  //---------------------------------------------------------------------------
  // Export the type of thing we're referencing.

  typedef T Type_t;

  //---------------------------------------------------------------------------
  // Reference can be created from a const ref.

  inline
  Reference(const T &reference)
    : reference_m(reference)
  { }

  //---------------------------------------------------------------------------
  // Copy constructor

  inline
  Reference(const Reference<T> &model)
    : reference_m(model.reference())
  { }

  //---------------------------------------------------------------------------
  // Reference can be converted to a const ref

  inline
  const T &reference() const
  {
    return reference_m;
  }

  //---------------------------------------------------------------------------
  // Conversion operators.

  operator const T& () const
  {
    return reference_m;
  }
  operator T& () const
  {
    return const_cast<T&>(reference_m);
  }

  const T &reference_m;
};

//-----------------------------------------------------------------------------
//
// CLASS NAME
//   DeReference<T>
//
// DESCRIPTION
//   DeReference is a simple traits class that unwraps the Reference struct.
//   If T is not a reference object then DeReference gives (const T &).
//   If T is a Reference object, then DeReference gives a const ref to the
//   wrapped object.
//
//-----------------------------------------------------------------------------

template<class T>
struct DeReference
{
  typedef const T &Return_t;
  typedef T Type_t;
  static inline Return_t apply(const T &a)
  {
    return a;
  }
};

template<class T>
struct DeReference<Reference<T> >
{
  typedef const T &Return_t;
  typedef T Type_t;
  static inline Return_t apply(const Reference<T> &a)
  {
    return a.reference();
  }
};

//-----------------------------------------------------------------------------
//
// CLASS NAME
//   UnaryNode<Op, Child>
//
// DESCRIPTION
//   A tree node for representing unary expressions. The node holds a
//   child (of type Child), which is the type of the expression sub tree and a
//   an operation (of type Op), which is typically the operation applied to
//   the sub tree.
//
//-----------------------------------------------------------------------------

template<class Op, class Child>
class UnaryNode
{
public:

  //---------------------------------------------------------------------------
  // Accessors making the operation and child available to the outside.

  inline
  const Op &operation() const
  {
    return op_m;
  }

  inline
  typename DeReference<Child>::Return_t
  child() const
  {
    return DeReference<Child>::apply(child_m);
  }

  //---------------------------------------------------------------------------
  // Constructor using both a operation and the child.

  inline
  UnaryNode(const Op &o, const Child &c)
    : op_m(o), child_m(c) { }

  //---------------------------------------------------------------------------
  // Constructor using just the child.

  inline
  UnaryNode(const Child &c)
    : child_m(c) { }

  //---------------------------------------------------------------------------
  // Copy constructor.

  inline
  UnaryNode(const UnaryNode<Op, Child> &t)
    : op_m(t.operation()), child_m(t.child()) { }

  //---------------------------------------------------------------------------
  // Constructor using a UnaryNode with a different child and/or a different
  // storage tag. Note: for this to work, a Child must be constructable
  // from an OtherChild.

  template<class OtherChild>
  inline
  UnaryNode(const UnaryNode<Op, OtherChild> &t)
    : op_m(t.operation()), child_m(t.child()) { }

  //---------------------------------------------------------------------------
  // Constructor using a UnaryNode with a different child,
  // some arbitrary argument, and a different storage tag.
  // Note: for this to work, a Child must be constructable
  // from an OtherChild and an Arg.

  template<class OtherChild, class Arg>
  inline
  UnaryNode(const UnaryNode<Op, OtherChild> &t, const Arg &a)
    : op_m(t.operation()), child_m(t.child(), a) { }

  //---------------------------------------------------------------------------
  // Constructor using a BinaryNode with a different Child and
  // two arbitrary arguments.
  // Note: for this to work, a Child  must be constructable
  // from an OtherChild and an Arg1 & Arg2.

  template<class OtherChild, class Arg1, class Arg2>
  inline
  UnaryNode(const UnaryNode<Op, OtherChild> &t,
            const Arg1 &a1, const Arg2 &a2)
    : op_m(t.operation()), child_m(t.child(), a1, a2)
  { }

private:

  Op    op_m;
  Child child_m;

};


//-----------------------------------------------------------------------------
//
// CLASS NAME
//   BinaryNode<Op, Left, Right>
//
// DESCRIPTION
//   A tree node for representing binary expressions. The node holds a
//   left child (of type Left), which is the type of the LHS expression
//   sub tree, a right child (of type Right), which is the type of the RHS
//   expression sub tree, and an operation (of type OP), which is applied
//   to the two sub trees.
//
//-----------------------------------------------------------------------------

template<class Op, class Left, class Right>
class BinaryNode
{
public:

  //---------------------------------------------------------------------------
  // Accessors making the operation and children available to the outside.

  inline
  const Op &operation() const
  {
    return op_m;
  }

  inline
  typename DeReference<Left>::Return_t
  left() const
  {
    return DeReference<Left>::apply(left_m);
  }

  inline
  typename DeReference<Right>::Return_t
  right() const
  {
    return DeReference<Right>::apply(right_m);
  }

  //---------------------------------------------------------------------------
  // Constructor using both the operation and the two children.

  inline
  BinaryNode(const Op &o, const Left &l, const Right &r)
    : op_m(o), left_m(l), right_m(r)
  { }

  //---------------------------------------------------------------------------
  // Constructor using just the two children.

  inline
  BinaryNode(const Left &l, const Right &r)
    : left_m(l), right_m(r)
  { }

  //---------------------------------------------------------------------------
  // Copy constructor.

  inline
  BinaryNode(const BinaryNode<Op, Left, Right> &t)
    : op_m(t.operation()), left_m(t.left()), right_m(t.right())
  { }

  //---------------------------------------------------------------------------
  // Constructor using a BinaryNode with a different Left/Right.
  // Note: for this to work, the Left/Right must be constructable
  // from an OtherLeft/OtherRight.

  template<class OtherLeft, class OtherRight>
  inline
  BinaryNode(const BinaryNode<Op, OtherLeft, OtherRight> &t)
    : op_m(t.operation()), left_m(t.left()), right_m(t.right())
  { }

  //---------------------------------------------------------------------------
  // Constructor using a BinaryNode with a different Left/Right and
  // some arbitrary argument.
  // Note: for this to work, a Left/Right must be constructable
  // from an OtherLeft/OtherRight and an Arg.

  template<class OtherLeft, class OtherRight, class Arg>
  inline
  BinaryNode(const BinaryNode<Op, OtherLeft, OtherRight> &t,
             const Arg &a)
    : op_m(t.operation()), left_m(t.left(), a), right_m(t.right(), a)
  { }

  //---------------------------------------------------------------------------
  // Constructor using a BinaryNode with a different Left/Right and
  // two arbitrary arguments.
  // Note: for this to work, a Left/Right must be constructable
  // from an OtherLeft/OtherRight and an Arg1 & Arg2.

  template<class OtherLeft, class OtherRight, class Arg1, class Arg2>
  inline
  BinaryNode(const BinaryNode<Op, OtherLeft, OtherRight> &t,
             const Arg1 &a1, const Arg2 &a2)
    : op_m(t.operation()),
      left_m(t.left(), a1, a2), right_m(t.right(), a1, a2)
  { }

private:

  //---------------------------------------------------------------------------
  // The operation and left/right sub expressions stored in this node of the
  // tree.

  Op    op_m;
  Left  left_m;
  Right right_m;

};


//-----------------------------------------------------------------------------
//
// CLASS NAME
//   TrinaryNode<Op, Left, Middle, Right>
//
// DESCRIPTION
//   A tree node for representing trinary expressions. The node holds a
//   Left child (of type Left), which is the type of the LHS expression
//   sub tree (typically a comparison operation); a Middle child (of type
//   Middle), which is the type of the middle (true branch) expression
//   sub tree; a Right child (of type Right), which is the type of
//   the expression (false branch) sub tree; and an operation (of type Op),
//   which is applied to the three sub trees.
//
//-----------------------------------------------------------------------------

template< class Op, class Left, class Middle, class Right>
class TrinaryNode
{
public:

  //---------------------------------------------------------------------------
  // Accessors making the operation and children available to the outside.

  inline
  const Op &operation() const
  {
    return op_m;
  }

  inline
  typename DeReference<Left>::Return_t
  left() const
  {
    return DeReference<Left>::apply(left_m);
  }

  inline
  typename DeReference<Right>::Return_t
  right() const
  {
    return DeReference<Right>::apply(right_m);
  }

  inline
  typename DeReference<Middle>::Return_t
  middle() const
  {
    return DeReference<Middle>::apply(middle_m);
  }

  //---------------------------------------------------------------------------
  // Constructor using the operation and three children.

  inline
  TrinaryNode(const Op &o, const Left &l, const Middle &m, const Right &r)
    : op_m(o), left_m(l), middle_m(m), right_m(r)
  { }

  //---------------------------------------------------------------------------
  // Constructor with just the three children.

  inline
  TrinaryNode(const Left &l, const Middle &m, const Right &r)
    : left_m(l), middle_m(m), right_m(r)
  { }

  //---------------------------------------------------------------------------
  // Copy constructor.

  inline
  TrinaryNode(const TrinaryNode<Op, Left, Middle, Right> &t)
    : op_m(t.operation()), left_m(t.left()), middle_m(t.middle()),
      right_m(t.right())
  { }

  //---------------------------------------------------------------------------
  // Constructor using a TrinaryNode with a different Left/Middle/Right.
  // Note: for this to work, the Left/Middle/Right must be constructable
  // from an OtherLeft/OtherMiddle/OtherRight.

  template<class OtherLeft, class OtherMiddle, class OtherRight>
  inline
  TrinaryNode(const TrinaryNode<Op, OtherLeft, OtherMiddle, OtherRight> & t)
    : op_m(t.operation()), left_m(t.left()), middle_m(t.middle()),
      right_m(t.right())
  { }

  //---------------------------------------------------------------------------
  // Constructor using a TrinaryNode with a different Left/Middle/Right and
  // some arbitrary argument.
  // Note: for this to work, a Left/Middle/Right must be constructable
  // from an OtherLeft/OtherMiddle/OtherRight and an Arg.

  template<class OtherLeft, class OtherMiddle, class OtherRight, class Arg>
  inline
  TrinaryNode(const TrinaryNode<Op, OtherLeft, OtherMiddle, OtherRight> &t,
              const Arg &a)
    : op_m(t.operation()), left_m(t.left(), a), middle_m(t.middle(), a),
      right_m(t.right(), a)
  { }

  //---------------------------------------------------------------------------
  // Constructor using a TrinaryNode with a different Left/Middle/Right and
  // two arbitrary arguments.
  // Note: for this to work, a Left/Middle/Right must be constructable
  // from an OtherLeft/OtherMiddle/OtherRight and an Arg1 & Arg2.

  template<class OtherLeft, class OtherMiddle, class OtherRight,
           class Arg1, class Arg2>
  inline
  TrinaryNode(const TrinaryNode<Op, OtherLeft, OtherMiddle, OtherRight> &t,
              const Arg1 &a1, const Arg2 &a2)
    : op_m(t.operation()), left_m(t.left(), a1, a2),
      middle_m(t.middle(), a1, a2) , right_m(t.right(), a1, a2)
  { }

private:

  //---------------------------------------------------------------------------
  // The operation and left, right, and middle sub trees stored at this node.

  Op     op_m;
  Left   left_m;
  Middle middle_m;
  Right  right_m;

};

}
#endif // PETE_PETE_TREENODES_H

// ACL:rcsinfo
