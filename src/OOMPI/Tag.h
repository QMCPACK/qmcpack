// -*- c++ -*-
//
// Copyright (c) 2002-2003 Indiana University.  All rights reserved.
// Copyright (c) 1996, 1997, 1998, 2000 University of Notre Dame.
//                         All rights reserved.
//
// This file is part of the OOMPI software package.  For license
// information, see the LICENSE file in the top level directory of the
// OOMPI source distribution.
//
// $Id$
//
// OOMPI Class library
// Tag base class
//

#ifndef _OOMPI_TAG_H_
#define _OOMPI_TAG_H_


class OOMPI_Tag
{
public:

  //
  // Note that copy constructor and assignment operator are
  // both deep copies
  //

  inline OOMPI_Tag(int tag = OOMPI_NO_TAG)
    : my_tag(tag)
  {};

  //
  // Access functions to get/set tag for this instance
  //

  inline int Get_tag(void) const
  {
    return my_tag;
  }
  inline void Set_tag(int tag)
  {
    my_tag = tag;
  }

protected:
  int my_tag;

private:
};

#endif




