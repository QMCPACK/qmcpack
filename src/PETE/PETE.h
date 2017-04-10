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

#ifndef PETE_PETE_PETE_H
#define PETE_PETE_PETE_H

///////////////////////////////////////////////////////////////////////////////
//
// This is the header file you should include if you want to use PETE, the
// Portable Expression Template Engine. You don't need add any member
// functions to make your container class PETE-aware, but you will need to
// create some traits classes and use the MakeOperator tool to create
// operator and math functions.
//
// See html/index.html for detailed instructions on using PETE.
//
///////////////////////////////////////////////////////////////////////////////

#if PETE_MAKE_EMPTY_CONSTRUCTORS

#define PETE_EMPTY_CONSTRUCTORS(CLASS)  \
  CLASS() { }   \
  CLASS(const CLASS &) { } \
  CLASS &operator=(const CLASS &) { return *this; }

#define PETE_EMPTY_CONSTRUCTORS_TEMPLATE(CLASS, ARG)  \
  CLASS() { }   \
  CLASS(const CLASS<ARG> &) { } \
  CLASS &operator=(const CLASS<ARG> &) { return *this; }

#else

#define PETE_EMPTY_CONSTRUCTORS(CLASS)
#define PETE_EMPTY_CONSTRUCTORS_TEMPLATE(CLASS, ARG)

#endif

#include "PETE/Scalar.h"
#include "PETE/TypeComputations.h"
#include "PETE/TreeNodes.h"
#include "PETE/OperatorTags.h"
#include "PETE/Functors.h"
#include "PETE/Combiners.h"
#include "PETE/ForEach.h"
#include "PETE/CreateLeaf.h"

// Some useful PETE definitions.

#define PETE_MAJOR_VERSION                 2
#define PETE_MINOR_VERSION                 1
#define PETE_PATCH_LEVEL                   1
#define PETE_VERSION_STRING                "PETE 2.1.1"
#define PETE_VERSION_NUM_STRING            "2.1.1"

#endif // PETE_PETE_PETE_H

// ACL:rcsinfo

