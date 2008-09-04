//////////////////////////////////////////////////////////////////
// (c) Copyright 1998-2002 by Jeongnim Kim
//
//   National Center for Supercomputing Applications &
//   Materials Computation Center
//   University of Illinois, Urbana-Champaign
//   Urbana, IL 61801
//   e-mail: jnkim@ncsa.uiuc.edu
//   Tel:    217-244-6319 (NCSA) 217-333-3324 (MCC)
//
// Supported by 
//   National Center for Supercomputing Applications, UIUC
//   Materials Computation Center, UIUC
//   Department of Physics, Ohio State University
//   Ohio Supercomputer Center
//////////////////////////////////////////////////////////////////
// -*- C++ -*-

#include "Utilities/OhmmsObject.h"
#include "Utilities/OhmmsInfo.h"

/**@file OhmmsObject.cpp
 *@brief Definitions of the member functions of OhmmsObject 
 */

// initialize the object counter
int OhmmsObject::ObjectCounter = 0;

/**default constructor
 *
 *Assign meaningless object and type names. Not givng a name to an object 
 *is OKAY as far as the object is not requested from a pool.
 */
OhmmsObject::OhmmsObject():
  OhmmsElementBase("none"), 
  TypeName("none"), 
  ElementByteSize(0)
{
  
  ObjectID = ObjectCounter;
  ObjectCounter++;
}

/**contructor
 *@param tname the name of the type
 *@param oname the name of the object
 *@brief Assign a unique ObjectID using the static data member ObjectCounter
 */
OhmmsObject::OhmmsObject(const std::string& tname, const std::string& oname):
  OhmmsElementBase(oname.c_str()), 
  TypeName(tname.c_str()), 
  ElementByteSize(0)
{
  
  ObjectID = ObjectCounter;
  ObjectCounter++;
}

OhmmsObject::~OhmmsObject() {
}

/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
