//////////////////////////////////////////////////////////////////
// (c) Copyright 1998- by Jeongnim Kim
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//   National Center for Supercomputing Applications &
//   Materials Computation Center
//   University of Illinois, Urbana-Champaign
//   Urbana, IL 61801
//   e-mail: jnkim@ncsa.uiuc.edu
//
// Supported by
//   National Center for Supercomputing Applications, UIUC
//   Materials Computation Center, UIUC
//   Department of Physics, Ohio State University
//   Ohio Supercomputer Center
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
#ifndef OHMMS_OHMMSOBJECT_H
#define OHMMS_OHMMSOBJECT_H
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#include "OhmmsData/OhmmsElementBase.h"
#include <string>

/**@file OhmmsObject.h
 *@brief Declaration of OhmmsObject
 */

/** An OhmmsElementBase with extra members to keep track the object/type.
 *
 *Additional data members of OhmmsObject are
 * - TypeName: the name of a type, e.g., PositionType
 * - ObjectID: an ID for this object assigned by a constructor
 * - ElementByteSize: the byte size of an object. A derived class specifies the value.
 * - ObjectCounter: static data member. Incremented whenever an OhmmsObject is created.
 *
 *It is common for class libraries or applications to have a base
 *class that has an object counter for garbage collections and other
 *functions.  This base class is not that smart but enables keeping
 *track the list unique for a process and makes it easy to retrieve an
 *object with a name.
 *
 *@todo This class should be integrated into OhmmsElementBase or
 *an equivalent base class.
 */
class OhmmsObject: public OhmmsElementBase
{
public:

  OhmmsObject();
  OhmmsObject(const std::string& tname, const std::string& oname);
  virtual ~OhmmsObject();

  ///return the ObjectID which is unique to each object
  int id() const
  {
    return ObjectID;
  }

  ///overwrite ObjectID
  void setID(int i)
  {
    ObjectID=i;
  }

  ///returns the byte size of the object
  int elementByteSize() const
  {
    return ElementByteSize;
  }

  ///set the type name
  void setTypeName(const std::string& tname)
  {
    TypeName = tname;
  }

  ///set the object name
  void setObjName(const std::string& oname)
  {
    myName = oname;
  }

  ///return the type name
  const std::string& typeName() const
  {
    return TypeName;
  }

  ///return the object name
  const std::string& objName() const
  {
    return myName;
  }

protected:

  ///the type name of this object
  std::string TypeName;

  ///the unique ID of this object
  int ObjectID;

  ///the byte size of this object
  int ElementByteSize;

private:
  static int ObjectCounter;
};
#endif

/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$
 ***************************************************************************/
