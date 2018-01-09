//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    


#include "Utilities/OhmmsObject.h"

/**@file OhmmsObject.cpp
 *@brief Definitions of the member functions of OhmmsObject
 */

// initialize the object counter
//int OhmmsObject::ObjectCounter = 0;

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
//  #pragma omp critical 
//  {
//    ObjectID = ObjectCounter;
//    ObjectCounter++;
//  }
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
//  #pragma omp critical 
//  {
//    ObjectID = ObjectCounter;
//    ObjectCounter++;
//  }
}

OhmmsObject::~OhmmsObject()
{
}

