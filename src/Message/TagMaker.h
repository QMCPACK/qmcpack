//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@intel.com, Intel Inc.
//                    Jeremy McMinnis, jmcminis@gmail.com, Navar Inc.
//
// File created by: Jeongnim Kim, jeongnim.kim@intel.com, Intel Inc.
//////////////////////////////////////////////////////////////////////////////////////
    
    




#ifndef OHMMS_TAG_MAKER_H
#define OHMMS_TAG_MAKER_H
/*!\class TagMaker
 * \brief Assign a unique tag whenver TagMaker::TagMaker() is called.
 */
class TagMaker
{
public:

  TagMaker()
  {
    MyTag = (++CurrentTag);
  }
  ~TagMaker() {}
  int operator()()const
  {
    return MyTag;
  }
private:
  int MyTag;
  static int CurrentTag;
};
#endif // OHMMS_TAG_MAKER_H

/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$
 ***************************************************************************/
